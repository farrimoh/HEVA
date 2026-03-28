#!/usr/bin/env python3
"""Optimize core epsilonLJ and sigmaLJ on top of a fixed capsid-only baseline."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent
if str(WORKFLOW_DIR) not in sys.path:
    sys.path.insert(0, str(WORKFLOW_DIR))

from assembly_opt.config import (
    CORE_PARAMETER_ORDER,
    DEFAULT_CORE_INITIAL_PARAMS,
    DEFAULT_CORE_PARAM_BOUNDS,
    build_config,
)
from assembly_opt.dependencies import require_runtime_dependencies
from assembly_opt.metrics import classify_outcome, summarize_outcomes, write_seed_outcomes, write_summary_json
from assembly_opt.trial import TrialExecution, TrialSpec, execute_trial
from assembly_opt.verify import run_static_verification


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Optimize core epsilonLJ and sigmaLJ using fixed capsid-only baseline parameters.")
    parser.add_argument("--skip-verify", action="store_true", help="Skip the static verification pass.")
    parser.add_argument("--verify-only", action="store_true", help="Run only the static verification pass.")
    parser.add_argument("--max-evals", type=int, default=8, help="Maximum scipy iterations.")
    parser.add_argument("--method", choices=("nelder-mead", "L-BFGS-B"), default="nelder-mead", help="scipy.optimize.minimize method.")
    parser.add_argument("--single-run", nargs=2, type=float, metavar=("EPSILONLJ", "SIGMALJ"), help="Evaluate one fixed core parameter set instead of running scipy minimization.")
    parser.add_argument("--runs-dir", default=None, help="Output directory for optimization runs.")
    parser.add_argument("--base-config", default=None, help="Override the workflow base config.")
    parser.add_argument("--assemble", default=None, help="Override the HEVA assemble binary path.")
    parser.add_argument("--init-mode", choices=("restart", "seed", "triangle", "pentamer", "hexamer"), default="pentamer", help="Initialization mode for each seed run.")
    parser.add_argument("--seed-shape", choices=("triangle", "pentamer", "hexamer"), default=None, help="Seed shape when --init-mode=seed.")
    parser.add_argument("--initial-structure", default=None, help="Restart file path when --init-mode=restart.")
    parser.add_argument("--max-sweeps", type=int, default=50000000, help="Maximum sweeps per seed run. Use 0 for no workflow-imposed limit.")
    parser.add_argument("--runs-per-eval", type=int, default=4, help="Number of HEVA seeds per core-parameter evaluation.")
    parser.add_argument("--max-parallel", type=int, default=4, help="Maximum concurrent local HEVA runs.")
    parser.add_argument("--base-seed", type=int, default=None, help="Override the first seed used by a new run directory.")
    parser.add_argument("--profile", choices=("test", "run", "extended"), default="run", help="HEVA engine profile.")
    parser.add_argument("--index-capacity", type=int, default=None, help="Required only for profile=extended.")
    parser.add_argument("--core-max-bonds", type=int, default=1000, help="Legacy maxbondsRNA control.")
    parser.add_argument("--dmud", "--drug-dmud", "--core-dmud", dest="dmud", type=float, default=0.0, help="Drug dmud value (legacy --core-dmud alias accepted).")
    parser.add_argument("--epsilonlj-init", type=float, default=DEFAULT_CORE_INITIAL_PARAMS["epsilonLJ"])
    parser.add_argument("--sigmalj-init", type=float, default=DEFAULT_CORE_INITIAL_PARAMS["sigmaLJ"])
    parser.add_argument("--epsilonlj-min", type=float, default=DEFAULT_CORE_PARAM_BOUNDS["epsilonLJ"][0])
    parser.add_argument("--epsilonlj-max", type=float, default=DEFAULT_CORE_PARAM_BOUNDS["epsilonLJ"][1])
    parser.add_argument("--sigmalj-min", type=float, default=DEFAULT_CORE_PARAM_BOUNDS["sigmaLJ"][0])
    parser.add_argument("--sigmalj-max", type=float, default=DEFAULT_CORE_PARAM_BOUNDS["sigmaLJ"][1])
    parser.add_argument("--target-yield", type=float, default=1.0)
    parser.add_argument("--target-ratio", type=float, default=0.9)
    parser.add_argument("--yield-weight", type=float, default=100.0)
    parser.add_argument("--ratio-weight", type=float, default=100.0)
    parser.add_argument("--offtarget-weight", type=float, default=10000.0)
    parser.add_argument("--failure-penalty", type=float, default=0.0)

    parser.add_argument("--epsilon0", type=float, default=2500.0)
    parser.add_argument("--kappa0", type=float, default=100.0)
    parser.add_argument("--kappaPhi0", type=float, default=1000.0)
    parser.add_argument("--theta0", type=float, default=0.240)
    parser.add_argument("--theta1", type=float, default=0.480)
    parser.add_argument("--ks0", type=float, default=0.002)
    parser.add_argument("--simulation-dg", type=float, default=0.100)
    parser.add_argument("--mudrug", type=float, default=0.0)
    parser.add_argument("--gdrug0", type=float, default=0.0)
    parser.add_argument("--kd0", type=float, default=0.0)

    parser.add_argument("--muCd", type=float, default=-7.0)
    parser.add_argument("--dmu", type=float, default=-3.0)
    parser.add_argument("--gb0", type=float, default=-6.8)
    parser.add_argument("--dg12", type=float, default=0.3)
    parser.add_argument("--dg01", type=float, default=0.1)
    parser.add_argument("--dg20", type=float, default=-0.1)
    parser.add_argument("--dg33", type=float, default=0.0)
    parser.add_argument("--dg00", type=float, default=-0.7)
    parser.add_argument("--dgother", type=float, default=-0.8)
    return parser.parse_args()


def log_verification(config) -> bool:
    ok = True
    for item in run_static_verification(config):
        log = getattr(logging, item.level, logging.info)
        log("Verification %s: %s", item.label, item.detail)
        if item.level == "error":
            ok = False
    return ok


def initialize_run_folder(run_dir: Path, base_seed: int) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    seed_path = run_dir / "next_seed.txt"
    if not seed_path.exists():
        seed_path.write_text(f"{base_seed}\n", encoding="utf-8")
    eval_path = run_dir / "next_eval.txt"
    if not eval_path.exists():
        eval_path.write_text("0\n", encoding="utf-8")
    record_path = run_dir / "submitted.csv"
    if not record_path.exists():
        with record_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "evaluation",
                    "epsilonLJ",
                    "sigmaLJ",
                    "objective",
                    "total_runs",
                    "failed_runs",
                    "closed_count",
                    "t4_count",
                    "t3_count",
                    "other_closed",
                    "target_yield",
                    "t4_ratio",
                    "evaluation_dir",
                ]
            )


def allocate_evaluation(run_dir: Path, runs_per_eval: int) -> tuple[int, list[int], Path]:
    eval_path = run_dir / "next_eval.txt"
    seed_path = run_dir / "next_seed.txt"
    eval_index = int(eval_path.read_text(encoding="utf-8").strip())
    next_seed = int(seed_path.read_text(encoding="utf-8").strip())
    seeds = list(range(next_seed, next_seed + runs_per_eval))
    evaluation_dir = run_dir / f"eval_{eval_index:04d}"
    eval_path.write_text(f"{eval_index + 1}\n", encoding="utf-8")
    seed_path.write_text(f"{next_seed + runs_per_eval}\n", encoding="utf-8")
    return eval_index, seeds, evaluation_dir


def write_commands(path: Path, executions: list[TrialExecution]) -> None:
    payload = {
        "commands": [
            {
                "run_dir": str(execution.run_dir),
                "config": str(execution.generated_config_path),
                "command": execution.command,
                "returncode": execution.returncode,
            }
            for execution in executions
        ]
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def append_evaluation_record(run_dir: Path, evaluation_index: int, summary, evaluation_dir: Path, epsilon_lj: float, sigma_lj: float) -> None:
    with (run_dir / "submitted.csv").open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                evaluation_index,
                f"{epsilon_lj:.6f}",
                f"{sigma_lj:.6f}",
                f"{summary.objective:.6f}",
                summary.total_runs,
                summary.failed_runs,
                summary.closed_count,
                summary.t4_count,
                summary.t3_count,
                summary.other_closed,
                f"{summary.target_yield:.6f}",
                f"{summary.t4_ratio:.6f}",
                evaluation_dir,
            ]
        )


def build_trial_spec(config, args, seed: int, output_dir: Path, epsilon_lj: float, sigma_lj: float) -> TrialSpec:
    return TrialSpec(
        workflow_dir=config.base_dir,
        assemble_path=config.resolve_assemble_path(),
        base_config_path=config.resolve_base_config_path(),
        seed=seed,
        init_mode=config.init_mode,
        seed_shape=config.seed_shape,
        initial_structure_path=config.resolve_initial_structure_path(),
        output_dir=output_dir,
        max_sweeps=config.max_sweeps,
        profile=config.profile,
        index_capacity=config.index_capacity,
        epsilon0=config.epsilon0,
        kappa0=config.kappa0,
        kappaPhi0=config.kappaPhi0,
        theta0=config.theta0,
        theta1=config.theta1,
        gb0=args.gb0,
        muCd=args.muCd,
        ks0=config.ks0,
        dmu=args.dmu,
        simulation_dg=config.simulation_dg,
        mudrug=config.mudrug,
        gdrug0=config.gdrug0,
        kd0=config.kd0,
        dmud=config.dmud,
        core_enabled=True,
        core_max_bonds=config.core_max_bonds,
        core_epsilon_lj=epsilon_lj,
        core_sigma_lj=sigma_lj,
        dg12=args.dg12,
        dg01=args.dg01,
        dg20=args.dg20,
        dg33=args.dg33,
        dg00=args.dg00,
        dgother=args.dgother,
    )


def evaluate_core_parameter_set(epsilon_lj: float, sigma_lj: float, config, args, run_dir: Path):
    eval_index, seeds, evaluation_dir = allocate_evaluation(run_dir, config.runs_per_eval)
    evaluation_dir.mkdir(parents=True, exist_ok=True)
    specs = [
        build_trial_spec(config, args, seed, evaluation_dir / f"seed-{seed}", epsilon_lj, sigma_lj)
        for seed in seeds
    ]

    executions: list[TrialExecution] = []
    outcomes = []
    if config.max_parallel > 1:
        with ThreadPoolExecutor(max_workers=config.max_parallel) as executor:
            futures = {executor.submit(execute_trial, spec, dry_run=False): spec.seed for spec in specs}
            for future in as_completed(futures):
                execution = future.result()
                outcome = classify_outcome(futures[future], execution.run_dir, execution.stdout, execution.returncode)
                executions.append(execution)
                outcomes.append(outcome)
    else:
        for spec in specs:
            execution = execute_trial(spec, dry_run=False)
            executions.append(execution)
            outcomes.append(classify_outcome(spec.seed, execution.run_dir, execution.stdout, execution.returncode))

    outcomes.sort(key=lambda item: item.seed)
    executions.sort(key=lambda item: item.run_dir.name)
    summary = summarize_outcomes(outcomes, config)
    params = {
        "epsilonLJ": float(epsilon_lj),
        "sigmaLJ": float(sigma_lj),
        "muCd": float(args.muCd),
        "dmu": float(args.dmu),
        "gb0": float(args.gb0),
        "dg12": float(args.dg12),
        "dg01": float(args.dg01),
        "dg20": float(args.dg20),
        "dg33": float(args.dg33),
        "dg00": float(args.dg00),
        "dgother": float(args.dgother),
    }
    write_seed_outcomes(evaluation_dir / "seed_outcomes.csv", outcomes)
    write_summary_json(evaluation_dir / "summary.json", summary, params, evaluation_dir)
    write_commands(evaluation_dir / "commands.json", executions)
    append_evaluation_record(run_dir, eval_index, summary, evaluation_dir, epsilon_lj, sigma_lj)
    logging.info(
        "Evaluation %s epsilonLJ %.4f sigmaLJ %.4f objective %.6f yield %.3f ratio %.3f closed %s/%s",
        eval_index,
        epsilon_lj,
        sigma_lj,
        summary.objective,
        summary.target_yield,
        summary.t4_ratio,
        summary.closed_count,
        summary.total_runs,
    )
    return summary


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    config = build_config(
        WORKFLOW_DIR,
        runs_dir=args.runs_dir,
        base_config=args.base_config,
        assemble=args.assemble,
        init_mode=args.init_mode,
        seed_shape=args.seed_shape,
        initial_structure=args.initial_structure,
        max_evals=args.max_evals,
        max_sweeps=args.max_sweeps,
        runs_per_eval=args.runs_per_eval,
        max_parallel=args.max_parallel,
        base_seed=args.base_seed if args.base_seed is not None else int(time.strftime("%H%M%S")),
        profile=args.profile,
        index_capacity=args.index_capacity,
        epsilon0=args.epsilon0,
        kappa0=args.kappa0,
        kappaPhi0=args.kappaPhi0,
        theta0=args.theta0,
        theta1=args.theta1,
        ks0=args.ks0,
        simulation_dg=args.simulation_dg,
        mudrug=args.mudrug,
        gdrug0=args.gdrug0,
        kd0=args.kd0,
        dmud=args.dmud,
        core_enabled=True,
        core_max_bonds=args.core_max_bonds,
        core_epsilon_lj=args.epsilonlj_init,
        core_sigma_lj=args.sigmalj_init,
        target_yield=args.target_yield,
        target_ratio=args.target_ratio,
        yield_weight=args.yield_weight,
        ratio_weight=args.ratio_weight,
        offtarget_weight=args.offtarget_weight,
        failure_penalty=args.failure_penalty,
    )

    if not args.skip_verify:
        if not log_verification(config):
            return 1
        if args.verify_only:
            return 0

    run_dir = config.resolve_runs_dir() / time.strftime("%Y%m%d-%H%M%S_core_opt")
    initialize_run_folder(run_dir, config.base_seed)

    def clamp(x):
        epsilon_lj = min(max(float(x[0]), args.epsilonlj_min), args.epsilonlj_max)
        sigma_lj = min(max(float(x[1]), args.sigmalj_min), args.sigmalj_max)
        return epsilon_lj, sigma_lj

    if args.single_run:
        epsilon_lj, sigma_lj = clamp(args.single_run)
        summary = evaluate_core_parameter_set(epsilon_lj, sigma_lj, config, args, run_dir)
        logging.info("Single-run objective: %.6f", summary.objective)
        return 0

    require_runtime_dependencies()
    from scipy.optimize import minimize

    x0 = [args.epsilonlj_init, args.sigmalj_init]
    bounds = [(args.epsilonlj_min, args.epsilonlj_max), (args.sigmalj_min, args.sigmalj_max)] if args.method == "L-BFGS-B" else None

    def objective(x):
        epsilon_lj, sigma_lj = clamp(x)
        summary = evaluate_core_parameter_set(epsilon_lj, sigma_lj, config, args, run_dir)
        return summary.objective

    result = minimize(
        objective,
        x0,
        method=args.method,
        bounds=bounds,
        options={"disp": True, "maxiter": args.max_evals},
    )

    epsilon_lj, sigma_lj = clamp(result.x)
    logging.info("Core optimization complete")
    logging.info("  epsilonLJ = %.6f", epsilon_lj)
    logging.info("  sigmaLJ = %.6f", sigma_lj)
    logging.info("Best objective: %.6f", float(result.fun))
    logging.info("Results directory: %s", run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
