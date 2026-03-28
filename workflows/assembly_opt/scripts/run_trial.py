#!/usr/bin/env python3
"""Run one fixed assembly-yield evaluation using HEVA."""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent
if str(WORKFLOW_DIR) not in sys.path:
    sys.path.insert(0, str(WORKFLOW_DIR))

from assembly_opt.config import DEFAULT_INITIAL_PARAMS, PARAMETER_ORDER, build_config
from assembly_opt.runner import evaluate_parameter_set, initialize_run_folder
from assembly_opt.verify import run_static_verification


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one fixed assembly optimization evaluation using HEVA.")
    parser.add_argument("--backend", choices=("local", "slurm"), default="local", help="Execution backend for per-seed trials.")
    parser.add_argument("--skip-verify", action="store_true", help="Skip the static verification pass.")
    parser.add_argument("--verify-only", action="store_true", help="Run only the static verification pass.")
    parser.add_argument("--runs-dir", default=None, help="Output directory for workflow evaluations.")
    parser.add_argument("--base-config", default=None, help="Override the workflow base config.")
    parser.add_argument("--assemble", default=None, help="Override the HEVA assemble binary path.")
    parser.add_argument("--init-mode", choices=("restart", "seed", "triangle", "pentamer", "hexamer"), default="pentamer", help="Initialization mode for each seed run.")
    parser.add_argument("--seed-shape", choices=("triangle", "pentamer", "hexamer"), default=None, help="Seed shape when --init-mode=seed.")
    parser.add_argument("--initial-structure", default=None, help="Restart file path when --init-mode=restart.")
    parser.add_argument("--max-sweeps", type=int, default=0, help="Maximum sweeps per seed run. Use 0 for no workflow-imposed limit.")
    parser.add_argument("--runs-per-eval", type=int, default=10, help="Number of HEVA seeds to evaluate for one parameter set.")
    parser.add_argument("--max-parallel", type=int, default=1, help="Maximum concurrent local HEVA runs.")
    parser.add_argument("--base-seed", type=int, default=None, help="Override the first seed used by a new run directory.")
    parser.add_argument("--profile", choices=("test", "run", "extended"), default="run", help="HEVA engine profile.")
    parser.add_argument("--index-capacity", type=int, default=None, help="Required only for profile=extended.")
    parser.add_argument("--slurm-account", default=None, help="Primary SLURM account when --backend=slurm.")
    parser.add_argument("--slurm-partition", default=None, help="Primary SLURM partition when --backend=slurm.")
    parser.add_argument("--slurm-qos", default=None, help="Optional SLURM QoS.")
    parser.add_argument("--slurm-time-limit", default="00:20:00", help="SLURM time limit.")
    parser.add_argument("--slurm-poll-interval", type=int, default=30, help="Seconds between SLURM queue polls.")
    parser.add_argument("--slurm-max-concurrent-jobs", type=int, default=200, help="Maximum active SLURM jobs for one evaluation.")
    parser.add_argument("--slurm-ntasks", type=int, default=1, help="SLURM ntasks value.")
    parser.add_argument("--slurm-cpus-per-task", type=int, default=1, help="SLURM cpus-per-task value.")
    parser.add_argument("--slurm-module-line", action="append", default=None, help="Shell line to emit before running HEVA inside the job script.")
    parser.add_argument("--slurm-extra-sbatch-line", action="append", default=None, help="Additional raw #SBATCH line content, for example --mem=4G.")
    parser.add_argument("--slurm-alternate-account", action="append", default=None, metavar="ACCOUNT:PARTITION", help="Alternate account/partition pair used round-robin across seeds.")
    parser.add_argument("--dry-run", action="store_true", help="Materialize configs and commands without running HEVA.")
    parser.add_argument("--epsilon0", type=float, default=2500.0)
    parser.add_argument("--kappa0", type=float, default=100.0)
    parser.add_argument("--kappaPhi0", type=float, default=1000.0)
    parser.add_argument("--theta0", type=float, default=0.240)
    parser.add_argument("--theta1", type=float, default=0.480)
    parser.add_argument("--ks0", type=float, default=0.002)
    parser.add_argument("--simulation-dg", type=float, default=0.100, help="Fixed simulation.dg value.")
    parser.add_argument("--mudrug", type=float, default=0.0)
    parser.add_argument("--gdrug0", type=float, default=0.0)
    parser.add_argument("--kd0", type=float, default=0.0)
    parser.add_argument("--core-enabled", action="store_true", help="Enable the core interaction model for each trial.")
    parser.add_argument("--core-max-bonds", type=int, default=1000, help="Legacy maxbondsRNA control when --core-enabled.")
    parser.add_argument("--core-epsilonlj", type=float, default=0.0, help="Core epsilonLJ value.")
    parser.add_argument("--core-sigmalj", type=float, default=0.0, help="Core sigmaLJ value.")
    parser.add_argument("--target-yield", type=float, default=1.0)
    parser.add_argument("--target-ratio", type=float, default=0.9)
    parser.add_argument("--yield-weight", type=float, default=100.0)
    parser.add_argument("--ratio-weight", type=float, default=100.0)
    parser.add_argument("--offtarget-weight", type=float, default=10000.0)
    parser.add_argument("--failure-penalty", type=float, default=0.0)
    for name in PARAMETER_ORDER:
        parser.add_argument(f"--{name}", type=float, default=DEFAULT_INITIAL_PARAMS[name])
    return parser.parse_args()


def log_verification(config) -> bool:
    ok = True
    for item in run_static_verification(config):
        log = getattr(logging, item.level, logging.info)
        log("Verification %s: %s", item.label, item.detail)
        if item.level == "error":
            ok = False
    return ok


def parse_alternate_accounts(values: list[str] | None) -> list[tuple[str, str]]:
    result: list[tuple[str, str]] = []
    for value in values or []:
        if ":" not in value:
            raise ValueError(f"invalid --slurm-alternate-account value: {value!r}")
        account, partition = value.split(":", 1)
        result.append((account, partition))
    return result


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    alternate_accounts = parse_alternate_accounts(args.slurm_alternate_account)

    config = build_config(
        WORKFLOW_DIR,
        runs_dir=args.runs_dir,
        base_config=args.base_config,
        assemble=args.assemble,
        init_mode=args.init_mode,
        seed_shape=args.seed_shape,
        initial_structure=args.initial_structure,
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
        core_enabled=args.core_enabled,
        core_max_bonds=args.core_max_bonds,
        core_epsilon_lj=args.core_epsilonlj,
        core_sigma_lj=args.core_sigmalj,
        target_yield=args.target_yield,
        target_ratio=args.target_ratio,
        yield_weight=args.yield_weight,
        ratio_weight=args.ratio_weight,
        offtarget_weight=args.offtarget_weight,
        failure_penalty=args.failure_penalty,
        backend=args.backend,
        slurm_account=args.slurm_account,
        slurm_partition=args.slurm_partition,
        slurm_qos=args.slurm_qos,
        slurm_time_limit=args.slurm_time_limit,
        slurm_poll_interval=args.slurm_poll_interval,
        slurm_max_concurrent_jobs=args.slurm_max_concurrent_jobs,
        slurm_ntasks=args.slurm_ntasks,
        slurm_cpus_per_task=args.slurm_cpus_per_task,
        slurm_module_lines=args.slurm_module_line,
        slurm_extra_sbatch_lines=args.slurm_extra_sbatch_line,
        slurm_alternate_accounts=alternate_accounts,
    )

    if not args.skip_verify:
        if not log_verification(config):
            return 1
        if args.verify_only:
            return 0

    run_dir = config.resolve_runs_dir() / time.strftime("%Y%m%d-%H%M%S_trial")
    initialize_run_folder(run_dir, config.base_seed)
    parameters = {name: getattr(args, name) for name in PARAMETER_ORDER}
    result = evaluate_parameter_set(parameters, config, run_dir, dry_run=args.dry_run)

    logging.info("Evaluation directory: %s", result.evaluation_dir)
    logging.info("Objective: %.6f", result.summary.objective)
    logging.info(
        "Counts: closed=%s t4=%s t3=%s failed=%s/%s",
        result.summary.closed_count,
        result.summary.t4_count,
        result.summary.t3_count,
        result.summary.failed_runs,
        result.summary.total_runs,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
