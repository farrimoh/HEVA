from __future__ import annotations

import csv
import json
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from .config import PARAMETER_ORDER, RuntimeConfig
from .metrics import EvaluationSummary, SeedOutcome, classify_outcome, summarize_outcomes, write_seed_outcomes, write_summary_json
from .slurm import execute_trial_batch
from .trial import TrialExecution, TrialSpec, execute_trial


@dataclass(frozen=True)
class EvaluationResult:
    parameters: dict[str, float]
    summary: EvaluationSummary
    evaluation_dir: Path
    command_paths: list[Path]


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
                    *PARAMETER_ORDER,
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


def _load_counter(path: Path) -> int:
    return int(path.read_text(encoding="utf-8").strip())


def _store_counter(path: Path, value: int) -> None:
    path.write_text(f"{value}\n", encoding="utf-8")


def allocate_evaluation(run_dir: Path, runs_per_eval: int) -> tuple[int, list[int], Path]:
    eval_path = run_dir / "next_eval.txt"
    seed_path = run_dir / "next_seed.txt"
    eval_index = _load_counter(eval_path)
    next_seed = _load_counter(seed_path)
    seeds = list(range(next_seed, next_seed + runs_per_eval))
    evaluation_dir = run_dir / f"eval_{eval_index:04d}"
    _store_counter(eval_path, eval_index + 1)
    _store_counter(seed_path, next_seed + runs_per_eval)
    return eval_index, seeds, evaluation_dir


def append_evaluation_record(run_dir: Path, evaluation_index: int, parameters: dict[str, float], result: EvaluationResult) -> None:
    with (run_dir / "submitted.csv").open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                evaluation_index,
                *[parameters[name] for name in PARAMETER_ORDER],
                f"{result.summary.objective:.6f}",
                result.summary.total_runs,
                result.summary.failed_runs,
                result.summary.closed_count,
                result.summary.t4_count,
                result.summary.t3_count,
                result.summary.other_closed,
                f"{result.summary.target_yield:.6f}",
                f"{result.summary.t4_ratio:.6f}",
                result.evaluation_dir,
            ]
        )


def build_trial_spec(config: RuntimeConfig, parameters: dict[str, float], seed: int, output_dir: Path) -> TrialSpec:
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
        gb0=parameters["gb0"],
        muCd=parameters["muCd"],
        ks0=config.ks0,
        dmu=parameters["dmu"],
        simulation_dg=config.simulation_dg,
        mudrug=config.mudrug,
        gdrug0=config.gdrug0,
        kd0=config.kd0,
        dmud=config.dmud,
        core_enabled=config.core_enabled,
        core_max_bonds=config.core_max_bonds,
        core_epsilon_lj=config.core_epsilon_lj,
        core_sigma_lj=config.core_sigma_lj,
        dg12=parameters["dg12"],
        dg01=parameters["dg01"],
        dg20=parameters["dg20"],
        dg33=parameters["dg33"],
        dg00=parameters["dg00"],
        dgother=parameters["dgother"],
    )


def _execute_seed(spec: TrialSpec, dry_run: bool) -> tuple[TrialExecution, SeedOutcome]:
    execution = execute_trial(spec, dry_run=dry_run)
    outcome = classify_outcome(spec.seed, execution.run_dir, execution.stdout, execution.returncode)
    return execution, outcome


def _classify_executions(executions: list[TrialExecution], seeds: list[int]) -> list[SeedOutcome]:
    return [
        classify_outcome(seed, execution.run_dir, execution.stdout, execution.returncode)
        for seed, execution in zip(seeds, executions, strict=True)
    ]


def _write_commands(path: Path, executions: list[TrialExecution]) -> None:
    payload = {
        "commands": [
            {
                "run_dir": str(execution.run_dir),
                "config": str(execution.generated_config_path),
                "command": execution.command,
                "returncode": execution.returncode,
                "launcher_path": None if execution.launcher_path is None else str(execution.launcher_path),
            }
            for execution in executions
        ]
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def evaluate_parameter_set(
    parameters: dict[str, float],
    config: RuntimeConfig,
    run_dir: Path,
    *,
    dry_run: bool = False,
) -> EvaluationResult:
    evaluation_index, seeds, evaluation_dir = allocate_evaluation(run_dir, config.runs_per_eval)
    evaluation_dir.mkdir(parents=True, exist_ok=True)

    executions: list[TrialExecution] = []
    outcomes: list[SeedOutcome] = []
    specs = [
        build_trial_spec(config, parameters, seed, evaluation_dir / f"seed-{seed}")
        for seed in seeds
    ]
    if dry_run and config.backend != "slurm":
        specs = specs[:1]

    if config.backend == "slurm":
        if config.slurm is None:
            raise ValueError("backend=slurm requires SLURM configuration")
        executions = execute_trial_batch(specs, config.slurm, dry_run=dry_run)
        outcomes = _classify_executions(executions, [spec.seed for spec in specs])
    elif config.max_parallel > 1 and not dry_run:
        with ThreadPoolExecutor(max_workers=config.max_parallel) as executor:
            futures = {
                executor.submit(_execute_seed, spec, False): spec.seed
                for spec in specs
            }
            for future in as_completed(futures):
                execution, outcome = future.result()
                executions.append(execution)
                outcomes.append(outcome)
    else:
        for spec in specs:
            execution, outcome = _execute_seed(spec, dry_run)
            executions.append(execution)
            outcomes.append(outcome)

    if dry_run:
        commands_path = evaluation_dir / "commands.json"
        _write_commands(commands_path, executions)
        summary = EvaluationSummary(
            total_runs=len(executions),
            failed_runs=0,
            closed_count=0,
            t4_count=0,
            t3_count=0,
            other_closed=0,
            target_yield=0.0,
            t4_ratio=0.0,
            objective=0.0,
        )
        result = EvaluationResult(parameters, summary, evaluation_dir, [execution.generated_config_path for execution in executions])
        append_evaluation_record(run_dir, evaluation_index, parameters, result)
        return result

    outcomes.sort(key=lambda item: item.seed)
    executions.sort(key=lambda item: item.run_dir.name)
    summary = summarize_outcomes(outcomes, config)
    write_seed_outcomes(evaluation_dir / "seed_outcomes.csv", outcomes)
    write_summary_json(evaluation_dir / "summary.json", summary, parameters, evaluation_dir)
    _write_commands(evaluation_dir / "commands.json", executions)

    result = EvaluationResult(parameters, summary, evaluation_dir, [execution.generated_config_path for execution in executions])
    append_evaluation_record(run_dir, evaluation_index, parameters, result)
    logging.info(
        "Evaluation %s objective %.6f yield %.3f ratio %.3f closed %s/%s",
        evaluation_index,
        summary.objective,
        summary.target_yield,
        summary.t4_ratio,
        summary.closed_count,
        summary.total_runs,
    )
    return result
