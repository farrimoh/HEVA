from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from .config import RuntimeConfig
from .dependencies import pd
from .metrics import combine_params, compute_error, compute_histogram
from .trial import TrialSpec, execute_trial


@dataclass(frozen=True)
class SimulationResult:
    error: float
    stdout: str
    command: str
    generated_config_path: Path
    trial_dir: Path


def initialize_run_folder(run_dir: Path, seed: int) -> None:
    logging.info("Initializing run folder: %s", run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "seed").write_text(f"{seed}\n", encoding="utf-8")
    record_path = run_dir / "minimize-record.csv"
    if not record_path.exists():
        record_path.write_text("kappa_l,kappa_theta,kappa_phi,error,status,trial_dir\n", encoding="utf-8")


def load_seed(run_dir: Path) -> int:
    seed = int((run_dir / "seed").read_text(encoding="utf-8").strip())
    (run_dir / "seed").write_text(str(seed + 1), encoding="utf-8")
    return seed


def append_trial_record(run_dir: Path, params: list[float], error: float, status: str, trial_dir: Path) -> None:
    with open(run_dir / "minimize-record.csv", "a", encoding="utf-8") as handle:
        handle.write(f"{params[0]:.0f},{params[1]:.0f},{params[2]:.0f},{error:.5f},{status},{trial_dir}\n")


def build_trial_spec(config: RuntimeConfig, seed: int, params: list[float], geometry_params: list[float], trial_dir: Path) -> TrialSpec:
    return TrialSpec(
        workflow_dir=config.base_dir,
        assemble_path=config.resolve_assemble_path(),
        prepare_initial_frame_path=config.resolve_prepare_initial_frame_path(),
        base_config_path=config.resolve_base_config_path(),
        initial_structure_path=config.resolve_initial_structure_path(),
        input_format=config.input_format,
        epsilon0=params[0],
        kappa0=params[1],
        kappaPhi0=params[2],
        theta0=geometry_params[2],
        theta1=geometry_params[3],
        theta2=0.2,
        theta3=0.2,
        l0=geometry_params[0],
        l1=geometry_params[1],
        phi33=geometry_params[4],
        phi12=geometry_params[5],
        phi01=geometry_params[6],
        phi20=geometry_params[7],
        seed=seed,
        relax_sweeps=config.relax_sweeps,
        frames=config.frames,
        sample_every=config.sample_every,
        profile=config.profile,
        resume=False,
        index_capacity=config.index_capacity,
        output_dir=trial_dir,
    )


def run_single_simulation(
    log_params: list[float],
    df_hist_aa,
    run_dir: Path,
    geometry_params: list[float],
    config: RuntimeConfig,
    *,
    dry_run: bool = False,
) -> SimulationResult:
    params = [10 ** log_params[0], 10 ** log_params[1], 10 ** log_params[2]]
    seed = load_seed(run_dir)
    trial_dir = run_dir / f"trial_seed_{seed}"
    spec = build_trial_spec(config, seed, params, geometry_params, trial_dir)
    execution = execute_trial(spec, dry_run=dry_run)
    command = " ".join(execution.command)

    if dry_run:
        append_trial_record(run_dir, params, 0.0, "DRY_RUN", execution.run_dir)
        return SimulationResult(0.0, execution.stdout, command, execution.generated_config_path, execution.run_dir)

    if execution.returncode != 0:
        logging.error("Simulation failed with return code %s", execution.returncode)
        if execution.stdout.strip():
            logging.error(execution.stdout.strip())
        append_trial_record(run_dir, params, config.penalty_error, "FAILED", execution.run_dir)
        return SimulationResult(config.penalty_error, execution.stdout, command, execution.generated_config_path, execution.run_dir)

    data_file = execution.run_dir / "data.dat"
    if not data_file.exists():
        logging.error("Output file not found: %s", data_file)
        append_trial_record(run_dir, params, config.penalty_error, "MISSING_DATA", execution.run_dir)
        return SimulationResult(config.penalty_error, execution.stdout, command, execution.generated_config_path, execution.run_dir)

    df_cg = pd.read_csv(data_file, header=0, names=["l0", "l1", "t0", "t1", "p33", "p12", "p01", "p20"]).astype(float)
    df_cg_combined = combine_params(df_cg)
    df_hist_cg = compute_histogram(df_cg_combined, df_hist_aa["bins"].tolist())
    error = compute_error(df_hist_cg, df_hist_aa, params, execution.run_dir)
    append_trial_record(run_dir, params, error, "OK", execution.run_dir)
    logging.info("Simulation error: %.3f", error)
    return SimulationResult(error, execution.stdout, command, execution.generated_config_path, execution.run_dir)
