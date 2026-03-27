from __future__ import annotations

import configparser
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from .prepare import prepare_initial_frame


@dataclass(frozen=True)
class TrialSpec:
    workflow_dir: Path
    assemble_path: Path
    prepare_initial_frame_path: Path
    base_config_path: Path
    initial_structure_path: Path
    input_format: str
    epsilon0: float
    kappa0: float
    kappaPhi0: float
    theta0: float | None
    theta1: float | None
    theta2: float | None
    theta3: float | None
    l0: float | None
    l1: float | None
    phi33: float | None
    phi12: float | None
    phi01: float | None
    phi20: float | None
    seed: int
    relax_sweeps: int
    frames: int
    sample_every: int
    profile: str
    resume: bool
    index_capacity: int | None
    output_dir: Path | None = None


@dataclass(frozen=True)
class TrialExecution:
    run_dir: Path
    prepared_restart_path: Path
    generated_config_path: Path
    prepare_command: list[str]
    command: list[str]
    prepare_stdout: str
    stdout: str
    returncode: int


def repo_root_from_workflow(workflow_dir: Path) -> Path:
    return workflow_dir.parents[1]


def load_config(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser()
    parser.optionxform = str
    with path.open("r", encoding="utf-8") as handle:
        parser.read_file(handle)
    return parser


def write_config(config: configparser.ConfigParser, path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        config.write(handle)


def format_float(value: float) -> str:
    return f"{value:.6f}"


def build_run_dir(workflow_dir: Path, explicit_output_dir: Path | None, suffix: str = "trial") -> Path:
    if explicit_output_dir is not None:
        return explicit_output_dir.resolve()
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    return (workflow_dir / "runs" / f"{stamp}_{suffix}").resolve()


def set_optional_float(config: configparser.ConfigParser, section: str, key: str, value: float | None) -> None:
    if value is not None:
        config[section][key] = format_float(value)


def materialize_trial(spec: TrialSpec, run_dir: Path, state_path: Path) -> tuple[Path, Path, list[str]]:
    config = load_config(spec.base_config_path)
    config["capsid_geometry"]["epsilon0"] = format_float(spec.epsilon0)
    config["capsid_geometry"]["kappa0"] = format_float(spec.kappa0)
    config["capsid_geometry"]["kappaPhi0"] = format_float(spec.kappaPhi0)
    set_optional_float(config, "capsid_geometry", "theta0", spec.theta0)
    set_optional_float(config, "capsid_geometry", "theta1", spec.theta1)
    set_optional_float(config, "capsid_geometry", "theta2", spec.theta2)
    set_optional_float(config, "capsid_geometry", "theta3", spec.theta3)
    set_optional_float(config, "capsid_geometry", "l0", spec.l0)
    set_optional_float(config, "capsid_geometry", "l1", spec.l1)
    set_optional_float(config, "capsid_geometry", "phi33", spec.phi33)
    set_optional_float(config, "capsid_geometry", "phi12", spec.phi12)
    set_optional_float(config, "capsid_geometry", "phi01", spec.phi01)
    set_optional_float(config, "capsid_geometry", "phi20", spec.phi20)

    total_sweeps = spec.relax_sweeps + spec.frames * spec.sample_every
    config["init"]["mode"] = "restart"
    config["init"]["path"] = str(state_path)
    config["runtime"]["seed"] = str(spec.seed)
    config["runtime"]["maxSweeps"] = str(total_sweeps)
    config["runtime"]["outputDir"] = str(run_dir)
    config["runtime"]["workflow"] = "relaxation"
    config["runtime"]["resume"] = "true" if spec.resume else "false"
    config["engine"]["profile"] = spec.profile
    config["cg_paramopt"]["sampleStartSweep"] = str(spec.relax_sweeps)
    config["cg_paramopt"]["sampleEvery"] = str(spec.sample_every)
    config["cg_paramopt"]["sampleOutputPath"] = "data.dat"

    if spec.profile == "extended":
        if spec.index_capacity is None:
            raise ValueError("--index-capacity is required when profile=extended")
        config["engine"]["indexCapacity"] = str(spec.index_capacity)
    else:
        config["engine"].pop("indexCapacity", None)

    generated_config_path = run_dir / "trial_config.in"
    write_config(config, generated_config_path)
    command = [str(spec.assemble_path), "--config", str(generated_config_path)]
    return run_dir, generated_config_path, command


def execute_trial(spec: TrialSpec, *, dry_run: bool = False) -> TrialExecution:
    run_dir = build_run_dir(spec.workflow_dir, spec.output_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    prepared_restart_path = spec.initial_structure_path
    prepare_command: list[str] = []
    prepare_stdout = ""
    if spec.input_format == "initial_frame":
        prepared = prepare_initial_frame(
            prepare_binary=spec.prepare_initial_frame_path,
            source_path=spec.initial_structure_path,
            destination_dir=run_dir / "prepared",
            index_capacity=spec.index_capacity,
            dry_run=dry_run,
        )
        prepared_restart_path = prepared.prepared_restart_path
        prepare_command = prepared.command
        prepare_stdout = prepared.stdout
        if prepared.returncode != 0:
            return TrialExecution(
                run_dir,
                prepared_restart_path,
                run_dir / "trial_config.in",
                prepare_command,
                [],
                prepare_stdout,
                "",
                prepared.returncode,
            )

    run_dir, generated_config_path, command = materialize_trial(spec, run_dir, prepared_restart_path)
    if dry_run:
        return TrialExecution(run_dir, prepared_restart_path, generated_config_path, prepare_command, command, prepare_stdout, "", 0)

    repo_root = repo_root_from_workflow(spec.workflow_dir)
    result = subprocess.run(
        command,
        cwd=repo_root / "src",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    return TrialExecution(run_dir, prepared_restart_path, generated_config_path, prepare_command, command, prepare_stdout, result.stdout, result.returncode)
