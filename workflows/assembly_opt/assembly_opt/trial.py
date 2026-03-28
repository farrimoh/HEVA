from __future__ import annotations

import configparser
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class TrialSpec:
    workflow_dir: Path
    assemble_path: Path
    base_config_path: Path
    seed: int
    init_mode: str
    seed_shape: str | None
    initial_structure_path: Path | None
    output_dir: Path
    max_sweeps: int
    profile: str
    index_capacity: int | None
    epsilon0: float
    kappa0: float
    kappaPhi0: float
    theta0: float
    theta1: float
    gb0: float
    muCd: float
    ks0: float
    dmu: float
    simulation_dg: float
    mudrug: float
    gdrug0: float
    kd0: float
    dg12: float
    dg01: float
    dg20: float
    dg33: float
    dg00: float
    dgother: float


@dataclass(frozen=True)
class TrialExecution:
    run_dir: Path
    generated_config_path: Path
    command: list[str]
    stdout: str
    returncode: int
    launcher_path: Path | None = None


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


def materialize_trial(spec: TrialSpec) -> tuple[Path, list[str]]:
    config = load_config(spec.base_config_path)

    config["capsid_geometry"]["epsilon0"] = format_float(spec.epsilon0)
    config["capsid_geometry"]["kappa0"] = format_float(spec.kappa0)
    config["capsid_geometry"]["kappaPhi0"] = format_float(spec.kappaPhi0)
    config["capsid_geometry"]["theta0"] = format_float(spec.theta0)
    config["capsid_geometry"]["theta1"] = format_float(spec.theta1)
    config["capsid_geometry"]["gb0"] = format_float(spec.gb0)
    config["capsid_geometry"]["dg12"] = format_float(spec.dg12)
    config["capsid_geometry"]["dg01"] = format_float(spec.dg01)
    config["capsid_geometry"]["dg20"] = format_float(spec.dg20)
    config["capsid_geometry"]["dg33"] = format_float(spec.dg33)
    config["capsid_geometry"]["dg00"] = format_float(spec.dg00)
    config["capsid_geometry"]["dgother"] = format_float(spec.dgother)

    config["simulation"]["muCd"] = format_float(spec.muCd)
    config["simulation"]["ks0"] = format_float(spec.ks0)
    config["simulation"]["dmu"] = format_float(spec.dmu)
    config["simulation"]["dg"] = format_float(spec.simulation_dg)

    config["drug"]["mudrug"] = format_float(spec.mudrug)
    config["drug"]["gdrug0"] = format_float(spec.gdrug0)
    config["drug"]["kd0"] = format_float(spec.kd0)

    config["init"]["mode"] = spec.init_mode
    if spec.init_mode == "restart":
        if spec.initial_structure_path is None:
            raise ValueError("restart mode requires an initial structure path")
        config["init"]["path"] = str(spec.initial_structure_path)
        config["init"].pop("seedShape", None)
    else:
        config["init"].pop("path", None)
        if spec.seed_shape:
            config["init"]["seedShape"] = spec.seed_shape
        else:
            config["init"].pop("seedShape", None)

    config["runtime"]["seed"] = str(spec.seed)
    config["runtime"]["maxSweeps"] = str(spec.max_sweeps)
    config["runtime"]["outputDir"] = str(spec.output_dir)
    config["runtime"]["workflow"] = "assembly"
    config["runtime"]["resume"] = "false"

    config["engine"]["profile"] = spec.profile
    if spec.profile == "extended":
        if spec.index_capacity is None:
            raise ValueError("--index-capacity is required when profile=extended")
        config["engine"]["indexCapacity"] = str(spec.index_capacity)
    else:
        config["engine"].pop("indexCapacity", None)

    generated_config_path = spec.output_dir / "trial_config.in"
    write_config(config, generated_config_path)
    command = [str(spec.assemble_path), "--config", str(generated_config_path)]
    return generated_config_path, command


def execute_trial(spec: TrialSpec, *, dry_run: bool = False) -> TrialExecution:
    spec.output_dir.mkdir(parents=True, exist_ok=True)
    generated_config_path, command = materialize_trial(spec)
    if dry_run:
        return TrialExecution(spec.output_dir, generated_config_path, command, "", 0, None)

    repo_root = repo_root_from_workflow(spec.workflow_dir)
    result = subprocess.run(
        command,
        cwd=repo_root / "src",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    return TrialExecution(spec.output_dir, generated_config_path, command, result.stdout, result.returncode, None)
