#!/usr/bin/env python3
"""Run one HEVA-backed CG-ParamOpt relaxation-sampling trial."""

from __future__ import annotations

import argparse
import configparser
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one HEVA CG-ParamOpt-style trial.")
    parser.add_argument("--initial-structure", type=Path, required=True, help="Initial-frame compatibility input or HEVA restart_lammps.dat.")
    parser.add_argument("--init-mode", choices=("initial_frame", "restart", "legacy_lammps"), default="initial_frame", help="How HEVA should read --initial-structure. Use initial_frame for CG-ParamOpt-style imports.")
    parser.add_argument("--epsilon0", type=float, required=True, help="Legacy kappa_l analog.")
    parser.add_argument("--kappa0", type=float, required=True, help="Legacy kappa_theta analog.")
    parser.add_argument("--kappaPhi0", type=float, required=True, help="Legacy kappa_phi analog.")
    parser.add_argument("--theta0", type=float, default=None, help="Equilibrium theta for type 0.")
    parser.add_argument("--theta1", type=float, default=None, help="Equilibrium theta for type 1.")
    parser.add_argument("--theta2", type=float, default=None, help="Equilibrium theta for type 2.")
    parser.add_argument("--theta3", type=float, default=None, help="Equilibrium theta for type 3.")
    parser.add_argument("--l0", type=float, default=None, help="Legacy AA-derived equilibrium length for type 0.")
    parser.add_argument("--l1", type=float, default=None, help="Legacy AA-derived equilibrium length for type 1.")
    parser.add_argument("--phi33", type=float, default=None, help="Legacy AA-derived phi33 target.")
    parser.add_argument("--phi12", type=float, default=None, help="Legacy AA-derived phi12 target.")
    parser.add_argument("--phi01", type=float, default=None, help="Legacy AA-derived phi01 target.")
    parser.add_argument("--phi20", type=float, default=None, help="Legacy AA-derived phi20 target.")
    parser.add_argument("--seed", type=int, default=521759, help="Random seed for the trial.")
    parser.add_argument("--relax-sweeps", type=int, default=50_000, help="Relaxation sweeps before sampling begins.")
    parser.add_argument("--frames", type=int, default=5_000, help="Number of sampled frames to write into data.dat.")
    parser.add_argument("--sample-every", type=int, default=1, help="Relaxation sweeps between sampled frames.")
    parser.add_argument("--profile", choices=("test", "run", "extended"), default="test", help="HEVA engine profile.")
    parser.add_argument("--index-capacity", type=int, default=None, help="Required only for profile=extended.")
    parser.add_argument("--assemble", type=Path, default=None, help="Path to the HEVA assemble binary.")
    parser.add_argument("--base-config", type=Path, default=None, help="Path to the workflow base config.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Explicit run output directory.")
    parser.add_argument("--dry-run", action="store_true", help="Write the generated config and print the command without running HEVA.")
    return parser.parse_args()


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


def build_run_dir(root: Path, explicit_output_dir: Path | None) -> Path:
    if explicit_output_dir is not None:
        return explicit_output_dir.resolve()
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    return (root / "workflows" / "cg_paramopt" / "runs" / f"{stamp}_trial").resolve()


def set_optional_float(config: configparser.ConfigParser, section: str, key: str, value: float | None) -> None:
    if value is not None:
        config[section][key] = format_float(value)


def main() -> int:
    args = parse_args()
    root = repo_root()

    assemble_path = args.assemble.resolve() if args.assemble else (root / "src" / "assemble").resolve()
    base_config_path = args.base_config.resolve() if args.base_config else (root / "workflows" / "cg_paramopt" / "config" / "base_config.in").resolve()
    initial_structure_path = args.initial_structure.resolve()
    run_dir = build_run_dir(root, args.output_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    if not assemble_path.exists():
        print(f"assemble binary not found: {assemble_path}", file=sys.stderr)
        print("Build HEVA first with `cd src && make`.", file=sys.stderr)
        return 1
    if not initial_structure_path.exists():
        print(f"initial structure not found: {initial_structure_path}", file=sys.stderr)
        return 1

    if args.relax_sweeps < 0 or args.frames <= 0 or args.sample_every <= 0:
        print("relax-sweeps must be >= 0, and frames/sample-every must be > 0.", file=sys.stderr)
        return 1

    config = load_config(base_config_path)
    config["capsid_geometry"]["epsilon0"] = format_float(args.epsilon0)
    config["capsid_geometry"]["kappa0"] = format_float(args.kappa0)
    config["capsid_geometry"]["kappaPhi0"] = format_float(args.kappaPhi0)
    set_optional_float(config, "capsid_geometry", "theta0", args.theta0)
    set_optional_float(config, "capsid_geometry", "theta1", args.theta1)
    set_optional_float(config, "capsid_geometry", "theta2", args.theta2)
    set_optional_float(config, "capsid_geometry", "theta3", args.theta3)
    set_optional_float(config, "capsid_geometry", "l0", args.l0)
    set_optional_float(config, "capsid_geometry", "l1", args.l1)
    set_optional_float(config, "capsid_geometry", "phi33", args.phi33)
    set_optional_float(config, "capsid_geometry", "phi12", args.phi12)
    set_optional_float(config, "capsid_geometry", "phi01", args.phi01)
    set_optional_float(config, "capsid_geometry", "phi20", args.phi20)

    total_sweeps = args.relax_sweeps + args.frames * args.sample_every
    config["init"]["mode"] = args.init_mode
    config["init"]["restartPath"] = str(initial_structure_path)
    config["runtime"]["seed"] = str(args.seed)
    config["runtime"]["maxSweeps"] = str(total_sweeps)
    config["runtime"]["outputDir"] = str(run_dir)
    config["runtime"]["workflow"] = "relaxation"
    config["engine"]["profile"] = args.profile
    config["cg_paramopt"]["sampleStartSweep"] = str(args.relax_sweeps)
    config["cg_paramopt"]["sampleEvery"] = str(args.sample_every)
    config["cg_paramopt"]["sampleOutputPath"] = "data.dat"

    if args.profile == "extended":
        if args.index_capacity is None:
            print("--index-capacity is required when --profile extended", file=sys.stderr)
            return 1
        config["engine"]["indexCapacity"] = str(args.index_capacity)
    else:
        config["engine"].pop("indexCapacity", None)

    generated_config_path = run_dir / "trial_config.in"
    write_config(config, generated_config_path)

    command = [str(assemble_path), "--config", str(generated_config_path)]

    print(f"Run directory: {run_dir}")
    print(f"Generated config: {generated_config_path}")
    print("Command:")
    print(" ".join(command))

    if args.dry_run:
        return 0

    subprocess.run(command, check=True, cwd=root / "src")
    print("Trial completed.")
    print(f"Sampled geometry: {run_dir / 'data.dat'}")
    print(f"Restart snapshot: {run_dir / 'restart_lammps.dat'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
