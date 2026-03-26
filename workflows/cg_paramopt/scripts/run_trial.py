#!/usr/bin/env python3
"""Run one HEVA-backed CG-ParamOpt relaxation-sampling trial."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent
if str(WORKFLOW_DIR) not in sys.path:
    sys.path.insert(0, str(WORKFLOW_DIR))

from cg_paramopt.trial import TrialSpec, execute_trial


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


def main() -> int:
    args = parse_args()

    workflow_dir = WORKFLOW_DIR.resolve()
    repo_root = workflow_dir.parents[1]
    assemble_path = args.assemble.resolve() if args.assemble else (repo_root / "src" / "assemble").resolve()
    base_config_path = args.base_config.resolve() if args.base_config else (workflow_dir / "config" / "base_config.in").resolve()
    initial_structure_path = args.initial_structure.resolve()

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

    spec = TrialSpec(
        workflow_dir=workflow_dir,
        assemble_path=assemble_path,
        base_config_path=base_config_path,
        initial_structure_path=initial_structure_path,
        init_mode=args.init_mode,
        epsilon0=args.epsilon0,
        kappa0=args.kappa0,
        kappaPhi0=args.kappaPhi0,
        theta0=args.theta0,
        theta1=args.theta1,
        theta2=args.theta2,
        theta3=args.theta3,
        l0=args.l0,
        l1=args.l1,
        phi33=args.phi33,
        phi12=args.phi12,
        phi01=args.phi01,
        phi20=args.phi20,
        seed=args.seed,
        relax_sweeps=args.relax_sweeps,
        frames=args.frames,
        sample_every=args.sample_every,
        profile=args.profile,
        index_capacity=args.index_capacity,
        output_dir=args.output_dir.resolve() if args.output_dir else None,
    )

    execution = execute_trial(spec, dry_run=args.dry_run)
    print(f"Run directory: {execution.run_dir}")
    print(f"Generated config: {execution.generated_config_path}")
    print("Command:")
    print(" ".join(execution.command))

    if args.dry_run:
        return 0

    if execution.stdout:
        print(execution.stdout, end="" if execution.stdout.endswith("\n") else "\n")

    if execution.returncode != 0:
        return execution.returncode

    print("Trial completed.")
    print(f"Sampled geometry: {execution.run_dir / 'data.dat'}")
    print(f"Restart snapshot: {execution.run_dir / 'restart_lammps.dat'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

