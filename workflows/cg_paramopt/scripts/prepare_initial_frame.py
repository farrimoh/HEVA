#!/usr/bin/env python3
"""Convert an AA-translated CG-ParamOpt Initial_frame.dat into a HEVA restart-style state file."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare a HEVA restart_lammps.dat from an AA-translated Initial_frame.dat input.")
    parser.add_argument("--input", type=Path, required=True, help="Path to the AA-translated Initial_frame.dat file.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory where restart_lammps.dat will be written.")
    parser.add_argument("--prepare-binary", type=Path, default=None, help="Path to the HEVA prepare_initial_frame binary.")
    parser.add_argument("--index-capacity", type=int, default=None, help="Optional index capacity override for the preparation step.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    workflow_dir = script_dir.parent
    repo_root = workflow_dir.parents[1]
    prepare_binary = args.prepare_binary.resolve() if args.prepare_binary else (repo_root / "src" / "prepare_initial_frame").resolve()
    input_path = args.input.resolve()
    output_dir = args.output_dir.resolve()

    if not prepare_binary.exists():
        print(f"prepare_initial_frame binary not found: {prepare_binary}", file=sys.stderr)
        print("Build HEVA first with `cd src && make`.", file=sys.stderr)
        return 1
    if not input_path.exists():
        print(f"input file not found: {input_path}", file=sys.stderr)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)
    command = [str(prepare_binary), str(input_path), str(output_dir)]
    if args.index_capacity is not None:
        command.append(str(args.index_capacity))

    print("Command:")
    print(" ".join(command))
    result = subprocess.run(command, check=False)
    if result.returncode == 0:
        print(output_dir / "restart_lammps.dat")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
