"""Shared data-path handling for the AA distribution plots.

Every script reads one data directory, given by --data-dir or the HEVA_AA_DATA environment variable:

    <data-dir>/observables/patches/   patch observable CSVs (capsid-calibrated junction estimator)
    <data-dir>/observables/capsid/     full T4 capsid observable CSV

Optional estimator-comparison folders (observables/calibrated, centroid, circumcenter) are used when present.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

PLOTS = Path(__file__).resolve().parent
SCRIPTS = PLOTS.parent
WORKFLOW = SCRIPTS.parent
CALIBRATION = WORKFLOW / "fixtures/aa_mapping/vertex_offset_calibration_T4.json"

LAYOUT = {
    "observables": "observables",
    "patches": "observables/patches",
    "capsid": "observables/capsid",
}


def add_arguments(parser: argparse.ArgumentParser, default_out: str) -> None:
    parser.add_argument("--data-dir", type=Path, default=os.environ.get("HEVA_AA_DATA"),
                        help="data directory (default: $HEVA_AA_DATA); layout in plots/README.md")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help=f"where figures are written (default: <data-dir>/{default_out})")
    parser.set_defaults(_default_out=default_out)


def resolve(args: argparse.Namespace, need: tuple[str, ...]) -> dict[str, Path]:
    if args.data_dir is None:
        sys.exit("error: give --data-dir or set HEVA_AA_DATA")
    root = Path(args.data_dir).expanduser().resolve()
    paths = {key: root / rel for key, rel in LAYOUT.items()}
    missing = [f"{key}: {paths[key]}" for key in need if not paths[key].exists()]
    if missing:
        sys.exit("error: missing data under --data-dir\n  " + "\n  ".join(missing))
    out = Path(args.out_dir).expanduser().resolve() if args.out_dir else root / args._default_out
    out.mkdir(parents=True, exist_ok=True)
    paths["root"], paths["out"] = root, out
    return paths
