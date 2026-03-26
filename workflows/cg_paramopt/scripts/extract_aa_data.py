#!/usr/bin/env python3
"""Extract AA geometry observables into the CG-ParamOpt CSV layout."""

from __future__ import annotations

import argparse
import csv
import glob
import statistics
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract AA geometry data for the HEVA CG-ParamOpt workflow.")
    parser.add_argument("--input-dir", type=Path, required=True, help="Directory containing T4-AA-*.dat files.")
    parser.add_argument("--output", type=Path, required=True, help="CSV file to write.")
    return parser.parse_args()


def extract_rows(path: Path) -> list[list[float]]:
    l0: list[float] = []
    l1: list[float] = []
    t0: list[float] = []
    t1: list[float] = []
    p33: list[float] = []
    p01: list[float] = []
    p12: list[float] = []
    p20: list[float] = []

    with path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()

    for line in lines[2:242]:
        kind, _, value = line.split()
        if int(float(kind)) == 0:
            l0.append(float(value))
        elif int(float(kind)) == 1:
            l1.append(float(value))

    for line in lines[243:483]:
        kind, _, value = line.split()
        if int(float(kind)) == 0:
            t0.append(float(value))
        elif int(float(kind)) == 1:
            t1.append(float(value))

    for line in lines[484:]:
        kind, _, _, value = line.split()
        phi_kind = int(float(kind))
        if phi_kind == 0:
            p33.append(float(value))
        elif phi_kind == 1:
            p12.append(float(value))
        elif phi_kind == 2:
            p01.append(float(value))
        elif phi_kind == 3:
            p20.append(float(value))

    avg_l = (statistics.mean(l0) + statistics.mean(l1)) / 2.0
    l0 = [value / avg_l for value in l0]
    l1 = [value / avg_l for value in l1]

    return [list(row) for row in zip(l0, l1, t0, t1, p33, p12, p01, p20)]


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_path = args.output.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data_files = [Path(path) for path in sorted(glob.glob(str(input_dir / "T4-AA-*.dat")))]
    if not data_files:
        raise SystemExit(f"No T4-AA-*.dat files found under {input_dir}")

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["l0", "l1", "t0", "t1", "p33", "p12", "p01", "p20"])
        for path in data_files:
            for row in extract_rows(path):
                writer.writerow(row)

    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
