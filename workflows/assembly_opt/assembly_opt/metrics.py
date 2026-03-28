from __future__ import annotations

import csv
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

from .config import RuntimeConfig


STOP_REASON_PATTERN = re.compile(r"STOP REASON\s+([A-Za-z0-9_]+)")
T4_SIGNATURE = (42, 120)
T3_SIGNATURE = (32, 90)


@dataclass(frozen=True)
class SeedOutcome:
    seed: int
    run_dir: Path
    returncode: int
    stop_reason: str
    closed: bool
    is_t4: bool
    is_t3: bool
    nv: int | None
    ne: int | None
    stdout_path: Path
    last_dat_path: Path | None


@dataclass(frozen=True)
class EvaluationSummary:
    total_runs: int
    failed_runs: int
    closed_count: int
    t4_count: int
    t3_count: int
    other_closed: int
    target_yield: float
    t4_ratio: float
    objective: float


def parse_stop_reason(stdout: str) -> str:
    match = STOP_REASON_PATTERN.search(stdout)
    if not match:
        return "unknown"
    return match.group(1)


def parse_last_dat(last_dat_path: Path) -> tuple[int | None, int | None]:
    try:
        line = last_dat_path.read_text(encoding="utf-8").strip()
    except OSError:
        return None, None
    if not line:
        return None, None
    fields = [item.strip() for item in line.split(",")]
    if len(fields) < 29:
        return None, None
    try:
        nv = int(float(fields[27]))
        ne = int(float(fields[28]))
    except ValueError:
        return None, None
    return nv, ne


def classify_outcome(seed: int, run_dir: Path, stdout: str, returncode: int) -> SeedOutcome:
    stdout_path = run_dir / "stdout.txt"
    stdout_path.write_text(stdout, encoding="utf-8")

    last_dat_path = run_dir / "last.dat"
    nv = None
    ne = None
    if last_dat_path.exists():
        nv, ne = parse_last_dat(last_dat_path)
    else:
        last_dat_path = None

    stop_reason = parse_stop_reason(stdout)
    closed = returncode == 0 and stop_reason == "closed"
    signature = (nv, ne)
    is_t4 = closed and signature == T4_SIGNATURE
    is_t3 = closed and signature == T3_SIGNATURE
    return SeedOutcome(
        seed=seed,
        run_dir=run_dir,
        returncode=returncode,
        stop_reason=stop_reason,
        closed=closed,
        is_t4=is_t4,
        is_t3=is_t3,
        nv=nv,
        ne=ne,
        stdout_path=stdout_path,
        last_dat_path=last_dat_path,
    )


def summarize_outcomes(outcomes: list[SeedOutcome], config: RuntimeConfig) -> EvaluationSummary:
    total_runs = len(outcomes)
    failed_runs = sum(1 for outcome in outcomes if outcome.returncode != 0)
    closed_count = sum(1 for outcome in outcomes if outcome.closed)
    t4_count = sum(1 for outcome in outcomes if outcome.is_t4)
    t3_count = sum(1 for outcome in outcomes if outcome.is_t3)
    target_count = t4_count + t3_count
    other_closed = closed_count - target_count
    target_yield = target_count / total_runs if total_runs > 0 else 0.0
    t4_ratio = t4_count / target_count if target_count > 0 else 0.0
    objective = (
        config.offtarget_weight * (other_closed ** 2)
        + config.yield_weight * ((config.target_yield - target_yield) ** 2)
        + config.ratio_weight * ((config.target_ratio - t4_ratio) ** 2)
        + config.failure_penalty * failed_runs
    )
    return EvaluationSummary(
        total_runs=total_runs,
        failed_runs=failed_runs,
        closed_count=closed_count,
        t4_count=t4_count,
        t3_count=t3_count,
        other_closed=other_closed,
        target_yield=target_yield,
        t4_ratio=t4_ratio,
        objective=objective,
    )


def write_seed_outcomes(path: Path, outcomes: list[SeedOutcome]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["seed", "returncode", "stop_reason", "closed", "is_t4", "is_t3", "nv", "ne", "run_dir"])
        for outcome in outcomes:
            writer.writerow([
                outcome.seed,
                outcome.returncode,
                outcome.stop_reason,
                int(outcome.closed),
                int(outcome.is_t4),
                int(outcome.is_t3),
                "" if outcome.nv is None else outcome.nv,
                "" if outcome.ne is None else outcome.ne,
                outcome.run_dir,
            ])


def write_summary_json(path: Path, summary: EvaluationSummary, params: dict[str, float], evaluation_dir: Path) -> None:
    payload = {
        "parameters": params,
        "summary": asdict(summary),
        "evaluation_dir": str(evaluation_dir),
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
