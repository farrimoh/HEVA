from __future__ import annotations

import os
import re
import shlex
import subprocess
import time
from pathlib import Path

from .config import SlurmConfig
from .trial import TrialExecution, TrialSpec, materialize_trial


JOB_ID_PATTERN = re.compile(r"Submitted batch job (\d+)")


def is_slurm_available() -> bool:
    try:
        result = subprocess.run(["squeue", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    except FileNotFoundError:
        return False
    return result.returncode == 0


def _account_partition(slurm: SlurmConfig, index: int) -> tuple[str, str]:
    options = [(slurm.account, slurm.partition)]
    options.extend((item.account, item.partition) for item in slurm.alternate_accounts)
    return options[index % len(options)]


def _render_job_script(spec: TrialSpec, command: list[str], slurm: SlurmConfig, job_index: int) -> str:
    account, partition = _account_partition(slurm, job_index)
    run_dir = spec.output_dir
    quoted_command = " ".join(shlex.quote(part) for part in command)

    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=heva-opt-{spec.seed}",
        "#SBATCH --output=slurm-%j.out",
        "#SBATCH --error=slurm-%j.err",
        f"#SBATCH --account={account}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --time={slurm.time_limit}",
        f"#SBATCH --ntasks={slurm.ntasks}",
        f"#SBATCH --cpus-per-task={slurm.cpus_per_task}",
    ]
    if slurm.qos:
        lines.append(f"#SBATCH --qos={slurm.qos}")
    lines.extend(f"#SBATCH {entry}" for entry in slurm.extra_sbatch_lines)
    lines.extend(
        [
            "",
            f"cd {shlex.quote(str(run_dir))}",
        ]
    )
    lines.extend(slurm.module_lines)
    lines.extend(
        [
            f"{quoted_command} > stdout.txt 2>&1",
            "status=$?",
            "printf '%s\\n' \"$status\" > returncode.txt",
            "exit \"$status\"",
        ]
    )
    return "\n".join(lines) + "\n"


def _submit_job(script_path: Path) -> str:
    result = subprocess.run(["sbatch", str(script_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or "sbatch failed")
    match = JOB_ID_PATTERN.search(result.stdout)
    if not match:
        raise RuntimeError(f"Could not parse SLURM job id from: {result.stdout.strip()}")
    return match.group(1)


def _running_jobs() -> set[str]:
    username = os.environ.get("USER", "")
    result = subprocess.run(["squeue", "-u", username, "-h", "-o", "%i"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or "squeue failed")
    jobs = set(result.stdout.splitlines())
    jobs.discard("")
    return jobs


def _collect_execution(spec: TrialSpec, generated_config_path: Path, command: list[str], script_path: Path) -> TrialExecution:
    stdout_path = spec.output_dir / "stdout.txt"
    stdout = stdout_path.read_text(encoding="utf-8") if stdout_path.exists() else ""
    returncode_path = spec.output_dir / "returncode.txt"
    if returncode_path.exists():
        try:
            returncode = int(returncode_path.read_text(encoding="utf-8").strip())
        except ValueError:
            returncode = 1
    else:
        returncode = 1
    return TrialExecution(spec.output_dir, generated_config_path, command, stdout, returncode, script_path)


def execute_trial_batch(specs: list[TrialSpec], slurm: SlurmConfig, *, dry_run: bool = False) -> list[TrialExecution]:
    executions: list[TrialExecution] = []
    materialized: list[tuple[TrialSpec, Path, list[str], Path]] = []

    for index, spec in enumerate(specs):
        spec.output_dir.mkdir(parents=True, exist_ok=True)
        generated_config_path, command = materialize_trial(spec)
        script_path = spec.output_dir / "job.sh"
        script_path.write_text(_render_job_script(spec, command, slurm, index), encoding="utf-8")
        script_path.chmod(0o755)
        materialized.append((spec, generated_config_path, command, script_path))
        if dry_run:
            executions.append(TrialExecution(spec.output_dir, generated_config_path, command, "", 0, script_path))

    if dry_run:
        return executions

    pending = list(materialized)
    active: dict[str, tuple[TrialSpec, Path, list[str], Path]] = {}
    completed: list[TrialExecution] = []

    while pending or active:
        while pending and len(active) < slurm.max_concurrent_jobs:
            spec, generated_config_path, command, script_path = pending.pop(0)
            job_id = _submit_job(script_path)
            active[job_id] = (spec, generated_config_path, command, script_path)

        if not active:
            continue

        running = _running_jobs()
        finished = [job_id for job_id in active if job_id not in running]
        if not finished:
            time.sleep(slurm.poll_interval)
            continue

        for job_id in finished:
            spec, generated_config_path, command, script_path = active.pop(job_id)
            completed.append(_collect_execution(spec, generated_config_path, command, script_path))

    completed.sort(key=lambda item: item.run_dir.name)
    return completed
