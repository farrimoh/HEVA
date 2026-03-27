from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PreparedStructure:
    source_path: Path
    prepared_restart_path: Path
    command: list[str]
    stdout: str
    returncode: int


def prepare_initial_frame(
    *,
    prepare_binary: Path,
    source_path: Path,
    destination_dir: Path,
    index_capacity: int | None = None,
    dry_run: bool = False,
) -> PreparedStructure:
    destination_dir.mkdir(parents=True, exist_ok=True)
    command = [str(prepare_binary), str(source_path), str(destination_dir)]
    if index_capacity is not None:
        command.append(str(index_capacity))

    prepared_restart_path = destination_dir / "restart_lammps.dat"
    if dry_run:
        return PreparedStructure(source_path, prepared_restart_path, command, "", 0)

    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    return PreparedStructure(source_path, prepared_restart_path, command, result.stdout, result.returncode)
