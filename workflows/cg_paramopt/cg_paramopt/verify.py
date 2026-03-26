from __future__ import annotations

import configparser
from dataclasses import dataclass
from pathlib import Path

from .config import RuntimeConfig


@dataclass(frozen=True)
class VerificationItem:
    level: str
    label: str
    detail: str


def load_config(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser()
    parser.optionxform = str
    with path.open("r", encoding="utf-8") as handle:
        parser.read_file(handle)
    return parser


def run_static_verification(config: RuntimeConfig) -> list[VerificationItem]:
    items: list[VerificationItem] = []

    assemble_path = config.resolve_assemble_path()
    items.append(
        VerificationItem(
            "info" if assemble_path.exists() else "error",
            "assemble binary",
            str(assemble_path),
        )
    )

    base_config_path = config.resolve_base_config_path()
    if not base_config_path.exists():
        items.append(VerificationItem("error", "base config", f"missing at {base_config_path}"))
        return items

    parser = load_config(base_config_path)
    init_mode = parser.get("init", "mode", fallback="")
    workflow = parser.get("runtime", "workflow", fallback="")
    has_cg_section = parser.has_section("cg_paramopt")
    items.append(VerificationItem("info" if init_mode in {"initial_frame", "legacy_lammps"} else "warning", "base init.mode", init_mode or "missing"))
    items.append(VerificationItem("info" if workflow == "relaxation" else "warning", "base runtime.workflow", workflow or "missing"))
    items.append(VerificationItem("info" if has_cg_section else "warning", "base [cg_paramopt]", "present" if has_cg_section else "missing"))

    initial_structure = config.resolve_initial_structure_path()
    items.append(
        VerificationItem(
            "info" if initial_structure.exists() else "warning",
            "initial structure",
            str(initial_structure),
        )
    )

    aa_data_path = config.resolve_aa_data_path()
    items.append(
        VerificationItem(
            "info" if aa_data_path.exists() else "warning",
            "AA data",
            str(aa_data_path),
        )
    )

    return items

