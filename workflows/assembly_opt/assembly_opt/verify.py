from __future__ import annotations

from dataclasses import dataclass

from .config import RuntimeConfig
from .slurm import is_slurm_available


@dataclass(frozen=True)
class VerificationItem:
    label: str
    detail: str
    level: str


def run_static_verification(config: RuntimeConfig) -> list[VerificationItem]:
    items: list[VerificationItem] = []

    assemble_path = config.resolve_assemble_path()
    items.append(
        VerificationItem(
            label="assemble",
            detail=f"assemble binary {'found' if assemble_path.exists() else 'missing'} at {assemble_path}",
            level="info" if assemble_path.exists() else "error",
        )
    )

    base_config_path = config.resolve_base_config_path()
    items.append(
        VerificationItem(
            label="base_config",
            detail=f"base config {'found' if base_config_path.exists() else 'missing'} at {base_config_path}",
            level="info" if base_config_path.exists() else "error",
        )
    )

    if config.init_mode == "restart":
        restart_path = config.resolve_initial_structure_path()
        exists = restart_path is not None and restart_path.exists()
        items.append(
            VerificationItem(
                label="initial_structure",
                detail=(
                    f"restart input {'found' if exists else 'missing'} at {restart_path}"
                    if restart_path is not None
                    else "restart mode selected but no initial structure was provided"
                ),
                level="info" if exists else "error",
            )
        )
    elif config.init_mode == "seed" and config.seed_shape is None:
        items.append(
            VerificationItem(
                label="seed_shape",
                detail="init_mode=seed requires seed_shape",
                level="error",
            )
        )

    if config.profile == "extended" and config.index_capacity is None:
        items.append(
            VerificationItem(
                label="index_capacity",
                detail="profile=extended requires index_capacity",
                level="error",
            )
        )

    if config.runs_per_eval <= 0:
        items.append(
            VerificationItem(
                label="runs_per_eval",
                detail="runs_per_eval must be positive",
                level="error",
            )
        )

    if config.max_parallel <= 0:
        items.append(
            VerificationItem(
                label="max_parallel",
                detail="max_parallel must be positive",
                level="error",
            )
        )

    if config.backend == "slurm":
        if config.slurm is None:
            items.append(
                VerificationItem(
                    label="slurm",
                    detail="backend=slurm requires account and partition configuration",
                    level="error",
                )
            )
        else:
            level = "info" if is_slurm_available() else "warning"
            items.append(
                VerificationItem(
                    label="slurm_runtime",
                    detail="SLURM commands are available" if level == "info" else "SLURM commands not found in this environment",
                    level=level,
                )
            )

    return items
