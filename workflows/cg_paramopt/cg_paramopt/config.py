from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


DEFAULT_PARAM_BOUNDS = {
    "kappa_l": (2, 3.5),
    "kappa_theta": (1, 3),
    "kappa_phi": (1, 3.2),
}


@dataclass(frozen=True)
class RuntimeConfig:
    base_dir: Path
    runs_dir: Path
    base_config_path: Path
    assemble_path: Path
    initial_structure_path: Path
    aa_data_path: Path
    max_evals: int
    frames: int
    relax_sweeps: int
    sample_every: int
    init_mode: str
    profile: str
    index_capacity: int | None
    penalty_error: float = 1e9

    def resolve_runs_dir(self) -> Path:
        if self.runs_dir.is_absolute():
            return self.runs_dir
        return (self.base_dir / self.runs_dir).resolve()

    def resolve_base_config_path(self) -> Path:
        if self.base_config_path.is_absolute():
            return self.base_config_path
        return (self.base_dir / self.base_config_path).resolve()

    def resolve_assemble_path(self) -> Path:
        if self.assemble_path.is_absolute():
            return self.assemble_path
        return (self.base_dir / self.assemble_path).resolve()

    def resolve_initial_structure_path(self) -> Path:
        if self.initial_structure_path.is_absolute():
            return self.initial_structure_path
        return (self.base_dir / self.initial_structure_path).resolve()

    def resolve_aa_data_path(self) -> Path:
        if self.aa_data_path.is_absolute():
            return self.aa_data_path
        return (self.base_dir / self.aa_data_path).resolve()


def build_config(
    base_dir: Path,
    *,
    runs_dir: str | None = None,
    base_config: str | None = None,
    assemble: str | None = None,
    initial_structure: str | None = None,
    aa_data: str | None = None,
    max_evals: int = 10,
    frames: int = 5000,
    relax_sweeps: int = 50000,
    sample_every: int = 1,
    init_mode: str = "initial_frame",
    profile: str = "test",
    index_capacity: int | None = None,
) -> RuntimeConfig:
    return RuntimeConfig(
        base_dir=base_dir.resolve(),
        runs_dir=Path(runs_dir) if runs_dir else Path("./runs"),
        base_config_path=Path(base_config) if base_config else Path("./config/base_config.in"),
        assemble_path=Path(assemble) if assemble else Path("../../src/assemble"),
        initial_structure_path=Path(initial_structure) if initial_structure else Path("./data/Initial_frame.dat"),
        aa_data_path=Path(aa_data) if aa_data else Path("./data/AAdata.csv"),
        max_evals=max_evals,
        frames=frames,
        relax_sweeps=relax_sweeps,
        sample_every=sample_every,
        init_mode=init_mode,
        profile=profile,
        index_capacity=index_capacity,
    )
