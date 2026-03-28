from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


PARAMETER_ORDER = [
    "muCd",
    "dmu",
    "gb0",
    "dg12",
    "dg01",
    "dg20",
    "dg33",
    "dg00",
    "dgother",
]

DEFAULT_INITIAL_PARAMS = {
    "muCd": -7.0,
    "dmu": -3.0,
    "gb0": -6.8,
    "dg12": 0.3,
    "dg01": 0.1,
    "dg20": -0.1,
    "dg33": 0.0,
    "dg00": -0.7,
    "dgother": -0.8,
}

DEFAULT_PARAM_BOUNDS = {
    "muCd": (-8.0, -5.5),
    "dmu": (-5.0, -2.0),
    "gb0": (-7.5, -6.2),
    "dg12": (0.0, 0.5),
    "dg01": (-0.3, 0.3),
    "dg20": (-0.3, 0.3),
    "dg33": (-0.3, 0.3),
    "dg00": (-1.0, 0.0),
    "dgother": (-1.0, 0.0),
}

CORE_PARAMETER_ORDER = [
    "epsilonLJ",
    "sigmaLJ",
]

DEFAULT_CORE_INITIAL_PARAMS = {
    "epsilonLJ": 0.05,
    "sigmaLJ": 0.30,
}

DEFAULT_CORE_PARAM_BOUNDS = {
    "epsilonLJ": (0.01, 0.50),
    "sigmaLJ": (0.20, 0.80),
}


@dataclass(frozen=True)
class SlurmAccountPair:
    account: str
    partition: str


@dataclass(frozen=True)
class SlurmConfig:
    account: str
    partition: str
    qos: str | None
    time_limit: str
    poll_interval: int
    max_concurrent_jobs: int
    ntasks: int
    cpus_per_task: int
    module_lines: tuple[str, ...]
    extra_sbatch_lines: tuple[str, ...]
    alternate_accounts: tuple[SlurmAccountPair, ...]


@dataclass(frozen=True)
class RuntimeConfig:
    base_dir: Path
    runs_dir: Path
    base_config_path: Path
    assemble_path: Path
    init_mode: str
    seed_shape: str | None
    initial_structure_path: Path | None
    max_evals: int
    max_sweeps: int
    runs_per_eval: int
    max_parallel: int
    base_seed: int
    profile: str
    index_capacity: int | None
    epsilon0: float
    kappa0: float
    kappaPhi0: float
    theta0: float
    theta1: float
    ks0: float
    simulation_dg: float
    mudrug: float
    gdrug0: float
    kd0: float
    core_enabled: bool
    core_max_bonds: int
    core_epsilon_lj: float
    core_sigma_lj: float
    target_yield: float
    target_ratio: float
    yield_weight: float
    ratio_weight: float
    offtarget_weight: float
    failure_penalty: float = 0.0
    backend: str = "local"
    slurm: SlurmConfig | None = None

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

    def resolve_initial_structure_path(self) -> Path | None:
        if self.initial_structure_path is None:
            return None
        if self.initial_structure_path.is_absolute():
            return self.initial_structure_path
        return (self.base_dir / self.initial_structure_path).resolve()


def build_config(
    base_dir: Path,
    *,
    runs_dir: str | None = None,
    base_config: str | None = None,
    assemble: str | None = None,
    init_mode: str = "pentamer",
    seed_shape: str | None = None,
    initial_structure: str | None = None,
    max_evals: int = 10,
    max_sweeps: int = 0,
    runs_per_eval: int = 10,
    max_parallel: int = 1,
    base_seed: int = 1000,
    profile: str = "run",
    index_capacity: int | None = None,
    epsilon0: float = 2500.0,
    kappa0: float = 100.0,
    kappaPhi0: float = 1000.0,
    theta0: float = 0.240,
    theta1: float = 0.480,
    ks0: float = 0.002,
    simulation_dg: float = 0.100,
    mudrug: float = 0.0,
    gdrug0: float = 0.0,
    kd0: float = 0.0,
    core_enabled: bool = False,
    core_max_bonds: int = 1000,
    core_epsilon_lj: float = 0.0,
    core_sigma_lj: float = 0.0,
    target_yield: float = 1.0,
    target_ratio: float = 0.9,
    yield_weight: float = 100.0,
    ratio_weight: float = 100.0,
    offtarget_weight: float = 10000.0,
    failure_penalty: float = 0.0,
    backend: str = "local",
    slurm_account: str | None = None,
    slurm_partition: str | None = None,
    slurm_qos: str | None = None,
    slurm_time_limit: str = "00:20:00",
    slurm_poll_interval: int = 30,
    slurm_max_concurrent_jobs: int = 200,
    slurm_ntasks: int = 1,
    slurm_cpus_per_task: int = 1,
    slurm_module_lines: list[str] | None = None,
    slurm_extra_sbatch_lines: list[str] | None = None,
    slurm_alternate_accounts: list[tuple[str, str]] | None = None,
) -> RuntimeConfig:
    def user_path(value: str) -> Path:
        return Path(value).expanduser().resolve()

    slurm = None
    if backend == "slurm" and slurm_account and slurm_partition:
        slurm = SlurmConfig(
            account=slurm_account,
            partition=slurm_partition,
            qos=slurm_qos,
            time_limit=slurm_time_limit,
            poll_interval=slurm_poll_interval,
            max_concurrent_jobs=slurm_max_concurrent_jobs,
            ntasks=slurm_ntasks,
            cpus_per_task=slurm_cpus_per_task,
            module_lines=tuple(slurm_module_lines or ()),
            extra_sbatch_lines=tuple(slurm_extra_sbatch_lines or ()),
            alternate_accounts=tuple(
                SlurmAccountPair(account=account, partition=partition)
                for account, partition in (slurm_alternate_accounts or ())
            ),
        )

    return RuntimeConfig(
        base_dir=base_dir.resolve(),
        runs_dir=user_path(runs_dir) if runs_dir else Path("./runs"),
        base_config_path=user_path(base_config) if base_config else Path("./config/base_config.in"),
        assemble_path=user_path(assemble) if assemble else Path("../../src/assemble"),
        init_mode=init_mode,
        seed_shape=seed_shape,
        initial_structure_path=user_path(initial_structure) if initial_structure else None,
        max_evals=max_evals,
        max_sweeps=max_sweeps,
        runs_per_eval=runs_per_eval,
        max_parallel=max_parallel,
        base_seed=base_seed,
        profile=profile,
        index_capacity=index_capacity,
        epsilon0=epsilon0,
        kappa0=kappa0,
        kappaPhi0=kappaPhi0,
        theta0=theta0,
        theta1=theta1,
        ks0=ks0,
        simulation_dg=simulation_dg,
        mudrug=mudrug,
        gdrug0=gdrug0,
        kd0=kd0,
        core_enabled=core_enabled,
        core_max_bonds=core_max_bonds,
        core_epsilon_lj=core_epsilon_lj,
        core_sigma_lj=core_sigma_lj,
        target_yield=target_yield,
        target_ratio=target_ratio,
        yield_weight=yield_weight,
        ratio_weight=ratio_weight,
        offtarget_weight=offtarget_weight,
        failure_penalty=failure_penalty,
        backend=backend,
        slurm=slurm,
    )
