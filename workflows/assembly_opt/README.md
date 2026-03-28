# HEVA Assembly Opt Workflow

This workflow searches HEVA assembly parameters against a yield-based objective.

Use `--max-sweeps 0` for production assembly optimization runs.

The default parameter set in this workflow is:

- `muCd = -7.0`
- `dmu = -3.0`
- `gb0 = -6.8`
- `dg12 = 0.3`
- `dg01 = 0.1`
- `dg20 = -0.1`
- `dg33 = 0.0`
- `dg00 = -0.7`
- `dgother = -0.8`

## What It Does

For each parameter set, the workflow:

- launches multiple stochastic HEVA assembly runs
- classifies the final state of each run from `last.dat`
- measures:
  - total closed shells
  - T4 count
  - T3 count
  - off-target closed shells
- scores the parameter set with:

```text
error = 10000 * other_closed^2
      +   100 * (1 - target_yield)^2
      +   100 * (0.9 - t4_ratio)^2
```

The optimized parameters are:

- `muCd`
- `dmu`
- `gb0`
- `dg12`
- `dg01`
- `dg20`
- `dg33`
- `dg00`
- `dgother`

## Fixed Model Settings

By default the workflow runs with:

- `epsilon0 = 2500`
- `kappa0 = 100`
- `kappaPhi0 = 1000`
- `theta0 = 0.240`
- `theta1 = 0.480`
- `ks0 = 0.002`
- `simulation.dg = 0.100`
- `mudrug = 0`
- `gdrug0 = 0`
- `kd0 = 0`

These can be overridden from the CLI when needed.

## Python Environment

For fixed trial runs, the standard system Python is enough.

For optimization with `scipy.optimize.minimize`:

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install scipy
```

## Run One Evaluation

From a deterministic seed start:

```bash
python3 workflows/assembly_opt/scripts/run_trial.py \
  --init-mode pentamer \
  --runs-per-eval 5 \
  --max-sweeps 0 \
  --muCd -7.0 \
  --dmu -3.0 \
  --gb0 -6.8 \
  --dg12 0.3 \
  --dg01 0.1 \
  --dg20 -0.1 \
  --dg33 0.0 \
  --dg00 -0.7 \
  --dgother -0.8
```

From a restart file:

```bash
python3 workflows/assembly_opt/scripts/run_trial.py \
  --init-mode restart \
  --initial-structure /path/to/restart_lammps.dat \
  --runs-per-eval 5 \
  --max-sweeps 0
```

Useful controls:

- `--runs-per-eval`
- `--max-parallel`
- `--max-sweeps`
- `--base-seed`
- `--profile`
- `--index-capacity`
- `--target-yield`
- `--target-ratio`
- `--yield-weight`
- `--ratio-weight`
- `--offtarget-weight`

## Local Parallelism

Local execution runs seed trials in parallel with `--max-parallel`.

```bash
python3 workflows/assembly_opt/scripts/run_trial.py \
  --init-mode pentamer \
  --runs-per-eval 10 \
  --max-parallel 4 \
  --max-sweeps 0
```

## Run The Optimizer

Evaluate one fixed parameter vector through the optimizer entrypoint:

```bash
python3 workflows/assembly_opt/scripts/optimize.py \
  --init-mode pentamer \
  --runs-per-eval 5 \
  --max-sweeps 0 \
  --single-run \
  -7.0 -3.0 -6.8 0.3 0.1 -0.1 0.0 -0.7 -0.8
```

Run a short search:

```bash
python3 workflows/assembly_opt/scripts/optimize.py \
  --init-mode pentamer \
  --runs-per-eval 5 \
  --max-sweeps 0 \
  --max-evals 10
```

Validated local sample:

```bash
.venv/bin/python workflows/assembly_opt/scripts/optimize.py \
  --init-mode pentamer \
  --runs-per-eval 2 \
  --max-parallel 2 \
  --max-sweeps 50000000 \
  --profile run \
  --method L-BFGS-B \
  --max-evals 1
```

## HPC / SLURM

The workflow supports both local execution and SLURM-backed per-seed fanout.

Use the SLURM backend with:

- `--backend slurm`
- `--slurm-account`
- `--slurm-partition`

Optional controls:

- `--slurm-qos`
- `--slurm-time-limit`
- `--slurm-poll-interval`
- `--slurm-max-concurrent-jobs`
- `--slurm-ntasks`
- `--slurm-cpus-per-task`
- repeated `--slurm-module-line`
- repeated `--slurm-extra-sbatch-line`
- repeated `--slurm-alternate-account ACCOUNT:PARTITION`

Dry-run job-script generation:

```bash
python3 workflows/assembly_opt/scripts/run_trial.py \
  --backend slurm \
  --slurm-account mrsec \
  --slurm-partition mrsec-compute \
  --slurm-alternate-account hagan-lab:hagan-compute \
  --slurm-module-line "module unload share_modules/ANACONDA/5.3_py3" \
  --slurm-module-line "module load share_modules/GSL/2.5" \
  --runs-per-eval 2 \
  --max-sweeps 0 \
  --dry-run
```

Live submission uses the same command without `--dry-run`.

## Outputs

Each evaluation creates:

- a workflow-level `submitted.csv`
- one `eval_####/` directory
- one per-seed run directory under that evaluation
- `seed_outcomes.csv`
- `summary.json`
- `commands.json`

Each seed directory contains the normal HEVA run artifacts, including:

- `trial_config.in`
- `run_config.out`
- `energy.dat`
- `last.dat`
- `restart_lammps.dat`

When using `--backend slurm`, each seed directory also contains:

- `job.sh`
- `slurm-*.out` and `slurm-*.err` after submission
- `returncode.txt`

## Classification

Runs are classified from the HEVA outputs as follows:

- `closed` means HEVA reported `STOP REASON closed`
- T4 means `Nv = 42` and `NE = 120`
- T3 means `Nv = 32` and `NE = 90`

## Files

```text
workflows/assembly_opt/
|-- README.md
|-- assembly_opt/
|   |-- config.py
|   |-- dependencies.py
|   |-- metrics.py
|   |-- runner.py
|   |-- slurm.py
|   |-- trial.py
|   `-- verify.py
|-- config/
|   `-- base_config.in
`-- scripts/
    |-- optimize.py
    `-- run_trial.py
```
