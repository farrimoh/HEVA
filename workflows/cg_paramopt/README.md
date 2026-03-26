# HEVA CG ParamOpt Workflow

This workflow reproduces the **legacy CG-ParamOpt run shape** on top of the maintained HEVA backend instead of copying the old C++ source tree into HEVA.

## Scope

This workflow is for:

- relaxing a fixed initialized structure
- sampling geometry without assembly moves
- writing legacy-style `data.dat` output for downstream optimization

It is not an assembly optimization workflow.

## Backend Split

The split between C++ and Python is intentional:

- C++ keeps geometry import, relaxation-only Monte Carlo, and per-sample observable export
- Python keeps trial-directory setup, config generation, seed/output management, and AA-data preparation

This keeps HEVA as the simulation backend while preserving the old workflow contract where it matters.

## What Was Added To HEVA

The HEVA backend now supports:

- `init.mode = initial_frame` for scaled import of CG-ParamOpt-style `Initial_frame.dat` inputs
- `init.mode = legacy_lammps` as a backward-compatible alias for the same compatibility path
- `runtime.workflow = relaxation` for move-vertex-only evolution
- optional `[cg_paramopt]` sampling controls that append legacy-style geometry rows into `data.dat`
- optional legacy geometry-reference parameters in `[capsid_geometry]`

## Files

```text
workflows/cg_paramopt/
|-- README.md
|-- config/
|   `-- base_config.in
`-- scripts/
    |-- extract_aa_data.py
    `-- run_trial.py
```

## Trial Run

From the HEVA repository root:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure /mnt/d/Coding/SimOpt/CG-ParamOpt/data/HE-AA-lammps1.dat \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800
```

By default this mirrors the legacy cadence:

- `50000` relaxation sweeps
- `5000` sampled frames
- one sample per sweep into `data.dat`

Useful overrides:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure /path/to/Initial_frame.dat \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800 \
  --l0 1.061 \
  --l1 0.974 \
  --theta0 0.219 \
  --theta1 0.467 \
  --phi33 1.047 \
  --phi12 1.012 \
  --phi01 1.168 \
  --phi20 0.952 \
  --frames 1000 \
  --relax-sweeps 10000
```

## AA Reference Extraction

To convert AA graph dumps into the CSV layout expected by legacy optimization code:

```bash
python3 workflows/cg_paramopt/scripts/extract_aa_data.py \
  --input-dir /path/to/AA-graphdata \
  --output workflows/cg_paramopt/data/AAdata.csv
```

## Outputs

Each run directory contains:

- `trial_config.in`
- `data.dat`
- `energy.dat`
- `restart_lammps.dat`
- `run_config.out`
- `parameters_run.out`

The key optimization-facing artifact is `data.dat`, which uses the legacy column layout:

```text
L0,L1,Theta0,Theta1,Phi33,Phi01,Phi12,Phi20
```
