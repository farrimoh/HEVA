# HEVA CG ParamOpt Workflow

This workflow reproduces the historical CG-ParamOpt run shape on top of the maintained HEVA backend instead of copying the old C++ source tree into HEVA.

## Scope

This workflow is for:

- relaxing a fixed initialized structure
- sampling geometry without assembly moves
- writing CG-ParamOpt-compatible `data.dat` output for downstream optimization

It is not an assembly optimization workflow.

## Backend Split

The split between C++ and Python is intentional:

- C++ keeps restart loading, relaxation-only Monte Carlo, per-sample observable export, and a small AA-translated `Initial_frame.dat` preparation binary that applies the CG-ParamOpt scaling/normalization step before writing a HEVA restart file
- Python keeps trial-directory setup, AA-translated initial-frame preparation orchestration, config generation, seed/output management, and AA-data preparation

This keeps HEVA as the simulation backend while preserving the historical workflow contract where it matters.

## What Was Added To HEVA

The HEVA backend now supports:

- `init.mode = restart | seed`
- `init.path` as the generic path for restart-style state files
- `runtime.resume = true | false` to distinguish continuation from fresh-start semantics
- `runtime.workflow = relaxation` for move-vertex-only evolution
- optional `[cg_paramopt]` sampling controls that append CG-ParamOpt-compatible geometry rows into `data.dat`
- optional CG-ParamOpt geometry-reference parameters in `[capsid_geometry]`
- `src/prepare_initial_frame` to apply the CG-ParamOpt AA-frame scaling/normalization and write a standard HEVA `restart_lammps.dat`

## Files

```text
workflows/cg_paramopt/
|-- README.md
|-- cg_paramopt/
|   |-- config.py
|   |-- metrics.py
|   |-- prepare.py
|   |-- runner.py
|   |-- trial.py
|   `-- verify.py
|-- config/
|   `-- base_config.in
`-- scripts/
    |-- extract_aa_data.py
    |-- optimize.py
    |-- prepare_initial_frame.py
    `-- run_trial.py
```

## Trial Run

From the HEVA repository root, using a CG-ParamOpt `Initial_frame.dat` input:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure /path/to/Initial_frame.dat \
  --input-format initial_frame \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800
```

In this workflow, `initial_frame` means the CG-ParamOpt `Initial_frame.dat` source that still needs the preparation step, including the AA-derived scaling/normalization needed before it becomes a normal HEVA restart file.

By default this mirrors the historical cadence:

- `50000` relaxation sweeps
- `5000` sampled frames
- one sample per sweep into `data.dat`

Useful overrides:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure /path/to/Initial_frame.dat \
  --input-format initial_frame \
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

If you already have a HEVA restart-style file, skip the AA-translated initial-frame prep step and start fresh from it:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure /path/to/restart_lammps.dat \
  --input-format restart \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800
```

To apply that AA-frame scaling/normalization step and write a restart-style file without running a trial:

```bash
python3 workflows/cg_paramopt/scripts/prepare_initial_frame.py \
  --input /path/to/Initial_frame.dat \
  --output-dir /path/to/prepared-frame
```

## AA Reference Extraction

To convert AA graph dumps into the CSV layout expected by the CG-ParamOpt optimizer:

```bash
python3 workflows/cg_paramopt/scripts/extract_aa_data.py \
  --input-dir /path/to/AA-graphdata \
  --output workflows/cg_paramopt/data/AAdata.csv
```

## Optimization Driver

The optimizer path expects a Python environment with:

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install "setuptools<81" numpy pandas matplotlib scipy hyperopt
```

`hyperopt 0.2.7` still imports `pkg_resources`, so the tested environment pins `setuptools` below `81`.

For one scored HEVA-backed trial:

```bash
python3 workflows/cg_paramopt/scripts/optimize.py \
  --aa-data /path/to/AAdata.csv \
  --initial-structure /path/to/Initial_frame.dat \
  --input-format initial_frame \
  --single-run 2.15 1.25 1.80
```

For a Hyperopt search:

```bash
python3 workflows/cg_paramopt/scripts/optimize.py \
  --aa-data /path/to/AAdata.csv \
  --initial-structure /path/to/Initial_frame.dat \
  --input-format initial_frame \
  --max-evals 25
```

This keeps the optimizer in Python while leaving preparation, restart loading, relaxation, and observable sampling in HEVA.

## Validated Scenarios
The following command shapes were run successfully during the cleanup pass. Replace the placeholder paths with your own inputs and output directory.

AA-translated initial-frame preparation only:

```bash
python3 workflows/cg_paramopt/scripts/prepare_initial_frame.py \
  --input /path/to/Initial_frame.dat \
  --output-dir /path/to/prepared-frame
```

AA-translated initial-frame input -> prepared restart -> fresh relaxation trial:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure /path/to/Initial_frame.dat \
  --input-format initial_frame \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800 \
  --frames 1 \
  --relax-sweeps 2 \
  --sample-every 1 \
  --output-dir /path/to/run-initial-frame
```

Restart-style input -> fresh relaxation trial:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure checks/fixtures/restart_lammps.dat \
  --input-format restart \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800 \
  --frames 1 \
  --relax-sweeps 2 \
  --sample-every 1 \
  --output-dir /path/to/run-restart-fresh
```

Restart-style input -> continuation:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure checks/fixtures/restart_lammps.dat \
  --input-format restart \
  --resume \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800 \
  --frames 1 \
  --relax-sweeps 80391 \
  --sample-every 1 \
  --output-dir /path/to/run-restart-resume
```

Optimizer single-run checks:

```bash
python3 workflows/cg_paramopt/scripts/optimize.py \
  --aa-data /path/to/AAdata.csv \
  --initial-structure /path/to/Initial_frame.dat \
  --input-format initial_frame \
  --single-run 2.15 1.25 1.80 \
  --frames 1 \
  --relax-sweeps 2 \
  --sample-every 1
```

```bash
python3 workflows/cg_paramopt/scripts/optimize.py \
  --aa-data /path/to/AAdata.csv \
  --initial-structure checks/fixtures/restart_lammps.dat \
  --input-format restart \
  --single-run 2.15 1.25 1.80 \
  --frames 1 \
  --relax-sweeps 2 \
  --sample-every 1
```

Short search checks:

```bash
python3 workflows/cg_paramopt/scripts/optimize.py \
  --aa-data /path/to/AAdata.csv \
  --initial-structure /path/to/Initial_frame.dat \
  --input-format initial_frame \
  --max-evals 2 \
  --frames 1 \
  --relax-sweeps 2 \
  --sample-every 1
```

```bash
python3 workflows/cg_paramopt/scripts/optimize.py \
  --aa-data /path/to/AAdata.csv \
  --initial-structure checks/fixtures/restart_lammps.dat \
  --input-format restart \
  --max-evals 2 \
  --frames 1 \
  --relax-sweeps 2 \
  --sample-every 1
```

Observed validation signals:

- all 1-frame trials wrote `61` lines to `data.dat`
- initial-frame input now goes through a prep step and then the standard restart path
- `resume = false` and `resume = true` both behaved as intended
- both 2-eval search runs completed with `OK` records in `minimize-record.csv`

## Current Status

- init semantics are now centered on `seed` versus `restart`
- continuation semantics are controlled by `runtime.resume`
- CG-ParamOpt AA-translated initial-frame compatibility is isolated behind `prepare_initial_frame`
- remaining cleanup is mostly polish: reducing noisy initial-frame preparation console prints and tightening the docs further if needed

## Outputs

Each run directory contains:

- `prepared/restart_lammps.dat` when `--input-format initial_frame` on an AA-translated input
- `trial_config.in`
- `data.dat`
- `energy.dat`
- `restart_lammps.dat`
- `run_config.out`
- `parameters_run.out`

The key optimization-facing artifact is `data.dat`, which uses the historical CG-ParamOpt column layout:

```text
L0,L1,Theta0,Theta1,Phi33,Phi01,Phi12,Phi20
```
