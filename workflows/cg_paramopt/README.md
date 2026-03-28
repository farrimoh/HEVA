# HEVA CG ParamOpt Workflow

This workflow fits coarse-grained geometry parameters against all-atom reference data using HEVA as the simulation backend.

## What It Does

For each trial parameter set, the workflow:

- starts from either an `Initial_frame.dat` input or a HEVA restart file
- runs HEVA in `workflow = relaxation` mode
- samples geometry into `data.dat`
- compares the sampled `l`, `theta`, and `phi` distributions to AA reference data
- scores the trial by KL divergence

The optimized parameters are:

- `kappa_l`
- `kappa_theta`
- `kappa_phi`

## Inputs

You need:

- a built HEVA tree with `src/assemble`
- for `initial_frame` input, a built `src/prepare_initial_frame`
- an initial structure:
  - `Initial_frame.dat` with `--input-format initial_frame`
  - or `restart_lammps.dat` with `--input-format restart`
- AA reference data in `AAdata.csv`

If you need to build `AAdata.csv` from AA graph dumps:

```bash
python3 workflows/cg_paramopt/scripts/extract_aa_data.py \
  --input-dir /path/to/AA-graphdata \
  --output /path/to/AAdata.csv
```

If you need to convert `Initial_frame.dat` into a HEVA restart file without running a trial:

```bash
python3 workflows/cg_paramopt/scripts/prepare_initial_frame.py \
  --input /path/to/Initial_frame.dat \
  --output-dir /path/to/prepared-frame
```

## Python Environment

The optimizer scripts expect:

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install "setuptools<81" numpy pandas matplotlib scipy hyperopt
```

`hyperopt 0.2.7` still imports `pkg_resources`, so `setuptools<81` is the safest tested choice here.

## Run One Trial

Using an `Initial_frame.dat` input:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure /path/to/Initial_frame.dat \
  --input-format initial_frame \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800
```

Using a HEVA restart file:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure /path/to/restart_lammps.dat \
  --input-format restart \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800
```

Useful trial controls:

- `--relax-sweeps`
- `--frames`
- `--sample-every`
- `--resume`
- `--output-dir`
- geometry overrides such as `--l0`, `--l1`, `--theta0`, `--theta1`, `--phi33`, `--phi12`, `--phi01`, `--phi20`

## Run The Optimizer

Evaluate one fixed parameter set:

```bash
python3 workflows/cg_paramopt/scripts/optimize.py \
  --aa-data /path/to/AAdata.csv \
  --initial-structure /path/to/Initial_frame.dat \
  --input-format initial_frame \
  --single-run 2.15 1.25 1.80
```

Run a short search:

```bash
python3 workflows/cg_paramopt/scripts/optimize.py \
  --aa-data /path/to/AAdata.csv \
  --initial-structure /path/to/Initial_frame.dat \
  --input-format initial_frame \
  --max-evals 25
```

Useful optimizer controls:

- `--max-evals`
- `--frames`
- `--relax-sweeps`
- `--sample-every`
- `--runs-dir`
- `--profile`
- `--index-capacity`

## Outputs

Each optimizer run creates a timestamped directory under `workflows/cg_paramopt/runs/` unless you override `--runs-dir`.

Typical outputs include:

- `minimize-record.csv`
- `data.dat`
- trial-specific HEVA output files
- comparison plots

Each single trial directory also contains:

- `trial_config.in`
- `run_config.out`
- `restart_lammps.dat`
- `energy.dat`
- `last.dat`

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
