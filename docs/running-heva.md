# Running HEVA

## Build
From the repository root:

```bash
cd src
make clean
make
```

The executable is:

```bash
src/assemble
```

## Main Entry Point
HEVA now runs through sectioned config files:

```bash
./assemble --config PATH
```

Run from `src/` or use the full path to `src/assemble`.

## Example Trial
The repository includes a sample config at:

```text
.local/docs/fast_capsid_trial.in
```

Run it with:

```bash
cd src
./assemble --config ../.local/docs/fast_capsid_trial.in
```

## Initialization Modes
The `[init]` section controls how a run starts.

- `mode = restart`
  Start from `restartPath`

- `mode = seed`
  Use `seedShape = triangle | pentamer | hexamer`

- `mode = initial_frame`
  Import a CG-ParamOpt-style `Initial_frame.dat` through the compatibility reader

- `mode = legacy_lammps`
  Backward-compatible alias for `initial_frame`

- `mode = triangle | pentamer | hexamer`
  Shorthand for deterministic seed initialization with that shape

`pentamer` and `hexamer` currently use deterministic startup paths.

`initial_frame` is useful when you want HEVA to relax and sample a prebuilt structure rather than grow one through the assembly workflow.

## Runtime Workflows
The `[runtime]` section also controls which Monte Carlo loop HEVA uses.

- `workflow = assembly`
  Full assembly moves, including growth and removal attempts

- `workflow = relaxation`
  Relaxation-only moves via `move_vertex`, which is the mode used by `workflows/cg_paramopt`

## CG ParamOpt Workflow
The HEVA-native CG parameter workflow lives under:

```text
workflows/cg_paramopt/
```

Typical usage from the repository root:

```bash
python3 workflows/cg_paramopt/scripts/run_trial.py \
  --initial-structure /path/to/Initial_frame.dat \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800
```

That workflow generates a sectioned HEVA config, runs `workflow = relaxation`, and writes legacy-style sampled geometry into `data.dat` for downstream optimization.

## Useful Overrides
Config files are the main interface, but a few run controls can still be overridden on the command line:

- `--max-sweeps`
- `--run-mode`
- `--index-capacity`
- `--init`
- `--seed-config`
- `--restart`
- `--output-dir`

Example:

```bash
cd src
./assemble --config ../checks/fixtures/smoke_config.in --max-sweeps 5
```

## Development Checks
From the repository root:

```bash
cd checks
make test_all
```

This runs config, invariant, smoke, regression, capacity-policy, open-shell, and deterministic-seed checks.
