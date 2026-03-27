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
  Start from `path` or from the default `restart_lammps.dat`

- `mode = seed`
  Use `seedShape = triangle | pentamer | hexamer`

- `mode = triangle | pentamer | hexamer`
  Shorthand for deterministic seed initialization with that shape

`pentamer` and `hexamer` currently use deterministic startup paths.

## Resume Behavior
The `[runtime]` section also controls whether a restart-source is treated as a fresh start or a continuation.

- `resume = false`
  Load the structure from `init.path`, reset sweep/output state, and start a new run

- `resume = true`
  Continue from the sweep/output context stored in the restart file

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
  --input-format initial_frame \
  --epsilon0 4200 \
  --kappa0 40 \
  --kappaPhi0 800
```

That workflow first converts the AA-translated initial frame into a standard HEVA `restart_lammps.dat`, then runs `workflow = relaxation` with `init.mode = restart` and `runtime.resume = false`. The sampled geometry is written into `data.dat` for downstream optimization.

## Validated Examples
The following concrete commands were exercised during the init/resume cleanup pass.

Seeded assembly starts:

```bash
cd /mnt/d/Coding/HEVA/HEVA/src
./assemble --config /tmp/heva-seed-triangle.in
./assemble --config /tmp/heva-seed-pentamer.in
```

Restart-source fresh start versus continuation:

```bash
cd /mnt/d/Coding/HEVA/HEVA/src
./assemble --config /tmp/heva-restart-fresh-from-seed.in
./assemble --config /tmp/heva-restart-resume-from-seed.in
```

These scenarios verified:

- `mode = seed` and shorthand seed modes still start normal HEVA runs
- `mode = restart` with `resume = false` resets sweep/output state and starts fresh
- `mode = restart` with `resume = true` continues from the stored sweep state

## Useful Overrides
Config files are the main interface, but a few run controls can still be overridden on the command line:

- `--max-sweeps`
- `--run-mode`
- `--index-capacity`
- `--init`
- `--init-path`
- `--seed-config`
- `--resume`
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
