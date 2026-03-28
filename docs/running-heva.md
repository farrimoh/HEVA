# Running HEVA

## Requirements
HEVA is currently tested as a Linux/WSL command-line build with:

- `g++`
- `make`
- GNU Scientific Library (GSL) development headers and libraries

If `pkg-config` can locate `gsl`, `src/Makefile` will use that automatically.
Otherwise the default fallback link flags are `-lgsl -lgslcblas -lm`.

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
The repository includes a committed smoke-test config at:

```text
checks/fixtures/smoke_config.in
```

Run it with:

```bash
cd src
./assemble --config ../checks/fixtures/smoke_config.in
```

Paths inside the config file are resolved relative to the config file location,
not relative to the shell's current working directory.

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

For CG-ParamOpt usage and examples, see:

- `workflows/cg_paramopt/README.md`

## Assembly Opt Workflow
The HEVA-native assembly optimization workflow lives under:

```text
workflows/assembly_opt/
```

For usage and examples, see:

- `workflows/assembly_opt/README.md`

## Validated Examples
The following concrete commands were exercised during the init/resume cleanup pass.

Seeded assembly starts:

```bash
cd src
./assemble --config /path/to/heva-seed-triangle.in
./assemble --config /path/to/heva-seed-pentamer.in
```

Restart-source fresh start versus continuation:

```bash
cd src
./assemble --config /path/to/heva-restart-fresh-from-seed.in
./assemble --config /path/to/heva-restart-resume-from-seed.in
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
./assemble --config ../checks/fixtures/smoke_config.in --max-sweeps 500
```

## Development Checks
From the repository root:

```bash
cd checks
make test_all
```

This runs config, invariant, smoke, regression, capacity-policy, open-shell, and deterministic-seed checks.
