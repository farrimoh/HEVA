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

- `mode = triangle | pentamer | hexamer`
  Shorthand for deterministic seed initialization with that shape

`pentamer` and `hexamer` currently use deterministic startup paths.

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
