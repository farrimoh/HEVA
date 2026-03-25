# HEVA Architecture

## Core Modules

### `src/geometry.cpp` / `src/geometry.hpp`
Owns:

- half-edge data structure
- geometry updates
- restart I/O
- deterministic initial builders such as triangle, pentamer, and hexamer

### `src/montecarlo.cpp` / `src/montecarlo.hpp`
Owns:

- Monte Carlo moves
- acceptance logic
- seed dispatch (`make_seed`, `make_seed_T3`)

### `src/simulation_config.cpp` / `src/simulation_config.hpp`
Owns:

- sectioned config parsing
- normalized config rendering
- runtime/engine validation

### `src/simulation_setup.cpp` / `src/simulation_setup.hpp`
Maps parsed config into geometry/model state.

### `src/simulation_runner.cpp` / `src/simulation_runner.hpp`
Owns:

- restart and seed initialization flow
- simulation loop settings
- run loop and finalization

### `src/run_assembly.cpp`
Main executable entry point and run-directory setup.

### `src/tri_tri_intersect.cpp` / `src/tri_tri_intersect.hpp`
Triangle intersection checks.

## Checks

The `checks/` directory contains:

- config parser checks
- geometry invariant checks
- index-capacity policy checks
- open-shell max-sweeps checks
- deterministic seed-init checks
- smoke restart runs
- regression reference runs

Run them with:

```bash
cd checks
make test_all
```
