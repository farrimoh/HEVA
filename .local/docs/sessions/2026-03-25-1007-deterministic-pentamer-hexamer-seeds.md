# 2026-03-25 10:07 Deterministic Pentamer And Hexamer Seeds

## Summary

Restored deterministic public seed initialization for `pentamer` and `hexamer`, and added explicit regression coverage so both modes build, run one sweep, and write a valid `last.dat` without crashing.

## Changes

- `src/seed_initialization.cpp`
  - Routed `make_seed()` to `make_initial_pentamer(g)` for `pentamer`.
  - Routed `make_seed()` to `make_initial_hexamer(g)` for `hexamer`.
  - Kept the RNG argument unused for these deterministic seed modes.
- `src/geometry.cpp`
  - Re-enabled the deterministic `geometry::make_hexamer()` builder.
  - Updated `force_add_dimer()` to carry forward boundary indices and boundary links before `update_boundary()`.
  - Updated `force_add_monomer()` to mirror the live boundary bookkeeping used by `add_monomer()`.
- `src/simulation_runner.cpp`
  - Changed failed seed init reporting to `Failed to initialize seed configuration: ...` so valid names are not mislabeled as unknown configs.
- `checks/Makefile`
  - Added `seed_init_run` to the check targets.
- `checks/seed_init_modes.sh`
  - Added deterministic startup coverage for `init=pentamer` and `init=hexamer`.
  - Verified the emitted `last.dat` records the expected vertex counts for both seed modes.

## Validation

- `make` in `src/`
- `make seed_init_run` in `checks/`
- `make test_all` in `checks/`

## Notes

The deterministic builders were present in older code paths, but the helper functions they relied on had drifted behind the newer boundary-topology rules. The fix was to bring `force_add_dimer()` and `force_add_monomer()` up to the same boundary bookkeeping contract used by the current accepted add/remove paths.
