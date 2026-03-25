# 2026-03-25 18:56 Seed Dispatch Folded Into Monte Carlo

## Summary

Removed the tiny standalone seed-initialization module and moved `make_seed()` / `make_seed_T3()` back into the Monte Carlo implementation.

## Changes

- Moved `make_seed()` and `make_seed_T3()` into `src/montecarlo.cpp`.
- Declared the seed entry points in `src/montecarlo.hpp`.
- Moved `make_initial_triangle()`, `make_initial_pentamer()`, and `make_initial_hexamer()` declarations into `src/geometry.hpp`, where the builders logically belong.
- Removed `src/seed_initialization.cpp` and `src/seed_initialization.hpp`.
- Removed the extra source file from `src/Makefile`.
- Updated `src/simulation_runner.cpp` and `src/geometry.cpp` includes accordingly.

## Validation

- `make clean && make` in `src/`
- `make test_all` in `checks/`

## Notes

This keeps the deterministic seed behavior from the earlier fix, but trims back a module that had become only a thin dispatch wrapper.
