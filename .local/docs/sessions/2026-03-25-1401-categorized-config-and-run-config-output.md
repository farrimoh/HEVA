# 2026-03-25 14:01 Categorized Config And Run-Config Output

## Summary

Replaced the old flat HEVA config model with a sectioned configuration format and made each run write a canonical `run_config.out` into its output folder.

## New Config Structure

- `[capsid_geometry]`
  Mechanical and interface terms:
  `epsilon0`, `kappa0`, `kappaPhi0`, `theta0`, `theta1`, `gb0`, `dg12`, `dg01`, `dg20`, `dg33`, `dg00`, `dgother`
- `[simulation]`
  Sweep-relevant assembly controls:
  `muCd`, `ks0`, `dmu`, `dg`
- `[drug]`
  Drug-related parameters:
  `mudrug`, `gdrug0`, `kd0`
- `[init]`
  Initialization choice and inputs:
  `mode`, `seedShape`, `restartPath`
- `[runtime]`
  Run-instance controls:
  `seed`, `maxSweeps`, `outputDir`
- `[engine]`
  Advanced execution controls:
  `profile`, `indexCapacity`

## Code Changes

- Refactored `SimulationConfig` into nested categories in `src/simulation_config.hpp`.
- Rewrote config parsing in `src/simulation_config.cpp` to require `--config PATH` and the sectioned format.
- Removed the legacy positional/flat config parsing path.
- Rewired `src/simulation_setup.cpp`, `src/simulation_runner.cpp`, and `src/run_assembly.cpp` to use the categorized config fields.
- Replaced the fragile `allparam.dat` branch with a canonical `run_config.out` written into the run directory on every run.
- Kept `parameters_run.out` as the verbose human-readable run report.

## Fixture And Test Updates

- Updated all check configs and generated check inputs to the new sectioned format.
- Updated `checks/test_simulation_config.cpp` to validate nested parsing and the new required `--config` flow.
- Updated `.local/docs/fast_capsid_trial.in` to the new categorized format.

## Validation

- `make` in `src/`
- `make test_simulation_config` in `checks/`
- `make test_all` in `checks/`

## Notes

This is an intentional format break. The new design treats reproducible sweep parameters, runtime controls, and low-level engine controls as different classes of input instead of saving them in one flat parameter bag.
