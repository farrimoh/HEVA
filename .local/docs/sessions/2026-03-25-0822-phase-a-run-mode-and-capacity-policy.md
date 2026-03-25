# Phase A: Run Mode and Capacity Policy

## Summary

Phase A was narrowed to a stability and UX fix for index-map capacity handling.

- Kept the memory policy out of the core simulation parameters by introducing `runMode=test|run|extended`.
- Left `indexCapacity` available only for advanced `extended` runs.
- Preserved backward compatibility for existing configs that already set `indexCapacity` by treating them as `extended`.
- Fixed the unsafe index-map initialization loop so the allocated capacity is fully cleared.
- Added fail-fast guards for new vertex and half-edge ids if they ever exceed the selected capacity.

## Effective Policy

- `run`: fixed preset of `2000000`
- `test`: fixed preset of `1000000`
- `extended`: user-supplied `indexCapacity`, validated to be at least `1000000`

This keeps the normal path simple while still allowing larger manual runs when needed.

## Files Updated

- Parser and usage: [simulation_config.cpp](/mnt/d/Coding/HEVA/HEVA/src/simulation_config.cpp)
- Config shape: [simulation_config.hpp](/mnt/d/Coding/HEVA/HEVA/src/simulation_config.hpp)
- Geometry safety: [geometry.cpp](/mnt/d/Coding/HEVA/HEVA/src/geometry.cpp), [geometry.hpp](/mnt/d/Coding/HEVA/HEVA/src/geometry.hpp)
- Invocation logging: [run_assembly.cpp](/mnt/d/Coding/HEVA/HEVA/src/run_assembly.cpp)
- Checks: [test_simulation_config.cpp](/mnt/d/Coding/HEVA/HEVA/checks/test_simulation_config.cpp), [index_capacity_policy.sh](/mnt/d/Coding/HEVA/HEVA/checks/index_capacity_policy.sh), [Makefile](/mnt/d/Coding/HEVA/HEVA/checks/Makefile)

## Verification

Completed successfully:

- `make` in [src](/mnt/d/Coding/HEVA/HEVA/src)
- `./test_simulation_config` in [checks](/mnt/d/Coding/HEVA/HEVA/checks)
- `make capacity_policy_run` in [checks](/mnt/d/Coding/HEVA/HEVA/checks)
- `make smoke_run` in [checks](/mnt/d/Coding/HEVA/HEVA/checks)
- `make regression_run` in [checks](/mnt/d/Coding/HEVA/HEVA/checks)
- `make invariants_run` in [checks](/mnt/d/Coding/HEVA/HEVA/checks)

The new policy regression confirms that a too-small manual capacity fails with a validation error instead of a segfault.

## Phase Run Snapshot

Run artifacts:

- Config: [phase_a_run.in](/mnt/d/Coding/HEVA/HEVA/.local/runs/2026-03-25-0822-phase-a-run/phase_a_run.in)
- Summary log: [run.stdout](/mnt/d/Coding/HEVA/HEVA/.local/runs/2026-03-25-0822-phase-a-run/run.stdout)
- Error log: [run.stderr](/mnt/d/Coding/HEVA/HEVA/.local/runs/2026-03-25-0822-phase-a-run/run.stderr)
- Energy tail: [last.dat](/mnt/d/Coding/HEVA/HEVA/.local/runs/2026-03-25-0822-phase-a-run/last.dat)
- Final restart: [restart_lammps.dat](/mnt/d/Coding/HEVA/HEVA/.local/runs/2026-03-25-0822-phase-a-run/restart_lammps.dat)
- Visualization: [heva_run_summary.svg](/mnt/d/Coding/HEVA/HEVA/.local/runs/2026-03-25-0822-phase-a-run/heva_run_summary.svg)

Observed result:

- restart-based verification run completed cleanly with `runMode=run`
- final `last.dat` row remained consistent with the reference morphology
- renderer produced a summary SVG from the phase run directory without extra post-processing

## Next

- Keep Phase B focused on coverage rather than more user-facing parameter surface area.
- Prefer the `runMode` abstraction in docs and examples; mention `indexCapacity` only as an advanced extended-run override.
