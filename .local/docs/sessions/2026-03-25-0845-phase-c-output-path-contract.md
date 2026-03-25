# Phase C: Output Path Contract Without `chdir`

## Summary

Phase C removed the remaining compatibility dependency on `chdir(output_dir)` for normal runs.

## What Changed

- Added a configured output-directory helper in [geometry.cpp](/mnt/d/Coding/HEVA/HEVA/src/geometry.cpp)
- Exposed the helper in [geometry.hpp](/mnt/d/Coding/HEVA/HEVA/src/geometry.hpp)
- Updated [run_assembly.cpp](/mnt/d/Coding/HEVA/HEVA/src/run_assembly.cpp) to open `energy.dat`, `allparam.dat`, and `parameters_run.out` by explicit path
- Updated [simulation_runner.cpp](/mnt/d/Coding/HEVA/HEVA/src/simulation_runner.cpp) so `last.dat` is also written through the configured output path

The existing dump helpers now honor the configured output directory directly, so the run no longer depends on changing the process working directory first.

## Practical Result

- wrappers can launch HEVA from any directory and still get outputs in the requested location
- restart input remains explicit
- `last.dat`, restart dumps, snapshot dumps, and parameter files all land under the chosen output directory
- the open-shell and closed-shell regression scripts continue to work without needing cwd-specific behavior

## Verification

Completed successfully:

- `make` in [src](/mnt/d/Coding/HEVA/HEVA/src)
- `make test_all` in [checks](/mnt/d/Coding/HEVA/HEVA/checks)

That verification included:

- parser checks
- invariant checks
- capacity-policy regression
- open-shell `max_sweeps` regression
- smoke run
- reference regression run

## Notes

- This keeps the external run contract cleaner without broadening the simulation parameter surface.
- The output-directory helper is still a compatibility layer around legacy dump functions, but it removes the wrapper-visible cwd coupling that Phase 3 had left in place.
