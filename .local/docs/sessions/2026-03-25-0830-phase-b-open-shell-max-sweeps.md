# Phase B: Open-Shell Max-Sweeps Coverage

## Summary

Phase B focused on closing the specific automated coverage gap that had been called out in the older roadmap notes: an explicitly open-shell path that reaches `STOP REASON max_sweeps`.

## What Landed

- Added a committed open-shell restart fixture: [open_shell_restart_lammps.dat](/mnt/d/Coding/HEVA/HEVA/checks/fixtures/open_shell_restart_lammps.dat)
- Added a dedicated regression script: [open_shell_max_sweeps.sh](/mnt/d/Coding/HEVA/HEVA/checks/open_shell_max_sweeps.sh)
- Updated the check harness: [Makefile](/mnt/d/Coding/HEVA/HEVA/checks/Makefile)

The fixture was generated from a small deterministic `triangle` seed run and captures a minimal shell that remains open across a short continuation run.

## What The New Check Proves

- restart-based file-driven runs can load a small open shell cleanly
- the simulation loop stays open long enough to stop on `max_sweeps`
- `last.dat` is written with a sweep count consistent with the requested limit
- the check suite now covers both closed-shell and open-shell short-run paths

## Verification

Completed successfully:

- `make open_shell_run` in [checks](/mnt/d/Coding/HEVA/HEVA/checks)
- `make test_all` in [checks](/mnt/d/Coding/HEVA/HEVA/checks)

That `test_all` pass included:

- config parser checks
- geometry invariants
- capacity policy regression
- open-shell `max_sweeps` regression
- restart smoke run
- reference regression run

## Notes

- The new fixture keeps Phase B focused on behavior coverage instead of adding more user-facing options.
- This also gives later wrapper/output work a stable open-shell example for short controlled runs.
