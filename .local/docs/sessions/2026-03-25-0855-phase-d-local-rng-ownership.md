# Phase D: Local RNG Ownership in `main`

## Summary

This Phase D step removes one small but real piece of driver-level global state.

## What Changed

- Removed the file-scope `gsl_rng *r` from [run_assembly.cpp](/mnt/d/Coding/HEVA/HEVA/src/run_assembly.cpp)
- Kept the RNG owned locally inside `main`
- Passed that local RNG through the existing setup/run/finalization calls

## Why It Matters

- the driver now has clearer ownership of the RNG lifecycle
- there is less hidden mutable state at file scope
- this is a safer base for future refactors around wrapper execution and testability

## Verification

Completed successfully:

- `make` in [src](/mnt/d/Coding/HEVA/HEVA/src)
- `make test_all` in [checks](/mnt/d/Coding/HEVA/HEVA/checks)

## Notes

- This is intentionally a narrow cleanup step, not a broad Monte Carlo API rewrite.
- The deeper Phase D/E cleanup work is still open, but the driver is now a little less coupled than before.
