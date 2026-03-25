# Phase E: Geometry Memory Cleanup

## Summary

This Phase E step tightens teardown ownership inside `geometry`.

## What Changed

- Initialized more pointer members explicitly in the constructor
- Freed nested `gb` rows in the destructor
- Freed nested `gdrug` rows in the destructor
- Freed `dist_points` rows in the destructor

All of this landed in [geometry.cpp](/mnt/d/Coding/HEVA/HEVA/src/geometry.cpp).

## Why It Matters

- the destructor now matches the actual allocation shape more closely
- repeated runs and test processes leave behind less leaked heap state
- this gives the geometry module a cleaner baseline for future internal cleanup

## Verification

Completed successfully:

- `make` in [src](/mnt/d/Coding/HEVA/HEVA/src)
- `make test_all` in [checks](/mnt/d/Coding/HEVA/HEVA/checks)

## Notes

- This does not change the scientific model or run contract.
- There is still deeper modernization work available later, but the current raw-array lifetime handling is now less incomplete than before.
