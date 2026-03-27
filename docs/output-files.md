# HEVA Output Files

Each run writes artifacts into the configured `outputDir`.

## Main Files

### `run_config.out`
Canonical normalized config snapshot for that run directory.

This is the main reproducibility artifact for rerunning or auditing a specific output folder.
Paths are written relative to the run directory when that is possible, so a run
folder can usually be moved or archived without rewriting the config snapshot.

### `energy.dat`
Append-style energy and progress log written during the run.

When a run starts from sweep `0`, the first line is a CSV header with these fields:

```text
sweep,seed,seconds,epsilon,kappa,kappaPhi,theta0,theta1,gb0,mu,dmu,dg,theta2,energy,binding_energy,Nv5,Nv6,NAB,NAB_in,NCD_Hex,NCD_other,NVin,Nhein,NCD_T4_in,NCD_T3_in,NCD_T4,NCD_T3,Nv,NE,Nsurf,Nboundary
```

Each later line is one sampled state at the current sweep.

### `last.dat`
Final one-line summary of the ending run state.

It uses the same comma-separated field layout as `energy.dat`, but contains only
the final record for the run. This is the file used by the regression checks.

### `restart_lammps.dat`
Restart snapshot used for continuation runs.

### `snap_*.dat`
Intermediate geometry snapshots.

### `parameters_run.out`
Human-readable run report and derived parameter dump.

This is useful for inspection, but `run_config.out` is the authoritative
machine-readable record of the normalized run configuration.

### `geo_param.dat`
Geometry-parameter side output from the geometry layer.

## Typical Use

- rerun or audit a folder with `run_config.out`
- continue a run from `restart_lammps.dat`
- inspect progress through `energy.dat`
- inspect the final state through `last.dat`
- inspect geometry snapshots through `snap_*.dat`
