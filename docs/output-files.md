# HEVA Output Files

Each run writes artifacts into the configured `outputDir`.

## Main Files

### `run_config.out`
Canonical normalized config snapshot for that run directory.

This is the main reproducibility artifact for rerunning or auditing a specific output folder.

### `energy.dat`
Append-style energy and progress log written during the run.

### `last.dat`
Final one-line summary of the ending run state.

### `restart_lammps.dat`
Restart snapshot used for continuation runs.

### `snap_*.dat`
Intermediate geometry snapshots.

### `parameters_run.out`
Human-readable run report and derived parameter dump.

### `geo_param.dat`
Geometry-parameter side output from the geometry layer.

## Typical Use

- rerun or audit a folder with `run_config.out`
- continue a run from `restart_lammps.dat`
- inspect progress through `energy.dat`
- inspect the final state through `last.dat`
- inspect geometry snapshots through `snap_*.dat`
