# HEVA Config Format

HEVA uses sectioned config files passed with:

```bash
./assemble --config PATH
```

## Sections

### `[capsid_geometry]`
Mechanical and interface parameters for the capsid model:

- `epsilon0`
- `kappa0`
- `kappaPhi0`
- `theta0`
- `theta1`
- `theta2` optional, defaults to `theta0`
- `theta3` optional, defaults to `theta0`
- `l0` optional, defaults to `1.05`
- `l1` optional, defaults to `0.95`
- `phi33` optional, defaults to `1.05`
- `phi12` optional, defaults to `1.17`
- `phi01` optional, defaults to `0.98`
- `phi20` optional, defaults to `1.05`
- `gb0`
- `dg12`
- `dg01`
- `dg20`
- `dg33`
- `dg00`
- `dgother`

### `[simulation]`
Assembly and sweep-relevant parameters:

- `muCd`
- `ks0`
- `dmu`
- `dg`

### `[drug]`
Drug-related parameters:

- `mudrug`
- `gdrug0`
- `kd0`

### `[init]`
Initialization behavior:

- `mode`
- `seedShape`
- `restartPath`

Supported `mode` values:

- `restart`
- `seed`
- `initial_frame`
- `legacy_lammps` alias for `initial_frame`
- `triangle`
- `pentamer`
- `hexamer`

`initial_frame` is the compatibility path for CG-ParamOpt-style `Initial_frame.dat` inputs. It reads the imported structure, rescales the coordinates the same way as the historical workflow, reconstructs the half-edge topology, and assigns the four edge types before the run begins.

### `[runtime]`
Per-run controls:

- `seed`
- `maxSweeps`
- `outputDir`
- `workflow`

### `[engine]`
Advanced implementation/runtime controls:

- `profile`
- `indexCapacity`

Supported `profile` values:

- `test`
- `run`
- `extended`

`indexCapacity` is an advanced engine control, not a scientific sweep parameter.

### `[cg_paramopt]`
Optional controls for relaxation-only geometry sampling used by the CG parameter workflow:

- `sampleStartSweep`
- `sampleEvery`
- `sampleOutputPath`

When present, this section is intended to be combined with `runtime.workflow = relaxation`.

## Example

```ini
[capsid_geometry]
epsilon0 = 4200.000
kappa0 = 40.000
kappaPhi0 = 800.000
theta0 = 0.240
theta1 = 0.480
theta2 = 0.200
theta3 = 0.200
l0 = 1.050
l1 = 0.950
phi33 = 1.050
phi12 = 1.170
phi01 = 0.980
phi20 = 1.050
gb0 = -9.480
dg12 = 0.300
dg01 = 0.100
dg20 = -0.100
dg33 = 0.000
dg00 = -0.800
dgother = -0.950

[simulation]
muCd = -10.800
ks0 = 0.020000
dmu = -3.500
dg = 0.100

[drug]
mudrug = 0.000
gdrug0 = 0.000
kd0 = 0.000000

[init]
mode = initial_frame
restartPath = /path/to/Initial_frame.dat

[runtime]
seed = 521759
maxSweeps = 55000
outputDir = /tmp/heva-run
workflow = relaxation

[engine]
profile = extended
indexCapacity = 1000000

[cg_paramopt]
sampleStartSweep = 50000
sampleEvery = 1
sampleOutputPath = data.dat
```
