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
- `triangle`
- `pentamer`
- `hexamer`

### `[runtime]`
Per-run controls:

- `seed`
- `maxSweeps`
- `outputDir`

### `[engine]`
Advanced implementation/runtime controls:

- `profile`
- `indexCapacity`

Supported `profile` values:

- `test`
- `run`
- `extended`

`indexCapacity` is an advanced engine control, not a scientific sweep parameter.

## Example

```ini
[capsid_geometry]
epsilon0 = 4200.000
kappa0 = 40.000
kappaPhi0 = 800.000
theta0 = 0.240
theta1 = 0.480
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
mode = triangle

[runtime]
seed = 521759
maxSweeps = 1000
outputDir = /tmp/heva-run

[engine]
profile = extended
indexCapacity = 1000000
```
