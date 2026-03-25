#!/usr/bin/env bash
set -euo pipefail

binary_path="$(readlink -f "${1:?missing assemble path}")"
restart_path="$(readlink -f "${2:?missing restart path}")"

workdir="$(mktemp -d /tmp/heva-regression.XXXXXX)"
output_dir="${workdir}/out"
config_path="${workdir}/regression.in"
trap 'rm -rf "${workdir}"' EXIT

mkdir -p "${output_dir}"

cat > "${config_path}" <<EOF
[capsid_geometry]
epsilon0=4200.000
kappa0=40.000
kappaPhi0=800.000
theta0=0.240
theta1=0.480
gb0=-9.480
dg12=0.300
dg01=0.100
dg20=-0.100
dg33=0.000
dg00=-0.800
dgother=-0.950

[simulation]
muCd=-10.800
ks0=0.020000
dmu=-3.500
dg=0.100

[drug]
mudrug=0.000
gdrug0=0.000
kd0=0.000000

[init]
mode=restart
restartPath=${restart_path}

[runtime]
seed=521759
outputDir=${output_dir}

[engine]
profile=extended
indexCapacity=1000000
EOF

"${binary_path}" --config "${config_path}" --max-sweeps 5 > "${workdir}/stdout.txt" 2> "${workdir}/stderr.txt"

test -s "${output_dir}/last.dat"

awk -F, '
NR != 1 { exit 1 }
$1 != 81389 { exit 1 }
$2 != 521759 { exit 1 }
$16 != 12 { exit 1 }
$17 != 30 { exit 1 }
$28 != 42 { exit 1 }
$29 != 120 { exit 1 }
$30 != 0 { exit 1 }
$31 != 1 { exit 1 }
($14 < 77.08168 || $14 > 77.08170) { exit 1 }
' "${output_dir}/last.dat"

grep -q "STOP REASON closed" "${workdir}/stdout.txt"

echo "heva regression reference passed"
