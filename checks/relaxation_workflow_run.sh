#!/usr/bin/env bash
set -euo pipefail

binary_path="$(readlink -f "${1:?missing assemble path}")"
restart_path="$(readlink -f "${2:?missing restart path}")"

workdir="$(mktemp -d /tmp/heva-relaxation-check.XXXXXX)"
output_dir="${workdir}/out"
config_path="${workdir}/relaxation.in"
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
ks0=0.500000
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
maxSweeps=25
outputDir=${output_dir}
workflow=relaxation

[engine]
profile=test
EOF

"${binary_path}" \
  --config "${config_path}" \
  > "${workdir}/stdout.txt" 2> "${workdir}/stderr.txt"

test -s "${output_dir}/energy.dat"
test -s "${output_dir}/last.dat"
test -s "${output_dir}/run_config.out"
grep -q "workflow = relaxation" "${output_dir}/run_config.out"
grep -q "# WORKFLOW relaxation" "${workdir}/stdout.txt"
grep -q "STOP REASON max_sweeps" "${workdir}/stdout.txt"
grep -q "MONOMER ADDED 0" "${workdir}/stdout.txt"
grep -q "DIMER ADDED 0" "${workdir}/stdout.txt"
grep -q "TYPE CHANGED 0" "${workdir}/stdout.txt"

awk -F, 'NR == 1 { if ($1 < 25) exit 1; seen = 1 } END { exit seen ? 0 : 1 }' "${output_dir}/last.dat"

echo "heva relaxation workflow run passed"
