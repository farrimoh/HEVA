#!/usr/bin/env bash
set -euo pipefail

binary_path="$(readlink -f "${1:?missing assemble path}")"
restart_path="$(readlink -f "${2:?missing open-shell restart path}")"

workdir="$(mktemp -d /tmp/heva-open-shell-check.XXXXXX)"
output_dir="${workdir}/out"
config_path="${workdir}/open_shell.in"
trap 'rm -rf "${workdir}"' EXIT

mkdir -p "${output_dir}"

cat > "${config_path}" <<EOF
seed=521759
epsilon0=4200.000
kappa0=40.000
kappaPhi0=800.000
theta0=0.240
theta1=0.480
gb0=-9.480
muCd=-10.800
ks0=0.020000
dmu=-3.500
dg=0.100
mudrug=0.000
gdrug0=0.000
kd0=0.000000
dg12=0.300
dg01=0.100
dg20=-0.100
dg33=0.000
dg00=-0.800
dgother=-0.950
runMode=test
restart=${restart_path}
outputDir=${output_dir}
EOF

"${binary_path}" \
  --config "${config_path}" \
  --max-sweeps 5 \
  > "${workdir}/stdout.txt" 2> "${workdir}/stderr.txt"

test -s "${output_dir}/energy.dat"
test -s "${output_dir}/last.dat"
grep -q "STOP REASON max_sweeps" "${workdir}/stdout.txt"

awk -F, 'NR == 1 { if ($1 < 5) exit 1; seen = 1 } END { exit seen ? 0 : 1 }' "${output_dir}/last.dat"

echo "heva open-shell max-sweeps run passed"
