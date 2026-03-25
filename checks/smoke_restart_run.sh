#!/usr/bin/env bash
set -euo pipefail

binary_path="$(readlink -f "${1:?missing assemble path}")"
restart_path="$(readlink -f "${2:?missing restart file path}")"

workdir="$(mktemp -d /tmp/heva-smoke.XXXXXX)"
output_dir="${workdir}/out"
config_path="${workdir}/smoke.in"
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

pushd /tmp >/dev/null
"${binary_path}" \
  --config "${config_path}" \
  --max-sweeps 5 \
  > "${workdir}/smoke.stdout" 2> "${workdir}/smoke.stderr"
popd >/dev/null

test -s "${output_dir}/energy.dat"
test -s "${output_dir}/last.dat"
grep -Eq "STOP REASON (closed|max_sweeps)" "${workdir}/smoke.stdout"

echo "heva smoke run passed"
