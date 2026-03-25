#!/usr/bin/env bash
set -euo pipefail

binary_path="$(readlink -f "${1:?missing assemble path}")"

run_seed_mode() {
  local mode="$1"
  local expected_nv="$2"
  local workdir
  workdir="$(mktemp -d "/tmp/heva-seed-${mode}.XXXXXX")"
  local output_dir="${workdir}/out"
  local config_path="${workdir}/${mode}.in"
  trap 'rm -rf "${workdir}"' RETURN

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
mode=${mode}

[runtime]
seed=521759
outputDir=${output_dir}

[engine]
profile=test
EOF

  "${binary_path}" --config "${config_path}" --max-sweeps 1 > "${workdir}/stdout.txt" 2> "${workdir}/stderr.txt"

  test -s "${output_dir}/last.dat"
  awk -F, -v expected_nv="${expected_nv}" 'NR == 1 { if ($30 != expected_nv) exit 1 }' "${output_dir}/last.dat"
}

run_seed_mode pentamer 5
run_seed_mode hexamer 6

echo "heva deterministic seed init passed"
