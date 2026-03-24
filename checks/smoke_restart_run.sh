#!/usr/bin/env bash
set -euo pipefail

binary_path="$(readlink -f "${1:?missing assemble path}")"
restart_path="$(readlink -f "${2:?missing restart file path}")"

workdir="$(mktemp -d /tmp/heva-smoke.XXXXXX)"
output_dir="${workdir}/out"
trap 'rm -rf "${workdir}"' EXIT

mkdir -p "${output_dir}"

pushd /tmp >/dev/null
"${binary_path}" \
  521759 4200.000 40.000 800.000 0.240 0.480 -9.480 -10.800 0.020000 -3.500 0.100 0.000 0.000 0.000000 0.300 0.100 -0.100 0.000 -0.800 -0.950 \
  --index-capacity 1000000 \
  --max-sweeps 5 \
  --restart "${restart_path}" \
  --output-dir "${output_dir}" \
  > "${workdir}/smoke.stdout" 2> "${workdir}/smoke.stderr"
popd >/dev/null

test -s "${output_dir}/energy.dat"
test -s "${output_dir}/last.dat"
grep -Eq "STOP REASON (closed|max_sweeps)" "${workdir}/smoke.stdout"

echo "heva smoke run passed"
