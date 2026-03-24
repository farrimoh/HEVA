#!/usr/bin/env bash
set -euo pipefail

binary_path="$(readlink -f "${1:?missing assemble path}")"
restart_path="$(readlink -f "${2:?missing restart file path}")"

workdir="$(mktemp -d /tmp/heva-smoke.XXXXXX)"
trap 'rm -rf "${workdir}"' EXIT

cp "${restart_path}" "${workdir}/restart_lammps.dat"

pushd "${workdir}" >/dev/null
"${binary_path}" \
  521759 4200.000 40.000 800.000 0.240 0.480 -9.480 -10.800 0.020000 -3.500 0.100 0.000 0.000 0.000000 0.300 0.100 -0.100 0.000 -0.800 -0.950 \
  1000000 5 \
  > smoke.stdout 2> smoke.stderr
popd >/dev/null

test -s "${workdir}/energy.dat"
test -s "${workdir}/last.dat"
grep -Eq "STOP REASON (closed|max_sweeps)" "${workdir}/smoke.stdout"

echo "heva smoke run passed"
