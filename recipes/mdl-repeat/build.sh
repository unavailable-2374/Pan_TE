#!/bin/bash
set -ex

# PORTABLE=1 drops -march=native so the binary is safe on bioconda CI and any consumer CPU.
make PORTABLE=1 -j "${CPU_COUNT}"

mkdir -p "$PREFIX/bin"
[ -x bin/mdl-repeat ] || { echo "mdl-repeat binary not produced" >&2; exit 1; }
install -m 0755 bin/mdl-repeat "$PREFIX/bin/mdl-repeat"
