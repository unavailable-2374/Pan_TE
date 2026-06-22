#!/bin/bash
set -ex

# The cargo workspace/crate lives under core/ (package `te-core`).
cd core

export CARGO_NET_GIT_FETCH_WITH_CLI=true
cargo build --release -j "${CPU_COUNT}"

mkdir -p "$PREFIX/bin"
n=0
for b in dtr te-discover te-refine te-seed; do
    if [ -x "target/release/$b" ]; then
        install -m 0755 "target/release/$b" "$PREFIX/bin/$b"
        n=$((n + 1))
    fi
done

# dtr is the Step-4 entry point the Pan_TE pipeline invokes; its absence is a hard error.
[ -x "$PREFIX/bin/dtr" ] || { echo "te-looker: dtr binary not produced" >&2; exit 1; }
[ "$n" -ge 1 ] || { echo "te-looker: no binaries produced" >&2; exit 1; }
