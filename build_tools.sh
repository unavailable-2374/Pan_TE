#!/usr/bin/env bash
# Build the 4 custom (non-bioconda) Pan_TE tools from source into the active conda env's
# bin/, completing the conda closed loop: environment.yml provides all deps + compilers
# (cmake / rust / c-compiler / make), this script compiles the bespoke binaries.
#
#   conda activate pan_te
#   ./build_tools.sh
#
# Each tool's source is resolved as: repo-local submodule/<tool>  ->  $PAN_TE_<TOOL>_SRC
# env override  ->  git clone from the known URL. A missing source is a hard error
# (never a silent skip / fake binary).
#
# Tools built:
#   look4ltrs   (CMake, C++17 + OpenMP)            -> look4ltrs
#   FastGA      (Make)                             -> FastGA GIXmake FAtoGDB GDBtoFA ... (suite)
#   mdl-repeat  (Make, PORTABLE=1 = no -march=native, container-safe) -> mdl-repeat
#   te-looker   (cargo --release)                  -> dtr te-discover te-refine te-seed
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMOD="$HERE/submodule"
PREFIX="${CONDA_PREFIX:?activate the pan_te env first (CONDA_PREFIX unset)}"
BIN="$PREFIX/bin"
JOBS="${JOBS:-$(nproc)}"
mkdir -p "$BIN"

log()  { echo "[build_tools] $*"; }
die()  { echo "[build_tools] ERROR: $*" >&2; exit 1; }

# Resolve a tool's source dir: submodule -> env override -> clone. Echoes the path.
resolve_src() {
    local name="$1" sub="$2" envvar="$3" url="$4"
    if [ -n "${!envvar:-}" ] && [ -d "${!envvar}" ]; then echo "${!envvar}"; return; fi
    if [ -d "$SUBMOD/$sub" ] && [ -n "$(ls -A "$SUBMOD/$sub" 2>/dev/null)" ]; then
        echo "$SUBMOD/$sub"; return
    fi
    [ -n "$url" ] || die "$name source not found (submodule/$sub empty, \$$envvar unset, no clone URL). Provide it via $envvar=/path."
    log "$name: cloning $url -> submodule/$sub"
    git clone --depth 1 "$url" "$SUBMOD/$sub" >&2 || die "$name clone failed"
    echo "$SUBMOD/$sub"
}

need() { command -v "$1" >/dev/null 2>&1 || die "build tool '$1' not on PATH (is the pan_te env active with compilers installed?)"; }
need cmake; need make; need cargo; need cc

# ── look4ltrs (CMake) ───────────────────────────────────────────────────────
L4="$(resolve_src look4ltrs Look4LTRs PAN_TE_LOOK4LTRS_SRC https://github.com/unavailable-2374/Look4LTRs.git)"
log "look4ltrs: cmake build in $L4"
( cd "$L4" && cmake -S . -B build_pante -DCMAKE_BUILD_TYPE=Release >/dev/null && cmake --build build_pante --target look4ltrs -j "$JOBS" )
LB="$(find "$L4" -name look4ltrs -type f -perm -u+x | head -1)"
[ -n "$LB" ] || die "look4ltrs binary not produced"
install -m755 "$LB" "$BIN/look4ltrs"; log "look4ltrs -> $BIN/look4ltrs"

# ── FastGA (Make; install the whole suite — FastGA needs GIXmake/GDB tools at runtime) ──
FG="$(resolve_src FastGA FastGA PAN_TE_FASTGA_SRC https://github.com/thegenemyers/FASTGA.git)"
log "FastGA: make in $FG"
( cd "$FG" && make -j "$JOBS" all )
fg_n=0
for b in FastGA GIXmake GIXshow GIXrm GIXmv GIXcp FAtoGDB GDBtoFA GDBstat GDBshow ALNtoPAF ALNshow ALNreset; do
    [ -x "$FG/$b" ] && install -m755 "$FG/$b" "$BIN/$b" && fg_n=$((fg_n+1))
done
[ "$fg_n" -gt 0 ] || die "no FastGA binaries produced"; log "FastGA: installed $fg_n binaries"

# ── mdl-repeat (Make, PORTABLE=1 → no -march=native, safe across CPUs/containers) ──
MD="$(resolve_src mdl-repeat mdl-repeat PAN_TE_MDL_SRC https://github.com/unavailable-2374/mdl-repeat.git)"
log "mdl-repeat: make PORTABLE=1 in $MD"
( cd "$MD" && make PORTABLE=1 -j "$JOBS" )
[ -x "$MD/bin/mdl-repeat" ] || die "mdl-repeat binary not produced"
install -m755 "$MD/bin/mdl-repeat" "$BIN/mdl-repeat"; log "mdl-repeat -> $BIN/mdl-repeat"

# ── te-looker / dtr (cargo). No public clone URL: must be a submodule or $PAN_TE_TE_LOOKER_SRC ──
TL="$(resolve_src te-looker te-looker PAN_TE_TE_LOOKER_SRC '')"
log "te-looker: cargo build --release in $TL/core"
( cd "$TL/core" && cargo build --release )
tl_n=0
for b in dtr te-discover te-refine te-seed; do
    [ -x "$TL/core/target/release/$b" ] && install -m755 "$TL/core/target/release/$b" "$BIN/$b" && tl_n=$((tl_n+1))
done
[ "$tl_n" -ge 1 ] || die "no te-looker binaries produced"; log "te-looker: installed $tl_n binaries (dtr + discovery core)"

log "DONE — custom tools installed into $BIN"
log "verify: for t in look4ltrs FastGA mdl-repeat dtr; do command -v \$t; done"
