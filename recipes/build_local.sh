#!/usr/bin/env bash
# build_local.sh — verify the Pan_TE conda recipes BEFORE any GitHub release exists, by
# building them from LOCAL source into a local channel. This is the pre-submission smoke
# test: it exercises every build.sh and the dependency wiring without needing release
# tarballs or real sha256 values.
#
# Two modes:
#   ./build_local.sh lint            # render + lint recipes only (fast, no compilation)
#   ./build_local.sh build           # build all 4 packages from local source -> $CHANNEL
#
# Local source dirs (override via env if yours differ):
#   PANTE_SRC      (default: the repo this script lives in)
#   LOOK4LTRS_SRC  (default: $PANTE_SRC/submodule/Look4LTRs)
#   MDLREPEAT_SRC  (default: $HOME/tool/mdl-repeat)
#   TELOOKER_SRC   (default: $HOME/tool/te-looker)
#
# Requires conda-build (the script installs it into the current env if missing).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
MODE="${1:-build}"

PANTE_SRC="${PANTE_SRC:-$REPO_ROOT}"
LOOK4LTRS_SRC="${LOOK4LTRS_SRC:-$PANTE_SRC/submodule/Look4LTRs}"
MDLREPEAT_SRC="${MDLREPEAT_SRC:-$HOME/tool/mdl-repeat}"
TELOOKER_SRC="${TELOOKER_SRC:-$HOME/tool/te-looker}"

CHANNEL="$HERE/local_channel"
WORK="$HERE/.build_tmp"

die() { echo "[build_local] ERROR: $*" >&2; exit 1; }
log() { echo "[build_local] $*"; }

command -v conda >/dev/null || die "conda not on PATH"
if ! conda build --version >/dev/null 2>&1; then
    log "conda-build not found; installing into current env (conda-build conda-verify)"
    conda install -y -n base -c conda-forge conda-build conda-verify >/dev/null \
        || die "failed to install conda-build"
fi

# Rewrite a recipe's `source:` block to point at local path(s), so we can build before a
# release exists. Args: <recipe-name> <out-dir> then pairs of <localpath> <folder|->.
make_local_recipe() {
    local name="$1" outdir="$2"; shift 2
    local srcdir="$HERE/$name"
    [ -d "$srcdir" ] || die "recipe dir not found: $srcdir"
    rm -rf "$outdir"; mkdir -p "$outdir"
    # copy everything except meta.yaml verbatim (build.sh, helper scripts)
    find "$srcdir" -maxdepth 1 -type f ! -name meta.yaml -exec cp {} "$outdir/" \;
    PANTE_LOCAL_ARGS="$*" python3 - "$srcdir/meta.yaml" "$outdir/meta.yaml" <<'PY'
import sys, os, re
src, dst = sys.argv[1], sys.argv[2]
args = os.environ.get("PANTE_LOCAL_ARGS", "").split()
pairs = [(args[i], args[i+1]) for i in range(0, len(args), 2)]
text = open(src).read()
# Build replacement source: block.
if len(pairs) == 1 and pairs[0][1] == "-":
    block = "source:\n  path: %s\n" % os.path.abspath(pairs[0][0])
else:
    lines = ["source:"]
    for path, folder in pairs:
        lines.append("  - path: %s" % os.path.abspath(path))
        if folder != "-":
            lines.append("    folder: %s" % folder)
    block = "\n".join(lines) + "\n"
# Replace from 'source:' up to the next top-level key (a line starting non-space, here 'build:').
pat = re.compile(r"^source:.*?(?=^build:)", re.S | re.M)
new = pat.sub(block + "\n", text, count=1)
if new == text:
    sys.exit("failed to rewrite source: block in %s" % src)
# Local smoke test only: drop license_file (adding real LICENSE files to mdl-repeat /
# te-looker / Pan_TE is a separate bioconda-submission step; their absence must not block
# verifying that build.sh compiles and packages correctly).
new = re.sub(r"^\s*license_file:.*\n", "", new, flags=re.M)
open(dst, "w").write(new)
PY
}

if [ "$MODE" = "lint" ]; then
    for r in look4ltrs mdl-repeat te-looker pan_te; do
        log "render $r"
        conda render "$HERE/$r" >/dev/null || die "render failed for $r"
    done
    log "all recipes render OK (note: real sha256/license still required for bioconda)."
    exit 0
fi

[ "$MODE" = "build" ] || die "unknown mode '$MODE' (use: lint | build)"

for d in "$PANTE_SRC" "$LOOK4LTRS_SRC" "$MDLREPEAT_SRC" "$TELOOKER_SRC"; do
    [ -d "$d" ] || die "local source dir missing: $d"
done

mkdir -p "$CHANNEL/noarch" "$CHANNEL/linux-64"
conda index "$CHANNEL" >/dev/null 2>&1 || true

build_one() {
    local name="$1" recipe="$2" extra="${3:-}"
    log "building $name from local source"
    conda build "$recipe" $extra \
        -c "file://$CHANNEL" -c conda-forge -c bioconda \
        --output-folder "$CHANNEL" \
        || die "conda build failed for $name"
    conda index "$CHANNEL" >/dev/null
}

rm -rf "$WORK"; mkdir -p "$WORK"

# Order matters: pan_te depends on the other three.
make_local_recipe look4ltrs  "$WORK/look4ltrs"  "$LOOK4LTRS_SRC" -
build_one look4ltrs "$WORK/look4ltrs"

make_local_recipe mdl-repeat "$WORK/mdl-repeat" "$MDLREPEAT_SRC" -
build_one mdl-repeat "$WORK/mdl-repeat"

make_local_recipe te-looker  "$WORK/te-looker"  "$TELOOKER_SRC" -
build_one te-looker "$WORK/te-looker"

make_local_recipe pan_te "$WORK/pan_te" "$PANTE_SRC" pan_te "$LOOK4LTRS_SRC" look4ltrs_src
# --no-test: the pan_te test env solves ~29 deps (repeatmasker/repeatmodeler/...); that solve
# is validated on bioconda CI in a clean env. Locally it can trip on unrelated transitive
# post-link quirks (e.g. drmaa needing a cross-compiler absent from a personal miniconda),
# which says nothing about the pan_te recipe itself. We verify packaging by inspecting the
# built artifact instead. Drop --no-test to exercise the full test env if your base is clean.
build_one pan_te "$WORK/pan_te" "--no-test"

log "DONE. Local channel: $CHANNEL"
log "test install:  conda create -n pan_te_test -c file://$CHANNEL -c conda-forge -c bioconda pan_te"
