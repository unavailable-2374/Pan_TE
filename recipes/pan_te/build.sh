#!/bin/bash
set -ex

# Two sources are unpacked side by side: $SRC_DIR/pan_te (main tree) and
# $SRC_DIR/look4ltrs_src (only for build_ltr_library.py).
SRC="$SRC_DIR/pan_te"
PKG="$PREFIX/share/pan_te"

mkdir -p "$PKG" "$PREFIX/bin"

# --- 1. install the pipeline tree (scripts + python subpackages) into share/pan_te ---
cp -r "$SRC/bin/." "$PKG/"

# strip caches / logs / local tooling config / dev-only scaffolding that should never ship
find "$PKG" -type d \( -name '__pycache__' -o -name '.claude' -o -name 'tests' \) \
    -prune -exec rm -rf {} + 2>/dev/null || true
find "$PKG" -type f \( -name '*.pyc' -o -name '*.log' -o -name '*.DONE' \) -delete 2>/dev/null || true
rm -rf "$PKG/Refiner/cache" "$PKG/Refiner/checkpoints" "$PKG/Refiner_mdl/checkpoints" 2>/dev/null || true

# ensure the script entry points are executable (cp usually preserves this, be explicit)
for f in Pan_TE clean_seq_fast index LTR_detect build_mdl run_Classifier renameTE \
         Combine_for_Two process_for_classify.py decode_gfa.pl; do
    [ -e "$PKG/$f" ] && chmod 0755 "$PKG/$f" || true
done

# --- 2. vendor build_ltr_library.py at the path Pan_TE resolves it from ---
# Pan_TE uses BIN_DIR/../submodule/Look4LTRs/build_ltr_library.py, and BIN_DIR is the
# real script dir ($PREFIX/share/pan_te), so the target is $PREFIX/share/submodule/...
mkdir -p "$PREFIX/share/submodule/Look4LTRs"
BLB="$(find "$SRC_DIR/look4ltrs_src" -name build_ltr_library.py -type f | head -1)"
[ -n "$BLB" ] || { echo "build_ltr_library.py not found in look4ltrs source" >&2; exit 1; }
install -m 0644 "$BLB" "$PREFIX/share/submodule/Look4LTRs/build_ltr_library.py"

# --- 3. PATH wrapper. A plain symlink would break Pan_TE: it uses
#        os.path.dirname(os.path.abspath(__file__)) (NOT realpath), so the symlink target
#        dir would be wrong. A wrapper that execs the real script keeps __file__ correct. ---
cat > "$PREFIX/bin/Pan_TE" <<'EOF'
#!/usr/bin/env bash
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec "$DIR/share/pan_te/Pan_TE" "$@"
EOF
chmod 0755 "$PREFIX/bin/Pan_TE"

# --- 4. post-install data fetcher ---
install -m 0755 "$RECIPE_DIR/pan_te-setup-data" "$PREFIX/bin/pan_te-setup-data"
