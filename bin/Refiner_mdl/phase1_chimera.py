"""
Phase 1 (rewrite) — Module M3: MSA-based, BLAST-free chimera detection.

mdl-repeat occasionally merges two distinct families into one consensus (an A+B
chimera). The old detector (phase1_consensus.detect_chimera) read this off blastn
qstart/qend coverage. This rewrite is fully decoupled from blastn: it works only on
the trimmed per-family MSA produced by M2 (phase1_boundary.refine_family), reading
the chimera signature directly off the alignment geometry.

Detection is a recursive **AND-gate** of two independent signals on the kept span:

  Signal 1 — occupancy valley: a sustained interior run of columns whose occupancy
             drops below `chimera_occupancy_ratio * median(occupancy)`. This marks
             where the recruited copies stop overlapping (the junction between the
             two merged families' instance pools).

  Signal 2 — copy-span bimodality coinciding with the valley: the copies fall into a
             left-reaching group and a right-reaching group, with few copies spanning
             across. A valley *without* two disjoint span groups (e.g. an AT-rich /
             indel-rich interior dip where every copy still spans the whole element)
             is NOT a chimera and is left intact.

Both must agree before a split is made. This is deliberate high-precision policy: a
false split permanently fractures a real family, so we favor under-splitting over
over-splitting (strategy §6 / plan §1.4, §6). The split point is mapped from the
valley column to a consensus-base coordinate by **re-walking the exact M2
build_consensus emit rules** (we call build_consensus itself, never a re-derived
copy), so the fragment sequences and the consensus coordinate can never drift.

The transient `_aln = (rows, metas, left, right)` handle attached by refine_family is
required input; `strip_aln` removes it before any record leaves Phase 1 (it carries
the full MSA and must never be pickled into a checkpoint).

See REFINE_IMPLEMENTATION_PLAN.md §1.4, §6, §9 (M3 config), §11, §12.M3.
"""

import logging
from typing import Dict, List, Tuple

import numpy as np

from phase1_align import occupancy_profile
from phase1_boundary import build_consensus

logger = logging.getLogger(__name__)

_GAP = frozenset('-.')


# ═══════════════════════════════════════════════════════════════════════
# Transient-handle / fragment helpers
# ═══════════════════════════════════════════════════════════════════════

def strip_aln(rec: Dict) -> Dict:
    """Return a copy of `rec` without the transient '_aln' MSA handle.

    Called on every record before it leaves Phase 1: '_aln' holds the full alignment
    (rows + metas) and must NOT reach Phase 2 or be pickled into a checkpoint. A
    record without '_aln' is returned unchanged (still copied, to avoid mutating the
    caller's dict)."""
    if '_aln' in rec:
        out = dict(rec)
        del out['_aln']
        return out
    return rec


def make_fragment(rec: Dict, sub_seq: str, sub_rows_slice: List[str],
                  suffix: str) -> Dict:
    """Build a child chimera-fragment record from a column sub-slice of the parent.

    The child carries the same per-copy `metas` (rows are column-sliced, not
    dropped, so every copy is still represented) and a fresh '_aln' spanning the full
    sub-slice ([0, len-1]); the recursion re-checks the sub-alignment without any
    re-extraction (strictly cheaper and consistent). `sub_seq` is the parent
    consensus restricted to this column block — it was produced by the same
    build_consensus emit rules, so it stays in lock-step with '_aln'."""
    metas = rec['_aln'][1]
    out = dict(rec)
    out.pop('_aln', None)
    out['id'] = rec['id'] + suffix
    out['sequence'] = sub_seq
    out['actual_length'] = len(sub_seq)
    out['is_chimera_fragment'] = True
    sub_len = len(sub_rows_slice[0]) if sub_rows_slice else 0
    out['_aln'] = (sub_rows_slice, metas, 0, sub_len - 1)
    return out


# ═══════════════════════════════════════════════════════════════════════
# Signal 1 — interior occupancy valley
# ═══════════════════════════════════════════════════════════════════════

def _longest_interior_valley(occ: np.ndarray, ratio: float,
                             margin: int) -> Tuple[int, int]:
    """Find the longest interior run of low-occupancy columns.

    A column is in a valley when occ < ratio * median(occ). The first/last `margin`
    columns are excluded (a true chimera junction is interior; terminal low-occupancy
    is the boundary that M2 trim already governs). Returns (run_length, run_start) of
    the longest qualifying interior run, or (0, -1) if none. `occ` is indexed in the
    kept-span coordinate frame (column 0 == span left)."""
    w = len(occ)
    if w == 0:
        return 0, -1
    med = float(np.median(occ))
    threshold = ratio * med
    valley = occ < threshold
    # Interior only: suppress the terminal `margin` columns on each side.
    lo = margin
    hi = w - margin
    best_len, best_start = 0, -1
    i = lo
    while i < hi:
        if valley[i]:
            j = i
            while j < hi and valley[j]:
                j += 1
            run_len = j - i
            if run_len > best_len:
                best_len, best_start = run_len, i
            i = j
        else:
            i += 1
    return best_len, best_start


# ═══════════════════════════════════════════════════════════════════════
# Signal 2 — copy-span bimodality at the valley
# ═══════════════════════════════════════════════════════════════════════

def _span_bounds(row: str) -> Tuple[int, int]:
    """First and last non-gap column of one MSA row, or (-1, -1) if all gaps.

    Non-gap (covers '-' and '.') includes N: a copy "reaches" a column when it has
    any aligned residue there, ambiguous or not — this is span coverage, not base
    quality."""
    first = -1
    last = -1
    for i, b in enumerate(row):
        if b not in _GAP:
            if first < 0:
                first = i
            last = i
    return first, last


def _bimodal_at(rows: List[str], v: int, margin: int,
                group_min_frac: float, max_span_frac: float) -> bool:
    """Do the copies split into a left group and a right group across column `v`?

    With `v` the valley midpoint (kept-span coordinate) and `margin` the guard band:
      left_only  : last_col  <  v + margin  AND  first_col < v
      right_only : first_col >  v - margin  AND  last_col  > v
      spanning   : first_col <= v - margin  AND  last_col >= v + margin
    Bimodal requires each group to clear `group_min_frac` and the spanning fraction
    to stay at/below `max_span_frac`. A valley with most copies spanning across it is
    an interior dip (AT-rich / indel), not a chimera."""
    n = len(rows)
    if n == 0:
        return False
    left_only = right_only = spanning = 0
    for r in rows:
        first, last = _span_bounds(r)
        if first < 0:
            continue
        if last < v + margin and first < v:
            left_only += 1
        if first > v - margin and last > v:
            right_only += 1
        if first <= v - margin and last >= v + margin:
            spanning += 1
    return (left_only / n >= group_min_frac
            and right_only / n >= group_min_frac
            and spanning / n <= max_span_frac)


# ═══════════════════════════════════════════════════════════════════════
# Recursive AND-gate chimera splitter
# ═══════════════════════════════════════════════════════════════════════

def detect_chimera(rec: Dict, config, depth: int = 0) -> List[Dict]:
    """MSA-based, BLAST-free chimera splitter — recursive two-signal AND-gate.

    Returns a list of records (one element = no split; >1 = split fragments). Every
    returned record has its transient '_aln' stripped. Conservative by construction:
    on any missing prerequisite, depth limit, or single-signal evidence, the family
    is returned intact — a real family is never fractured on weak evidence.

    Prerequisites and short-circuits (each -> keep the family intact):
      * no '_aln' (original-seed fallback families have no alignment) -> cannot split
      * depth >= chimera_max_depth
      * n < chimera_min_hits  or  span width < 2 * min_chimera_fragment
      * no qualifying interior occupancy valley (Signal 1)
      * valley present but copies do not split into two span groups (Signal 2)
    """
    # No alignment handle -> cannot reason about geometry; never force a split.
    if '_aln' not in rec:
        return [strip_aln(rec)]
    if depth >= config.chimera_max_depth:
        return [strip_aln(rec)]

    rows_full, metas, left, right = rec['_aln']
    # Work on the kept span only (column 0 here == full-coordinate column `left`).
    span = [r[left:right + 1] for r in rows_full]
    w = right - left + 1
    n = len(span)
    if n < config.chimera_min_hits or w < 2 * config.min_chimera_fragment:
        return [strip_aln(rec)]

    margin = config.chimera_min_valley_cols
    occ = occupancy_profile(span)

    # --- Signal 1: interior occupancy valley -------------------------------
    run_len, run_start = _longest_interior_valley(
        occ, config.chimera_occupancy_ratio, margin)
    if run_len < config.chimera_min_valley_cols:
        return [strip_aln(rec)]
    v = run_start + run_len // 2          # valley midpoint, kept-span coordinate

    # --- Signal 2: copy-span bimodality at the valley ----------------------
    if not _bimodal_at(span, v, margin,
                       config.chimera_span_group_min_frac,
                       config.chimera_max_span_frac):
        # Valley without two disjoint span groups = AT-rich / indel dip, not a
        # chimera. AND-gate precision: do not split.
        return [strip_aln(rec)]

    # --- Optional Signal 3: BED independence corroboration (default off) ----
    # enable_bed_independence_chimera would verify each half occurs as an
    # independent genomic insertion (inverse of Phase 0 co-occurrence). It only
    # *raises* confidence and is never required for a split; left unimplemented by
    # design (default False) so the common path stays blastn/BED-free.

    # --- Split: map valley column v to a consensus-base coordinate ----------
    # split_col is the full-coordinate boundary: frag1 = columns [left, split_col),
    # frag2 = columns [split_col, right]. split_base is the number of consensus bases
    # build_consensus emits over frag1's columns — computed by CALLING build_consensus
    # over [left, split_col - 1] so the coordinate is identical to the parent
    # consensus by construction (no re-derived emit logic to drift from M2).
    split_col = left + v
    seq = rec['sequence']
    split_base = len(build_consensus(rows_full, left, split_col - 1, metas, config))

    frag_rows_l = [r[left:split_col] for r in rows_full]
    frag_rows_r = [r[split_col:right + 1] for r in rows_full]
    frag1 = make_fragment(rec, seq[:split_base], frag_rows_l, '_chimfrag1')
    frag2 = make_fragment(rec, seq[split_base:], frag_rows_r, '_chimfrag2')

    out: List[Dict] = []
    for frag in (frag1, frag2):
        if len(frag['sequence']) < config.min_chimera_fragment:
            continue                      # discard a sub-min fragment, never emit it
        out.extend(detect_chimera(frag, config, depth + 1))
    if not out:
        # Both fragments fell below the floor (degenerate split) -> keep the family.
        return [strip_aln(rec)]
    logger.debug("detect_chimera: split %s at consensus base %d (valley col %d, "
                 "depth %d) -> %d fragment(s)",
                 rec.get('id', '?'), split_base, split_col, depth, len(out))
    return out
