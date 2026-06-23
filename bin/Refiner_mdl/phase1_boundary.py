"""
Phase 1 (rewrite) — Module M2b: extend-align-trim boundary refinement.

This is where a padded per-family MSA becomes a single refined consensus. The
defining invariant, preserved verbatim from the old code, is **never regress below
the mdl-repeat seed**: every degenerate or low-confidence path returns the original
sequence with consensus_source='original'. Refinement may only *replace* a family's
consensus when it can demonstrably improve it.

Two behaviors fall out of one mechanism (occupancy-trimming the padded MSA):
  * EXTENSION — pad columns beyond the original BED edge that stay solid across
    copies are kept, so a 5'-truncated terminus is recovered.
  * TRIM       — per-copy-unique flank collapses below the occupancy floor and is
    cut. Only termini move inward; interior low-occupancy columns are NEVER trimmed
    (guards AT-rich / indel-rich interiors).

See REFINE_IMPLEMENTATION_PLAN.md §1.3, §4, §5, §12.M2.
"""

import logging
from statistics import mean
from typing import Dict, List, Tuple

from phase1_align import (
    Alignment,
    build_msa,
    column_majority,
    occupancy_profile,
    select_regime,
)

logger = logging.getLogger(__name__)

_ACGT = frozenset('ACGT')


# ═══════════════════════════════════════════════════════════════════════
# Two-pointer occupancy trim
# ═══════════════════════════════════════════════════════════════════════

def trim_termini(rows: Alignment, config) -> Tuple[int, int]:
    """Two-pointer occupancy trim. Returns (left, right) inclusive column bounds.

    A column is `solid` if its present-count clears a floor that scales with family
    size — `max(trim_min_copies, trim_occupancy_floor * n)` — AND its majority base
    holds at least half of the present bases. The absolute floor of `trim_min_copies`
    stops a single chance-similar flank from masquerading as an extension; the
    fractional floor scales with family size (CIAlign / trimAl style).

    ONLY the two termini move inward to the first solid column. Interior
    low-occupancy columns between left..right are never trimmed. Shared pad columns
    beyond the original BED edge that stay solid produce extension; per-copy-unique
    flank collapses below the floor and is trimmed — both from this one scan."""
    n = len(rows)
    if n == 0:
        return 0, -1
    aln_len = len(rows[0])
    floor = max(config.trim_min_copies, config.trim_occupancy_floor * n)

    def solid(col: int) -> bool:
        base, cnt, npres = column_majority(rows, col)
        return npres >= floor and base is not None and cnt / max(npres, 1) >= 0.5

    left = 0
    while left < aln_len and not solid(left):
        left += 1
    right = aln_len - 1
    while right >= left and not solid(right):
        right -= 1
    return left, right


# ═══════════════════════════════════════════════════════════════════════
# Occupancy-aware consensus over the kept span
# ═══════════════════════════════════════════════════════════════════════

def build_consensus(rows: Alignment, left: int, right: int,
                    metas: List, config) -> str:
    """Occupancy-aware majority consensus over columns [left, right] (inclusive).

    Reuses the proven emit rules of the old build_majority_consensus, restricted to
    the kept span and operating on List[str]:
      * pure-insertion columns (no definite base) are dropped;
      * near-empty columns (< max(2, 0.1*n) present) are dropped as insertions;
      * otherwise emit the majority base if it holds >= consensus_min_occupancy of
        the present bases, else 'N'.
    Only ACGTN is emitted — no IUPAC ambiguity codes (downstream RepeatClassifier /
    hmmbuild prefer ACGTN). Optional divergence weighting (config-gated, default off)
    lets cleaner (low-divergence) copies dominate the vote."""
    if not rows or right < left:
        return ""
    n = len(rows)
    weighted = getattr(config, 'consensus_divergence_weighted', False)
    weights = None
    if weighted and metas:
        weights = [max(0.0, 1.0 - m.divergence) for m in metas]

    min_real = max(2, int(n * 0.1))
    out: List[str] = []
    for col in range(left, right + 1):
        counts: Dict[str, float] = {}
        total = 0.0
        present = 0
        for i, r in enumerate(rows):
            b = r[col]
            if b in _ACGT:
                w = weights[i] if weights is not None else 1.0
                counts[b] = counts.get(b, 0.0) + w
                total += w
                present += 1
        if present == 0:
            continue                    # pure insertion column -> drop
        if present < min_real:
            continue                    # near-empty insertion -> drop
        base = max(counts, key=lambda k: (counts[k], k))
        if counts[base] / total >= config.consensus_min_occupancy:
            out.append(base)
        else:
            out.append('N')
    return ''.join(out)


def alignment_qc_ok(rows: Alignment, config) -> bool:
    """List[str] equivalent of the old check_alignment_quality.

    Rejects an alignment whose fraction of low-occupancy columns (occupancy <
    consensus_min_occupancy) exceeds consensus_occupancy_fail_ratio."""
    if not rows:
        return False
    occ = occupancy_profile(rows)
    aln_len = len(occ)
    if aln_len == 0:
        return False
    low = int((occ < config.consensus_min_occupancy).sum())
    return (low / aln_len) <= config.consensus_occupancy_fail_ratio


# ═══════════════════════════════════════════════════════════════════════
# Per-family refinement (never regresses below the mdl-repeat seed)
# ═══════════════════════════════════════════════════════════════════════

def refine_family(rec: Dict, copies: List, config) -> Dict:
    """Produce the refined record for ONE family. Never worse than mdl-repeat.

    Fallback chain (any rung -> keep original seed, consensus_source='original'):
      * fewer than min_copies_for_msa copies          -> original seed
      * build_msa returned None (MAFFT failed / abPOA deferred) -> original seed
      * alignment_qc_ok is False                       -> original seed
      * trimmed span < min_chimera_fragment            -> original seed
      * consensus shorter than min_chimera_fragment    -> original seed
    On success: write the refined sequence, actual_length, consensus_source=<regime>,
    and a transient '_aln' handle (rows, metas, left, right) for the later chimera
    step. '_aln' MUST be stripped before the record leaves Phase 1."""
    out = dict(rec)
    out['sequence'] = rec['sequence']
    out['consensus_source'] = 'original'

    if len(copies) < config.min_copies_for_msa:
        return out

    mean_div = mean(c.divergence for c in copies) if copies else 0.0
    avg_len = int(mean(len(c.sequence) for c in copies)) if copies else 0
    regime = select_regime(len(copies), mean_div, config, avg_len=avg_len)
    fam = build_msa(copies, regime, config)
    if fam is None:
        return out
    rows, metas = fam
    if not alignment_qc_ok(rows, config):
        return out

    left, right = trim_termini(rows, config)
    if right - left + 1 < config.min_chimera_fragment:
        return out

    cons = build_consensus(rows, left, right, metas, config)
    if len(cons) < config.min_chimera_fragment:
        return out

    out['sequence'] = cons
    out['actual_length'] = len(cons)
    out['consensus_source'] = regime
    out['_aln'] = (rows, metas, left, right)   # transient; stripped before P1 exit
    return out
