"""
Phase 1 — per-copy structural boundary trim (TSD-anchored).

Occupancy trimming is blind to a host-gene flank that is SHARED across copies (it stays
solid → kept → the consensus over-extends into CDS). This trims each copy to its true TE
boundary BEFORE the MSA, using the Target Site Duplication (TSD) — the short direct repeat
a TE creates in the host on insertion, flanking the element on both sides:

    [ ...host... TSD ][ TE ][ TSD ...host... ]

Each padded copy carries the instance offset (te_start/te_end). We search the 5' and 3'
insertion junctions (±slack) for a TSD: a k-mer ending at a 5' boundary b5 that re-occurs
starting at a 3' boundary b3. The TE is [b5:b3]; everything outside (the padding / shared
flank) is cut. This both removes over-extension AND recovers boundaries mdl under-called.
Trimming a subset of copies lowers the flank columns' occupancy, so the downstream
occupancy trim then contracts the consensus too — partial per-copy success still helps.

Conservative: acts ONLY on positive TSD evidence (else the copy is left padded for the
occupancy trim to decide); low-complexity / homopolymer TSDs are rejected.
"""

import logging
from typing import List

logger = logging.getLogger(__name__)


def _low_complexity(s: str) -> bool:
    if not s:
        return True
    if 'N' in s:
        return True
    return len(set(s)) <= 1            # homopolymer -> chance TSD, reject


def _find_tsd(seq: str, ts: int, te: int, config):
    """Return (b5, b3) TE boundaries if a flanking TSD is found, else None."""
    L = len(seq)
    slack = getattr(config, 'tsd_slack', 15)
    kmin = getattr(config, 'tsd_min', 4)
    kmax = getattr(config, 'tsd_max', 20)
    min_te = getattr(config, 'struct_trim_min_te', 50)
    lo5, hi5 = max(kmax, ts - slack), min(L, ts + slack)
    lo3, hi3 = max(0, te - slack), min(L - kmin, te + slack)
    if hi5 <= lo5 or hi3 <= lo3:
        return None
    # Prefer a longer TSD and a boundary closer to the instance edge.
    best = None
    for k in range(kmax, kmin - 1, -1):
        # index 3'-window k-mers -> start positions
        idx = {}
        for p3 in range(lo3, hi3 + 1):
            if p3 + k > L:
                continue
            idx.setdefault(seq[p3:p3 + k], []).append(p3)
        if not idx:
            continue
        for b5 in range(lo5, hi5 + 1):
            tsd = seq[b5 - k:b5]
            if _low_complexity(tsd) or tsd not in idx:
                continue
            for b3 in idx[tsd]:
                if b3 - b5 < min_te:
                    continue
                # score: prefer boundaries near the instance edges + longer TSD
                score = k - 0.05 * (abs(b5 - ts) + abs(b3 - te))
                if best is None or score > best[0]:
                    best = (score, b5, b3)
        if best is not None:
            break                       # take the longest-k TSD found
    return (best[1], best[2]) if best else None


def trim_copies_structural(copies: List, config) -> List:
    """Per-copy TSD trim with FAMILY-CONSENSUS gating.

    A single copy's TSD search is chance-prone (a short k-mer matches random flanks), so
    we never trust one copy alone. Each copy votes its best TSD boundary OFFSET (d5 = b5 −
    te_start, d3 = b3 − te_end); a real family TSD makes many copies agree on the SAME
    offset, while chance matches scatter. Only when an offset is supported by enough copies
    do we trim EVERY copy to that consensus boundary (removing the shared flank uniformly).
    """
    if not getattr(config, 'enable_structural_trim', False):
        return copies
    min_cn = getattr(config, 'struct_trim_min_copies', 4)
    if len(copies) < min_cn:
        return copies

    from collections import Counter
    # Real TSDs differ in SEQUENCE per copy (duplicated host target site), so we vote the
    # boundary OFFSET (binned to tolerate small instance-edge jitter), not the TSD bases.
    # The 5' and 3' boundaries are INDEPENDENT, so vote each separately (a joint vote
    # fragments support across the d5×d3 grid and never reaches a consensus).
    binbp = max(1, getattr(config, 'struct_trim_offset_bin', 3))
    v5, v3 = Counter(), Counter()
    for c in copies:
        seq = c.sequence
        ts, te = getattr(c, 'te_start', 0), getattr(c, 'te_end', 0)
        if te <= ts or te > len(seq) or (ts <= 0 and te >= len(seq)):
            continue
        hit = _find_tsd(seq, ts, te, config)
        if hit:
            v5[int(round((hit[0] - ts) / binbp)) * binbp] += 1
            v3[int(round((hit[1] - te) / binbp)) * binbp] += 1
    if not v5 or not v3:
        return copies
    d5, c5 = v5.most_common(1)[0]
    d3, c3 = v3.most_common(1)[0]
    frac = getattr(config, 'struct_trim_min_frac', 0.3)
    need = max(min_cn, int(frac * len(copies)))
    if c5 < need or c3 < need:
        return copies                                   # no consensus TSD -> leave to occupancy

    # Apply the consensus boundary offset to every copy (uniform flank removal).
    min_te = getattr(config, 'struct_trim_min_te', 50)
    n_trim = trimmed_bp = 0
    for c in copies:
        ts, te = getattr(c, 'te_start', 0), getattr(c, 'te_end', 0)
        if te <= ts:
            continue
        b5 = max(0, ts + d5)
        b3 = min(len(c.sequence), te + d3)
        if b3 - b5 < min_te or (b5 <= 0 and b3 >= len(c.sequence)):
            continue
        removed = len(c.sequence) - (b3 - b5)
        if removed <= 0:
            continue
        c.sequence = c.sequence[b5:b3]
        c.te_start = 0
        c.te_end = len(c.sequence)
        n_trim += 1
        trimmed_bp += removed
    if n_trim:
        logger.debug("Per-copy structural trim: consensus TSD offset d5=%d(%d) d3=%d(%d) "
                     "of %d copies -> trimmed %d (%d bp)", d5, c5, d3, c3, len(copies),
                     n_trim, trimmed_bp)
    return copies
