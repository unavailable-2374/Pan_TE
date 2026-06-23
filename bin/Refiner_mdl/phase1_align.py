"""
Phase 1 (rewrite) — Module M2a: MSA construction + alignment-row utilities.

This module owns the boundary between the MSA *engine* (MAFFT for small/medium
families, abPOA for large ones) and every downstream routine (occupancy, trim,
consensus, chimera). MAFFT and abPOA emit different containers, so both are
normalized to a single representation at the `build_msa` boundary:

    Alignment = List[str]   # one equal-length uppercase string per copy, '-' = gap

Nothing downstream ever touches a Bio.Align object; that decoupling is the whole
point of this module.

Key behavioral change from the old phase1_consensus.py: copies arriving here are
ALREADY strand-co-oriented by phase1_extract (minus-strand instances were reverse-
complemented at extraction). MAFFT is therefore called with adjustdirection=False —
the old code relied on `--adjustdirection` to *guess* orientation, which could
mis-flip short or divergent copies. run_mafft (Refiner/utils/alignment_utils.py)
exposes `adjustdirection` as a real parameter, so we pass it through directly.

See REFINE_IMPLEMENTATION_PLAN.md §1.2, §4, §5, §12.M2.
"""

import logging
import os
import random
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# ── Reuse Refiner MSA driver (real run_mafft signature, verified) ─────────
# run_mafft(sequences: List[Dict], algorithm='localpair', maxiterate=1000,
#           adjustdirection=True, thread=None, config=None) -> MultipleSeqAlignment
# It accepts a real `adjustdirection` flag, so the "copies are pre-oriented ->
# disable MAFFT auto-orientation" semantics is honored by passing False.
_refiner_dir = os.path.join(os.path.dirname(__file__), '..', 'Refiner')
if _refiner_dir not in sys.path:
    sys.path.insert(0, _refiner_dir)
from utils.alignment_utils import run_mafft  # noqa: E402


# ═══════════════════════════════════════════════════════════════════════
# Shared types (canonical home — phase1_extract imports CopyMeta from here)
# ═══════════════════════════════════════════════════════════════════════

# An MSA as a list of equal-length uppercase strings, one per copy, '-' = gap.
Alignment = List[str]
# rows + parallel per-copy metadata (one CopyMeta per row, same order).
AlignedFamily = Tuple[Alignment, List["CopyMeta"]]


@dataclass
class CopyMeta:
    """Per-row metadata carried parallel to MSA rows (chimera span analysis).

    Canonical definition lives here (the align module) per plan §1.1; phase1_extract
    imports it so M1's TODO(M2) placeholder is retired."""
    id: str
    divergence: float
    strand: str


# ═══════════════════════════════════════════════════════════════════════
# Regime selection
# ═══════════════════════════════════════════════════════════════════════

def select_regime(n_copies: int, mean_div: float, config, avg_len: int = 0) -> str:
    """Return 'mafft_linsi' | 'mafft_auto' | 'abpoa'.

    L-INS-i is reserved for small, low-divergence families where its accuracy is
    affordable; medium families use MAFFT --auto; large families defer to abPOA
    (POA scales where progressive MSA does not).

        n_copies × avg_len > msa_abpoa_bp_threshold                       -> abpoa
        n <= small_family_threshold and mean_div <= small_family_max_div -> mafft_linsi
        n <= large_family_threshold                                      -> mafft_auto
        else                                                             -> abpoa

    Cost is driven by TOTAL MSA work (n_copies × sequence length), not copy count
    alone: MAFFT on ~100 sequences of 8.5 kb is minutes (measured: R=1567 refine >300 s),
    so any cluster whose n_copies × avg_len exceeds `msa_abpoa_bp_threshold` is routed to
    abPOA regardless of count — POA is far faster on many/long sequences. Short-sequence
    clusters keep MAFFT's accuracy.
    """
    if avg_len and n_copies * avg_len > getattr(config, 'msa_abpoa_bp_threshold', 150_000):
        return 'abpoa'
    if (n_copies <= config.small_family_threshold
            and mean_div <= config.small_family_max_div):
        return 'mafft_linsi'
    if n_copies <= config.large_family_threshold:
        return 'mafft_auto'
    return 'abpoa'


# ═══════════════════════════════════════════════════════════════════════
# MSA construction
# ═══════════════════════════════════════════════════════════════════════

def _subsample_copies(copies: List, cap: int, seed: int, max_total_bp: int = 0) -> List:
    """Deterministic down-sample of co-oriented copies to <= cap.

    Stratified by divergence (ascending) then sequence length (descending) so the
    consensus stays representative of the full age spectrum rather than being
    dominated by the youngest/most-numerous subfamily. Mirrors
    phase1_extract.stratified_subsample but operates on Copy objects (which carry
    `divergence` and a `sequence`, not BED start/end).

    MSA cost scales with n_copies × sequence length, so when `max_total_bp` is set the
    count cap is tightened to `max_total_bp / avg_len` for long-sequence clusters: a
    114×8.5 kb cluster otherwise takes ~80 s and (for these large divergent families)
    falls back to the original seed anyway — a representative ~30-copy subset gives the
    same consensus at a fraction of the cost."""
    if max_total_bp and copies:
        avg_len = max(1, sum(len(c.sequence) for c in copies) // len(copies))
        bp_cap = max(2, max_total_bp // avg_len)
        cap = min(cap, bp_cap)
    if len(copies) <= cap:
        return list(copies)
    rng = random.Random(seed)
    ordered = sorted(copies, key=lambda c: c.divergence)
    n = len(ordered)
    n_strata = min(10, max(1, n // cap + 1))
    strata: List[List] = []
    base, rem = divmod(n, n_strata)
    idx = 0
    for s in range(n_strata):
        size = base + (1 if s < rem else 0)
        chunk = ordered[idx:idx + size]
        idx += size
        rng.shuffle(chunk)
        chunk.sort(key=lambda c: len(c.sequence), reverse=True)
        strata.append(chunk)
    selection: List = []
    pos = 0
    while len(selection) < cap:
        progressed = False
        for chunk in strata:
            if pos < len(chunk):
                selection.append(chunk[pos])
                progressed = True
                if len(selection) >= cap:
                    break
        if not progressed:
            break
        pos += 1
    return selection


def _normalize_row(s: str) -> str:
    """Uppercase a row and normalize every gap glyph ('.', ' ') to '-'."""
    return s.upper().replace('.', '-').replace(' ', '-')


def _bio_alignment_to_rows(alignment, copies) -> Optional[AlignedFamily]:
    """Convert a Bio MultipleSeqAlignment to (rows, metas) in *input* copy order.

    MAFFT may reorder records, so we index aligned rows by sequence id and re-emit
    them in the order of `copies`. Duplicate ids (two instances at identical coords)
    are drained in encounter order via a per-id list cursor. Returns None if any
    input copy is missing from the alignment or rows are not equal length (never
    silently truncates)."""
    from collections import defaultdict, deque

    by_id: "defaultdict[str, deque]" = defaultdict(deque)
    for rec in alignment:
        by_id[rec.id].append(_normalize_row(str(rec.seq)))

    rows: List[str] = []
    metas: List[CopyMeta] = []
    for c in copies:
        bucket = by_id.get(c.id)
        if not bucket:
            logger.warning("build_msa: copy id %s missing from MAFFT output", c.id)
            return None
        rows.append(bucket.popleft())
        metas.append(CopyMeta(id=c.id, divergence=c.divergence, strand=c.strand))

    aln_len = len(rows[0])
    if any(len(r) != aln_len for r in rows):
        logger.warning("build_msa: unequal MAFFT row lengths (%s)",
                       sorted({len(r) for r in rows}))
        return None
    return rows, metas


def _parse_fasta_with_gaps(path: str) -> List[Tuple[str, str]]:
    """Tolerant FASTA-with-gaps reader for abPOA MSA output (file order preserved).

    abPOA's `-r 1` mode writes an aligned FASTA (it calls the format PIR, but the
    bytes are `>name` headers followed by gapped rows). A strict PIR parser would
    choke on it, so we read it as plain FASTA: each '>' starts a record (name = first
    whitespace-delimited token), subsequent lines are concatenated as the row. Gap
    glyphs are NOT normalized here — the caller does that via _normalize_row.
    Returns [(name, raw_seq), ...] in the order they appear in the file."""
    records: List[Tuple[str, List[str]]] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line[0] == '>':
                records.append((line[1:].split()[0] if line[1:].strip() else '',
                                []))
            elif records:
                records[-1][1].append(line.strip())
    return [(name, ''.join(parts)) for name, parts in records]


def _oversized_for_msa(work: List, config, engine: str) -> bool:
    """Return True when the (already-subsampled) work set is too big for an MSA engine.

    The danger case for BOTH engines is a small number of VERY LONG copies — a giant
    LTR/retrotransposon family (chr4 R=2: 65 copies up to ~400 kb with padding, ~3.6 Mbp
    total). abPOA's banded DP OOM-kills on it; MAFFT --auto does not OOM but grinds for
    its full wall-clock timeout (minutes per family) before falling back — the dominant
    cost on a full run. So we apply the SAME copy-set / total-bp guard to both paths:
    an oversized family skips MSA entirely and the caller keeps the original mdl-repeat
    seed up front (never-regress, never faked, no wasted compute). 0 disables a guard."""
    max_seq_bp = getattr(config, 'abpoa_max_seq_bp', 0) or 0
    max_total_bp = getattr(config, 'abpoa_max_total_bp', 0) or 0
    if max_seq_bp > 0:
        longest = max(len(c.sequence) for c in work)
        if longest > max_seq_bp:
            logger.warning("build_msa(%s): skipping — longest copy %d bp > "
                           "abpoa_max_seq_bp %d (giant element; keeping seed)",
                           engine, longest, max_seq_bp)
            return True
    if max_total_bp > 0:
        total_bp = sum(len(c.sequence) for c in work)
        if total_bp > max_total_bp:
            logger.warning("build_msa(%s): skipping — total %d bp over %d copies > "
                           "abpoa_max_total_bp %d (keeping seed)",
                           engine, total_bp, len(work), max_total_bp)
            return True
    return False


def _run_abpoa(copies: List, config) -> Optional[AlignedFamily]:
    """Large-family partial-order MSA via the abPOA binary (regime == 'abpoa').

    Copies are subsampled to msa_subsample_cap_abpoa (already strand-co-oriented).
    They are written with positional names (s0, s1, ...) so the parsed rows map back
    to input order by name regardless of any reordering. On any failure path
    (binary missing, non-zero exit, empty output, wrong record count, or unequal row
    lengths) we return None and the caller keeps the original mdl-repeat seed — a
    fabricated alignment is never produced."""
    work = _subsample_copies(copies, config.msa_subsample_cap_abpoa,
                             config.subsample_seed,
                             max_total_bp=getattr(config, 'msa_max_total_bp', 0))
    if len(work) < 2:
        return None

    # Pre-launch size guard (§N6): a giant-element work set OOM-kills abPOA before the
    # wall clock can fire, so we skip it up front and keep the seed (never-regress).
    if _oversized_for_msa(work, config, 'abpoa'):
        return None

    abpoa_exe = getattr(config, 'abpoa_exe', 'abpoa') or 'abpoa'
    if os.path.sep in abpoa_exe and not os.path.exists(abpoa_exe):
        abpoa_exe = 'abpoa'  # PATH fallback when the pinned path is absent

    tmp_dir = tempfile.mkdtemp(prefix='abpoa_msa_')
    in_path = os.path.join(tmp_dir, 'in.fa')
    out_path = os.path.join(tmp_dir, 'out.msa')
    name_to_copy: Dict[str, object] = {}
    try:
        with open(in_path, 'w') as fh:
            for i, c in enumerate(work):
                tag = f"s{i}"
                name_to_copy[tag] = c
                fh.write(f">{tag}\n{c.sequence}\n")

        cmd = [abpoa_exe, '-r', '1', '-o', out_path, in_path]
        wall_s = getattr(config, 'msa_wall_s', 300) or 300
        try:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    timeout=wall_s)
        except FileNotFoundError:
            logger.warning("build_msa(abpoa): binary not found (%s)", abpoa_exe)
            return None
        except subprocess.TimeoutExpired:
            logger.warning("build_msa(abpoa): wall-clock timeout (%ds) over %d copies "
                           "(%d bp) — keeping seed (never-regress)",
                           wall_s, len(work), sum(len(c.sequence) for c in work))
            return None
        except Exception as e:  # noqa: BLE001 — report + fall back, never fake an MSA
            logger.warning("build_msa(abpoa): execution failed (%s)", e)
            return None
        if result.returncode != 0:
            logger.warning("build_msa(abpoa): exit %d: %s",
                           result.returncode, (result.stderr or '')[:200])
            return None
        if not os.path.exists(out_path):
            logger.warning("build_msa(abpoa): no output file produced")
            return None

        parsed = _parse_fasta_with_gaps(out_path)
        if not parsed:
            logger.warning("build_msa(abpoa): empty MSA output")
            return None

        by_name = {name: seq for name, seq in parsed}
        rows: List[str] = []
        metas: List[CopyMeta] = []
        for i, c in enumerate(work):
            raw = by_name.get(f"s{i}")
            if raw is None:
                logger.warning("build_msa(abpoa): copy s%d missing from MSA output", i)
                return None
            rows.append(_normalize_row(raw))
            metas.append(CopyMeta(id=c.id, divergence=c.divergence, strand=c.strand))

        aln_len = len(rows[0])
        if aln_len == 0 or any(len(r) != aln_len for r in rows):
            logger.warning("build_msa(abpoa): unequal/empty row lengths (%s)",
                           sorted({len(r) for r in rows}))
            return None
        return rows, metas
    finally:
        import shutil
        shutil.rmtree(tmp_dir, ignore_errors=True)


def build_msa(copies: List, regime: str, config) -> Optional[AlignedFamily]:
    """Align co-oriented copies; return (rows, metas) or None on any failure.

    Copies are ALREADY strand-co-oriented, so MAFFT is invoked with
    adjustdirection=False (the key change from the old code). Small/medium families
    take the MAFFT path; large families (regime 'abpoa') take partial-order alignment
    via the abPOA binary. Every failure path returns None and the caller falls back
    to the original mdl-repeat seed — safe, never fabricated."""
    if regime == 'abpoa':
        return _run_abpoa(copies, config)
    if regime not in ('mafft_linsi', 'mafft_auto'):
        logger.warning("build_msa: unknown regime %r", regime)
        return None
    if len(copies) < 2:
        return None

    work = _subsample_copies(copies, config.msa_subsample_cap_mafft,
                             config.subsample_seed,
                             max_total_bp=getattr(config, 'msa_max_total_bp', 0))
    # Same giant-element guard as the abPOA path (§N6): a few very long copies make
    # MAFFT --auto grind for its full 600 s wall-clock before falling back. Skip up
    # front and keep the seed — never-regress, no wasted minutes per giant family.
    if _oversized_for_msa(work, config, regime):
        return None
    seqs = [{'id': c.id, 'sequence': c.sequence} for c in work]
    algorithm = 'localpair' if regime == 'mafft_linsi' else 'auto'

    alignment = run_mafft(seqs, algorithm=algorithm, adjustdirection=False,
                          thread=1, config=config)
    if alignment is None or len(alignment) == 0:
        return None
    return _bio_alignment_to_rows(alignment, work)


# ═══════════════════════════════════════════════════════════════════════
# Alignment-row utilities (operate on Alignment = List[str])
# ═══════════════════════════════════════════════════════════════════════

_ACGT = frozenset('ACGT')


def occupancy_profile(rows: Alignment) -> np.ndarray:
    """occ[col] = (# rows with a definite base in {A,C,G,T}) / n_rows.

    N, '-' and '.' all count as *absent*. Returned as a float ndarray of length
    aln_len. Empty input -> empty array."""
    n = len(rows)
    if n == 0:
        return np.zeros(0, dtype=float)
    aln_len = len(rows[0])
    occ = np.zeros(aln_len, dtype=float)
    for col in range(aln_len):
        present = sum(1 for r in rows if r[col] in _ACGT)
        occ[col] = present / n
    return occ


def column_majority(rows: Alignment, col: int) -> Tuple[Optional[str], int, int]:
    """Return (majority_base, majority_count, n_present) over ACGT at one column.

    n_present excludes N/'-'/'.'. With no definite base present, returns
    (None, 0, 0)."""
    counts: Dict[str, int] = {}
    for r in rows:
        b = r[col]
        if b in _ACGT:
            counts[b] = counts.get(b, 0) + 1
    if not counts:
        return None, 0, 0
    base = max(counts, key=lambda k: (counts[k], k))
    return base, counts[base], sum(counts.values())
