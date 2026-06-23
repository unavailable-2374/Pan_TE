"""
Phase 1 (rewrite) — Module M1: coordinate math + BED-seeded copy extraction.

This is the geometric foundation of the Phase 1 rewrite (whole-genome BLAST
recruitment → BED-seeded extend–align–trim). Every downstream module trusts the
coordinates and strand orientation produced here, so an off-by-one or a missed
reverse-complement would silently corrupt the entire TE library. The pure
coordinate functions are unit-tested in both directions (tests/test_m1_extract.py).

System constraints honored exactly (verified, not assumed):
  * samtools 1.4.1 has NO `-r` region-file flag → regions are passed positionally,
    batched at config.extract_batch_size.
  * BLAST `-parse_seqids` subject ids carry an NCBI `ref|...|` wrapper → stripped by
    strip_seqid_prefix (a no-op on the clean `chr1..chrN` BED path, load-bearing only
    on the blastn fallback in a later milestone).

See REFINE_IMPLEMENTATION_PLAN.md §1.1, §2.2, §3, §11, §12.M1.
"""

import logging
import os
import random
import subprocess
from collections import OrderedDict, deque
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# Data structures (REFINE_IMPLEMENTATION_PLAN.md §2.2)
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class Instance:
    """One genomic copy locus of a family (from mdl_repeat.instances.bed)."""
    chrom: str
    start: int          # 0-based half-open BED start
    end: int            # 0-based half-open BED end
    strand: str         # '+' / '-'
    divergence: float   # 1 - score/1000


@dataclass
class Copy:
    """An extracted, strand-co-oriented genomic copy ready for MSA."""
    id: str
    sequence: str
    divergence: float
    strand: str
    # Offset of the instance (TE core, pre-pad) within `sequence`, strand-oriented.
    # Lets per-copy structural trim anchor a TSD search at the 5'/3' insertion junctions.
    te_start: int = 0
    te_end: int = 0


# CopyMeta's canonical definition lives in phase1_align (the module that produces MSA
# rows it parallels), per plan §1.1. M1 defined it locally with a TODO(M2); that TODO
# is now retired — we re-export the align-module definition for any M1-era importer.
from phase1_align import CopyMeta  # noqa: E402,F401


# ═══════════════════════════════════════════════════════════════════════
# Coordinate math (pure functions — unit-tested in both directions)
# ═══════════════════════════════════════════════════════════════════════

def bed_to_samtools_region(chrom: str, b_start: int, b_end: int,
                           pad: int, chrom_len: int) -> Tuple[str, int, int]:
    """Convert a BED interval (+pad) to a 1-based-inclusive samtools region.

    BED is 0-based half-open [b_start, b_end). samtools faidx 'chrom:X-Y' is
    1-based inclusive. Core (no pad) -> chrom:(b_start+1)-(b_end). With pad and
    clamping to [1, chrom_len]:

        s1 = max(1, b_start + 1 - pad)
        e1 = min(chrom_len, b_end + pad)

    Returns (region_str, s1, e1) where region_str == f"{chrom}:{s1}-{e1}".
    Raises ValueError if e1 < s1 (corrupt coords or locus past chromosome end);
    the caller drops + counts that instance, never fabricating a copy.
    """
    s1 = max(1, b_start + 1 - pad)
    e1 = min(chrom_len, b_end + pad)
    if e1 < s1:
        raise ValueError(
            f"empty/corrupt region for {chrom}: b_start={b_start} b_end={b_end} "
            f"pad={pad} chrom_len={chrom_len} -> s1={s1} e1={e1}")
    return f"{chrom}:{s1}-{e1}", s1, e1


def strip_seqid_prefix(sid: str) -> str:
    """Strip an NCBI-style 'ref|NC_x.y|' wrapper -> 'NC_x.y'.

    No-op for clean chr names (the BED path). Retained for the blastn fallback whose
    subject ids inherit makeblastdb -parse_seqids pipes (MEMORY: BLAST DB seqid
    format mismatch with samtools). Identical logic to the original _extract_hits
    prefix handling: first non-empty token after the first '|', else the whole id.
    """
    if '|' not in sid:
        return sid
    parts = sid.split('|')
    return next((p for p in parts[1:] if p), parts[0])


def load_fai_lengths(fai_path: str) -> Dict[str, int]:
    """Parse genome.fa.fai -> {chrom: length} (col 0 = name, col 1 = length).

    Loaded once per worker process (one row per chromosome, small)."""
    lengths: Dict[str, int] = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            lengths[parts[0]] = int(parts[1])
    return lengths


# Full IUPAC complement table; everything outside it (and any unknown char) -> 'N'.
_IUPAC_COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
    'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N',
}


def reverse_complement(seq: str) -> str:
    """ACGTN-safe reverse complement with full IUPAC table; unknown chars -> 'N'."""
    comp = (_IUPAC_COMPLEMENT.get(base, 'N') for base in reversed(seq.upper()))
    return ''.join(comp)


def compute_pad(family_len: int, config) -> int:
    """pad = int(clamp(config.pad_fraction * family_len, config.pad_min,
                       config.pad_cap)). Default 0.2*L clamped to [50, 500] bp.

    Extension is intrinsically capped at `pad` per side (only pad bp of flank is
    ever extracted), so no separate extension cap is needed."""
    raw = config.pad_fraction * family_len
    return int(max(config.pad_min, min(config.pad_cap, raw)))


def stratified_subsample(instances: List[Instance], cap: int,
                         seed: int) -> List[Instance]:
    """Down-sample to <= cap, stratified by divergence then length.

    Uniform random over-represents the youngest/most-numerous subfamily; stratifying
    by divergence keeps the consensus representative of the full age spectrum
    (strategy §4.1). Deterministic for a fixed seed (G10 posture)."""
    if len(instances) <= cap:
        return list(instances)

    rng = random.Random(seed)
    ordered = sorted(instances, key=lambda i: i.divergence)
    n = len(ordered)
    n_strata = min(10, max(1, n // cap + 1))

    # Contiguous equal-count divergence bins.
    strata: List[List[Instance]] = []
    base, rem = divmod(n, n_strata)
    idx = 0
    for s in range(n_strata):
        size = base + (1 if s < rem else 0)
        chunk = ordered[idx:idx + size]
        idx += size
        # Within a bin: length descending, ties broken by a deterministic shuffle.
        rng.shuffle(chunk)
        chunk.sort(key=lambda i: i.end - i.start, reverse=True)
        strata.append(chunk)

    # Round-robin one per bin until cap reached.
    selection: List[Instance] = []
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


# ═══════════════════════════════════════════════════════════════════════
# samtools faidx batch extraction
# ═══════════════════════════════════════════════════════════════════════

def _ensure_fai(genome_file: str, samtools_exe: str) -> None:
    """Build genome_file + '.fai' once if missing; raise loudly on failure.

    We cannot proceed without random access — never silently skip extraction."""
    fai_path = genome_file + '.fai'
    if os.path.exists(fai_path):
        return
    result = subprocess.run([samtools_exe, 'faidx', genome_file],
                            capture_output=True, text=True, timeout=600)
    if result.returncode != 0 or not os.path.exists(fai_path):
        raise RuntimeError(
            f"samtools faidx failed to index {genome_file}: {result.stderr[:300]}")


def _parse_faidx_fasta(text: str) -> Dict[str, str]:
    """Parse samtools faidx multi-region FASTA output -> {header: uppercased seq}.

    The header is the region string exactly as requested, e.g. '>testA:101-300'."""
    seqs: Dict[str, List[str]] = OrderedDict()
    cur: Optional[str] = None
    for line in text.split('\n'):
        if not line:
            continue
        if line[0] == '>':
            cur = line[1:].split()[0]
            seqs.setdefault(cur, [])
        elif cur is not None:
            seqs[cur].append(line.strip())
    return {h: ''.join(parts).upper() for h, parts in seqs.items()}


def extract_padded_copies(instances: List[Instance], genome_file: str,
                          fai_lengths: Dict[str, int], pad_fn: Callable[[int], int],
                          config) -> List[Copy]:
    """Batch-extract co-oriented, padded genomic copies for one family.

    Algorithm (plan §1.1, steps 1–7):
      1. For each instance compute a region; skip + count any that raise on corrupt
         or out-of-range coords — never fabricate a copy.
      2. Build region_to_meta : region_str -> deque[Instance] (a region string may
         repeat if two instances share coords).
      3. Sort unique regions by (chrom, start) so faidx seeks are monotonic.
      4. Ensure genome_file + '.fai' exists; fail loudly if indexing errors.
      5. Batch positional samtools faidx calls (config.extract_batch_size).
      6. Parse FASTA output; for each region emit one Copy per Instance sharing it
         (draining its deque deterministically), reverse-complementing '-' instances.
      7. Return List[Copy].

    pad_fn is applied per instance using the instance's own length as the element-length
    proxy, so a short truncated copy and a full-length copy get appropriately scaled
    flanks. Resident memory is O(sum of extracted bp for this family) only.
    """
    region_to_meta: "OrderedDict[str, deque]" = OrderedDict()
    region_sortkey: Dict[str, Tuple[str, int]] = {}
    skipped_bad_coords = 0
    skipped_no_chrom = 0

    for inst in instances:
        chrom_len = fai_lengths.get(inst.chrom)
        if chrom_len is None:
            skipped_no_chrom += 1
            continue
        pad = pad_fn(inst.end - inst.start)
        try:
            region_str, s1, _e1 = bed_to_samtools_region(
                inst.chrom, inst.start, inst.end, pad, chrom_len)
        except ValueError:
            skipped_bad_coords += 1
            continue
        if region_str not in region_to_meta:
            region_to_meta[region_str] = deque()
            region_sortkey[region_str] = (inst.chrom, s1)
        region_to_meta[region_str].append(inst)

    if skipped_bad_coords or skipped_no_chrom:
        logger.debug("extract_padded_copies: skipped %d bad-coord + %d unknown-chrom "
                     "instances", skipped_bad_coords, skipped_no_chrom)

    if not region_to_meta:
        return []

    # Step 3: monotonic seek order.
    unique_regions = sorted(region_to_meta.keys(),
                            key=lambda r: region_sortkey[r])

    # Step 4: random access requires the index.
    _ensure_fai(genome_file, config.samtools_exe)

    # Step 5: positional batches (samtools 1.4.1 has no -r flag).
    seq_by_region: Dict[str, str] = {}
    batch_size = config.extract_batch_size
    for batch_start in range(0, len(unique_regions), batch_size):
        batch = unique_regions[batch_start:batch_start + batch_size]
        cmd = [config.samtools_exe, 'faidx', genome_file] + batch
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        except Exception as e:  # noqa: BLE001 - report + continue, never fake a copy
            logger.warning("samtools faidx extraction failed for a batch: %s", e)
            continue
        if result.returncode != 0:
            logger.warning("samtools faidx returned %d: %s",
                           result.returncode, result.stderr[:200])
            continue
        seq_by_region.update(_parse_faidx_fasta(result.stdout))

    # Step 6–7: emit copies in monotonic region order, draining each deque.
    copies: List[Copy] = []
    for region_str in unique_regions:
        seq = seq_by_region.get(region_str)
        if not seq:
            # Batch failed or empty output for this region; those copies are lost
            # but never fabricated. Family may fall back to its seed downstream.
            continue
        rc_seq = None
        s1 = region_sortkey[region_str][1]          # 1-based extract start
        ext_len = len(seq)
        for inst in region_to_meta[region_str]:
            if inst.strand == '-':
                if rc_seq is None:
                    rc_seq = reverse_complement(seq)
                oriented = rc_seq
            else:
                oriented = seq
            # Instance offset within the (strand-oriented) extract: extract covers
            # genome [s1-1 : e1) 0-based, so the instance [inst.start:inst.end) sits at
            # [inst.start-(s1-1) : inst.end-(s1-1)); reverse for '-' (RC flips positions).
            ts = inst.start - (s1 - 1)
            te = inst.end - (s1 - 1)
            if inst.strand == '-':
                ts, te = ext_len - te, ext_len - ts
            cid = f"{inst.chrom}:{inst.start}-{inst.end}({inst.strand})"
            copies.append(Copy(id=cid, sequence=oriented,
                               divergence=inst.divergence, strand=inst.strand,
                               te_start=max(0, ts), te_end=min(ext_len, te)))

    return copies
