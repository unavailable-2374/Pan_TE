"""
Phase 0: Metadata parsing + hard filtering + adaptive tiering +
         fragment assembly + two-round deduplication.

Input:  mdl-repeat FASTA (>R=N length=L copies=C mdl=M)
Output: Tiered, deduplicated sequences ready for Phase 1
"""

import logging
import os
import re
import shutil
import subprocess
import tempfile
import numpy as np
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

logger = logging.getLogger(__name__)

# ── Reuse Refiner utilities ──────────────────────────────────────────
import sys
_refiner_dir = os.path.join(os.path.dirname(__file__), '..', 'Refiner')
_bin_dir = os.path.join(os.path.dirname(__file__), '..')
for _d in (_refiner_dir, _bin_dir):
    if _d not in sys.path:
        sys.path.insert(0, _d)
from utils.complexity_utils import (calculate_shannon_entropy, calculate_dust_score,
                                     calculate_lowcomplexity_fraction)
# Shared chain-stitching primitives (pure; safe to import here — see fragment_stitcher).
import fragment_stitcher


# ═══════════════════════════════════════════════════════════════════════
# Step 0.1  Metadata Parsing
# ═══════════════════════════════════════════════════════════════════════

# Required core fields written by mdl-repeat output.c since v0
_HEADER_RE = re.compile(
    r'>R=(\d+)\s+length=(\d+)\s+copies=(\d+)\s+mdl=([\d.]+)')

# Optional fields emitted since mdl-repeat v6.1 — match independently so
# headers that lack them (older mdl-repeat builds) still parse cleanly.
_DIV_RE = re.compile(r'\bdiv=([\d.]+)')
_TOPO_RE = re.compile(r'\btopo=(linear|cyclic|complex)')
# mdl-repeat's NATIVE per-family confidence (header: `accept=... tier=...`). These are
# mdl's own verdict on each family; consumed by the Phase-2 quality gate. `tier` here is
# mdl's core/warn/reject — distinct from Refiner's processing tier (T1/T2/T3 in rec['tier']),
# so they are stored under mdl_accept / mdl_tier to avoid collision.
_ACCEPT_RE = re.compile(r'\baccept=(\w+)')
_MDLTIER_RE = re.compile(r'\btier=(\w+)')


def parse_mdl_fasta(fasta_path: str) -> List[Dict]:
    """Parse mdl-repeat FASTA.  O(n) single pass.

    Returns list of dicts with keys:
        id, R, length, copies, mdl, sequence, divergence (opt), topology (opt)

    `divergence` and `topology` are populated when the FASTA was produced by
    mdl-repeat >= v6.1 (header carries `div=` / `topo=`).  Older headers
    leave these fields absent; downstream code must not assume them.
    """
    records = []
    current = None
    seq_parts = []

    def _flush():
        if current is not None:
            current['sequence'] = ''.join(seq_parts).upper()
            current['actual_length'] = len(current['sequence'])
            records.append(current)

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                _flush()
                seq_parts = []
                m = _HEADER_RE.match(line)
                if m:
                    r_id = int(m.group(1))
                    current = {
                        'id': f'mdl_R{r_id}',
                        'R': r_id,
                        'length': int(m.group(2)),
                        'copies': int(m.group(3)),
                        'mdl': float(m.group(4)),
                    }
                    div_m = _DIV_RE.search(line)
                    if div_m:
                        current['divergence'] = float(div_m.group(1))
                    topo_m = _TOPO_RE.search(line)
                    if topo_m:
                        current['topology'] = topo_m.group(1)
                    acc_m = _ACCEPT_RE.search(line)
                    if acc_m:
                        current['mdl_accept'] = acc_m.group(1)
                    mt_m = _MDLTIER_RE.search(line)
                    if mt_m:
                        current['mdl_tier'] = mt_m.group(1)
                else:
                    # Fallback: non-standard header
                    hdr = line[1:].split()[0]
                    current = {
                        'id': hdr,
                        'R': 0,
                        'length': 0,
                        'copies': 1,
                        'mdl': 1.0,
                    }
            else:
                seq_parts.append(line)
    _flush()

    n_with_topo = sum(1 for r in records if 'topology' in r)
    n_with_div = sum(1 for r in records if 'divergence' in r)
    n_with_tier = sum(1 for r in records if 'mdl_tier' in r)
    logger.info(
        f"Parsed {len(records)} sequences from {fasta_path} "
        f"(topology: {n_with_topo}, divergence: {n_with_div}, mdl_tier: {n_with_tier})"
    )
    return records


# ═══════════════════════════════════════════════════════════════════════
# Step 0.2  Hard Filtering
# ═══════════════════════════════════════════════════════════════════════

def hard_filter(records: List[Dict], config) -> Tuple[List[Dict], Dict]:
    """Apply hard filters.  Returns passing records."""
    passed = []
    stats = {'total': len(records), 'too_short': 0, 'single_copy': 0,
             'low_entropy': 0, 'high_n': 0, 'low_mdl': 0, 'ssr_dust': 0,
             'cyclic_topology': 0}

    for rec in records:
        seq = rec['sequence']
        slen = len(seq)

        # Cyclic topology — mdl-repeat flags this when discovered instances
        # form a closed circle, which is overwhelmingly plasmid / mitochondrial
        # / chloroplast contamination rather than a true TE family.  Drop
        # before any length / entropy check (saves DUST work).
        if rec.get('topology') == 'cyclic':
            stats['cyclic_topology'] += 1
            continue

        if slen < config.min_length:
            stats['too_short'] += 1
            continue
        if rec['copies'] <= 1:
            # mdl-repeat already filters mdl_score <= 0 at output.c — single-copy
            # families have negative MDL, so this branch should be ~0 in
            # practice.  A non-zero count signals an mdl-repeat MDL anomaly.
            stats['single_copy'] += 1
            continue
        if rec['mdl'] <= config.min_mdl_score:
            stats['low_mdl'] += 1
            continue

        # N content
        n_frac = seq.count('N') / slen if slen > 0 else 0
        if n_frac > config.max_n_percent:
            stats['high_n'] += 1
            continue

        # Shannon entropy
        ent = calculate_shannon_entropy(seq, k=1)
        thresh = (config.entropy_threshold_long if slen >= 100
                  else config.entropy_threshold_short)
        if ent < thresh:
            stats['low_entropy'] += 1
            continue

        # DUST score: primary SSR/simple-repeat filter
        # SSR (AT)n, (ATC)n etc. have DUST ~1.0; real TEs ~0.1-0.4
        dust = calculate_dust_score(seq)
        if dust > config.dust_score_threshold:
            stats['ssr_dust'] += 1
            continue

        # k=2 (dinucleotide) entropy: catches (CA)n-style low-complexity that k=1 passes
        ent2 = calculate_shannon_entropy(seq, k=2)
        thresh2 = (config.entropy2_threshold_long if slen >= 100
                   else config.entropy2_threshold_short)
        if ent2 < thresh2:
            stats['low_entropy2'] = stats.get('low_entropy2', 0) + 1
            continue

        # Low-complexity fraction: mostly-simple consensi that slip the average DUST gate
        lc_frac = calculate_lowcomplexity_fraction(seq, dust_cut=config.dust_score_threshold)
        if lc_frac > config.lowcomplexity_fraction_max:
            stats['high_lc_fraction'] = stats.get('high_lc_fraction', 0) + 1
            continue

        rec['entropy'] = ent
        rec['dust_score'] = dust
        rec['n_frac'] = n_frac
        passed.append(rec)

    stats['passed'] = len(passed)
    logger.info(f"Hard filter: {stats['total']} -> {stats['passed']} "
                f"(cyclic={stats['cyclic_topology']}, short={stats['too_short']}, "
                f"single={stats['single_copy']}, entropy={stats['low_entropy']}, "
                f"SSR/dust={stats['ssr_dust']}, N%={stats['high_n']}, "
                f"mdl={stats['low_mdl']})")
    return passed, stats


# ═══════════════════════════════════════════════════════════════════════
# Step 0.3  Adaptive Tiering
# ═══════════════════════════════════════════════════════════════════════

def assign_tiers(records: List[Dict], config) -> List[Dict]:
    """Compute mdl_per_copy and assign T1/T2/T3 tiers."""
    if not records:
        return records

    for rec in records:
        rec['mdl_per_copy'] = rec['mdl'] / max(rec['copies'], 1)

    mpc_values = np.array([r['mdl_per_copy'] for r in records])
    n = len(records)

    if n < config.small_input_threshold:
        # Small input: absolute thresholds
        p75 = config.tier_absolute_high
        p25 = config.tier_absolute_low
        logger.info(f"Small input ({n} seq): using absolute tier thresholds "
                    f"T1>{p75}, T3<{p25}")
    else:
        p25 = float(np.percentile(mpc_values, 25))
        p75 = float(np.percentile(mpc_values, 75))
        logger.info(f"Adaptive tier thresholds: P25={p25:.2f}, P75={p75:.2f}")

    # Divergence ceiling for T1 promotion: a family with mdl_per_copy in the
    # top quartile but high mean instance divergence has noisy copies, so its
    # consensus is less trustworthy.  Cap T1 at div <= 0.20 when the field is
    # available; older mdl-repeat headers without div= are unaffected.
    DIV_T1_CEILING = 0.20

    t1_count = t2_count = t3_count = 0
    t1_demoted_high_div = 0
    for rec in records:
        slen = len(rec['sequence'])
        mpc = rec['mdl_per_copy']
        div = rec.get('divergence')

        is_t1_candidate = (mpc >= p75 or (rec['copies'] >= 50 and slen >= 80))
        if is_t1_candidate and div is not None and div > DIV_T1_CEILING:
            # Demote: high divergence means consensus polish must work harder;
            # treat as T2 so it gets full BLASTN+MAFFT instead of T1 fast-track.
            is_t1_candidate = False
            t1_demoted_high_div += 1

        if is_t1_candidate:
            rec['tier'] = 'T1'
            t1_count += 1
        elif mpc < p25:
            rec['tier'] = 'T3'
            t3_count += 1
        else:
            rec['tier'] = 'T2'
            t2_count += 1

    logger.info(f"Tiering: T1={t1_count}, T2={t2_count}, T3={t3_count}")
    if t1_demoted_high_div:
        logger.info(f"  T1→T2 demotions due to div>{DIV_T1_CEILING}: "
                    f"{t1_demoted_high_div}")

    # T2 downgrade: if too many T2 sequences, downgrade short ones to T1
    if t2_count > config.t2_downgrade_count:
        downgraded = 0
        for rec in records:
            if rec['tier'] == 'T2' and len(rec['sequence']) < config.t2_downgrade_length:
                rec['tier'] = 'T1'
                downgraded += 1
        if downgraded:
            logger.info(f"T2 downgrade: {downgraded} short sequences (<{config.t2_downgrade_length}bp) "
                        f"moved to T1 fast track")

    return records


# ═══════════════════════════════════════════════════════════════════════
# Step 0.4  Fragment Assembly — Directed Greedy Chain Algorithm
#
# Why not Union-Find?  Union-Find produces connected components via
# transitive closure: A–B, B–C, C–D → {A,B,C,D}.  In TE-rich genomes
# different families nest inside each other, creating long adjacency
# chains that merge thousands of unrelated fragments into one scaffold.
#
# Directed chains model TEs correctly: a TE is a *linear* sequence of
# domains.  Each fragment has at most 1 left neighbor and 1 right
# neighbor.  This naturally caps scaffold size without arbitrary limits.
#
# Multi-dimensional edge filtering:
#   1. Co-occurrence ratio   — must exceed min fraction of the rarer family
#   2. Direction consistency — same order in ≥80% of co-occurrences
#   3. Copy-number ratio     — true fragments have similar copy counts
#   4. Overlap sanity        — excessive overlap = redundancy, not fragments
# ═══════════════════════════════════════════════════════════════════════


def _bed_genome_coord_check(bed_path: str, genome_file: str) -> None:
    """Fail loud if the instances BED chrom names don't intersect the genome's sequence
    names. A mismatch (e.g. a full-genome BED chr1..chrN reused against a SAMPLE whose
    sequences are windowed names like chr6_w6713774_7713773) makes every Phase-1 copy
    extraction silently miss — Phase 1 becomes a no-op and raw mdl fragments are shipped
    as the library. Stopping here is strictly better than emitting a fragmented library."""
    if not bed_path or not os.path.isfile(bed_path) or not genome_file:
        return
    names = set()
    fai = genome_file + '.fai'
    if os.path.isfile(fai):
        with open(fai) as fh:
            for line in fh:
                if line.strip():
                    names.add(line.split('\t', 1)[0])
    elif os.path.isfile(genome_file):
        with open(genome_file) as fh:
            for line in fh:
                if line.startswith('>'):
                    names.add(line[1:].split()[0])
    if not names:
        return  # cannot determine genome names; don't block
    bed_chroms = set()
    with open(bed_path) as fh:
        for i, line in enumerate(fh):
            if line.strip():
                bed_chroms.add(line.split('\t', 1)[0])
            if i >= 200000:  # a prefix is enough to detect a total coordinate-space mismatch
                break
    if bed_chroms and names.isdisjoint(bed_chroms):
        raise RuntimeError(
            "BED/genome coordinate-space mismatch: instances-BED chroms "
            f"{sorted(bed_chroms)[:3]}… share NO sequence name with the refine genome "
            f"{sorted(names)[:3]}… ({genome_file}). Phase 1 copy extraction would recruit "
            "nothing (silent no-op) and ship raw mdl-repeat fragments. mdl-repeat must run "
            "on the SAME genome Refiner extracts from; in sampled mode the cached "
            "full-genome BED must be regenerated on the sample (delete the stale "
            "mdl_repeat_complete checkpoint)."
        )


def _parse_bed_instances(bed_path: str) -> Dict[str, List[Tuple]]:
    """Parse mdl-repeat BED instances file.

    Expected format: chr  start  end  R=N  score  strand
    Returns: {family_id: [(chr, start, end, strand), ...]}
    """
    instances = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            fam_id = parts[3]  # e.g. "R=42"
            strand = parts[5] if len(parts) > 5 else '+'
            instances[fam_id].append((chrom, start, end, strand))
    return instances


def _score_and_filter_edges(directed_pairs: Dict, copies: Dict,
                            config, check_direction: bool = True) -> List[Dict]:
    """Score + filter adjacency edges. Delegates to the shared implementation in
    fragment_stitcher (single source of truth, shared with the te-looker track)."""
    return fragment_stitcher.score_filter_edges(
        directed_pairs, copies, config, check_direction=check_direction)


def _build_chains_and_scaffolds(scored_edges: List[Dict], id_to_rec: Dict,
                                records: List[Dict], config) -> List[Dict]:
    """Build directed chains from scored edges and construct scaffold records.

    Each node can have at most one left neighbor and one right neighbor,
    preventing the transitive-closure over-assembly of Union-Find.
    """
    if not scored_edges:
        logger.info("Fragment assembly: no edges passed multi-dimensional filtering")
        return records

    # Greedy one-left/one-right directed-DAG chain extraction (shared implementation).
    chains = fragment_stitcher.extract_chains(scored_edges)

    # Adaptive max scaffold length: max(P95 × 2, stitch_max_len), capped at 50kb.
    # The floor is config.stitch_max_len (default 25 kb) so a full-length LTR
    # retrotransposon fragmented across families can be reassembled here rather than
    # being truncated — mdl-repeat does fragment such elements (they are NOT always
    # caught intact by Look4LTRs). The 50 kb hard cap still guards against runaway.
    stitch_max = getattr(config, 'stitch_max_len', 25000)
    max_len = config.max_scaffold_length
    if max_len <= 0:
        input_lens = sorted(len(r['sequence']) for r in records)
        p95_idx = max(0, int(len(input_lens) * 0.95) - 1)
        p95 = input_lens[p95_idx] if input_lens else 5000
        max_len = min(max(p95 * 2, stitch_max), 50000)
        logger.info(f"Adaptive max scaffold length: {max_len}bp "
                    f"(P95 input = {p95}bp)")

    # Build scaffolds from validated chains
    new_records = []
    merged_ids = set()
    scaffold_id = 0
    skipped_long = 0

    for chain in chains:
        member_recs = [id_to_rec[m] for m in chain if m in id_to_rec]
        if len(member_recs) <= 1:
            continue

        total_len = (sum(len(r['sequence']) for r in member_recs)
                     + 10 * (len(member_recs) - 1))
        if total_len > max_len:
            skipped_long += 1
            continue

        scaffold_id += 1
        # Concatenate in chain order (preserves genomic arrangement)
        scaffold_seq = ('N' * 10).join(r['sequence'] for r in member_recs)

        scaffold_rec = {
            'id': f'scaffold_{scaffold_id}',
            'R': member_recs[0]['R'],
            'length': len(scaffold_seq),
            'copies': min(r['copies'] for r in member_recs),
            'mdl': min(r['mdl'] for r in member_recs),
            'sequence': scaffold_seq,
            'actual_length': len(scaffold_seq),
            'mdl_per_copy': min(r.get('mdl_per_copy', 0) for r in member_recs),
            'tier': max((r['tier'] for r in member_recs),
                        key=lambda t: {'T1': 3, 'T2': 2, 'T3': 1}[t]),
            'is_scaffold': True,
            'scaffold_members': [r['id'] for r in member_recs],
        }
        new_records.append(scaffold_rec)
        merged_ids.update(m for m in chain if m in id_to_rec)

    kept = [r for r in records if r['id'] not in merged_ids]
    kept.extend(new_records)
    logger.info(
        f"Directed chain assembly: {len(chains)} chains, "
        f"{len(merged_ids)} fragments -> {len(new_records)} scaffolds"
        f"{f', {skipped_long} chains exceeded {max_len}bp' if skipped_long else ''}; "
        f"{len(kept)} total sequences")
    return kept


def assemble_fragments_bed(records: List[Dict], bed_path: str,
                           config) -> List[Dict]:
    """Fragment assembly using BED instance co-occurrence.

    Uses directed greedy chain algorithm: each family can have at most
    one left neighbor and one right neighbor, preventing transitive
    over-assembly.
    """
    instances = _parse_bed_instances(bed_path)
    if not instances:
        logger.warning("No instances parsed from BED; skipping assembly")
        return records

    id_to_rec = {rec['id']: rec for rec in records}
    r_to_id = {}
    for rec in records:
        r_to_id[f"R={rec['R']}"] = rec['id']

    copies = {rec['id']: rec.get('copies', 1) for rec in records}
    seq_lens = {rec['id']: len(rec['sequence']) for rec in records}

    # Collect all instances sorted by (chrom, start)
    all_inst = []
    for fam_key, locs in instances.items():
        mid = r_to_id.get(fam_key)
        if mid is None:
            continue
        for chrom, start, end, strand in locs:
            all_inst.append((chrom, start, end, strand, mid))
    all_inst.sort(key=lambda x: (x[0], x[1]))

    # Sliding-window cross-family adjacency edges (shared implementation): looks ahead
    # up to stitch_window_k instances within stitch_gap, so a third family interleaved
    # between two fragments of one element no longer breaks the chain (the old scan only
    # linked the immediate i/i+1 pair). stitch_gap (~300) also recovers co-linear
    # fragments separated by a short divergent linker that fragment_gap=50 excluded.
    stitch_gap = getattr(config, 'stitch_gap', config.fragment_gap)
    window_k = getattr(config, 'stitch_window_k', 4)
    directed_pairs = fragment_stitcher.build_edges(
        all_inst, seq_lens, stitch_gap, window_k, config.max_overlap_fraction)

    scored = _score_and_filter_edges(directed_pairs, copies, config,
                                     check_direction=True)
    return _build_chains_and_scaffolds(scored, id_to_rec, records, config)


def assemble_fragments_blastn(records: List[Dict], config) -> List[Dict]:
    """Fragment assembly via BLASTN-to-genome coordinate co-occurrence.

    Fallback when BED instances are not available.  Uses directed greedy
    chain algorithm (without direction-consistency filter since BLASTN
    hits lack strand aggregation from the BED path).
    """
    if not config.genome_blast_db:
        logger.info("No BLAST DB for fragment assembly; skipping")
        return records

    # Write query FASTA
    tmp_dir = tempfile.mkdtemp(prefix='frag_asm_')
    query_path = os.path.join(tmp_dir, 'query.fa')
    with open(query_path, 'w') as fh:
        for rec in records:
            fh.write(f">{rec['id']}\n{rec['sequence']}\n")

    out_path = os.path.join(tmp_dir, 'blast.out')
    cmd = [
        config.blastn_exe,
        '-query', query_path,
        '-db', config.genome_blast_db,
        '-outfmt', '6 qseqid sseqid sstart send qlen pident',
        '-max_target_seqs', '500',
        '-max_hsps', '3',
        '-evalue', '1e-10',
        '-num_threads', str(config.threads),
        '-dust', 'no',
        '-out', out_path,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if result.returncode != 0:
            logger.warning(f"BLASTN for fragment assembly failed: {result.stderr[:200]}")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return records
    except subprocess.TimeoutExpired:
        logger.warning("BLASTN for fragment assembly timed out")
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return records

    # Parse hits: family -> [(chrom, start, end)]
    hits_by_fam = defaultdict(list)
    with open(out_path) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            qid = parts[0]
            sid = parts[1]
            sstart, send = int(parts[2]), int(parts[3])
            if sstart > send:
                sstart, send = send, sstart
            hits_by_fam[qid].append((sid, sstart, send))

    shutil.rmtree(tmp_dir, ignore_errors=True)

    # Flatten + sort all hits by position
    all_hits = []
    for fam_id, locs in hits_by_fam.items():
        for sid, s, e in locs:
            all_hits.append((sid, s, e, fam_id))
    all_hits.sort(key=lambda x: (x[0], x[1]))

    # Scan adjacent pairs → directed edges
    id_to_rec = {r['id']: r for r in records}
    copies = {r['id']: r.get('copies', 1) for r in records}
    seq_lens = {r['id']: len(r['sequence']) for r in records}

    directed_pairs = defaultdict(lambda: {'count': 0})
    gap_limit = config.fragment_gap
    max_overlap_frac = config.max_overlap_fraction

    for i in range(len(all_hits) - 1):
        chr1, start1, end1, id1 = all_hits[i]
        chr2, start2, end2, id2 = all_hits[i + 1]
        if chr1 != chr2 or id1 == id2:
            continue
        gap = start2 - end1
        if gap > gap_limit:
            continue
        if gap < 0:
            shorter = min(seq_lens.get(id1, 0), seq_lens.get(id2, 0))
            if shorter > 0 and (-gap) > shorter * max_overlap_frac:
                continue
        directed_pairs[(id1, id2)]['count'] += 1

    # No direction check: BLASTN hits lose strand info after coordinate normalization
    scored = _score_and_filter_edges(directed_pairs, copies, config,
                                     check_direction=False)
    return _build_chains_and_scaffolds(scored, id_to_rec, records, config)


# ═══════════════════════════════════════════════════════════════════════
# Step 0.5  Two-Round Deduplication
# ═══════════════════════════════════════════════════════════════════════

def _run_cdhit(input_path: str, output_path: str, identity: float,
               word_size: int, config) -> Optional[str]:
    """Run CD-HIT-EST and return output path or None on failure."""
    cmd = [
        config.cdhit_exe,
        '-i', input_path,
        '-o', output_path,
        '-c', str(identity),
        '-n', str(word_size),
        '-aS', '0.8',
        '-g', '1',
        '-r', '1',
        '-d', '0',
        '-T', str(config.threads),
        '-M', '0',
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode == 0 and os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            return output_path
        logger.warning(f"CD-HIT-EST failed (rc={result.returncode}): {result.stderr[:200]}")
    except subprocess.TimeoutExpired:
        logger.warning("CD-HIT-EST timed out")
    except Exception as e:
        logger.warning(f"CD-HIT-EST error: {e}")
    return None


def deduplicate(records: List[Dict], config) -> List[Dict]:
    """Single-round CD-HIT-EST deduplication at near-exact identity.

    Pre-refine dedup collapses only duplicate *descriptions* of the same family
    (`dedup_round1_identity`, default 0.98).  The looser 0.90 second round was
    removed from Phase 0: boundary-convergent duplicates are only exposed after
    Phase 1's extend-align-trim, so that merge now lives in Phase 1's post-refine
    pass (`merge_post_refine` @ `post_refine_merge_identity`).  `dedup_round2_identity`
    is retained as a back-compat config field but is no longer referenced here.
    """
    if len(records) < 2:
        return records

    tmp_dir = tempfile.mkdtemp(prefix='dedup_')
    in_path = os.path.join(tmp_dir, 'input.fa')

    # Write input sorted by (mdl_score desc, length desc) so cd-hit's first-comer
    # cluster representative aligns with mdl-repeat's MDL ranking, not raw length.
    # Length-only sort would let a long but low-mdl family overshadow a shorter
    # but high-copy/high-mdl family in the same cluster.
    records_sorted = sorted(
        records,
        key=lambda r: (r.get('mdl', 0.0), len(r['sequence'])),
        reverse=True,
    )
    with open(in_path, 'w') as fh:
        for rec in records_sorted:
            fh.write(f">{rec['id']}\n{rec['sequence']}\n")

    before = len(records)

    # Single round at near-exact identity (word size 11 for 0.95-1.0).
    r1_out = os.path.join(tmp_dir, 'round1')
    r1_path = _run_cdhit(in_path, r1_out, config.dedup_round1_identity, 11, config)
    if not r1_path:
        logger.warning("Dedup CD-HIT-EST failed; skipping deduplication")
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return records

    # Read surviving IDs
    surviving_ids = set()
    with open(r1_path) as fh:
        for line in fh:
            if line.startswith('>'):
                surviving_ids.add(line[1:].strip().split()[0])

    result = [r for r in records if r['id'] in surviving_ids]

    shutil.rmtree(tmp_dir, ignore_errors=True)

    logger.info(f"Deduplication: {before} -> {len(result)} sequences "
                f"(single round @ {config.dedup_round1_identity})")
    return result


# ═══════════════════════════════════════════════════════════════════════
# Step 0.6  Per-family instance index (Phase 0 → Phase 1)
# ═══════════════════════════════════════════════════════════════════════

_MEMBER_R_RE = re.compile(r'mdl_R(\d+)$')


def _length_class_expected_copies(length: int, table) -> int:
    """Expected minimum BED copies for a family of this consensus length.

    `table` is config.length_class_copy_expectation: an iterable of
    (min_length_inclusive, expected_min_copies) entries.  Walk it ordered by
    length floor descending and return the first class whose floor <= length.
    The default table's `(0, ...)` catch-all guarantees a hit.  Returns 0 if the
    table is empty (signal effectively disabled).
    """
    best_floor = -1
    best_expect = 0
    for floor, expect in table:
        if length >= floor and floor > best_floor:
            best_floor = floor
            best_expect = expect
    return best_expect


def write_instance_index(records: List[Dict], bed_path: str,
                         out_path: str,
                         presignal_path: Optional[str] = None,
                         config=None) -> Optional[str]:
    """Stream mdl_repeat.instances.bed once → a per-record instance index TSV.

    Output columns (plan §2.1), one row per surviving BED instance:
        0 record_id    owner record (scaffold members carry the scaffold id)
        1 member_R     original family number from BED col 3 (R=N)
        2 chrom        BED col 0 (clean chr1..chrN)
        3 start        BED col 1 (0-based half-open)
        4 end          BED col 2 (0-based half-open)
        5 strand       BED col 5 (+/-)
        6 divergence   1 - int(BED col 4)/1000

    The `R=N → record_id` routing maps each family's BED rows to the record that
    owns them: a normal record owns its own `R`, a scaffold record owns every member
    family's `R` (member ids `mdl_R<N>` are parsed back to N).  Instances of families
    that were hard-filtered out in Phase 0 have no owner and are dropped (their
    families are not refined).  The file is finally `sort -k1,1` so Phase 1 can route
    by a single streaming group-by with no in-RAM index — what keeps 10-30 Gb
    tractable.  Returns out_path, or None if the BED is unreadable.

    N1 extension (REFINE_STRATEGY_DESIGN_v2.md §2.4): when `presignal_path` and
    `config` are supplied, the SAME single BED pass also accumulates a per-record
    completeness pre-signal (instance count, BED-divergence mean and spread, and
    a length-class copy-expectation flag) and writes it to `presignal_path`.  This
    adds only O(1) per-instance accumulation and one O(records) write at the end —
    no extra BED pass, no genome touch.  The instances_index output and its sort
    are unchanged whether or not the pre-signal is requested.
    """
    if not bed_path or not os.path.isfile(bed_path):
        logger.warning("write_instance_index: BED file missing (%s); no index",
                       bed_path)
        return None

    # Build R=N -> owner record_id (scaffold members remap to the scaffold id).
    r_to_recid: Dict[int, str] = {}
    for rec in records:
        if rec.get('is_scaffold'):
            for member_id in rec.get('scaffold_members', []):
                m = _MEMBER_R_RE.match(member_id)
                if m:
                    r_to_recid[int(m.group(1))] = rec['id']
        else:
            r_to_recid[int(rec['R'])] = rec['id']

    # N1: per-record streaming accumulators for the completeness pre-signal.
    # Welford-style sum/sumsq give mean + population std in one pass (O(1)/inst).
    emit_presignals = bool(presignal_path) and config is not None
    acc_n: Dict[str, int] = defaultdict(int)
    acc_sum: Dict[str, float] = defaultdict(float)
    acc_sumsq: Dict[str, float] = defaultdict(float)

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    tmp_path = out_path + '.unsorted'
    n_written = 0
    n_no_owner = 0
    with open(bed_path) as fh, open(tmp_path, 'w') as out:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start = parts[1]
            end = parts[2]
            fam = parts[3]                      # "R=N"
            score = parts[4]
            strand = parts[5]
            if not fam.startswith('R='):
                continue
            try:
                r_num = int(fam[2:])
            except ValueError:
                continue
            recid = r_to_recid.get(r_num)
            if recid is None:
                n_no_owner += 1
                continue
            try:
                divergence = 1.0 - int(score) / 1000.0
            except ValueError:
                divergence = 0.0
            out.write(f"{recid}\t{r_num}\t{chrom}\t{start}\t{end}\t{strand}\t"
                      f"{divergence:.4f}\n")
            n_written += 1
            if emit_presignals:
                acc_n[recid] += 1
                acc_sum[recid] += divergence
                acc_sumsq[recid] += divergence * divergence

    # Sort by record_id (col 1) so Phase 1 routing is a streaming group-by.
    try:
        with open(out_path, 'w') as out_fh:
            rc = subprocess.run(['sort', '-k1,1', tmp_path], stdout=out_fh,
                                stderr=subprocess.PIPE, text=True)
        if rc.returncode != 0:
            logger.warning("instance index sort failed (%s); using unsorted order",
                           rc.stderr[:200])
            shutil.move(tmp_path, out_path)
        else:
            os.remove(tmp_path)
    except Exception as e:  # noqa: BLE001
        logger.warning("instance index sort error: %s; using unsorted order", e)
        if os.path.exists(tmp_path):
            shutil.move(tmp_path, out_path)

    logger.info("Instance index: %d rows written to %s (%d instances had no "
                "surviving owner family)", n_written, out_path, n_no_owner)

    if emit_presignals:
        _write_presignals(records, acc_n, acc_sum, acc_sumsq,
                          presignal_path, config)

    return out_path


def _write_presignals(records: List[Dict], acc_n: Dict[str, int],
                      acc_sum: Dict[str, float], acc_sumsq: Dict[str, float],
                      presignal_path: str, config) -> None:
    """Write the per-record completeness pre-signal TSV (N1, §2.4).

    One row per surviving record (every record gets a row even if it has 0 BED
    instances, so Phase 1 completeness assessment sees the whole family set).

    Columns:
        0 record_id
        1 n_instances              BED instance count owned by this record
        2 bed_divergence_mean      mean per-instance divergence (NaN-safe: '' if n=0)
        3 bed_divergence_spread    population std of per-instance divergence
        4 length                   record consensus length
        5 length_class_expectation expected min copies for this length class
        6 length_class_flag        1 if n_instances < expectation else 0
    """
    table = getattr(config, 'length_class_copy_expectation', ())
    os.makedirs(os.path.dirname(os.path.abspath(presignal_path)), exist_ok=True)
    n_rows = 0
    with open(presignal_path, 'w') as out:
        out.write("record_id\tn_instances\tbed_divergence_mean\t"
                  "bed_divergence_spread\tlength\tlength_class_expectation\t"
                  "length_class_flag\n")
        for rec in records:
            rid = rec['id']
            n = acc_n.get(rid, 0)
            length = rec.get('actual_length', len(rec.get('sequence', '')))
            expect = _length_class_expected_copies(length, table)
            if n > 0:
                mean = acc_sum[rid] / n
                var = max(0.0, acc_sumsq[rid] / n - mean * mean)
                spread = var ** 0.5
                mean_s = f"{mean:.4f}"
                spread_s = f"{spread:.4f}"
            else:
                mean_s = ""
                spread_s = ""
            flag = 1 if n < expect else 0
            out.write(f"{rid}\t{n}\t{mean_s}\t{spread_s}\t{length}\t"
                      f"{expect}\t{flag}\n")
            n_rows += 1
    logger.info("Completeness pre-signals: %d records -> %s", n_rows,
                presignal_path)


# ═══════════════════════════════════════════════════════════════════════
# Phase 0 Main Entry
# ═══════════════════════════════════════════════════════════════════════

def run_phase0(config) -> Dict:
    """Execute Phase 0: triage pipeline.

    Returns dict with keys:
        records: List[Dict]  — tiered, deduplicated records
        stats: Dict          — per-step statistics
    """
    logger.info("=" * 60)
    logger.info("Phase 0: Metadata Parsing + Quality Triage")
    logger.info("=" * 60)

    stats = {}

    # Step 0.1: Parse
    records = parse_mdl_fasta(config.input_file)
    stats['input_count'] = len(records)

    # Step 0.2: Hard filter
    records, filter_stats = hard_filter(records, config)
    stats['filter'] = filter_stats

    # Step 0.3: Tiering
    records = assign_tiers(records, config)
    tier_counts = defaultdict(int)
    for r in records:
        tier_counts[r['tier']] += 1
    stats['tiers_after_filter'] = dict(tier_counts)

    # Fail loud on a BED/genome coordinate-space mismatch BEFORE it silently no-ops
    # Phase 1 (full-genome BED chr1..chrN vs a sample with windowed names chr*_w*):
    # extraction would recruit nothing and raw mdl fragments would be shipped.
    _bed_genome_coord_check(config.bed_file, config.genome_file)

    # Step 0.4: Fragment assembly
    if config.bed_file and os.path.isfile(config.bed_file):
        logger.info(f"Using BED instances for fragment assembly: {config.bed_file}")
        records = assemble_fragments_bed(records, config.bed_file, config)
    elif config.genome_blast_db:
        logger.info("Using BLASTN for fragment assembly")
        records = assemble_fragments_blastn(records, config)
    else:
        logger.info("No BED or BLAST DB; skipping fragment assembly")
    stats['after_assembly'] = len(records)

    # Step 0.5: Deduplication
    records = deduplicate(records, config)
    stats['after_dedup'] = len(records)

    # Re-tier after assembly/dedup may have changed composition
    # (scaffolds get the best tier of their members, which is correct)
    tier_counts_final = defaultdict(int)
    for r in records:
        tier_counts_final[r['tier']] += 1
    stats['tiers_final'] = dict(tier_counts_final)

    # Step 0.6: Per-family instance index for Phase 1 BED-seeded extraction.
    # N1: the same single BED pass also emits per-family completeness pre-signals
    # (§2.4) when a BED is present — no extra pass, no genome touch.
    instance_index = None
    completeness_presignals = None
    if config.bed_file and os.path.isfile(config.bed_file):
        index_path = os.path.join(config.temp_dir, 'instances_index.tsv')
        presignal_path = os.path.join(config.temp_dir,
                                      'completeness_presignals.tsv')
        instance_index = write_instance_index(
            records, config.bed_file, index_path,
            presignal_path=presignal_path, config=config)
        if instance_index is not None:
            completeness_presignals = presignal_path
    else:
        logger.info("No BED file; Phase 1 will run without a BED instance index "
                    "(families fall back to original seeds / blastn recruitment)")
    stats['instance_index'] = instance_index
    stats['completeness_presignals'] = completeness_presignals

    logger.info(f"Phase 0 complete: {stats['input_count']} -> {len(records)} sequences")
    logger.info(f"  Final tiers: {dict(tier_counts_final)}")

    return {'records': records, 'stats': stats,
            'instance_index': instance_index,
            'completeness_presignals': completeness_presignals}
