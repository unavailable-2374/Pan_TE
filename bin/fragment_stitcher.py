"""
fragment_stitcher — shared fragment-chain stitching for both de novo TE tracks.

A de novo repeat finder emits one real element as several *fragment* families
(conserved sub-blocks separated by divergent linker) that recur in a fixed co-linear
order across genomic copies. This module detects those co-linear chains from a BED of
member instances and (for tracks that need it) re-builds one full-length consensus per
chain by re-extracting the spanning genomic copies and aligning them.

Design notes
------------
* The chain-DETECTION functions (`build_edges`, `score_filter_edges`, `extract_chains`)
  are PURE (stdlib only) and live at module top level. `Refiner_mdl/phase0_triage.py`
  imports them, so keeping them dependency-free guarantees no import cycle.
* The re-CONSENSUS functions (`restitch_consensus`, `stitch_te_looker`) reuse the
  Refiner_mdl Phase-1 machinery (copy extraction + MSA + consensus). Those heavy imports
  are LAZY (inside the functions) so importing this module never drags in Phase 1.

Used by:
  - mdl track: phase0_triage.assemble_fragments_bed (detection only; re-consensus is the
    existing Phase 1, fed by the scaffold spanning instances).
  - te-looker track: Pan_TE.process_te_looker via `stitch_te_looker` (detection + re-consensus).
"""

import os
import sys
import logging
import subprocess
from collections import defaultdict

logger = logging.getLogger(__name__)

# Make the Refiner_mdl / Refiner dirs importable for the lazy Phase-1 reuse below.
_BIN = os.path.dirname(os.path.abspath(__file__))
for _d in (_BIN, os.path.join(_BIN, 'Refiner_mdl'), os.path.join(_BIN, 'Refiner')):
    if _d not in sys.path:
        sys.path.insert(0, _d)


# ════════════════════════════════════════════════════════════════════════
# Chain detection (PURE — no Refiner imports; safe for phase0_triage to import)
# ════════════════════════════════════════════════════════════════════════

def build_edges(all_inst, seq_lens, gap, window_k, max_overlap_frac):
    """Sliding-window cross-family adjacency edges.

    `all_inst` : list of (chrom, start, end, strand, famid) SORTED by (chrom, start).
    For each instance i, look ahead up to `window_k` instances on the same chrom+strand
    while the inter-instance gap <= `gap`, recording a directed cross-family edge
    famid(i) -> famid(j) for each DISTINCT different family seen (deduped per i). Unlike
    the old i/i+1-only scan, a third family interleaved between two fragments of one
    element no longer breaks the adjacency.

    Returns {(idA, idB): {'count': n}}.
    """
    directed = defaultdict(lambda: {'count': 0})
    n = len(all_inst)
    for i in range(n):
        c1, _s1, e1, st1, id1 = all_inst[i]
        seen = set()
        for j in range(i + 1, min(i + 1 + window_k, n)):
            c2, s2, _e2, st2, id2 = all_inst[j]
            if c2 != c1:
                break                       # sorted by chrom: nothing further on c1
            gapv = s2 - e1
            if gapv > gap:
                break                       # sorted by start: gap only grows from here
            if st2 != st1 or id2 == id1 or id2 in seen:
                continue
            if gapv < 0:                    # overlap guard (redundant descriptions)
                shorter = min(seq_lens.get(id1, 0), seq_lens.get(id2, 0))
                if shorter > 0 and (-gapv) > shorter * max_overlap_frac:
                    continue
            # Strand-normalize: a '-'-strand copy of an element has its fragments in
            # REVERSED genomic order, so the genomic edge id1->id2 corresponds to the
            # element's id2->id1. Flip on '-' so both strands vote for the SAME chain
            # direction; otherwise a bidirectional element splits its votes fwd/rev and
            # the direction-consistency gate wrongly kills the edge.
            pair = (id1, id2) if st1 == '+' else (id2, id1)
            directed[pair]['count'] += 1
            seen.add(id2)
    return dict(directed)


def score_filter_edges(directed_pairs, copies, config, check_direction=True,
                       cooccur_thr=None, cooccur_ratio=None, dir_con_min=None,
                       copy_ratio_min=None):
    """Aggregate directed adjacency pairs into undirected edges and apply the
    multi-dimensional coherence gates (co-occurrence count + ratio, direction
    consistency, copy-number consistency). Returns [{'left','right','ratio'}] sorted by
    ratio desc. (Generic; phase0_triage delegates to this so there is one implementation.)

    The four thresholds default to the strict mdl-tuned config values; callers may
    override them (the te-looker track relaxes co-occurrence/copy ratios because its
    discovery is noisier — a real chain still records far fewer adjacencies per copy than
    mdl-repeat — while keeping direction consistency as the chimera guard).
    """
    thr = config.cooccurrence_threshold if cooccur_thr is None else cooccur_thr
    ratio_min = config.cooccurrence_ratio if cooccur_ratio is None else cooccur_ratio
    dcon = config.min_direction_consistency if dir_con_min is None else dir_con_min
    cratio = config.min_copy_ratio if copy_ratio_min is None else copy_ratio_min

    undirected = {}
    for (a, b), stats in directed_pairs.items():
        key = (min(a, b), max(a, b))
        if key not in undirected:
            undirected[key] = {'fwd': 0, 'rev': 0}
        if a <= b:
            undirected[key]['fwd'] += stats['count']
        else:
            undirected[key]['rev'] += stats['count']

    scored = []
    for (a, b), stats in undirected.items():
        total = stats['fwd'] + stats['rev']
        if total < thr:
            continue
        min_cop = min(copies.get(a, 1), copies.get(b, 1))
        ratio = total / max(min_cop, 1)
        if ratio < ratio_min:
            continue
        if check_direction:
            dir_con = max(stats['fwd'], stats['rev']) / total
            if dir_con < dcon:
                continue
        cop_a, cop_b = copies.get(a, 1), copies.get(b, 1)
        copy_ratio = min(cop_a, cop_b) / max(cop_a, cop_b)
        if copy_ratio < cratio:
            continue
        left, right = (a, b) if stats['fwd'] >= stats['rev'] else (b, a)
        scored.append({'left': left, 'right': right, 'ratio': ratio})

    scored.sort(key=lambda e: e['ratio'], reverse=True)
    return scored


def extract_chains(scored_edges):
    """Greedy one-left/one-right directed-DAG chain extraction. Each node gets at most one
    left and one right neighbor (strongest-evidence first), preventing transitive
    over-assembly. Returns a list of ordered famid chains (length >= 2)."""
    right_of, left_of = {}, {}
    for e in scored_edges:
        l, r = e['left'], e['right']
        if l in right_of or r in left_of:
            continue
        right_of[l] = r
        left_of[r] = l

    visited, chains = set(), []
    for node in right_of:
        if node in visited or node in left_of:
            continue                        # not a chain head
        chain = [node]
        visited.add(node)
        cur = node
        while cur in right_of:
            nxt = right_of[cur]
            if nxt in visited:
                break                       # cycle guard
            chain.append(nxt)
            visited.add(nxt)
            cur = nxt
        if len(chain) > 1:
            chains.append(chain)
    return chains


# ════════════════════════════════════════════════════════════════════════
# Re-consensus (LAZY Refiner_mdl Phase-1 reuse; used by the te-looker track)
# ════════════════════════════════════════════════════════════════════════

def restitch_consensus(chain_famids, inst_by_fam, genome_file, fai_lengths,
                       config, stitch_gap, stitch_max_len):
    """Re-build ONE full-length consensus for a co-linear chain.

    inst_by_fam: {famid: [(chrom, start, end, strand), ...]}. For each (chrom, strand),
    merge the chain members' instances within `stitch_gap` into spanning intervals; keep
    intervals that cover >= 2 distinct chain members and whose length is in
    [80, stitch_max_len] (real co-occurrence, not a lone fragment). Re-extract those
    genomic spans and run the Phase-1 extract+MSA+consensus to produce a full-length
    consensus. Returns (sequence, n_copies) or None (never fabricates on failure).
    """
    try:
        from phase1_extract import Instance, extract_padded_copies
        from phase1_boundary import refine_family

        per = defaultdict(list)             # (chrom, strand) -> [(start, end, fam)]
        for fam in chain_famids:
            for (chrom, start, end, strand) in inst_by_fam.get(fam, []):
                per[(chrom, strand)].append((start, end, fam))

        spans = []
        for (chrom, strand), locs in per.items():
            locs.sort()
            cs = ce = None
            members = set()
            for (s, e, fam) in locs:
                if cs is None:
                    cs, ce, members = s, e, {fam}
                elif s - ce <= stitch_gap:
                    ce = max(ce, e)
                    members.add(fam)
                else:
                    if len(members) >= 2 and 80 <= (ce - cs) <= stitch_max_len:
                        spans.append(Instance(chrom, cs, ce, strand, 0.0))
                    cs, ce, members = s, e, {fam}
            if cs is not None and len(members) >= 2 and 80 <= (ce - cs) <= stitch_max_len:
                spans.append(Instance(chrom, cs, ce, strand, 0.0))

        if len(spans) < config.min_copies_for_msa:
            return None

        def pad_fn(elem_len):
            return min(max(int(elem_len * config.pad_fraction), config.pad_min),
                       config.pad_cap)

        copies = extract_padded_copies(spans, genome_file, fai_lengths, pad_fn, config)
        if len(copies) < config.min_copies_for_msa:
            return None

        seed = max((c.sequence for c in copies), key=len)
        rec = {'id': 'stitch', 'R': 0, 'copies': len(copies), 'sequence': seed,
               'length': len(seed), 'actual_length': len(seed), 'mdl': 1.0, 'tier': 'T2'}
        out = refine_family(rec, copies, config)
        seq = (out.get('sequence') or '').replace('-', '')
        if len(seq) < 80:
            return None
        return seq, len(copies)
    except Exception as e:                   # noqa: BLE001 - report + None, never fake
        logger.warning("restitch_consensus failed for chain %s: %s",
                       chain_famids[:3], e)
        return None


def _restitch_worker(args):
    """ProcessPoolExecutor worker: re-stitch ONE chain. Top-level (picklable).
    Re-asserts sys.path so the lazy phase1_* imports resolve under spawn as well as fork."""
    import os as _os
    _bin = _os.path.dirname(_os.path.abspath(__file__))
    for _p in (_bin, _os.path.join(_bin, 'Refiner_mdl')):
        if _p not in sys.path:
            sys.path.insert(0, _p)
    (chain, sub_inst, genome_file, fai_lengths, config, gap, maxlen) = args
    return restitch_consensus(chain, sub_inst, genome_file, fai_lengths,
                              config, gap, maxlen)


def _read_fasta_records(path):
    """Yield (header_line_without_'>', sequence) preserving header verbatim."""
    header, parts = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header = line[1:].rstrip('\n')
                parts = []
            else:
                parts.append(line.strip())
    if header is not None:
        yield header, ''.join(parts)


def stitch_te_looker(families_fasta, members_bed, genome_file, out_fasta,
                     config=None, stitch_gap=300, window_k=4, stitch_max_len=25000,
                     cooccur_ratio=0.1, copy_ratio=0.1, dir_con=0.7):
    """Stitch co-linear fragment chains in the te-looker discovery output.

    Detects chains from `members_bed` (6-col: chrom start end fam_N . strand), re-builds a
    full-length consensus per chain from the genome, and writes `out_fasta` = stitched
    consensi + every UNCHAINED families.fasta record verbatim. Returns out_fasta.
    """
    from phase1_extract import load_fai_lengths
    if config is None:
        from Refiner_mdl.config import RefinerMdlConfig  # package import (unambiguous)
        config = RefinerMdlConfig()
        config.genome_file = genome_file

    fai = genome_file + '.fai'
    if not os.path.exists(fai):
        subprocess.run([config.samtools_exe, 'faidx', genome_file], check=False)
    if not os.path.exists(fai):
        raise RuntimeError(f"genome index missing and could not be built: {fai}")
    fai_lengths = load_fai_lengths(fai)

    # families.fasta: '>fam_N_ref members=M len=L' -> famid 'fam_N'
    fam_records = {}                          # famid -> (header, seq)
    seq_lens = {}
    for header, seq in _read_fasta_records(families_fasta):
        famtok = header.split()[0]
        famid = famtok[:-4] if famtok.endswith('_ref') else famtok
        fam_records[famid] = (header, seq)
        seq_lens[famid] = len(seq)

    # disc.members.bed: famid = col 4 verbatim ('fam_N')
    inst_by_fam = defaultdict(list)
    all_inst = []
    with open(members_bed) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 6:
                continue
            chrom, start, end, famid, strand = p[0], int(p[1]), int(p[2]), p[3], p[5]
            inst_by_fam[famid].append((chrom, start, end, strand))
            all_inst.append((chrom, start, end, strand, famid))
    all_inst.sort(key=lambda x: (x[0], x[1]))

    copies = {fam: max(len(inst_by_fam.get(fam, [])), 1) for fam in fam_records}
    # seq_lens for any famid present only in the BED but not the FASTA: default 0 (skip).
    for fam in inst_by_fam:
        seq_lens.setdefault(fam, 0)

    directed = build_edges(all_inst, seq_lens, stitch_gap, window_k,
                           config.max_overlap_fraction)
    # te-looker: relax co-occurrence/copy ratios (noisier discovery), keep direction
    # consistency + co-occurrence count as the chimera guards (a real chain recurs in a
    # consistent strand-normalized order across >=cooccurrence_threshold loci).
    scored = score_filter_edges(directed, copies, config,
                                cooccur_ratio=cooccur_ratio,
                                copy_ratio_min=copy_ratio,
                                dir_con_min=dir_con)
    chains = extract_chains(scored)
    logger.info("te-looker stitch: %d families, %d edges -> %d scored -> %d chains",
                len(fam_records), len(directed), len(scored), len(chains))

    # Keep only chains whose families have a consensus record (>=2).
    valid = [(ci, [f for f in chain if f in fam_records])
             for ci, chain in enumerate(chains)]
    valid = [(ci, c) for ci, c in valid if len(c) >= 2]

    # Re-stitch each chain in PARALLEL — each is independent and MSA-bound. Pass only the
    # chain's own instances (small) per task, not the whole 2M-row index.
    stitched = []
    consumed = set()
    max_workers = max(1, int(getattr(config, 'threads', 1) or 1))
    tasks = [(c, {f: inst_by_fam.get(f, []) for f in c},
              genome_file, fai_lengths, config, stitch_gap, stitch_max_len)
             for _ci, c in valid]
    if max_workers > 1 and len(tasks) > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            futs = {ex.submit(_restitch_worker, t): (ci, c)
                    for (ci, c), t in zip(valid, tasks)}
            for fut in as_completed(futs):
                ci, chain = futs[fut]
                try:
                    res = fut.result()
                except Exception as e:                       # noqa: BLE001
                    logger.warning("restitch worker failed for chain %s: %s", chain[:3], e)
                    res = None
                if res:
                    seq, ncopies = res
                    stitched.append((f"stitched_{ci + 1}", ncopies, seq))
                    consumed.update(chain)
    else:
        for ci, chain in valid:
            res = restitch_consensus(chain, inst_by_fam, genome_file, fai_lengths,
                                     config, stitch_gap, stitch_max_len)
            if res:
                seq, ncopies = res
                stitched.append((f"stitched_{ci + 1}", ncopies, seq))
                consumed.update(chain)

    n_in = len(fam_records)
    with open(out_fasta, 'w') as out:
        for name, ncopies, seq in stitched:
            out.write(f">{name} members={ncopies} len={len(seq)}\n")
            for k in range(0, len(seq), 60):
                out.write(seq[k:k + 60] + "\n")
        for famid, (header, seq) in fam_records.items():
            if famid in consumed:
                continue
            out.write(f">{header}\n")
            for k in range(0, len(seq), 60):
                out.write(seq[k:k + 60] + "\n")

    n_out = len(stitched) + (n_in - len(consumed))
    logger.info("te-looker stitch: %d chains found, %d stitched (%d fragment families "
                "consumed); families %d -> %d",
                len(chains), len(stitched), len(consumed), n_in, n_out)
    return out_fasta
