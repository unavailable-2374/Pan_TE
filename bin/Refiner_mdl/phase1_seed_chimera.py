"""
Phase 1 — N7: rmblastn copy recruitment + nested/chimeric SEED resolution.

Two coupled jobs, both off ONE per-family rmblastn of the seed vs the genome:

  (#1) COPY SOURCE — the recruited genomic loci (merged HSPs → Instance objects) replace
       mdl's BED instances as the copy set fed to the MSA. rmblastn (25p divergence
       matrix) finds the diverged members the BED misses, so the consensus is more
       representative. The Instance objects feed the EXISTING extract_padded_copies.

  (#3 / N7) DE-NESTING — the mdl SEED itself may be a nested/chimeric assembly of two TEs
       (A+B). Detected from the HSP QUERY-COVERAGE profile: a query position p where few
       loci span across p AND the loci covering the left are a DISJOINT genomic population
       from those covering the right (A-copies vs B-copies). The seed is cut at p and each
       component re-recruited independently. Composes with N2 clustering and M3 (mosaic
       safety net) — N7 runs upstream, on the raw HSP geometry, before recruitment/MSA.

Both default OFF (ship dark). Skips cleanly (returns no instances) if rmblastn / matrix /
genome DB are unavailable, so the caller falls back to the BED path.
"""

import logging
import os
import subprocess
import tempfile
import concurrent.futures
from collections import defaultdict
from typing import Dict, List, Tuple

from phase1_extract import Instance
from phase2_copy_recruit import _find_matrix

logger = logging.getLogger(__name__)


def _rmblastn_hsps(seq: str, cfg) -> List[Tuple]:
    """rmblastn one sequence vs the genome DB. Returns HSPs as
    (sid, gmin, gmax, strand, qa, qb, pident). Empty on any failure."""
    db = getattr(cfg, 'genome_blast_db', '') or ''
    if not db:
        return []
    matrix = _find_matrix(cfg)
    if not matrix:
        return []
    exe = getattr(cfg, 'rmblastn_exe', 'rmblastn')
    work = tempfile.mkdtemp(prefix='seed_rmb_')
    try:
        q = os.path.join(work, 'seed.fa')
        with open(q, 'w') as fh:
            fh.write(f">seed\n{seq}\n")
        out = os.path.join(work, 'h.tsv')
        cmd = [exe, '-query', q, '-db', db, '-matrix', matrix,
               '-gapopen', str(getattr(cfg, 'copy_recruit_gapopen', 30)),
               '-gapextend', str(getattr(cfg, 'copy_recruit_gapextend', 6)),
               '-complexity_adjust', '-dust', 'no',
               '-evalue', str(getattr(cfg, 'copy_recruit_evalue', 1e-5)),
               '-num_threads', '1',
               '-outfmt', '6 sseqid sstart send qstart qend pident',
               '-out', out]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=getattr(cfg, 'seed_chimera_blast_timeout_s', 600))
            if r.returncode != 0:
                return []
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return []
        hsps = []
        if os.path.exists(out):
            for ln in open(out):
                p = ln.rstrip('\n').split('\t')
                if len(p) < 6:
                    continue
                try:
                    ss, se, qs, qe = int(p[1]), int(p[2]), int(p[3]), int(p[4])
                    pid = float(p[5])
                except ValueError:
                    continue
                strand = '+' if ss <= se else '-'
                gmin, gmax = (ss, se) if ss <= se else (se, ss)
                qa, qb = (qs, qe) if qs <= qe else (qe, qs)
                hsps.append((p[0], gmin, gmax, strand, qa, qb, pid))
        return hsps
    finally:
        import shutil
        shutil.rmtree(work, ignore_errors=True)


def _union_len(spans):
    iv = sorted(spans)
    tot = 0
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            tot += ce - cs + 1
            cs, ce = s, e
    return tot + ce - cs + 1


def hsps_to_instances_loci(hsps: List[Tuple], gap: int,
                           qlen: int = 0) -> Tuple[List[Instance], List[Dict]]:
    """Merge HSPs per (scaffold,strand) within `gap` into genomic loci. Returns
    (instances, loci) in parallel. Each locus dict carries the query footprint
    (qstart,qend), a locus key (for N7 breakpoint detection), and `cov` = the fraction of
    the consensus its HSPs cover (union/qlen) — used to drop partial gene-paralog matches
    from COPY EXTRACTION while keeping every locus for the N7 coverage profile."""
    by = defaultdict(list)
    for sid, gmin, gmax, strand, qa, qb, pid in hsps:
        by[(sid, strand)].append((gmin, gmax, qa, qb, pid))
    instances, loci = [], []
    for (sid, strand), hs in by.items():
        hs.sort()
        cur_g = [hs[0][0], hs[0][1]]
        cur_q = [(hs[0][2], hs[0][3])]
        cur_div = [hs[0][4]]
        def _flush(g, qspans, divs):
            qstart = min(a for a, _ in qspans); qend = max(b for _, b in qspans)
            div = 1.0 - (sum(divs) / len(divs)) / 100.0
            cov = (_union_len(qspans) / qlen) if qlen else 1.0
            instances.append(Instance(chrom=sid, start=g[0] - 1, end=g[1],
                                      strand=strand, divergence=max(0.0, div)))
            loci.append({'key': (sid, strand, g[0], g[1]),
                         'qstart': qstart, 'qend': qend, 'cov': cov})
        for gmin, gmax, qa, qb, pid in hs[1:]:
            if gmin <= cur_g[1] + gap:
                cur_g[1] = max(cur_g[1], gmax); cur_q.append((qa, qb)); cur_div.append(pid)
            else:
                _flush(cur_g, cur_q, cur_div)
                cur_g = [gmin, gmax]; cur_q = [(qa, qb)]; cur_div = [pid]
        _flush(cur_g, cur_q, cur_div)
    return instances, loci


def detect_seed_breakpoints(loci: List[Dict], qlen: int, cfg) -> List[int]:
    """N7 breakpoint detection from the locus query-footprint geometry."""
    min_comp = getattr(cfg, 'chimera_seed_min_component_bp', 80)
    min_side = getattr(cfg, 'chimera_seed_min_side_copies', 4)
    max_span = getattr(cfg, 'chimera_seed_max_span_frac', 0.15)
    min_disj = getattr(cfg, 'chimera_seed_min_disjoint', 0.8)
    margin = getattr(cfg, 'chimera_seed_margin_bp', 30)
    binbp = max(1, getattr(cfg, 'chimera_seed_bin_bp', 20))
    if len(loci) < 2 * min_side or qlen < 2 * min_comp:
        return []

    def depth(j):
        return sum(1 for L in loci if L['qstart'] <= j <= L['qend'])

    cands = []
    for p in range(min_comp, qlen - min_comp, binbp):
        span = sum(1 for L in loci if L['qstart'] <= p - margin and L['qend'] >= p + margin)
        left = [L for L in loci if L['qend'] < p + margin and L['qstart'] < p]
        right = [L for L in loci if L['qstart'] > p - margin and L['qend'] > p]
        if len(left) < min_side or len(right) < min_side:
            continue
        dloc = max(1, (depth(p - margin) + depth(p + margin)) // 2)
        if span / dloc > max_span:
            continue
        lset = {L['key'] for L in left}; rset = {L['key'] for L in right}
        inter = len(lset & rset)
        disj = 1.0 - inter / max(1, min(len(lset), len(rset)))
        if disj < min_disj:
            continue
        score = (1 - span / dloc) * disj * min(len(left), len(right)) / dloc
        cands.append((score, p))

    chosen = []
    for s, p in sorted(cands, reverse=True):
        if all(abs(p - c) >= min_comp for c in chosen):
            chosen.append(p)
        if len(chosen) >= getattr(cfg, 'chimera_seed_max_breakpoints', 3):
            break
    return sorted(chosen)


def split_seed(rec: Dict, bps: List[int], cfg) -> List[Dict]:
    """Cut the seed at breakpoints into component sub-seeds (>=2 or [] -> caller keeps)."""
    seq = rec['sequence']; Q = len(seq)
    min_comp = getattr(cfg, 'chimera_seed_min_component_bp', 80)
    cuts = [0] + list(bps) + [Q]
    out = []
    for k, (a, b) in enumerate(zip(cuts, cuts[1:]), 1):
        if b - a < min_comp:
            continue
        sub = dict(rec)
        sub['id'] = f"{rec['id']}_seedpart{k}"
        sub['sequence'] = seq[a:b]
        sub['actual_length'] = b - a
        sub['is_seed_chimera_part'] = True
        sub['consensus_source'] = 'original'
        sub.pop('_aln', None)
        out.append(sub)
    return out if len(out) >= 2 else []
