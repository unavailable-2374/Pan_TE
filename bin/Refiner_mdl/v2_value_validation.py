#!/usr/bin/env python3
"""
DECISIVE value-validation experiment for mdl-repeat refine v2 "subfamily multi-consensus".

Question: For REAL TEs (families carrying a TE structural signal), does v2's
multi-consensus (one consensus per qualifying subfamily cluster) improve member
recovery over the single-consensus baseline?

Pipeline (all on REAL full-Arabidopsis data; no fabricated numbers):
  1. From mdl_repeat.raw.fa select candidates: copies >= 10 AND div >= 0.2.
  2. Run filter_by_te_structure on candidate consensus -> TE_DIV (TE-signal subset).
  3. For each TE_DIV family: extract padded BED copies, run
       MULTI  = cluster_copies + qualifying_clusters -> refine_family per cluster
       SINGLE = all copies as one cluster -> refine_family once
  4. blastn each consensus set back to genome; per-family recovery = fraction of
     that family's BED instances hit by ANY consensus (pident>=70, >=50% instance covered).
  5. Aggregate: TE_DIV size, fraction splitting into >=2 subfamilies, multi-vs-single
     per-family recovery distribution, divergent-member recovery, flagship examples.

Outputs a JSON of all real counts and distributions.
"""
import json
import logging
import os
import re
import subprocess
import sys
import tempfile
import time
from statistics import mean, median

# Run from inside Refiner_mdl so sibling imports resolve.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from config import RefinerMdlConfig
from phase1_extract import (Instance, compute_pad, extract_padded_copies,
                            load_fai_lengths)
from phase1_cluster import cluster_copies, qualifying_clusters, Cluster
from phase1_boundary import refine_family
from te_structure_filter import filter_by_te_structure

logging.basicConfig(level=logging.WARNING,
                    format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger('v2val')
log.setLevel(logging.INFO)

# ── REAL data paths (confirmed) ──────────────────────────────────────────
BASE = '/scratch/shuoc/TE/Arabidopsis_thaliana/demo'
RAW_FA = f'{BASE}/mdl-repeat/tmp/mdl_repeat.raw.fa'
BED = f'{BASE}/mdl-repeat/tmp/mdl_repeat.instances.bed'
GENOME = f'{BASE}/genome/genome.fa'
FAI = f'{BASE}/genome/genome.fa.fai'
GENOME_DB = GENOME  # makeblastdb output prefix == genome.fa

OUT_JSON = f'{HERE}/v2_value_validation.json'
DONE = f'{HERE}/v2_value_validation.DONE'

# Candidate selection thresholds.
MIN_COPIES = 10
MIN_DIV = 0.2

# Recovery blastn acceptance (aligned with mdl-repeat's identity>=70% acceptance).
REC_PIDENT = 70.0
REC_INST_COVFRAC = 0.50   # >= 50% of the BED instance length covered by an HSP
# Divergent-member band: BED instance divergence in [0.12, 0.50].
DIV_MEMBER_LO = 0.12
DIV_MEMBER_HI = 0.50

HDR_RE = re.compile(
    r'>R=(\d+)\s+length=(\d+)\s+copies=(\d+)\s+mdl=([\d.]+)\s+div=([\d.]+)')


def parse_raw_fa(path):
    """Parse mdl_repeat.raw.fa -> {R: rec dict} mirroring phase0 record shape."""
    recs = {}
    cur = None
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if cur is not None:
                    s = ''.join(seq).upper()
                    cur['sequence'] = s
                    cur['actual_length'] = len(s)
                    recs[cur['R']] = cur
                m = HDR_RE.match(line)
                if not m:
                    cur = None
                    seq = []
                    continue
                r = int(m.group(1))
                cur = {
                    'id': f'mdl_R{r}',
                    'R': r,
                    'length': int(m.group(2)),
                    'copies': int(m.group(3)),
                    'mdl': float(m.group(4)),
                    'div': float(m.group(5)),
                }
                seq = []
            elif cur is not None:
                seq.append(line.strip())
    if cur is not None:
        s = ''.join(seq).upper()
        cur['sequence'] = s
        cur['actual_length'] = len(s)
        recs[cur['R']] = cur
    return recs


def load_bed_by_R(path):
    """Parse instances BED -> {R: [Instance,...]}. score=int(1000*(1-div))."""
    by_r = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            chrom, start, end, name, score, strand = parts[:6]
            m = re.match(r'R=(\d+)', name)
            if not m:
                continue
            r = int(m.group(1))
            try:
                s = int(start); e = int(end); sc = int(score)
            except ValueError:
                continue
            div = 1.0 - sc / 1000.0
            by_r.setdefault(r, []).append(
                Instance(chrom=chrom, start=s, end=e, strand=strand,
                         divergence=div))
    return by_r


def make_config():
    cfg = RefinerMdlConfig(
        input_file=RAW_FA,
        genome_file=GENOME,
        bed_file=BED,
        genome_blast_db=GENOME_DB,
        threads=8,
    )
    cfg.samtools_exe = 'samtools'   # from PATH (conda)
    cfg.abpoa_exe = 'abpoa'         # from PATH (conda)
    return cfg


def refine_multi(rec, copies, cfg):
    """v2 MULTI path: cluster -> qualify -> refine each qualifying cluster.

    Returns (list_of_consensus_seqs, n_subfamilies, audit)."""
    clusters = cluster_copies(copies, cfg)
    qcl, audit = qualifying_clusters(clusters, copies, cfg)
    seqs = []
    for cl in qcl:
        out = refine_family(dict(rec), cl.members, cfg)
        s = out.get('sequence', '')
        if s:
            seqs.append(s)
    return seqs, len(qcl), audit


def refine_single(rec, copies, cfg):
    """SINGLE baseline: all copies as ONE cluster -> one refined consensus."""
    out = refine_family(dict(rec), copies, cfg)
    s = out.get('sequence', '')
    return ([s] if s else []), out.get('consensus_source', 'original')


def blastn_consensus_to_genome(seqs, tag, cfg, work):
    """blastn a set of consensus seqs vs the genome. Returns list of HSP dicts
    with chrom/sstart/send/pident (1-based, normalized so sstart<=send)."""
    if not seqs:
        return []
    qfa = os.path.join(work, f'{tag}.fa')
    with open(qfa, 'w') as fh:
        for i, s in enumerate(seqs):
            fh.write(f'>{tag}_c{i}\n{s}\n')
    out = os.path.join(work, f'{tag}.tsv')
    cmd = [cfg.blastn_exe, '-query', qfa, '-db', GENOME_DB,
           '-outfmt', '6 qseqid sseqid pident length sstart send slen',
           '-evalue', '1e-5', '-dust', 'no',
           '-num_threads', str(cfg.threads), '-out', out]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        log.warning('blastn %s exit %d: %s', tag, r.returncode, r.stderr[:200])
        return []
    hsps = []
    if os.path.exists(out):
        with open(out) as fh:
            for line in fh:
                p = line.rstrip('\n').split('\t')
                if len(p) < 7:
                    continue
                try:
                    pident = float(p[2])
                    ss = int(p[4]); se = int(p[5])
                except ValueError:
                    continue
                chrom = p[1]
                lo, hi = (ss, se) if ss <= se else (se, ss)
                hsps.append({'chrom': chrom, 'lo': lo, 'hi': hi,
                             'pident': pident})
    return hsps


def recovery_for_family(instances, hsps):
    """Fraction of BED instances (0-based half-open) recovered by >=1 HSP.

    An instance is recovered if some HSP on the same chrom has pident>=REC_PIDENT
    and the overlap covers >= REC_INST_COVFRAC of the instance length.
    Returns (recovered_count, total, recovered_flags_list)."""
    # Index HSPs by chrom for speed.
    by_chrom = {}
    for h in hsps:
        if h['pident'] < REC_PIDENT:
            continue
        by_chrom.setdefault(h['chrom'], []).append((h['lo'], h['hi']))
    flags = []
    for inst in instances:
        ilen = inst.end - inst.start
        if ilen <= 0:
            flags.append(False)
            continue
        # BED 0-based half-open -> 1-based inclusive for overlap math.
        i_lo = inst.start + 1
        i_hi = inst.end
        rec = False
        for (lo, hi) in by_chrom.get(inst.chrom, ()):
            ov = min(i_hi, hi) - max(i_lo, lo) + 1
            if ov > 0 and ov >= REC_INST_COVFRAC * ilen:
                rec = True
                break
        flags.append(rec)
    return sum(flags), len(flags), flags


def main():
    t0 = time.time()
    log.info('Parsing raw.fa ...')
    recs = parse_raw_fa(RAW_FA)
    log.info('  total families parsed: %d', len(recs))

    # Step 1: candidate selection.
    cands = [r for r in recs.values()
             if r['copies'] >= MIN_COPIES and r['div'] >= MIN_DIV]
    log.info('Candidates (copies>=%d & div>=%.2f): %d',
             MIN_COPIES, MIN_DIV, len(cands))

    cfg = make_config()
    fai_lengths = load_fai_lengths(FAI)
    bed_by_r = load_bed_by_R(BED)

    # Step 2: TE structural filter on candidate consensus.
    log.info('Running filter_by_te_structure on %d candidates ...', len(cands))
    cand_records = [{'id': r['id'], 'sequence': r['sequence']} for r in cands]
    passed, te_stats = filter_by_te_structure(cand_records, cfg)
    te_ids = {p['id'] for p in passed}
    te_div = [r for r in cands if r['id'] in te_ids]
    log.info('TE_DIV (candidates with TE signal): %d / %d',
             len(te_div), len(cands))
    log.info('  TE-signal check hits: %s', te_stats)

    # Step 3-4: per-family multi vs single + recovery.
    per_family = []
    work_root = tempfile.mkdtemp(prefix='v2val_')
    n = len(te_div)
    for idx, r in enumerate(te_div, 1):
        R = r['R']
        instances = bed_by_r.get(R, [])
        if not instances:
            continue
        fam_work = os.path.join(work_root, f'R{R}')
        os.makedirs(fam_work, exist_ok=True)
        try:
            copies = extract_padded_copies(
                instances, GENOME, fai_lengths,
                lambda L: compute_pad(L, cfg), cfg)
            if len(copies) < cfg.min_copies_for_msa:
                continue

            multi_seqs, n_sf, audit = refine_multi(r, copies, cfg)
            single_seqs, single_src = refine_single(r, copies, cfg)

            m_hsps = blastn_consensus_to_genome(multi_seqs, 'multi', cfg, fam_work)
            s_hsps = blastn_consensus_to_genome(single_seqs, 'single', cfg, fam_work)

            m_rec, total, m_flags = recovery_for_family(instances, m_hsps)
            s_rec, _, s_flags = recovery_for_family(instances, s_hsps)

            # Divergent-member band recovery.
            dm_idx = [i for i, inst in enumerate(instances)
                      if DIV_MEMBER_LO <= inst.divergence <= DIV_MEMBER_HI]
            dm_total = len(dm_idx)
            dm_m = sum(1 for i in dm_idx if m_flags[i])
            dm_s = sum(1 for i in dm_idx if s_flags[i])

            per_family.append({
                'R': R, 'id': r['id'],
                'length': r['length'], 'copies_hdr': r['copies'],
                'div_hdr': r['div'],
                'n_instances': total,
                'n_copies_extracted': len(copies),
                'n_subfamilies': n_sf,
                'n_multi_consensus': len(multi_seqs),
                'n_single_consensus': len(single_seqs),
                'single_source': single_src,
                'fallback_single_cluster': audit.get('fallback_single_cluster'),
                'recovery_multi': m_rec / total if total else 0.0,
                'recovery_single': s_rec / total if total else 0.0,
                'recovered_multi': m_rec, 'recovered_single': s_rec,
                'div_member_total': dm_total,
                'div_member_rec_multi': dm_m,
                'div_member_rec_single': dm_s,
            })
        except Exception as e:  # noqa: BLE001
            log.warning('R=%d failed: %s', R, e)
            continue
        if idx % 25 == 0:
            log.info('  processed %d/%d families (%.0fs)',
                     idx, n, time.time() - t0)

    # Step 5: aggregate.
    split_fams = [f for f in per_family if f['n_subfamilies'] >= 2]
    n_te_div_processed = len(per_family)

    def dist(xs):
        if not xs:
            return {'n': 0}
        return {'n': len(xs), 'mean': mean(xs), 'median': median(xs),
                'min': min(xs), 'max': max(xs)}

    # multi vs single on the families that ACTUALLY split into >=2 subfamilies.
    split_delta = [f['recovery_multi'] - f['recovery_single'] for f in split_fams]
    split_multi = [f['recovery_multi'] for f in split_fams]
    split_single = [f['recovery_single'] for f in split_fams]

    # divergent-member recovery on split families (where members exist).
    dm_split = [f for f in split_fams if f['div_member_total'] > 0]
    dm_multi = [f['div_member_rec_multi'] / f['div_member_total'] for f in dm_split]
    dm_single = [f['div_member_rec_single'] / f['div_member_total'] for f in dm_split]

    # Flagship: split families with largest multi-over-single recovery gain.
    flagships = sorted(split_fams,
                       key=lambda f: f['recovery_multi'] - f['recovery_single'],
                       reverse=True)[:5]

    summary = {
        'data': {
            'raw_fa': RAW_FA, 'bed': BED, 'genome': GENOME,
            'total_families': len(recs),
            'note': 'REAL full-Arabidopsis data; all numbers from real runs.',
        },
        'thresholds': {
            'min_copies': MIN_COPIES, 'min_div': MIN_DIV,
            'rec_pident': REC_PIDENT, 'rec_inst_covfrac': REC_INST_COVFRAC,
            'div_member_band': [DIV_MEMBER_LO, DIV_MEMBER_HI],
            'subfamily_divergence_cut': cfg.subfamily_divergence_cut,
            'subfamily_min_members': cfg.subfamily_min_members,
            'subfamily_intra_homogeneity_floor': cfg.subfamily_intra_homogeneity_floor,
        },
        'step1_candidates': len(cands),
        'step2_te_structure': te_stats,
        'step2_te_div_count': len(te_div),
        'step3_te_div_processed': n_te_div_processed,
        'step3_split_families': len(split_fams),
        'step3_split_fraction': (len(split_fams) / n_te_div_processed
                                 if n_te_div_processed else 0.0),
        'recovery_all_te_div': {
            'multi': dist([f['recovery_multi'] for f in per_family]),
            'single': dist([f['recovery_single'] for f in per_family]),
            'delta_multi_minus_single': dist(
                [f['recovery_multi'] - f['recovery_single'] for f in per_family]),
        },
        'recovery_split_only': {
            'multi': dist(split_multi),
            'single': dist(split_single),
            'delta_multi_minus_single': dist(split_delta),
            'n_multi_gt_single': sum(1 for d in split_delta if d > 1e-9),
            'n_multi_lt_single': sum(1 for d in split_delta if d < -1e-9),
            'n_equal': sum(1 for d in split_delta if abs(d) <= 1e-9),
        },
        'divergent_member_recovery_split': {
            'n_families': len(dm_split),
            'multi': dist(dm_multi),
            'single': dist(dm_single),
            'delta': dist([m - s for m, s in zip(dm_multi, dm_single)]),
        },
        'flagships': flagships,
        'per_family': per_family,
        'elapsed_s': time.time() - t0,
    }

    with open(OUT_JSON, 'w') as fh:
        json.dump(summary, fh, indent=2)
    log.info('Wrote %s (%.0fs)', OUT_JSON, time.time() - t0)
    with open(DONE, 'w') as fh:
        fh.write('done\n')


if __name__ == '__main__':
    main()
