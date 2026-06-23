"""
Phase 2: genomic copy-number recruitment via rmblastn.

The `copies` field inherited from mdl-repeat is its MDL-based `accept=exclusive` instance
count, further lowered by fragment-assembly (min) and subfamily division, and NEVER
updated by N4 recall. Measured against RepeatMasker it disagrees with genomic reality in
both directions (~1/3 of the library gets a different verdict), so thresholding on it
wrongly drops real recurrent families and keeps non-families.

This module recruits the REAL genomic copy number: one rmblastn of every consensus
against the whole-genome BLAST DB (reusing phase1_fallback.ensure_genome_blast_db), then
merges HSPs into genomic INSTANCES (loci) per consensus and counts an instance as a COPY
when the union of consensus positions its HSPs cover >= copy_recruit_cov of the consensus
length. The result is written to rec['genomic_copies'], which the low-copy / hard-floor
filter reads in preference to the unreliable mdl `copies`.

rmblastn parameters (see config): RepeatMasker `25p<GC>g` matrix chosen by genome GC
(divergence-tolerant to capture diverged family members without re-creating the mdl
under-count), gapopen 30 / gapextend 6, -complexity_adjust, -dust no. Deterministic.
On any failure (no DB, rmblastn missing, no matrix) genomic_copies is left UNSET and the
filter falls back to mdl copies — never a silent zero.
"""

import logging
import os
import subprocess
import tempfile
from collections import defaultdict
from typing import Dict, List

logger = logging.getLogger(__name__)

_GC_BINS = [35, 37, 39, 41, 43, 45, 47, 49, 51, 53]


def _genome_gc_bin(genome_file: str, sample_bp: int = 20_000_000) -> int:
    """Pick the nearest RepeatMasker GC matrix bin from genome GC (sampled)."""
    gc = at = 0
    try:
        with open(genome_file) as fh:
            for line in fh:
                if line.startswith('>'):
                    continue
                s = line.strip().upper()
                gc += s.count('G') + s.count('C')
                at += s.count('A') + s.count('T')
                if gc + at >= sample_bp:
                    break
    except OSError:
        return 41  # neutral default
    if gc + at == 0:
        return 41
    pct = gc / (gc + at) * 100
    return min(_GC_BINS, key=lambda b: abs(b - pct))


def _find_matrix(config) -> str:
    """Locate the rmblastn ncbi/nt matrix 25p<GC>g.matrix for this genome's GC."""
    mdir = getattr(config, 'rm_matrix_dir', '') or os.path.join(
        os.environ.get('CONDA_PREFIX', ''), 'share', 'RepeatMasker',
        'Matrices', 'ncbi', 'nt')
    if not os.path.isdir(mdir):
        return ''
    gc_bin = _genome_gc_bin(config.genome_file)
    div = getattr(config, 'copy_recruit_divergence', 25)
    name = f"{div}p{gc_bin}g.matrix"
    path = os.path.join(mdir, name)
    return path if os.path.isfile(path) else ''


def _count_instances(hsps: List, qlen: int, cov_cut: float, gap: int) -> int:
    """Merge HSPs (per scaffold) into instances within `gap`; count those whose union
    consensus coverage >= cov_cut * qlen. hsps: list of (sid, gmin, gmax, qa, qb)."""
    by_scaf = defaultdict(list)
    for sid, gmin, gmax, qa, qb in hsps:
        by_scaf[sid].append((gmin, gmax, qa, qb))
    copies = 0
    for hs in by_scaf.values():
        hs.sort()
        cur_qs = [(hs[0][2], hs[0][3])]
        cur_end = hs[0][1]
        for gmin, gmax, qa, qb in hs[1:]:
            if gmin <= cur_end + gap:
                cur_qs.append((qa, qb)); cur_end = max(cur_end, gmax)
            else:
                if _union(cur_qs) >= cov_cut * qlen:
                    copies += 1
                cur_qs = [(qa, qb)]; cur_end = gmax
        if _union(cur_qs) >= cov_cut * qlen:
            copies += 1
    return copies


def _union(intervals) -> int:
    iv = sorted(intervals)
    tot = 0
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            tot += ce - cs + 1
            cs, ce = s, e
    return tot + ce - cs + 1


def recruit_genomic_copies(records: List[Dict], config) -> Dict:
    """rmblastn every consensus vs the genome; write rec['genomic_copies']. Returns stats.

    Leaves genomic_copies unset on any failure (caller falls back to mdl copies)."""
    stats = {'input': len(records), 'enabled': bool(getattr(config, 'enable_copy_recruit', False))}
    if not records or not getattr(config, 'enable_copy_recruit', False):
        stats['skipped_reason'] = 'disabled' if records else 'no_records'
        return stats

    # Lazy genome BLAST DB (shared with phase1 fallback recruiter). Bare import matches
    # the sibling convention (phase1_consensus); main.py puts Refiner_mdl on sys.path.
    from phase1_fallback import ensure_genome_blast_db
    db = ensure_genome_blast_db(config)
    if not db:
        stats['skipped_reason'] = 'genome_blast_db unavailable'
        logger.warning("Copy recruitment SKIPPED — no genome BLAST DB; falling back to mdl copies")
        return stats

    matrix = _find_matrix(config)
    if not matrix:
        stats['skipped_reason'] = 'rmblastn matrix not found'
        logger.warning("Copy recruitment SKIPPED — no RepeatMasker matrix; falling back to mdl copies")
        return stats

    exe = getattr(config, 'rmblastn_exe', 'rmblastn')
    cov_cut = getattr(config, 'copy_recruit_cov', 0.5)
    # rmblastn's -num_threads does not scale on a single query stream, so parallelise by
    # SPLITTING the query into n_jobs chunks run as concurrent processes (each 1 thread).
    n_jobs = getattr(config, 'copy_recruit_jobs', 0) or config.threads
    n_jobs = max(1, min(n_jobs, len(records)))
    work = tempfile.mkdtemp(prefix='copy_recruit_')
    try:
        qlen = {r['id']: len(r['sequence']) for r in records}
        # round-robin records into n_jobs chunk files (balances long/short across chunks)
        chunk_paths = [os.path.join(work, f'q{j}.fa') for j in range(n_jobs)]
        fhs = [open(p, 'w') for p in chunk_paths]
        for i, rec in enumerate(records):
            fhs[i % n_jobs].write(f">{rec['id']}\n{rec['sequence']}\n")
        for fh in fhs:
            fh.close()

        base = [exe, '-db', db, '-matrix', matrix,
                '-gapopen', str(getattr(config, 'copy_recruit_gapopen', 30)),
                '-gapextend', str(getattr(config, 'copy_recruit_gapextend', 6)),
                '-complexity_adjust', '-dust', 'no',
                '-evalue', str(getattr(config, 'copy_recruit_evalue', 1e-5)),
                '-num_threads', '1',
                '-outfmt', '6 qseqid sseqid pident length qstart qend qlen sstart send']

        def _run_chunk(j):
            outp = os.path.join(work, f'o{j}.tsv')
            try:
                r = subprocess.run(base + ['-query', chunk_paths[j], '-out', outp],
                                   capture_output=True, text=True,
                                   timeout=getattr(config, 'copy_recruit_timeout_s', 21600))
                return (j, r.returncode, r.stderr[:200], outp)
            except (subprocess.TimeoutExpired, FileNotFoundError) as e:
                return (j, -1, type(e).__name__, outp)

        logger.info(f"Copy recruitment: rmblastn {len(records)} consensi vs genome "
                    f"({os.path.basename(matrix)}, cov>={cov_cut}) — {n_jobs}-way parallel")
        import concurrent.futures
        out_files = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as ex:
            for j, rc, err, outp in ex.map(_run_chunk, range(n_jobs)):
                if rc != 0:
                    stats['skipped_reason'] = f'rmblastn chunk {j} rc={rc} ({err})'
                    logger.warning(f"Copy recruitment chunk {j} failed ({rc} {err}); "
                                   f"falling back to mdl copies")
                    return stats
                out_files.append(outp)

        hsp = defaultdict(list)
        for outp in out_files:
            if not os.path.exists(outp):
                continue
            for ln in open(outp):
                p = ln.rstrip('\n').split('\t')
                if len(p) < 9:
                    continue
                try:
                    qs, qe = int(p[4]), int(p[5])
                    ss, se = int(p[7]), int(p[8])
                except ValueError:
                    continue
                gmin, gmax = (ss, se) if ss <= se else (se, ss)
                qa, qb = (qs, qe) if qs <= qe else (qe, qs)
                hsp[p[0]].append((p[1], gmin, gmax, qa, qb))

        n_set = 0
        for rec in records:
            ql = qlen[rec['id']]
            gc = _count_instances(hsp.get(rec['id'], []), ql, cov_cut, ql) if rec['id'] in hsp else 0
            rec['genomic_copies'] = gc
            n_set += 1
        stats.update({'recruited': n_set, 'matrix': os.path.basename(matrix), 'cov_cut': cov_cut})
        logger.info(f"Copy recruitment: wrote genomic_copies for {n_set} consensi "
                    f"(no-hit -> 0); these now drive the low-copy / hard-floor filter")
        return stats
    finally:
        import shutil
        shutil.rmtree(work, ignore_errors=True)
