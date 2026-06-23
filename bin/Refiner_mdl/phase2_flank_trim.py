"""
Phase 2: protein-homology flank trimming.

Phase-1 occupancy trimming removes flanks that are NOT shared across copies, but it keeps
a flank that IS consistently present across copies — including a host-GENE flank that a
TE consensus over-extends into (the over-extension masks CDS at every copy). This step
removes that gene flank using PROTEIN homology, which the occupancy trim is blind to.

Driver = a CELLULAR-protein blastx hit (the cellular DB has TE proteins removed, so a hit
is specifically a host gene). For a consensus with a PARTIAL terminal gene hit (the cov
band between flank_trim_min_cov and gene_excl_drop_cov — the chimeric TE+gene case), the
gene-homologous TERMINAL span is cut and the TE core kept. (Near-full-length gene hits,
cov >= gene_excl_drop_cov, are domesticated/host genes dropped wholesale by gene exclusion,
not trimmed here.) TEsorter / TE structure is a GUARD, not the cut driver: we only trim a
TERMINAL gene span and require a substantial TE core to remain, so a real TE is never cut
down to its protein domain.

Conservative: only terminal gene spans are trimmed (internal gene insertions are left for
the chimera/de-nesting step); the retained core must stay >= flank_trim_retain_min bp.
"""

import logging
import os
import subprocess
import tempfile
import shutil
import concurrent.futures
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)


def _blastx_spans(records: List[Dict], db_path: str, config, work_dir: str) -> Dict[str, Tuple]:
    """Per record best cellular hit -> (cov, qstart, qend, qlen, pident). Chunked parallel."""
    n_jobs = getattr(config, 'flank_trim_jobs', 0) or getattr(config, 'threads', 1)
    n_jobs = max(1, min(n_jobs, len(records)))
    chunk_paths, fhs = [], []
    for j in range(n_jobs):
        p = os.path.join(work_dir, f'ft_q{j}.fa'); chunk_paths.append(p); fhs.append(open(p, 'w'))
    for i, rec in enumerate(records):
        fhs[i % n_jobs].write(f">{rec['id']}\n{rec['sequence']}\n")
    for fh in fhs:
        fh.close()
    base = [getattr(config, 'blastx_exe', 'blastx'), '-db', db_path,
            '-outfmt', '6 qseqid pident length qlen qstart qend evalue',
            '-evalue', str(getattr(config, 'gene_excl_evalue', 1e-5)),
            '-max_target_seqs', '3', '-num_threads', '1', '-seg', 'yes']

    def _run(j):
        op = os.path.join(work_dir, f'ft_o{j}.out')
        try:
            r = subprocess.run(base + ['-query', chunk_paths[j], '-out', op],
                               capture_output=True, text=True, timeout=14400)
            return (r.returncode, op)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return (-1, op)

    outs = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as ex:
        for rc, op in ex.map(_run, range(n_jobs)):
            if rc != 0:
                logger.warning("flank-trim blastx chunk failed; skipping flank trim")
                return {}
            outs.append(op)

    best: Dict[str, Tuple] = {}
    pid_floor = getattr(config, 'gene_excl_min_pident', 35.0)
    for op in outs:
        if not os.path.exists(op):
            continue
        for line in open(op):
            p = line.rstrip('\n').split('\t')
            if len(p) < 6:
                continue
            qid = p[0]
            try:
                pid = float(p[1]); aln = int(p[2]); ql = int(p[3])
                qs, qe = int(p[4]), int(p[5])
            except ValueError:
                continue
            if ql <= 0 or pid < pid_floor:
                continue
            cov = aln * 3.0 / ql
            if qid not in best or cov > best[qid][0]:
                best[qid] = (cov, min(qs, qe), max(qs, qe), ql, pid)
    return best


def trim_gene_flanks(records: List[Dict], config) -> Tuple[List[Dict], Dict]:
    """Trim terminal host-gene flanks (cov band [min,drop)) off TE consensi. Mutates
    rec['sequence'] in place for trimmed records; returns (records, stats)."""
    stats = {'input': len(records), 'enabled': bool(getattr(config, 'enable_flank_trim', False)),
             'trimmed': 0}
    if not records or not getattr(config, 'enable_flank_trim', False):
        stats['skipped_reason'] = 'disabled' if records else 'no_records'
        return records, stats
    db = getattr(config, 'cellular_protein_db', '')
    if not (db and (os.path.isfile(db + '.pin') or os.path.isfile(db + '.pal'))):
        stats['skipped_reason'] = 'cellular_protein_db not ready'
        logger.warning("Flank trim SKIPPED — cellular_protein_db not ready")
        return records, stats

    lo = getattr(config, 'flank_trim_min_cov', 0.3)
    hi = getattr(config, 'gene_excl_drop_cov', 0.7)
    margin = getattr(config, 'flank_trim_margin', 30)
    retain_min = getattr(config, 'flank_trim_retain_min', 100)

    work = tempfile.mkdtemp(prefix='flank_trim_')
    try:
        hits = _blastx_spans(records, db, config, work)
        n_trim = trimmed_bp = 0
        for rec in records:
            hit = hits.get(rec['id'])
            if not hit:
                continue
            cov, qs, qe, ql, pid = hit
            if not (lo <= cov < hi):          # only the partial-terminal chimeric band
                continue
            seq = rec['sequence']; L = len(seq)
            # blastx qstart/qend are 1-based on the nucleotide query
            gene5 = qs <= margin                       # gene span touches the 5' terminus
            gene3 = qe >= L - margin                   # gene span touches the 3' terminus
            new = None; side = ''
            if gene5 and not gene3 and (L - qe) >= retain_min:
                new = seq[qe:]; side = "5'"            # cut [0:qe], keep TE core 3'
            elif gene3 and not gene5 and (qs - 1) >= retain_min:
                new = seq[:qs - 1]; side = "3'"        # cut [qs:], keep TE core 5'
            if new and len(new) >= retain_min:
                rec['flank_trim'] = {'side': side, 'removed_bp': L - len(new),
                                     'cov': round(cov, 3), 'pident': pid}
                trimmed_bp += L - len(new)
                rec['sequence'] = new
                n_trim += 1
        stats.update({'trimmed': n_trim, 'trimmed_bp': trimmed_bp})
        logger.info(f"Flank trim: {n_trim} consensi had a terminal host-gene flank cut "
                    f"({trimmed_bp} bp removed; cov band [{lo},{hi}))")
        return records, stats
    finally:
        shutil.rmtree(work, ignore_errors=True)
