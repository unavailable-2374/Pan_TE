"""
Phase 2: Protein-homology gene exclusion (negative filter).

Multi-copy host gene families — F-box (~700 genes in Arabidopsis), NBS-LRR disease
resistance, kinases, transporters, glycosyltransferases — are recurrent and dispersed,
so a de novo repeat finder (mdl-repeat) collects them as "repeats". In the masking
library they over-mask coding sequence (measured: ~19.5% of all CDS on Col-CEN T2T),
which is false-positive masking, not TE sensitivity.

This module removes that noise WITHOUT using the genome's own gene annotation (which
does not exist yet at TE-annotation time) and WITHOUT structural gating (which sheds
divergent / dark-matter TEs). It is a NEGATIVE filter: a consensus is dropped iff

    (a) it strongly matches a CELLULAR protein in an external, TE-removed protein DB
        (coverage >= gene_excl_min_cov AND pident >= gene_excl_min_pident), AND
    (b) it shows ZERO TE evidence — the te_signal list from annotate_te_structure is
        empty (no protein_domain / tir / ltr / known_te_homology).

Because ANY TE signal keeps the sequence, divergent, no-ORF, and structural TEs are
never dropped (validation on Col-CEN: only 1.1% of EDTA-confirmed TE masking touched,
vs 28.5% of CDS over-masking removed). Domesticated transposase genes (FAR1/FHY3/FRS)
typically also hit RepeatPeps and are spared by guard (b).

Data integrity: if cellular_protein_db is unset or its BLAST DB files are missing, the
step is SKIPPED (every record kept) and logged — it never drops a sequence without the
discriminating evidence in hand.
"""

import logging
import os
import subprocess
import tempfile
import shutil
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)


def _cellular_db_ready(db_path: str) -> bool:
    """A protein BLAST DB is usable if the .pin (or .pal multi-volume) index exists."""
    if not db_path:
        return False
    return os.path.isfile(db_path + '.pin') or os.path.isfile(db_path + '.pal')


def _blastx_cellular(fasta_path: str, db_path: str, config,
                     work_dir: str) -> Dict[str, Tuple[float, float, str]]:
    """blastx the consensus library vs the cellular protein DB.

    Returns {qid: (max_coverage, pident_at_max_cov, subject_title)} where coverage =
    aligned_aa * 3 / nucleotide_qlen, taken as the per-query maximum over HSPs that
    meet the identity floor. Empty dict on any failure (caller treats as no hit).
    """
    # blastx -num_threads scales poorly on one query stream; split the query into
    # n_jobs chunks run as concurrent processes (each 1 thread) to use the cores.
    n_jobs = getattr(config, 'gene_excl_jobs', 0) or getattr(config, 'threads', 1)
    recs = []
    rid, seq = None, []
    for ln in open(fasta_path):
        if ln.startswith('>'):
            if rid:
                recs.append((rid, ''.join(seq)))
            rid = ln[1:].rstrip('\n'); seq = []
        else:
            seq.append(ln.strip())
    if rid:
        recs.append((rid, ''.join(seq)))
    n_jobs = max(1, min(n_jobs, len(recs)))
    chunk_paths, fhs = [], []
    for j in range(n_jobs):
        p = os.path.join(work_dir, f'gx_q{j}.fa'); chunk_paths.append(p); fhs.append(open(p, 'w'))
    for i, (h, s) in enumerate(recs):
        fhs[i % n_jobs].write(f">{h}\n{s}\n")
    for fh in fhs:
        fh.close()

    base = [getattr(config, 'blastx_exe', 'blastx'), '-db', db_path,
            '-outfmt', '6 qseqid pident length qlen evalue stitle',
            '-evalue', str(config.gene_excl_evalue),
            '-max_target_seqs', str(config.gene_excl_max_target_seqs),
            '-num_threads', '1', '-seg', 'yes']

    def _run(j):
        op = os.path.join(work_dir, f'gx_o{j}.out')
        try:
            r = subprocess.run(base + ['-query', chunk_paths[j], '-out', op],
                               capture_output=True, text=True, timeout=14400)
            return (j, r.returncode, r.stderr[:200], op)
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            return (j, -1, type(e).__name__, op)

    import concurrent.futures
    out_files = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as ex:
        for j, rc, err, op in ex.map(_run, range(n_jobs)):
            if rc != 0:
                logger.warning(f"gene-exclusion blastx chunk {j} failed ({rc} {err}); "
                               f"skipping exclusion")
                return {}
            out_files.append(op)

    pid_floor = config.gene_excl_min_pident
    best: Dict[str, Tuple[float, float, str]] = {}
    for out_path in out_files:
        if not os.path.exists(out_path):
            continue
        with open(out_path) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 6:
                    continue
                qid = parts[0]
                try:
                    pident = float(parts[1])
                    aln_len = int(parts[2])
                    qlen = int(parts[3])
                except ValueError:
                    continue
                if qlen <= 0 or pident < pid_floor:
                    continue
                cov = aln_len * 3.0 / qlen
                if qid not in best or cov > best[qid][0]:
                    best[qid] = (cov, pident, parts[5])
    return best


def exclude_gene_derived(records: List[Dict], signal_map: Dict[str, List[str]],
                         config) -> Tuple[List[Dict], List[Dict], Dict]:
    """Drop consensi that are cellular-gene-derived and carry no TE evidence.

    Args:
        records:    QC-passed consensus records (each has 'id', 'sequence').
        signal_map: {id: [te signals]} from annotate_te_structure (TE evidence guard).
        config:     RefinerMdlConfig.

    Returns (kept, dropped, stats). `dropped` records each get a 'gene_exclusion'
    annotation (coverage, pident, matched protein) for the provenance TSV.
    """
    stats = {'input': len(records), 'enabled': bool(config.enable_gene_exclusion),
             'dropped': 0, 'skipped_reason': ''}
    if not records:
        return records, [], stats
    if not config.enable_gene_exclusion:
        stats['skipped_reason'] = 'disabled'
        return records, [], stats

    db_path = config.cellular_protein_db
    if not _cellular_db_ready(db_path):
        stats['skipped_reason'] = f'cellular_protein_db not ready ({db_path or "unset"})'
        logger.warning(f"Gene exclusion SKIPPED — {stats['skipped_reason']}; "
                       f"all {len(records)} records kept (never dropped without DB)")
        return records, [], stats

    work_dir = tempfile.mkdtemp(prefix='gene_excl_')
    try:
        fasta_path = os.path.join(work_dir, 'consensi.fa')
        with open(fasta_path, 'w') as fh:
            for rec in records:
                fh.write(f">{rec['id']}\n{rec['sequence']}\n")

        logger.info(f"Gene exclusion: blastx {len(records)} consensi vs cellular DB "
                    f"({os.path.basename(db_path)})")
        cell = _blastx_cellular(fasta_path, db_path, config, work_dir)

        cov_cut = config.gene_excl_min_cov
        drop_cov = getattr(config, 'gene_excl_drop_cov', 0.7)
        kept, dropped = [], []
        n_cell_strong = n_domesticated = 0
        for rec in records:
            rid = rec['id']
            te_evidence = bool(signal_map.get(rid))      # any TE signal
            hit = cell.get(rid)
            cov = hit[0] if hit is not None else 0.0
            strong = cov >= cov_cut
            if strong:
                n_cell_strong += 1
            # Drop a host gene that the cellular DB recognizes near-full-length EVEN IF it
            # carries a TE domain (domesticated transposase: FAR1/FHY3/MUSTANG …). The DB
            # has TE proteins removed, so a >= drop_cov hit is specifically a host gene, not
            # a real TE. Below that, keep anything with TE evidence (chimeric → flank-trim).
            unconditional = cov >= drop_cov
            drop = unconditional or (strong and not te_evidence)
            if drop:
                rec['gene_exclusion'] = {
                    'coverage': round(cov, 3), 'pident': hit[1], 'protein': hit[2],
                    'reason': 'domesticated_gene(TE+gene)' if (unconditional and te_evidence)
                              else 'gene_derived',
                }
                if unconditional and te_evidence:
                    n_domesticated += 1
                dropped.append(rec)
            else:
                kept.append(rec)

        stats.update({
            'cellular_strong_hits': n_cell_strong,
            'domesticated_dropped': n_domesticated,
            'dropped': len(dropped),
            'kept': len(kept),
        })
        logger.info(f"Gene exclusion: {n_cell_strong}/{len(records)} strong cellular hit; "
                    f"dropped {len(dropped)} (incl {n_domesticated} domesticated TE+gene at "
                    f"cov>={drop_cov}); kept {len(kept)}")
        return kept, dropped, stats
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
