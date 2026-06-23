"""
Phase 2: Non-TE repeat separation (classification / routing, not deletion).

A de novo repeat finder collects every recurrent sequence, so the consensus library
contains non-transposon repeats that should NOT live in a TE library:

  * organelle insertions (NUMT / NUPT) — chloroplast/mito DNA in the nuclear genome
  * ribosomal DNA arrays (45S = 18S/5.8S/28S, and 5S) — megabase-scale tandem rDNA
  * satellite / simple tandem repeats — centromeric / telomeric arrays

On Col-CEN these are few entries but carry enormous masking weight (one 45S rDNA
consensus alone masks ~12.6 Mb), so leaving them tagged "TE" overstates genome TE
content by ~10 percentage points. This module routes them OUT of the TE library into
dedicated tracks (nonte_organelle.fa / nonte_rrna.fa / nonte_satellite.fa) rather than
deleting them — they remain available for whole-genome masking, just not credited as TE.

Each detector independently uses an external tool / reference and is INDEPENDENTLY
skippable: a missing reference or binary disables only that class (logged), never a
silent drop. Priority on multi-class hits: organelle > rrna > satellite.

This is orthogonal to gene exclusion (phase2_gene_exclusion): that removes protein-coding
host-gene noise from CDS; this reclassifies non-TE structural/tandem repeats.
"""

import logging
import os
import subprocess
import tempfile
import shutil
from typing import Dict, List, Set, Tuple

logger = logging.getLogger(__name__)


def _write_fasta(records: List[Dict], path: str) -> None:
    with open(path, 'w') as fh:
        for rec in records:
            fh.write(f">{rec['id']}\n{rec['sequence']}\n")


def _detect_organelle(fasta_path: str, config, work_dir: str) -> Set[str]:
    """blastn vs the organelle reference (chloroplast + mito). Returns ids whose
    summed hit coverage >= organelle_min_cov at >= organelle_min_pid identity."""
    ref = getattr(config, 'organelle_ref', '')
    if not ref or not os.path.isfile(ref):
        logger.info("organelle_ref unset/missing; skipping organelle detection")
        return set()
    db = os.path.join(work_dir, 'organelle_db')
    try:
        subprocess.run([config.makeblastdb_exe, '-in', ref, '-dbtype', 'nucl',
                        '-out', db], capture_output=True, text=True, timeout=300, check=True)
    except Exception as e:
        logger.warning(f"organelle makeblastdb failed: {e}; skipping")
        return set()
    out = os.path.join(work_dir, 'organelle.out')
    try:
        subprocess.run([config.blastn_exe, '-query', fasta_path, '-db', db,
                        '-outfmt', '6 qseqid pident length qlen', '-evalue', '1e-10',
                        '-max_target_seqs', '1', '-max_hsps', '5', '-dust', 'no',
                        '-num_threads', str(config.threads), '-out', out],
                       capture_output=True, text=True, timeout=3600, check=True)
    except Exception as e:
        logger.warning(f"organelle blastn failed: {e}; skipping")
        return set()
    pid_floor = config.organelle_min_pid
    cov_sum: Dict[str, float] = {}
    qlen: Dict[str, int] = {}
    if os.path.exists(out):
        for ln in open(out):
            p = ln.split('\t')
            if len(p) < 4:
                continue
            try:
                pid = float(p[1]); aln = int(p[2]); ql = int(p[3])
            except ValueError:
                continue
            if pid < pid_floor or ql <= 0:
                continue
            cov_sum[p[0]] = cov_sum.get(p[0], 0) + aln
            qlen[p[0]] = ql
    return {q for q in cov_sum if cov_sum[q] / qlen[q] >= config.organelle_min_cov}


def _detect_rrna(fasta_path: str, config, work_dir: str) -> Set[str]:
    """barrnap (eukaryote) — any rRNA feature flags the sequence as rDNA."""
    exe = getattr(config, 'barrnap_exe', 'barrnap')
    gff = os.path.join(work_dir, 'barrnap.gff')
    try:
        with open(gff, 'w') as gfh:
            subprocess.run([exe, '--kingdom', 'euk', '--threads', str(config.threads),
                            fasta_path], stdout=gfh, stderr=subprocess.DEVNULL,
                           timeout=3600, check=True)
    except FileNotFoundError:
        logger.info("barrnap not found; skipping rRNA detection")
        return set()
    except Exception as e:
        logger.warning(f"barrnap failed: {e}; skipping rRNA detection")
        return set()
    hits = set()
    if os.path.exists(gff):
        for ln in open(gff):
            if ln.startswith('#') or not ln.strip():
                continue
            hits.add(ln.split('\t')[0])
    return hits


def _detect_satellite(fasta_path: str, config, work_dir: str,
                      seqlen: Dict[str, int]) -> Set[str]:
    """TRF — flag sequences whose tandem-repeat coverage >= satellite_tandem_frac."""
    exe = getattr(config, 'trf_exe', 'trf')
    out = os.path.join(work_dir, 'trf.ngs')
    # TRF -ngs streams results to stdout; run inside work_dir (it litters cwd).
    try:
        with open(out, 'w') as ofh:
            subprocess.run([exe, os.path.basename(fasta_path),
                            '2', '7', '7', '80', '10', '50', '500', '-h', '-ngs'],
                           stdout=ofh, stderr=subprocess.DEVNULL, cwd=work_dir,
                           timeout=3600, check=False)  # TRF exits non-zero on success
    except FileNotFoundError:
        logger.info("trf not found; skipping satellite detection")
        return set()
    except Exception as e:
        logger.warning(f"trf failed: {e}; skipping satellite detection")
        return set()
    cov: Dict[str, int] = {}
    cur = None
    if os.path.exists(out):
        for ln in open(out):
            if ln.startswith('@'):
                cur = ln[1:].split()[0]
                continue
            p = ln.split()
            if len(p) >= 15 and cur:
                try:
                    s, e = int(p[0]), int(p[1])
                except ValueError:
                    continue
                cov[cur] = cov.get(cur, 0) + (e - s + 1)
    frac = config.satellite_tandem_frac
    return {q for q in cov if seqlen.get(q, 0) > 0 and cov[q] / seqlen[q] >= frac}


def classify_nonte(records: List[Dict], config) -> Tuple[List[Dict], Dict[str, List[Dict]], Dict]:
    """Split records into TE-clean vs non-TE classes (organelle / rrna / satellite).

    Returns (te_records, nonte_by_class, stats). Each non-TE record gets a
    'nonte_class' annotation. TE-clean records are returned unchanged.
    """
    stats = {'input': len(records), 'enabled': bool(getattr(config, 'enable_nonte_separation', False))}
    if not records or not getattr(config, 'enable_nonte_separation', False):
        stats['skipped_reason'] = 'disabled' if records else 'no_records'
        return records, {}, stats

    work_dir = tempfile.mkdtemp(prefix='nonte_')
    try:
        fasta_path = os.path.join(work_dir, 'consensi.fa')
        _write_fasta(records, fasta_path)
        seqlen = {r['id']: len(r['sequence']) for r in records}

        organelle = _detect_organelle(fasta_path, config, work_dir)
        rrna = _detect_rrna(fasta_path, config, work_dir)
        satellite = _detect_satellite(fasta_path, config, work_dir, seqlen)

        # Priority: organelle > rrna > satellite
        cls: Dict[str, str] = {}
        for rid in organelle:
            cls[rid] = 'organelle'
        for rid in rrna:
            cls.setdefault(rid, 'rrna')
        for rid in satellite:
            cls.setdefault(rid, 'satellite')

        te_records, nonte_by_class = [], {'organelle': [], 'rrna': [], 'satellite': []}
        for rec in records:
            c = cls.get(rec['id'])
            if c:
                rec['nonte_class'] = c
                nonte_by_class[c].append(rec)
            else:
                te_records.append(rec)

        stats.update({
            'organelle': len(nonte_by_class['organelle']),
            'rrna': len(nonte_by_class['rrna']),
            'satellite': len(nonte_by_class['satellite']),
            'nonte_total': len(records) - len(te_records),
            'te_clean': len(te_records),
        })
        logger.info(f"Non-TE separation: {len(records)} -> TE-clean {len(te_records)} "
                    f"(organelle={stats['organelle']}, rrna={stats['rrna']}, "
                    f"satellite={stats['satellite']})")
        return te_records, nonte_by_class, stats
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
