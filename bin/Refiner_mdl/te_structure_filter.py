"""
TE structure verification filter.

Keeps sequences that show ANY evidence of TE origin:
  1. TE protein domain (blastx vs RepeatPeps.lib)
  2. Terminal Inverted Repeats (self reverse-complement alignment)
  3. Long Terminal Repeats (self direct-repeat alignment)
  4. Known TE homology (blastn vs RepeatMasker.lib)

Sequences with zero TE signals are excluded from both libraries.
"""

import logging
import os
import subprocess
import tempfile
import shutil
from typing import Dict, List, Set, Tuple

logger = logging.getLogger(__name__)

# Default paths (RepeatMasker Libraries)
_RM_LIBS = os.path.join(
    os.environ.get('CONDA_PREFIX', ''),
    'share', 'RepeatMasker', 'Libraries')


def _find_library(config, name: str) -> str:
    """Locate a RepeatMasker library file."""
    # Check config override first
    path = getattr(config, f'{name}_path', '')
    if path and os.path.isfile(path):
        return path
    # Default location
    path = os.path.join(_RM_LIBS, name)
    if os.path.isfile(path):
        return path
    return ''


def _build_self_blastdb(fasta_path: str, config, work_dir: str) -> str:
    """Build a nucleotide BLAST DB from the query sequences.

    Returns db_path on success, empty string on failure.
    """
    db_path = os.path.join(work_dir, 'self_db')
    cmd = [
        config.makeblastdb_exe,
        '-in', fasta_path,
        '-dbtype', 'nucl',
        '-out', db_path,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            logger.warning(f"makeblastdb for self-alignment failed: {result.stderr[:200]}")
            return ''
    except Exception as e:
        logger.warning(f"makeblastdb for self-alignment failed: {e}")
        return ''
    return db_path


def _detect_protein_domains(fasta_path: str, config, work_dir: str) -> Set[str]:
    """Detect TE protein domains via blastx vs RepeatPeps.lib.

    Returns set of sequence IDs with at least one protein domain hit.
    """
    peps_db = _find_library(config, 'RepeatPeps.lib')
    if not peps_db:
        logger.warning("RepeatPeps.lib not found; skipping protein domain check")
        return set()

    out_path = os.path.join(work_dir, 'blastx_peps.out')
    cmd = [
        getattr(config, 'blastx_exe', 'blastx'),
        '-query', fasta_path,
        '-db', peps_db,
        '-outfmt', '6 qseqid',
        '-evalue', '1e-5',
        '-max_target_seqs', '1',
        '-max_hsps', '1',
        '-num_threads', str(config.threads),
        '-seg', 'yes',
        '-out', out_path,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        if result.returncode != 0:
            logger.warning(f"blastx failed: {result.stderr[:200]}")
            return set()
    except subprocess.TimeoutExpired:
        logger.warning("blastx timed out (2h)")
        return set()
    except FileNotFoundError:
        logger.warning("blastx not found; skipping protein domain check")
        return set()

    hits = set()
    if os.path.exists(out_path):
        with open(out_path) as fh:
            for line in fh:
                qid = line.strip().split('\t')[0]
                if qid:
                    hits.add(qid)
    return hits


def _detect_terminal_repeats(fasta_path: str, db_path: str, config,
                             work_dir: str, strand: str, label: str,
                             min_repeat_len: int = 10,
                             max_repeat_fraction: float = 0.5) -> Set[str]:
    """Detect terminal repeats via self-alignment BLASTN.

    Args:
        strand: 'minus' for TIR (inverted repeats), 'plus' for LTR (direct repeats)
        label: name for logging and output files (e.g. 'tir', 'ltr')
        min_repeat_len: minimum alignment length to count as a repeat
        max_repeat_fraction: maximum fraction of sequence length for a repeat unit.
            Filters out trivial full-length self-hits on plus strand.

    Looks for repeats at sequence termini (within first/last 20% of the
    sequence). Returns set of IDs with terminal repeats detected.
    """
    if not db_path:
        return set()

    out_path = os.path.join(work_dir, f'{label}_self.out')
    cmd = [
        config.blastn_exe,
        '-query', fasta_path,
        '-db', db_path,
        '-outfmt', '6 qseqid sseqid qstart qend qlen sstart send slen pident length',
        '-evalue', '1e-3',
        '-max_target_seqs', '1',
        '-max_hsps', '5',
        '-strand', strand,
        '-dust', 'no',
        '-num_threads', str(config.threads),
        '-out', out_path,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            logger.warning(f"{label.upper()} BLASTN failed: {result.stderr[:200]}")
            return set()
    except subprocess.TimeoutExpired:
        logger.warning(f"{label.upper()} BLASTN timed out")
        return set()

    hits = set()
    if os.path.exists(out_path):
        with open(out_path) as fh:
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                qid, sid = parts[0], parts[1]
                if qid != sid:
                    continue  # Self-hit only
                qstart, qend = int(parts[2]), int(parts[3])
                qlen = int(parts[4])
                sstart, send = int(parts[5]), int(parts[6])
                aln_len = int(parts[9])

                if aln_len < min_repeat_len:
                    continue

                # Filter trivial full-length self-hits (especially on plus strand).
                # A real terminal repeat can't exceed half the sequence length.
                if aln_len > qlen * max_repeat_fraction:
                    continue

                # Check if hit is at opposite termini
                # One end in first 20%, other end in last 20%
                boundary = max(1, int(qlen * 0.2))
                q_at_start = min(qstart, qend) <= boundary
                q_at_end = max(qstart, qend) >= qlen - boundary
                s_at_start = min(sstart, send) <= boundary
                s_at_end = max(sstart, send) >= qlen - boundary

                if (q_at_start and s_at_end) or (q_at_end and s_at_start):
                    hits.add(qid)

    return hits


def _detect_known_te_homology(fasta_path: str, config, work_dir: str) -> Set[str]:
    """Detect homology to known TEs via blastn vs RepeatMasker.lib.

    Returns set of IDs with significant hits to known TE consensus sequences.
    """
    te_db = _find_library(config, 'RepeatMasker.lib')
    if not te_db:
        logger.warning("RepeatMasker.lib not found; skipping TE homology check")
        return set()

    out_path = os.path.join(work_dir, 'blastn_telib.out')
    cmd = [
        config.blastn_exe,
        '-task', 'blastn',
        '-word_size', '9',
        '-query', fasta_path,
        '-db', te_db,
        '-outfmt', '6 qseqid pident length qlen',
        '-evalue', '1e-5',
        '-max_target_seqs', '1',
        '-max_hsps', '1',
        '-dust', 'no',
        '-num_threads', str(config.threads),
        '-out', out_path,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            logger.warning(f"TE homology BLASTN failed: {result.stderr[:200]}")
            return set()
    except subprocess.TimeoutExpired:
        logger.warning("TE homology BLASTN timed out")
        return set()

    hits = set()
    if os.path.exists(out_path):
        with open(out_path) as fh:
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                qid = parts[0]
                aln_len = int(parts[2])
                # Minimum meaningful match: 30bp
                if aln_len >= 30:
                    hits.add(qid)
    return hits


def annotate_te_structure(records: List[Dict], config) -> Tuple[Dict[str, List[str]], Dict]:
    """Run the four TE-structure checks and RETURN per-record signals — never drop.

    Structure / homology signals CANNOT decide whether a sequence is a true repeat:
    divergent, decayed, novel, and "dark-matter" TEs (the very targets of refine +
    TE-looker) carry no detectable protein domain / TIR / LTR / known-TE homology, so
    using these as a hard gate systematically discards real elements (false negatives).
    This function therefore only TAGS: each record id maps to the list of structural
    signals it does show (possibly empty). Dropping by recurrence/coherence happens
    upstream (Phase 1 N5, structure-free); true/false-TE class is deferred downstream
    (RepeatClassifier / ClassifyTE / TE-looker).

    Returns (signal_map, stats) where signal_map[id] = sorted list of any of
    {'protein_domain','tir','ltr','known_te_homology'}.
    """
    if not records:
        return {}, {}

    work_dir = tempfile.mkdtemp(prefix='te_annot_')
    try:
        fasta_path = os.path.join(work_dir, 'query.fa')
        with open(fasta_path, 'w') as fh:
            for rec in records:
                fh.write(f">{rec['id']}\n{rec['sequence']}\n")

        logger.info(f"TE structure annotation (tag-only, no drop): {len(records)} sequences")
        db_path = _build_self_blastdb(fasta_path, config, work_dir)

        protein_hits = _detect_protein_domains(fasta_path, config, work_dir)
        tir_hits = _detect_terminal_repeats(fasta_path, db_path, config, work_dir,
                                            strand='minus', label='tir', min_repeat_len=10)
        ltr_hits = _detect_terminal_repeats(fasta_path, db_path, config, work_dir,
                                            strand='plus', label='ltr', min_repeat_len=50)
        te_homology_hits = _detect_known_te_homology(fasta_path, config, work_dir)

        signal_map: Dict[str, List[str]] = {}
        for rec in records:
            rid = rec['id']
            sigs = []
            if rid in protein_hits:
                sigs.append('protein_domain')
            if rid in tir_hits:
                sigs.append('tir')
            if rid in ltr_hits:
                sigs.append('ltr')
            if rid in te_homology_hits:
                sigs.append('known_te_homology')
            signal_map[rid] = sigs

        n_with = sum(1 for s in signal_map.values() if s)
        stats = {
            'input': len(records),
            'protein_domain': len(protein_hits),
            'tir': len(tir_hits),
            'ltr': len(ltr_hits),
            'known_te_homology': len(te_homology_hits),
            'any_te_signal': n_with,
            'no_te_signal': len(records) - n_with,
        }
        logger.info(f"  Tagged with >=1 TE signal: {n_with}/{len(records)} "
                    f"({n_with/len(records)*100:.1f}%) — NONE dropped on this basis")
        return signal_map, stats
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def filter_by_te_structure(records: List[Dict], config) -> Tuple[List[Dict], Dict]:
    """Filter sequences: keep only those with at least one TE structural signal.

    Runs four independent checks (protein domain, TIR, LTR, known TE homology)
    and keeps any sequence that passes at least one.

    NOTE: structure-based filtering has systematic false negatives on divergent /
    dark-matter TEs and is DISABLED by default (Phase 2 uses annotate_te_structure
    for tagging instead). Retained behind config.enable_te_structure_gate for
    explicit opt-in / backward comparison only.

    Returns (passed_records, stats_dict).
    """
    if not records:
        return records, {}

    work_dir = tempfile.mkdtemp(prefix='te_struct_')

    try:
        # Write all sequences to FASTA
        fasta_path = os.path.join(work_dir, 'query.fa')
        with open(fasta_path, 'w') as fh:
            for rec in records:
                fh.write(f">{rec['id']}\n{rec['sequence']}\n")

        logger.info(f"TE structure verification: {len(records)} sequences")

        # Build self-BLAST DB once, reuse for TIR and LTR detection
        db_path = _build_self_blastdb(fasta_path, config, work_dir)

        # Run four checks
        protein_hits = _detect_protein_domains(fasta_path, config, work_dir)
        logger.info(f"  Protein domains: {len(protein_hits)} sequences")

        tir_hits = _detect_terminal_repeats(
            fasta_path, db_path, config, work_dir,
            strand='minus', label='tir', min_repeat_len=10)
        logger.info(f"  TIR detected: {len(tir_hits)} sequences")

        ltr_hits = _detect_terminal_repeats(
            fasta_path, db_path, config, work_dir,
            strand='plus', label='ltr', min_repeat_len=50)
        logger.info(f"  LTR detected: {len(ltr_hits)} sequences")

        te_homology_hits = _detect_known_te_homology(fasta_path, config, work_dir)
        logger.info(f"  Known TE homology: {len(te_homology_hits)} sequences")

        # Union of all hits
        all_te_ids = protein_hits | tir_hits | ltr_hits | te_homology_hits
        logger.info(f"  Total with TE signal: {len(all_te_ids)} / {len(records)} "
                     f"({len(all_te_ids)/len(records)*100:.1f}%)")

        passed = [r for r in records if r['id'] in all_te_ids]

        stats = {
            'input': len(records),
            'protein_domain': len(protein_hits),
            'tir': len(tir_hits),
            'ltr': len(ltr_hits),
            'known_te_homology': len(te_homology_hits),
            'any_te_signal': len(all_te_ids),
            'no_te_signal': len(records) - len(all_te_ids),
        }

        logger.info(f"  Removed (no TE signal): {stats['no_te_signal']} sequences")

        return passed, stats

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
