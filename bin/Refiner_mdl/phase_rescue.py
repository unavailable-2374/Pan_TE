"""
Copy-gate rescue: recover low-copy consensi that carry TE structural evidence.

The recurrence hard-floor (>= hard_min_copies genomic copies) is a blunt instrument:
genuinely-TE families that happen to be low-copy (young low-copy insertions, or families
mdl split fine) would be dropped with the segmental-dup / artifact noise. This module
gives every dropped sequence a SECOND CHANCE on positive TE evidence only — it can rescue
but never drop, so it cannot create the divergent-TE false-negatives that a structural
GATE would (a no-evidence sequence stays dropped on copy count, not on structure).

Evidence (rescue if >= 1 fires):
  1. TE protein domain via HMM — TEsorter / REXdb hmmscan (GAG/PROT/RT/RH/INT/TPase ...),
     far more sensitive to diverged TE domains than blastx; the primary signal and the
     one robust to the ragged boundaries of an un-refined mdl consensus.
  2. LTR / TIR terminal structure + known-TE homology — reused from annotate_te_structure
     (self-alignment + nhmmer/blastn). Boundary-dependent, so weaker pre-refinement, but
     additive.

Skips cleanly (rescues nothing, drops nothing) if TEsorter is unavailable — the structural
signals still run.
"""

import logging
import os
import subprocess
import tempfile
import shutil
from typing import Dict, List, Set, Tuple

logger = logging.getLogger(__name__)


def _tesorter_domains(fasta_path: str, config, work_dir: str) -> Set[str]:
    """Run TEsorter (REXdb HMM) and return ids classified with >=1 TE protein domain."""
    exe = getattr(config, 'tesorter_exe', 'TEsorter')
    db = getattr(config, 'tesorter_db', 'rexdb')
    pre = os.path.join(work_dir, 'ts')
    cmd = [exe, fasta_path, '-db', db, '-p', str(config.threads), '-pre', pre]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=getattr(config, 'tesorter_timeout_s', 14400), cwd=work_dir)
        if r.returncode != 0:
            logger.warning(f"TEsorter rc={r.returncode}: {r.stderr[-200:]}; HMM rescue signal off")
            return set()
    except FileNotFoundError:
        logger.warning("TEsorter not found; HMM rescue signal off (structural signals still run)")
        return set()
    except subprocess.TimeoutExpired:
        logger.warning("TEsorter timed out; HMM rescue signal off")
        return set()
    hits = set()
    cls = pre + '.cls.tsv'
    if os.path.exists(cls):
        for ln in open(cls):
            if ln.startswith('#') or not ln.strip():
                continue
            hits.add(ln.split('\t')[0])
    return hits


def rescue_te_structured(dropped: List[Dict], config) -> Tuple[List[Dict], Dict]:
    """Return (rescued_records, stats). Each rescued record gets a 'rescued_by' tag."""
    stats = {'input_dropped': len(dropped), 'rescued': 0}
    if not dropped or not getattr(config, 'enable_rescue', False):
        stats['skipped_reason'] = 'disabled' if dropped else 'no_records'
        return [], stats

    work_dir = tempfile.mkdtemp(prefix='rescue_')
    try:
        fasta_path = os.path.join(work_dir, 'dropped.fa')
        with open(fasta_path, 'w') as fh:
            for rec in dropped:
                fh.write(f">{rec['id']}\n{rec['sequence']}\n")

        # Run the independent detectors CONCURRENTLY (each already multi-threaded), so the
        # provided cores stay busy instead of serialising TEsorter then 3 BLAST/HMM passes.
        # The blastx-vs-RepeatPeps protein check is intentionally DROPPED — TEsorter's REXdb
        # HMM is the same signal (TE protein domains) but more sensitive and hmmscan -cpu
        # scales, so it fully subsumes it.
        from te_structure_filter import (_build_self_blastdb, _detect_terminal_repeats,
                                          _detect_known_te_homology)
        import concurrent.futures

        db_path = _build_self_blastdb(fasta_path, config, work_dir)  # needed for TIR/LTR

        def _hmm():
            return ('hmm', _tesorter_domains(fasta_path, config, work_dir))

        def _tir():
            return ('tir', _detect_terminal_repeats(fasta_path, db_path, config, work_dir,
                                                    strand='minus', label='tir', min_repeat_len=10))

        def _ltr():
            return ('ltr', _detect_terminal_repeats(fasta_path, db_path, config, work_dir,
                                                    strand='plus', label='ltr', min_repeat_len=50))

        def _known():
            return ('known', _detect_known_te_homology(fasta_path, config, work_dir))

        sig: Dict[str, Set[str]] = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
            for key, ids in ex.map(lambda f: f(), [_hmm, _tir, _ltr, _known]):
                sig[key] = ids

        rescued = []
        counts = {'hmm': 0, 'tir': 0, 'ltr': 0, 'known': 0}
        for rec in dropped:
            rid = rec['id']
            by = [k for k in ('hmm', 'tir', 'ltr', 'known') if rid in sig.get(k, set())]
            if by:
                for k in by:
                    counts[k] += 1
                rec['rescued_by'] = ','.join(by)
                rescued.append(rec)

        stats.update({'rescued': len(rescued), 'by': counts})
        logger.info(f"Rescue: {len(dropped)} dropped -> rescued {len(rescued)} "
                    f"(HMM TE-domain={counts['hmm']}, TIR={counts['tir']}, "
                    f"LTR={counts['ltr']}, known-TE={counts['known']})")
        return rescued, stats
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
