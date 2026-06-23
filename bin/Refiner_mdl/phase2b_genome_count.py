"""
Phase 2b — genome-wide copy count + final recurrence floor (large-genome SAMPLED mode).

In sampled mode (genome > large-genome threshold) every upstream rmblastn — the copy
gate, Phase-1 per-family recruitment, recall, fallback — runs against a SAMPLE of the
genome (config.genome_file). Sample-based copy counts undercount real families, so the
copy gate DEFERS its hard floor (config.defer_copy_floor). This module restores it once,
at the end, on the small FINAL consensus set against the COMPLETE genome:

  1. read the final masking + analysis libraries (Phase 2 output) as records;
  2. ONE makeblastdb of the full genome + batched rmblastn of the consensi vs it
     (reusing phase2_copy_recruit.recruit_genomic_copies on a full-genome config) →
     genome-wide rec['genomic_copies'];
  3. apply the SAME filter_short_lowcopy (hard floor + genome-size-adaptive short-low-copy
     gate) used by the copy gate, but now on genome-wide counts;
  4. rescue dropped-but-TE-structured consensi (phase_rescue), identical to the copy gate;
  5. rewrite consensus_masking.fa + phase3_analysis_library.fa, preserving the original
     FASTA records verbatim (only dropping non-survivors), masking ⊆ analysis as before.

This is, deliberately, "run the copy gate again on the full genome after Phase 2". The one
genome-wide alignment in the whole pipeline lives here, on a few-thousand-sequence library
rather than on every family — the scaling win that makes sampling worthwhile.

Hard safety: if the genome-wide recruit fails for any reason (no DB / no matrix / rmblastn
error), genomic_copies is left unset and we KEEP every consensus (never drop on a failed
count → never a faked floor); the libraries are returned unchanged with a logged warning.
"""

import logging
import os
import shutil
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)


def _read_fasta_records(path: str) -> List[Dict]:
    """Parse a FASTA into ordered records preserving each entry's exact text.

    Each record: {'id', 'sequence' (concatenated, no newlines), 'raw' (verbatim
    header+sequence lines, newline-terminated)}. 'raw' is what gets re-emitted so the
    rewritten library is byte-identical to the input for every surviving sequence.
    """
    records: List[Dict] = []
    if not path or not os.path.exists(path):
        return records
    rid = None
    raw_lines: List[str] = []
    seq_parts: List[str] = []

    def _flush():
        if rid is not None:
            records.append({'id': rid, 'sequence': ''.join(seq_parts),
                            'raw': ''.join(raw_lines)})

    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                _flush()
                rid = line[1:].split()[0] if len(line) > 1 else ''
                raw_lines = [line if line.endswith('\n') else line + '\n']
                seq_parts = []
            else:
                raw_lines.append(line if line.endswith('\n') else line + '\n')
                seq_parts.append(line.strip())
    _flush()
    return records


def _ids_in(path: str) -> set:
    ids = set()
    if not path or not os.path.exists(path):
        return ids
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                ids.add(line[1:].split()[0] if len(line) > 1 else '')
    return ids


def _rewrite(path: str, records: List[Dict], keep_ids: set) -> int:
    """Rewrite `path` keeping only records whose id is in keep_ids, in original order,
    re-emitting each surviving record's verbatim text. Returns count written."""
    n = 0
    tmp = path + '.phase2b.tmp'
    with open(tmp, 'w') as out:
        for rec in records:
            if rec['id'] in keep_ids:
                out.write(rec['raw'])
                n += 1
    os.replace(tmp, path)
    return n


def run_genome_count(phase2_output: Dict, config) -> Dict:
    """Apply the genome-wide recurrence floor to the Phase 2 libraries. Returns an
    updated copy of phase2_output (paths unchanged; stats augmented with 'phase2b')."""
    from copy import copy as _shallow
    from phase2_copy_recruit import recruit_genomic_copies
    from phase2_lowcopy_filter import filter_short_lowcopy
    from phase_rescue import rescue_te_structured

    genome_full = getattr(config, 'genome_full_file', '') or ''
    masking_path = phase2_output.get('masking_library', '')
    analysis_path = phase2_output.get('analysis_library', '')

    out = dict(phase2_output)
    if not genome_full or not os.path.exists(genome_full):
        logger.info("phase2b genome count: no --genome-full set; skipping (not sampled mode)")
        return out
    if not analysis_path or not os.path.exists(analysis_path):
        logger.warning("phase2b genome count: analysis library missing (%s); skipping",
                       analysis_path)
        return out

    logger.info("=" * 60)
    logger.info("Phase 2b: genome-wide copy count + recurrence floor (sampled mode)")
    logger.info("  Full genome: %s (%.2f Gb)", genome_full,
                os.path.getsize(genome_full) / 1e9)
    logger.info("=" * 60)

    # Analysis library is the superset (masking ⊆ analysis); count once over it.
    records = _read_fasta_records(analysis_path)
    masking_ids = _ids_in(masking_path)
    if not records:
        logger.warning("phase2b genome count: analysis library empty; skipping")
        return out

    # Full-genome config: same parameters, but align against the COMPLETE genome with a
    # dedicated DB/temp dir, and never defer (this IS the deferred floor).
    cfg_full = _shallow(config)
    cfg_full.genome_file = genome_full
    cfg_full.genome_blast_db = ''
    cfg_full.temp_dir = os.path.join(config.temp_dir, 'phase2b_genome_count')
    cfg_full.defer_copy_floor = False
    cfg_full.enable_copy_recruit = True
    cfg_full.enable_lowcopy_filter = True
    os.makedirs(cfg_full.temp_dir, exist_ok=True)

    # 1) genome-wide genomic_copies (one makeblastdb of the full genome + batched rmblastn)
    recruit_stats = recruit_genomic_copies(records, cfg_full)
    if 'recruited' not in recruit_stats:
        # Recruit was skipped/failed → counts are absent. NEVER drop on an absent count.
        logger.error("phase2b genome count: genome-wide recruit FAILED (%s) — keeping ALL "
                     "%d consensi (no floor applied; libraries unchanged)",
                     recruit_stats.get('skipped_reason', 'unknown'), len(records))
        out.setdefault('stats', {})
        out['phase2b'] = {'status': 'recruit_failed', 'recruit': recruit_stats,
                          'kept': len(records)}
        return out

    # 2) recurrence floor on genome-wide counts (hard floor + short-low-copy joint gate)
    kept, dropped, lc_stats = filter_short_lowcopy(records, cfg_full)

    # 3) rescue dropped-but-TE-structured consensi (same engine as the copy gate)
    rescued, rescue_stats = rescue_te_structured(dropped, config)

    survivor_ids = {r['id'] for r in kept} | {r['id'] for r in rescued}

    # 4) rewrite both libraries, masking ⊆ analysis preserved
    analysis_n = _rewrite(analysis_path, records, survivor_ids)
    masking_n = _rewrite(masking_path, records, survivor_ids & masking_ids) \
        if masking_path and os.path.exists(masking_path) else 0

    if not config.keep_temp:
        shutil.rmtree(cfg_full.temp_dir, ignore_errors=True)

    phase2b_stats = {
        'status': 'ok', 'input': len(records),
        'kept_by_floor': len(kept), 'dropped': len(dropped), 'rescued': len(rescued),
        'survivors': len(survivor_ids),
        'analysis_after': analysis_n, 'masking_after': masking_n,
        'recruit': recruit_stats, 'lowcopy': lc_stats, 'rescue': rescue_stats,
    }
    out['phase2b'] = phase2b_stats
    logger.info("Phase 2b: %d consensi -> %d survivors (floor-kept %d + rescued %d of "
                "%d dropped); analysis %d, masking %d",
                len(records), len(survivor_ids), len(kept), len(rescued),
                len(dropped), analysis_n, masking_n)
    return out
