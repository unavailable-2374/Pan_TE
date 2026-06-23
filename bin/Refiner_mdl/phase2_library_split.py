"""
Phase 2: Final QC + library splitting.

Produces two FASTA libraries:
  - consensus_masking.fa       (T1+T2, length-dependent copies) -> genome masking -> RepGraph
  - phase3_analysis_library.fa (masking + qualified T3)          -> Combine classification

Includes TE structure verification: sequences without any TE signal
(protein domain, TIR, or known TE homology) are excluded.
"""

import logging
import json
import os
from typing import Dict, List
from pathlib import Path

logger = logging.getLogger(__name__)

# ── Reuse Refiner utilities ──────────────────────────────────────────
import sys
_refiner_dir = os.path.join(os.path.dirname(__file__), '..', 'Refiner')
if _refiner_dir not in sys.path:
    sys.path.insert(0, _refiner_dir)
from utils.complexity_utils import (calculate_shannon_entropy, calculate_dust_score,
                                     calculate_lowcomplexity_fraction)

from .te_structure_filter import filter_by_te_structure, annotate_te_structure
from .phase2_gene_exclusion import exclude_gene_derived
from .phase2_nonte_separation import classify_nonte


def _format_header(rec: Dict) -> str:
    """Format FASTA header with metadata."""
    tier = rec.get('tier', 'T2')
    copies = rec.get('copies', 0)
    mdl = rec.get('mdl', 0)
    slen = len(rec.get('sequence', ''))
    return (f">{rec['id']} length={slen} copies={copies} "
            f"mdl={mdl:.1f} tier={tier}")


def run_phase2(phase1_output: Dict, config) -> Dict:
    """Execute Phase 2: QC + library splitting.

    Returns dict with keys:
        masking_library: str   — path to consensus_masking.fa
        analysis_library: str  — path to phase3_analysis_library.fa
        statistics: str        — path to stats JSON
        stats: Dict
    """
    records = phase1_output['records']
    logger.info("=" * 60)
    logger.info(f"Phase 2: Final QC + Library Split ({len(records)} sequences)")
    logger.info("=" * 60)

    stats = {'input_count': len(records), 'qc_failed': 0,
             'too_short': 0, 'low_entropy': 0, 'high_n': 0, 'ssr_dust': 0,
             'low_entropy2': 0, 'high_lc_fraction': 0}

    # ── QC Filter ────────────────────────────────────────────────────
    passed = []
    for rec in records:
        seq = rec.get('sequence', '')
        slen = len(seq)

        if slen < config.min_length:
            stats['too_short'] += 1
            stats['qc_failed'] += 1
            continue

        ent = calculate_shannon_entropy(seq, k=1)
        thresh = (config.entropy_threshold_long if slen >= 100
                  else config.entropy_threshold_short)
        if ent < thresh:
            stats['low_entropy'] += 1
            stats['qc_failed'] += 1
            continue

        # DUST score: SSR/simple-repeat filter (consistent with Phase 0)
        dust = calculate_dust_score(seq)
        if dust > config.dust_score_threshold:
            stats['ssr_dust'] += 1
            stats['qc_failed'] += 1
            continue

        # k=2 (dinucleotide) entropy: dinucleotide-repeat low-complexity past the k=1 gate
        ent2 = calculate_shannon_entropy(seq, k=2)
        thresh2 = (config.entropy2_threshold_long if slen >= 100
                   else config.entropy2_threshold_short)
        if ent2 < thresh2:
            stats['low_entropy2'] += 1
            stats['qc_failed'] += 1
            continue

        # Low-complexity fraction: mostly-simple consensi that slip the average DUST gate
        if calculate_lowcomplexity_fraction(seq, dust_cut=config.dust_score_threshold) \
                > config.lowcomplexity_fraction_max:
            stats['high_lc_fraction'] += 1
            stats['qc_failed'] += 1
            continue

        n_frac = seq.count('N') / slen if slen > 0 else 0
        if n_frac > config.max_n_percent:
            stats['high_n'] += 1
            stats['qc_failed'] += 1
            continue

        passed.append(rec)

    logger.info(f"QC: {len(records)} -> {len(passed)} "
                f"(failed: short={stats['too_short']}, "
                f"entropy={stats['low_entropy']}, entropy2={stats['low_entropy2']}, "
                f"SSR/dust={stats['ssr_dust']}, LC-frac={stats['high_lc_fraction']}, "
                f"N%={stats['high_n']})")

    # ── Output dir ───────────────────────────────────────────────────
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # NOTE: genomic copy recruitment + recurrence hard-floor + TE-structure rescue now
    # run in the copy gate BEFORE Phase 1 (phase05_copy_gate), so the floor applies to
    # whole-family copy counts and Phase 1 only refines survivors. Phase 2 no longer
    # re-filters on copies.

    # ── TE Structure: TAG by default (do NOT drop) ─────────────────────
    # Structure / homology cannot decide true-vs-false repeat — divergent and
    # dark-matter TEs lack detectable structure, so a hard gate sheds real elements.
    # Default: annotate each record with its TE signals (sidecar TSV) and keep every
    # QC-passing repeat; real-vs-noise dropping is done upstream by Phase 1 N5
    # (recurrence/coherence, structure-free); class is deferred to downstream
    # RepeatClassifier / ClassifyTE / TE-looker. The legacy hard gate is opt-in.
    if getattr(config, 'enable_te_structure_gate', False):
        pre_filter_count = len(passed)
        passed, te_stats = filter_by_te_structure(passed, config)
        stats['te_structure_filter'] = te_stats
        stats['no_te_signal'] = pre_filter_count - len(passed)
        # Survivors all carry >=1 TE signal; protect them all from gene exclusion.
        signal_map = {r['id']: ['te_gate_passed'] for r in passed}
        logger.info(f"TE structure HARD GATE (opt-in): {pre_filter_count} -> {len(passed)} "
                    f"(removed {pre_filter_count - len(passed)} with no TE signal)")
    else:
        signal_map, te_stats = annotate_te_structure(passed, config)
        stats['te_structure_annotation'] = te_stats
        stats['no_te_signal'] = 0  # tag-only: nothing dropped on structure
        for rec in passed:
            rec['te_signal'] = signal_map.get(rec['id'], [])
        annot_path = str(output_dir / 'te_signal_annotation.tsv')
        with open(annot_path, 'w') as afh:
            afh.write("id\tte_signals\n")
            for rec in passed:
                sigs = rec.get('te_signal', [])
                afh.write(f"{rec['id']}\t{','.join(sigs) if sigs else 'none'}\n")
        stats['te_signal_annotation_path'] = annot_path
        logger.info(f"TE structure TAG-ONLY (default): {len(passed)} kept; "
                    f"{te_stats['any_te_signal']} carry >=1 signal "
                    f"(annotation -> {annot_path}); no repeat dropped on structure")

    # ── mdl-repeat NATIVE quality gate ──────────────────────────────────
    # mdl-repeat scores every family (header tier=core/warn/reject, accept=…) but
    # Refiner_mdl otherwise ignores it, passing the low-confidence bulk (warn / standalone /
    # few copies) straight through and inflating the library with redundant fragments.
    # Keep a family iff it is mdl-high-confidence (tier in keep_tiers) OR carries independent
    # TE structural signal OR has high copy support; drop only the no-signal ∩ low-tier ∩
    # low-copy intersection. SAFETY: skip the gate entirely when no structural signal could
    # be computed (missing RepeatPeps/RepeatMasker libs) — never drop without the rescue,
    # and never gate records whose header lacks tier= (older mdl-repeat builds).
    if getattr(config, 'enable_mdl_quality_gate', True):
        annotation_available = any(signal_map.get(r['id']) for r in passed)
        n_tiered = sum(1 for r in passed if r.get('mdl_tier') is not None)
        if not annotation_available:
            logger.warning("mdl quality gate SKIPPED: no TE structural signal available "
                            "(RepeatPeps/RepeatMasker libs missing?) — refusing to gate "
                            "without the signal rescue")
        elif n_tiered == 0:
            logger.info("mdl quality gate SKIPPED: no mdl tier= in headers (older build)")
        else:
            keep_tiers = {t.strip() for t in str(getattr(config, 'mdl_quality_keep_tiers', 'core')).split(',') if t.strip()}
            hc = int(getattr(config, 'mdl_quality_highcopy_keep', 50))
            kept, dropped_q = [], []
            for rec in passed:
                mt = rec.get('mdl_tier')
                if (mt is None or mt in keep_tiers
                        or signal_map.get(rec['id']) or rec.get('copies', 0) >= hc):
                    kept.append(rec)
                else:
                    dropped_q.append(rec)
            stats['mdl_quality_gate'] = {
                'input': len(passed), 'kept': len(kept), 'dropped': len(dropped_q),
                'keep_tiers': sorted(keep_tiers), 'highcopy_keep': hc,
            }
            if dropped_q:
                qdrop_path = str(output_dir / 'mdl_quality_dropped.tsv')
                with open(qdrop_path, 'w') as qfh:
                    qfh.write("id\tmdl_tier\tmdl_accept\tcopies\tlength\n")
                    for rec in dropped_q:
                        qfh.write(f"{rec['id']}\t{rec.get('mdl_tier')}\t{rec.get('mdl_accept')}\t"
                                  f"{rec.get('copies')}\t{rec.get('length')}\n")
                stats['mdl_quality_gate']['dropped_tsv'] = qdrop_path
            logger.info(f"mdl quality gate: {len(passed)} -> {len(kept)} "
                        f"(dropped {len(dropped_q)} low-confidence no-signal low-copy)")
            passed = kept

    # ── Protein-homology gene exclusion (negative filter) ──────────────
    # Drop consensi that strongly match a cellular protein (external TE-removed DB)
    # AND carry zero TE evidence (te_signal empty). Any TE signal protects the seq;
    # DB absent -> skipped. Removes multi-copy host-gene noise that over-masks CDS.
    passed, gene_dropped, gene_excl_stats = exclude_gene_derived(passed, signal_map, config)
    stats['gene_exclusion'] = gene_excl_stats
    if gene_dropped:
        drop_path = str(output_dir / 'gene_exclusion_dropped.tsv')
        with open(drop_path, 'w') as dfh:
            dfh.write("id\tcoverage\tpident\tmatched_protein\n")
            for rec in gene_dropped:
                ge = rec.get('gene_exclusion', {})
                dfh.write(f"{rec['id']}\t{ge.get('coverage','')}\t"
                          f"{ge.get('pident','')}\t{ge.get('protein','')}\n")
        stats['gene_exclusion']['dropped_tsv'] = drop_path
        logger.info(f"Gene exclusion: dropped {len(gene_dropped)} gene-derived consensi "
                    f"(provenance -> {drop_path})")

    # ── Protein-homology flank trim (cut terminal host-gene over-extension) ──
    # Complements Phase-1 occupancy trim: removes a terminal gene flank that is shared
    # across copies (so occupancy keeps it) but matches a host gene (partial cov band).
    from .phase2_flank_trim import trim_gene_flanks
    passed, flank_stats = trim_gene_flanks(passed, config)
    stats['flank_trim'] = flank_stats

    # ── Non-TE repeat separation (route organelle / rRNA / satellite out) ──
    # These are recurrent but not transposons; keep them out of the TE library
    # (written to dedicated nonte_*.fa tracks) so TE content is not overstated.
    passed, nonte_by_class, nonte_stats = classify_nonte(passed, config)
    stats['nonte_separation'] = nonte_stats
    for cls_name, recs in (nonte_by_class or {}).items():
        if not recs:
            continue
        track_path = str(output_dir / f'nonte_{cls_name}.fa')
        with open(track_path, 'w') as tfh:
            for rec in recs:
                tfh.write(f"{_format_header(rec)}\n{rec['sequence']}\n")
        logger.info(f"Non-TE track '{cls_name}': {len(recs)} consensi -> {track_path}")

    masking_path = str(output_dir / 'consensus_masking.fa')
    analysis_path = str(output_dir / 'phase3_analysis_library.fa')

    masking_count = 0
    analysis_count = 0

    with open(masking_path, 'w') as mask_fh, \
         open(analysis_path, 'w') as anal_fh:

        for rec in passed:
            tier = rec.get('tier', 'T2')
            header = _format_header(rec)
            seq = rec['sequence']
            slen = len(seq)
            copies = rec.get('copies', 0)
            entry = f"{header}\n{seq}\n"

            # Masking library: T1+T2 with length-dependent copies threshold
            # Longer seqs have higher specificity → need fewer copies.
            in_masking = False
            if tier in ('T1', 'T2'):
                if slen >= 500:
                    min_c = config.masking_copies_long
                elif slen >= 200:
                    min_c = config.masking_copies_medium
                elif slen >= 100:
                    min_c = config.masking_copies_short
                else:
                    min_c = config.masking_copies_vshort
                if copies >= min_c:
                    in_masking = True

            if in_masking:
                mask_fh.write(entry)
                masking_count += 1

            # Analysis library: masking seqs + T3 with sufficient copies
            if in_masking or (tier == 'T3' and copies >= config.analysis_t3_min_copies):
                anal_fh.write(entry)
                analysis_count += 1

    stats['masking_count'] = masking_count
    stats['analysis_count'] = analysis_count
    stats['qc_passed'] = len(passed)
    stats['qc_excluded'] = len(passed) - analysis_count

    # Count tiers in final output
    tier_counts = {}
    for rec in passed:
        t = rec.get('tier', 'unknown')
        tier_counts[t] = tier_counts.get(t, 0) + 1
    stats['tiers'] = tier_counts

    # Calculate total library sizes
    masking_size = os.path.getsize(masking_path) if os.path.exists(masking_path) else 0
    analysis_size = os.path.getsize(analysis_path) if os.path.exists(analysis_path) else 0
    stats['masking_library_bytes'] = masking_size
    stats['analysis_library_bytes'] = analysis_size

    # ── Statistics Report ────────────────────────────────────────────
    stats_path = str(output_dir / 'refiner_mdl_stats.json')

    # Merge all phase stats
    all_stats = {
        'phase0': phase1_output.get('phase0_stats', {}),
        'phase1': phase1_output.get('stats', {}),
        'phase2': stats,
    }
    with open(stats_path, 'w') as fh:
        json.dump(all_stats, fh, indent=2, default=str)

    logger.info(f"Phase 2 complete:")
    logger.info(f"  Masking library:  {masking_count} sequences -> {masking_path}")
    logger.info(f"  Analysis library: {analysis_count} sequences -> {analysis_path}")
    logger.info(f"  Excluded (QC+quality): {len(passed) - analysis_count} sequences")
    logger.info(f"  Tiers: {tier_counts}")
    logger.info(f"  Statistics: {stats_path}")

    return {
        'masking_library': masking_path,
        'analysis_library': analysis_path,
        'statistics': stats_path,
        'stats': stats,
        'summary': (f"{analysis_count} consensus sequences "
                    f"(masking: {masking_count}, T1={tier_counts.get('T1', 0)}, "
                    f"T2={tier_counts.get('T2', 0)}, T3={tier_counts.get('T3', 0)})"),
    }
