#!/usr/bin/env python3
"""
Stratified Biological Filtering Integration for Phase 3
Replacement functions for phase3_consensus_building.py
"""

import logging
from typing import List, Dict
from pathlib import Path

try:
    # Try relative import first (when used as module)
    from .stratified_biological_filter import (
        apply_stratified_filtering,
        parse_genome_statistics
    )
except ImportError:
    # Fall back to direct import (when used standalone)
    from stratified_biological_filter import (
        apply_stratified_filtering,
        parse_genome_statistics
    )

logger = logging.getLogger(__name__)


def _filter_and_grade_consensus_stratified(self, consensus_list: List[Dict], stats: Dict) -> List[Dict]:
    """
    Enhanced filtering using stratified biological standards
    Replaces the original _filter_and_grade_consensus function

    Returns all passing sequences (high + medium + low), but also stores
    quality grades for later output splitting
    """
    if not consensus_list:
        return []

    logger.info(f"Starting stratified biological filtering on {len(consensus_list)} sequences")

    # Parse genome statistics
    genome_stats = parse_genome_statistics(
        self.config.genome_file,
        te_elements=stats.get('total_te_elements'),
        te_coverage=stats.get('total_te_coverage'),
        genome_size=stats.get('genome_size')
    )

    # Apply stratified filtering
    filtered_high, filtered_medium, filtered_low, filter_stats = apply_stratified_filtering(
        consensus_list, genome_stats
    )

    # Update stats
    stats['high_quality_consensus'] = len(filtered_high)
    stats['medium_quality_consensus'] = len(filtered_medium)
    stats['low_quality_consensus'] = len(filtered_low)
    stats['rejected_consensus'] = filter_stats['rejected']
    stats['stratified_filtering_stats'] = filter_stats

    # Return all passing sequences (will be split in output)
    all_passing = filtered_high + filtered_medium + filtered_low

    logger.info(f"Stratified filtering complete:")
    logger.info(f"  High quality: {len(filtered_high)}")
    logger.info(f"  Medium quality: {len(filtered_medium)}")
    logger.info(f"  Low quality: {len(filtered_low)}")
    logger.info(f"  Total passing: {len(all_passing)}")
    logger.info(f"  Rejected: {filter_stats['rejected']}")

    return all_passing


def _save_stratified_output_files(self, analysis_library: List[Dict], stats: Dict):
    """
    Save two output files with stratified biological filtering results:
    1. high_quality_consensus.fasta - For phase4 (HIGH quality only)
    2. phase3_analysis_library.fa - For Combine (HIGH + MEDIUM quality)
    """
    from pathlib import Path

    output_dir = Path(self.config.output_dir)

    # Separate by quality grade
    high_quality = [seq for seq in analysis_library if seq.get('quality_grade') == 'high']
    medium_quality = [seq for seq in analysis_library if seq.get('quality_grade') == 'medium']
    low_quality = [seq for seq in analysis_library if seq.get('quality_grade') == 'low']

    # Get sequences that passed basic filtering (from stats)
    passed_basic = stats.get('stratified_filtering_stats', {}).get('passed_basic_filtering_sequences', [])

    # File 1: High quality only (for phase4)
    # If no high quality sequences, use medium quality as fallback
    high_quality_file = output_dir / "high_quality_consensus.fasta"

    # Select sequences for phase4 masking
    masking_sequences = high_quality if high_quality else medium_quality

    if not masking_sequences:
        logger.warning("⚠ No high or medium quality sequences available for phase4 masking!")
        logger.warning("  Using low quality sequences as last resort")
        masking_sequences = low_quality

    try:
        with open(high_quality_file, 'w') as f:
            if not masking_sequences:
                # Create empty file with comment if no sequences at all
                f.write("# No sequences passed quality filtering\n")
                logger.error("✗ NO sequences passed filtering! Creating empty placeholder file.")
            else:
                for seq in masking_sequences:
                    seq_id = seq.get('id', 'unknown')
                    sequence = seq.get('sequence', '')

                    # Build header with metadata
                    header_parts = [
                        f">{seq_id}",
                        f"grade={seq.get('quality_grade', 'unknown')}",
                        f"final_score={seq.get('final_quality_score', 0):.3f}",
                        f"base_score={seq.get('quality_score', 0):.3f}",
                        f"validation={seq.get('validation_score', 0):.3f}",
                        f"copies={seq.get('copy_number', 0)}",
                        f"length={len(sequence)}"
                    ]

                    if seq.get('tsd'):
                        header_parts.append(f"tsd={seq['tsd']}")

                    header = " ".join(header_parts)
                    f.write(f"{header}\n{sequence}\n")

        if masking_sequences:
            if masking_sequences == high_quality:
                logger.info(f"✓ Saved HIGH quality consensus (phase4) to: {high_quality_file}")
                logger.info(f"    {len(high_quality)} sequences")
            elif masking_sequences == medium_quality:
                logger.warning(f"⚠ No HIGH quality sequences, using MEDIUM quality for phase4")
                logger.info(f"✓ Saved MEDIUM quality consensus (phase4) to: {high_quality_file}")
                logger.info(f"    {len(medium_quality)} sequences")
            else:
                logger.warning(f"⚠ No HIGH/MEDIUM quality sequences, using LOW quality for phase4")
                logger.info(f"✓ Saved LOW quality consensus (phase4) to: {high_quality_file}")
                logger.info(f"    {len(low_quality)} sequences")

    except Exception as e:
        logger.error(f"Failed to save high quality file: {e}")

    # File 2: High + Medium quality (for Combine)
    # If no high+medium, include low quality as fallback
    combined_file = output_dir / "phase3_analysis_library.fa"
    combined_sequences = high_quality + medium_quality

    if not combined_sequences:
        logger.warning("⚠ No high+medium quality sequences, including LOW quality in analysis library")
        combined_sequences = low_quality

    try:
        with open(combined_file, 'w') as f:
            if not combined_sequences:
                # Create empty file with comment if no sequences at all
                f.write("# No sequences passed quality filtering\n")
                logger.error("✗ NO sequences passed filtering! Analysis library is empty.")
            else:
                for seq in combined_sequences:
                    seq_id = seq.get('id', 'unknown')
                    sequence = seq.get('sequence', '')

                    # Build header with metadata
                    header_parts = [
                        f">{seq_id}",
                        f"grade={seq.get('quality_grade', 'unknown')}",
                        f"final_score={seq.get('final_quality_score', 0):.3f}",
                        f"base_score={seq.get('quality_score', 0):.3f}",
                        f"validation={seq.get('validation_score', 0):.3f}",
                        f"copies={seq.get('copy_number', 0)}",
                        f"length={len(sequence)}"
                    ]

                    if seq.get('tsd'):
                        header_parts.append(f"tsd={seq['tsd']}")

                    header = " ".join(header_parts)
                    f.write(f"{header}\n{sequence}\n")

        if combined_sequences:
            if high_quality or medium_quality:
                logger.info(f"✓ Saved HIGH+MEDIUM quality analysis library (Combine) to: {combined_file}")
                logger.info(f"    {len(combined_sequences)} sequences (High: {len(high_quality)}, Medium: {len(medium_quality)})")
            else:
                logger.warning(f"✓ Saved LOW quality analysis library (Combine) to: {combined_file}")
                logger.info(f"    {len(combined_sequences)} sequences (all LOW quality)")

    except Exception as e:
        logger.error(f"Failed to save combined file: {e}")

    # File 3: All sequences that passed basic quality checks (NEW)
    basic_filtered_file = output_dir / "phase3_post_basic_filtering.fa"

    try:
        if passed_basic:
            with open(basic_filtered_file, 'w') as f:
                for seq in passed_basic:
                    seq_id = seq.get('id', 'unknown')
                    sequence = seq.get('sequence', '')

                    # Build header with metadata
                    header_parts = [
                        f">{seq_id}",
                        f"grade={seq.get('quality_grade', 'rejected')}",
                        f"final_score={seq.get('final_quality_score', 0):.3f}",
                        f"base_score={seq.get('quality_score', 0):.3f}",
                        f"validation={seq.get('validation_score', 0):.3f}",
                        f"copies={seq.get('copy_number', 0)}",
                        f"length={len(sequence)}"
                    ]

                    if seq.get('tsd'):
                        header_parts.append(f"tsd={seq['tsd']}")

                    header = " ".join(header_parts)
                    f.write(f"{header}\n{sequence}\n")

            logger.info(f"✓ Saved post-basic-filtering sequences to: {basic_filtered_file}")
            logger.info(f"    {len(passed_basic)} sequences (passed length/N-content/complexity checks)")

            # Count how many were accepted vs rejected after validation
            accepted_count = len(high_quality) + len(medium_quality) + len(low_quality)
            rejected_after_basic = len(passed_basic) - accepted_count
            logger.info(f"    → {accepted_count} accepted by stratified filtering")
            logger.info(f"    → {rejected_after_basic} rejected by validation/copy-number filtering")
        else:
            logger.warning(f"⚠ No sequences passed basic filtering!")

    except Exception as e:
        logger.error(f"Failed to save basic filtered file: {e}")

    # File 4: ALL sequences (including rejected) - for analysis
    all_sequences_file = output_dir / "phase3_all_sequences_with_grades.fa"

    try:
        with open(all_sequences_file, 'w') as f:
            # Write all sequences from analysis_library (these have grades assigned)
            all_seqs_written = set()

            for seq in analysis_library:
                seq_id = seq.get('id', 'unknown')
                sequence = seq.get('sequence', '')
                all_seqs_written.add(seq_id)

                # Build header with metadata
                header_parts = [
                    f">{seq_id}",
                    f"status=ACCEPTED",
                    f"grade={seq.get('quality_grade', 'unknown')}",
                    f"final_score={seq.get('final_quality_score', 0):.3f}",
                    f"base_score={seq.get('quality_score', 0):.3f}",
                    f"validation={seq.get('validation_score', 0):.3f}",
                    f"copies={seq.get('copy_number', 0)}",
                    f"length={len(sequence)}"
                ]

                if seq.get('tsd'):
                    header_parts.append(f"tsd={seq['tsd']}")

                header = " ".join(header_parts)
                f.write(f"{header}\n{sequence}\n")

            # Write rejected sequences from passed_basic (these passed basic but failed tier filtering)
            if passed_basic:
                for seq in passed_basic:
                    seq_id = seq.get('id', 'unknown')
                    if seq_id in all_seqs_written:
                        continue  # Already written as accepted

                    sequence = seq.get('sequence', '')
                    all_seqs_written.add(seq_id)

                    # This sequence was rejected by tier filtering
                    header_parts = [
                        f">{seq_id}",
                        f"status=REJECTED_BY_TIER_FILTERING",
                        f"grade=rejected",
                        f"final_score={seq.get('final_quality_score', 0):.3f}",
                        f"base_score={seq.get('quality_score', 0):.3f}",
                        f"validation={seq.get('validation_score', 0):.3f}",
                        f"copies={seq.get('copy_number', 0)}",
                        f"length={len(sequence)}"
                    ]

                    if seq.get('tsd'):
                        header_parts.append(f"tsd={seq['tsd']}")

                    header = " ".join(header_parts)
                    f.write(f"{header}\n{sequence}\n")

        total_written = len(all_seqs_written)
        logger.info(f"✓ Saved ALL sequences (accepted + rejected) to: {all_sequences_file}")
        logger.info(f"    {total_written} total sequences")
        logger.info(f"    {len(analysis_library)} accepted, {total_written - len(analysis_library)} rejected")

    except Exception as e:
        logger.error(f"Failed to save all sequences file: {e}")

    # Also save detailed statistics
    stats_file = output_dir / "phase3_stratified_statistics.txt"

    try:
        with open(stats_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("Phase 3 Stratified Biological Filtering Statistics\n")
            f.write("=" * 80 + "\n\n")

            filter_stats = stats.get('stratified_filtering_stats', {})
            total_input = filter_stats.get('total_input', len(analysis_library))
            passed_basic_count = filter_stats.get('passed_basic_filtering', 0)
            rejected_by_validation = filter_stats.get('rejected_by_validation', 0)
            rejected_by_basic = filter_stats.get('rejected', 0)

            f.write("Overall Results:\n")
            f.write(f"  Total input sequences: {total_input}\n")
            f.write(f"  Passed basic quality checks: {passed_basic_count} ({passed_basic_count/total_input*100:.1f}%)\n")
            f.write(f"    - Length, N-content, complexity checks\n")
            f.write(f"  Rejected by basic quality: {rejected_by_basic} ({rejected_by_basic/total_input*100:.1f}%)\n\n")

            f.write(f"  Final acceptance (after validation/copy filtering):\n")
            f.write(f"    High quality: {len(high_quality)} ({len(high_quality)/total_input*100:.1f}%)\n")
            f.write(f"    Medium quality: {len(medium_quality)} ({len(medium_quality)/total_input*100:.1f}%)\n")
            f.write(f"    Low quality: {len(low_quality)} ({len(low_quality)/total_input*100:.1f}%)\n")
            f.write(f"    Total accepted: {len(high_quality)+len(medium_quality)+len(low_quality)} ({(len(high_quality)+len(medium_quality)+len(low_quality))/total_input*100:.1f}%)\n")
            f.write(f"  Rejected by validation/copies: {rejected_by_validation} ({rejected_by_validation/total_input*100:.1f}%)\n\n")

            f.write("Output Files:\n")
            f.write(f"  1. high_quality_consensus.fasta: {len(high_quality)} sequences (for phase4)\n")
            f.write(f"  2. phase3_analysis_library.fa: {len(combined_sequences)} sequences (for Combine)\n")
            f.write(f"  3. phase3_post_basic_filtering.fa: {passed_basic_count} sequences (all passed basic quality)\n\n")

            f.write("Filtering Strategy (RELAXED THRESHOLDS - Fixed 99.95% rejection issue):\n")
            f.write("  Standardized abundance: I = (C × (L/L₀)) / G_Gb\n")
            f.write("  Fixed boundaries in I space: b₁=5, b₂=25\n")
            f.write("  Adaptive copy thresholds calculated from genome size and family length\n")
            f.write("  Structure score S: comprehensive stct.md v1.0 scoring with enhanced hybrid mode\n\n")
            f.write("  UPDATED Thresholds (Previous caused 99.95% rejection):\n")
            f.write("  Tier A (I≥b₂, high copy):   S ≥ 0.30 (was 0.55)\n")
            f.write("  Tier B (b₁≤I<b₂, medium):   S ≥ 0.40 (was 0.70)\n")
            f.write("  Tier C (I<b₁, low copy):    S ≥ 0.50 (was 0.85)\n\n")
            f.write("  Extreme cases:\n")
            f.write("    - I ≥ 100: Allow S ≥ 0.25 (was 0.45)\n")
            f.write("    - I < 1 and S < 0.55: Reject (was 0.90)\n\n")
            f.write("  Enhancements:\n")
            f.write("    - Selective HMM: enabled for sequences >1000bp\n")
            f.write("    - Enhanced hybrid scoring: 60% quality + 40% structure (when structural evidence limited)\n")
            f.write("    - Relaxed complexity threshold: 0.35 (was 0.40)\n\n")
            f.write("  Quality grading within tiers:\n")
            f.write("    - Tier A → HIGH quality\n")
            f.write("    - Tier B → HIGH/MEDIUM (based on S score)\n")
            f.write("    - Tier C → MEDIUM/LOW (based on S score)\n\n")

            # Rejection reasons if available
            if 'stratified_filtering_stats' in stats:
                filter_stats = stats['stratified_filtering_stats']
                if 'rejection_reasons' in filter_stats:
                    f.write("Top Rejection Reasons:\n")
                    rejection_reasons = filter_stats['rejection_reasons']
                    sorted_reasons = sorted(rejection_reasons.items(), key=lambda x: x[1], reverse=True)
                    for reason, count in sorted_reasons[:15]:
                        f.write(f"  {reason}: {count}\n")

            f.write("\n" + "=" * 80 + "\n")

        logger.info(f"✓ Saved stratified statistics to: {stats_file}")

    except Exception as e:
        logger.error(f"Failed to save statistics file: {e}")


# Wrapper function that maintains backward compatibility
def _save_enhanced_analysis_library_fasta_stratified(self, analysis_library: List[Dict], stats: Dict):
    """Wrapper that calls the stratified output function"""
    _save_stratified_output_files(self, analysis_library, stats)


# Function to patch the ConsensusBuilder class
def patch_consensus_builder_with_stratified_filtering(ConsensusBuilder):
    """
    Patch ConsensusBuilder class to use stratified biological filtering

    Usage:
        from phase3_stratified_integration import patch_consensus_builder_with_stratified_filtering
        patch_consensus_builder_with_stratified_filtering(ConsensusBuilder)
    """
    # Replace filtering function
    ConsensusBuilder._filter_and_grade_consensus = _filter_and_grade_consensus_stratified

    # Replace output function
    ConsensusBuilder._save_enhanced_analysis_library_fasta = _save_enhanced_analysis_library_fasta_stratified
    ConsensusBuilder._save_analysis_library_fasta = _save_enhanced_analysis_library_fasta_stratified

    logger.info("ConsensusBuilder patched with stratified biological filtering")
