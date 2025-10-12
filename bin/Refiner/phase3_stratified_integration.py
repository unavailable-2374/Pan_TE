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

    # File 1: High quality only (for phase4)
    high_quality_file = output_dir / "high_quality_consensus.fasta"

    try:
        with open(high_quality_file, 'w') as f:
            for seq in high_quality:
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

        logger.info(f"✓ Saved HIGH quality consensus (phase4) to: {high_quality_file}")
        logger.info(f"    {len(high_quality)} sequences")

    except Exception as e:
        logger.error(f"Failed to save high quality file: {e}")

    # File 2: High + Medium quality (for Combine)
    combined_file = output_dir / "phase3_analysis_library.fa"
    combined_sequences = high_quality + medium_quality

    try:
        with open(combined_file, 'w') as f:
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

        logger.info(f"✓ Saved HIGH+MEDIUM quality analysis library (Combine) to: {combined_file}")
        logger.info(f"    {len(combined_sequences)} sequences (High: {len(high_quality)}, Medium: {len(medium_quality)})")

    except Exception as e:
        logger.error(f"Failed to save combined file: {e}")

    # Also save detailed statistics
    stats_file = output_dir / "phase3_stratified_statistics.txt"

    try:
        with open(stats_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("Phase 3 Stratified Biological Filtering Statistics\n")
            f.write("=" * 80 + "\n\n")

            f.write("Overall Results:\n")
            f.write(f"  Total input sequences: {stats.get('total_input', len(analysis_library))}\n")
            f.write(f"  High quality: {len(high_quality)} ({len(high_quality)/len(analysis_library)*100:.1f}%)\n")
            f.write(f"  Medium quality: {len(medium_quality)} ({len(medium_quality)/len(analysis_library)*100:.1f}%)\n")
            f.write(f"  Low quality: {len(low_quality)} ({len(low_quality)/len(analysis_library)*100:.1f}%)\n")
            f.write(f"  Rejected: {stats.get('rejected_consensus', 0)}\n\n")

            f.write("Output Files:\n")
            f.write(f"  1. high_quality_consensus.fasta: {len(high_quality)} sequences (for phase4)\n")
            f.write(f"  2. phase3_analysis_library.fa: {len(combined_sequences)} sequences (for Combine)\n\n")

            f.write("Filtering Strategy:\n")
            f.write("  Tier 1 (≥5 copies): Standard quality requirements\n")
            f.write("  Tier 2 (4-5 copies): Requires TSD + quality≥0.75 + validation≥0.70\n")
            f.write("  Tier 3 (2-3 copies): Requires TSD + quality≥0.80 + validation≥0.80\n\n")

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
