#!/usr/bin/env python3
"""
Stratified Biological Filtering Module
Evidence-based tiered filtering with copy number-adjusted thresholds
"""

import re
from typing import List, Dict, Tuple
from collections import Counter
import logging

logger = logging.getLogger(__name__)


def calculate_n_content(sequence: str) -> float:
    """Calculate N content percentage"""
    if not sequence:
        return 0.0
    n_count = sequence.upper().count('N')
    return (n_count / len(sequence) * 100)


def calculate_sequence_complexity(sequence: str) -> float:
    """Calculate sequence complexity (entropy)"""
    if not sequence:
        return 0.0

    base_counts = Counter(sequence.upper())
    total = sum(base_counts.values())

    if total == 0:
        return 0.0

    entropy = 0.0
    for count in base_counts.values():
        if count > 0:
            freq = count / total
            import math
            entropy -= freq * math.log2(freq)

    # Normalize to 0-1 range (max entropy for DNA is 2 bits)
    return entropy / 2.0


def parse_genome_statistics(genome_file: str, te_elements: int = None, te_coverage: int = None,
                            genome_size: int = None) -> Dict:
    """
    Parse genome statistics
    If te_elements, te_coverage, or genome_size not provided, estimate from file
    """
    import os

    genome_stats = {
        'genome_size': genome_size or 0,
        'te_coverage': te_coverage or 0,
        'te_elements': te_elements or 0,
        'avg_te_length': 0,
        'te_percentage': 0
    }

    # Estimate genome size from file if not provided
    if not genome_stats['genome_size'] and os.path.exists(genome_file):
        try:
            genome_stats['genome_size'] = os.path.getsize(genome_file)
        except:
            genome_stats['genome_size'] = 100000000  # Default 100Mb

    # Calculate derived statistics
    if genome_stats['te_elements'] > 0:
        if genome_stats['te_coverage'] > 0:
            genome_stats['avg_te_length'] = genome_stats['te_coverage'] / genome_stats['te_elements']
        if genome_stats['genome_size'] > 0:
            genome_stats['te_percentage'] = (genome_stats['te_coverage'] / genome_stats['genome_size']) * 100

    return genome_stats


def determine_stratified_thresholds(genome_stats: Dict, num_families: int,
                                   copy_numbers: list = None) -> Dict:
    """
    Determine biologically-informed DYNAMIC thresholds
    STRATIFIED: Different requirements based on copy number

    Strategy: Hybrid approach using median + expected copies + biological minimums
    """
    import numpy as np

    thresholds = {}

    # Expected copies per family from genome statistics
    if genome_stats['te_elements'] > 0 and num_families > 0:
        expected_copies_per_family = genome_stats['te_elements'] / num_families
    else:
        expected_copies_per_family = 20.0

    # Calculate copy number distribution statistics if provided
    if copy_numbers:
        median_copies = np.median(copy_numbers)
        mean_copies = np.mean(copy_numbers)
        q25 = np.percentile(copy_numbers, 25)
        q75 = np.percentile(copy_numbers, 75)
    else:
        # Use expected copies as fallback
        median_copies = expected_copies_per_family * 0.5
        mean_copies = expected_copies_per_family
        q25 = expected_copies_per_family * 0.2
        q75 = expected_copies_per_family * 1.2

    # DYNAMIC THRESHOLDS using Hybrid Strategy
    # Strategy: median-based with expected copies constraint and biological minimums

    # High quality threshold
    # max(5, min(median * 1.2, expected * 0.4))
    # - Uses median * 1.2 as baseline (20% above median)
    # - Constrained by expected * 0.4 (genome-informed)
    # - Minimum 5 (biological floor)
    thresholds['min_copies_high_grade'] = max(5, min(
        int(median_copies * 1.2),
        int(expected_copies_per_family * 0.4)
    ))

    # Medium quality threshold
    # max(3, median * 0.4)
    # - 40% of median copies
    # - Minimum 3 (biological floor)
    thresholds['min_copies_confident'] = max(3, int(median_copies * 0.4))

    # Low quality threshold
    # max(2, median * 0.15)
    # - 15% of median copies
    # - Minimum 2 (biological floor)
    thresholds['min_copies_low'] = max(2, int(median_copies * 0.15))

    # Quality score thresholds for medium and low confidence sequences
    thresholds['medium_confidence_min_quality'] = 0.75
    thresholds['medium_confidence_min_validation'] = 0.70
    thresholds['low_confidence_min_quality'] = 0.80
    thresholds['low_confidence_min_validation'] = 0.80

    # Quality score thresholds - biological standards
    thresholds['min_quality_high'] = 0.80
    thresholds['min_quality_medium'] = 0.70
    thresholds['min_quality_low'] = 0.60

    # Base validation threshold
    thresholds['min_validation'] = 0.50

    # Length thresholds - based on known TE biology
    thresholds['min_length'] = 50
    thresholds['typical_min_length'] = 100
    thresholds['max_length'] = 20000

    # Quality metrics thresholds
    thresholds['max_n_content'] = 15.0
    thresholds['min_complexity'] = 0.40

    # Log threshold calculation details
    logger.info(f"Dynamic copy number thresholds calculated:")
    logger.info(f"  Genome: {genome_stats['te_elements']:,} TE elements")
    logger.info(f"  Library: {num_families} families")
    logger.info(f"  Expected copies/family: {expected_copies_per_family:.1f}")
    if copy_numbers:
        logger.info(f"  Median copies: {median_copies:.1f}")
        logger.info(f"  Copy range (Q25-Q75): {q25:.1f} - {q75:.1f}")
    logger.info(f"  → High grade: ≥{thresholds['min_copies_high_grade']} copies")
    logger.info(f"  → Medium grade: ≥{thresholds['min_copies_confident']} copies")
    logger.info(f"  → Low grade: ≥{thresholds['min_copies_low']} copies")

    return thresholds


def apply_stratified_filtering(consensus_list: List[Dict], genome_stats: Dict) -> Tuple[List[Dict], List[Dict], List[Dict], Dict]:
    """
    Apply stratified biological filtering with FIXED biological thresholds

    Strategy:
    - Fixed biological thresholds for accept/reject decisions
    - Dynamic thresholds for quality grading (high/medium/low)

    Biological thresholds (FIXED):
      ≥10 copies:  High confidence TE, ACCEPT UNCONDITIONALLY (no quality requirement)
      6-9 copies:  Medium confidence, require quality≥0.70 + validation≥0.60
      4-5 copies:  Require TSD + quality≥0.75 + validation≥0.70
      2-3 copies:  Require TSD + quality≥0.80 + validation≥0.80
      <2 copies:   Reject

    Returns:
        (filtered_high, filtered_medium, filtered_low, rejection_stats)
    """
    if not consensus_list:
        return [], [], [], {}

    # Fixed biological thresholds (for accept/reject decisions)
    BIOLOGICAL_HIGH_CONFIDENCE = 10  # ≥10 copies: high confidence TE
    BIOLOGICAL_MEDIUM_CONFIDENCE = 6  # 6-9 copies: medium confidence

    # Extract copy numbers for dynamic threshold calculation
    copy_numbers = [consensus.get('copy_number', 0) for consensus in consensus_list if consensus.get('copy_number', 0) > 0]

    # Determine DYNAMIC thresholds (for quality grading only)
    thresholds = determine_stratified_thresholds(genome_stats, len(consensus_list), copy_numbers)

    # Filter results
    filtered_high = []
    filtered_medium = []
    filtered_low = []
    rejected = []
    rejection_reasons = Counter()

    logger.info(f"Using fixed biological thresholds:")
    logger.info(f"  ≥{BIOLOGICAL_HIGH_CONFIDENCE} copies: ACCEPT UNCONDITIONALLY (no quality requirement)")
    logger.info(f"  {BIOLOGICAL_MEDIUM_CONFIDENCE}-{BIOLOGICAL_HIGH_CONFIDENCE-1} copies: Medium confidence (quality≥0.70)")
    logger.info(f"  4-5 copies: TSD + quality≥0.75 + validation≥0.70")
    logger.info(f"  2-3 copies: TSD + quality≥0.80 + validation≥0.80")

    for consensus in consensus_list:
        # Extract metadata
        sequence = consensus.get('sequence', '')
        final_score = consensus.get('final_quality_score', consensus.get('quality_score', 0))
        validation = consensus.get('validation_score', 0)
        copies = consensus.get('copy_number', 0)
        length = len(sequence)
        tsd = consensus.get('tsd', 'N/A')

        # Calculate additional quality metrics
        n_content = calculate_n_content(sequence)
        complexity = calculate_sequence_complexity(sequence)

        # Rejection criteria (applies to all tiers)
        reject = False
        reasons = []

        # 1. Absolute length limits
        if length < thresholds['min_length']:
            reject = True
            reasons.append(f"too_short({length}<{thresholds['min_length']})")
        elif length > thresholds['max_length']:
            reject = True
            reasons.append(f"too_long({length}>{thresholds['max_length']})")

        # 2. Sequence quality
        if n_content > thresholds['max_n_content']:
            reject = True
            reasons.append(f"high_N({n_content:.1f}%)")

        if complexity < thresholds['min_complexity']:
            reject = True
            reasons.append(f"low_complexity({complexity:.2f})")

        # 3. Minimum validation requirement (only for low copy number)
        if copies < BIOLOGICAL_HIGH_CONFIDENCE and validation < thresholds['min_validation']:
            reject = True
            reasons.append(f"low_validation({validation:.2f})")

        if reject:
            rejected.append(consensus)
            for reason in reasons:
                rejection_reasons[reason] += 1
            continue

        # Extract TSD information
        has_tsd = (tsd not in ['N/A', 'NNNN', None] and not str(tsd).startswith('NNNN'))

        # STRATIFIED FILTERING WITH FIXED BIOLOGICAL THRESHOLDS

        # Tier 1: ≥10 copies (HIGH CONFIDENCE TE)
        # ACCEPT UNCONDITIONALLY - no quality requirements
        # ≥10 independent copies is strong biological evidence
        if copies >= BIOLOGICAL_HIGH_CONFIDENCE:
            # Use DYNAMIC thresholds for quality grading only
            if copies >= thresholds['min_copies_high_grade'] and final_score >= 0.80:
                consensus['quality_grade'] = 'high'
                filtered_high.append(consensus)
            elif copies >= thresholds['min_copies_confident'] and final_score >= 0.70:
                consensus['quality_grade'] = 'medium'
                filtered_medium.append(consensus)
            else:
                # Accept even with low quality score
                consensus['quality_grade'] = 'low'
                filtered_low.append(consensus)

        # Tier 2: 6-9 copies (MEDIUM CONFIDENCE)
        # Require good quality evidence
        elif copies >= BIOLOGICAL_MEDIUM_CONFIDENCE:
            if final_score >= 0.70 and validation >= 0.60:
                consensus['quality_grade'] = 'medium'
                filtered_medium.append(consensus)
            else:
                rejected.append(consensus)
                reasons = []
                if final_score < 0.70:
                    reasons.append(f'low_quality({final_score:.2f})')
                if validation < 0.60:
                    reasons.append(f'low_validation({validation:.2f})')
                rejection_reasons[f'6-9_copies_insufficient_evidence:' + ','.join(reasons)] += 1

        # Tier 3: 4-5 copies (REQUIRE STRONG EVIDENCE)
        # Must have TSD + quality≥0.75 + validation≥0.70
        elif copies in [4, 5]:
            if (has_tsd and
                final_score >= 0.75 and
                validation >= 0.70):
                # Passed strict requirements for 4-5 copies
                if final_score >= 0.80:
                    consensus['quality_grade'] = 'high'
                    filtered_high.append(consensus)
                else:
                    consensus['quality_grade'] = 'medium'
                    filtered_medium.append(consensus)
            else:
                rejected.append(consensus)
                reasons = []
                if not has_tsd:
                    reasons.append('no_TSD')
                if final_score < 0.75:
                    reasons.append(f'low_quality({final_score:.2f})')
                if validation < 0.70:
                    reasons.append(f'low_validation({validation:.2f})')
                rejection_reasons[f'4-5_copies_insufficient_evidence:' + ','.join(reasons)] += 1

        # Tier 4: 2-3 copies (REQUIRE EXTREME EVIDENCE)
        # Must have TSD + quality≥0.80 + validation≥0.80
        elif copies in [2, 3]:
            if (has_tsd and
                final_score >= 0.80 and
                validation >= 0.80):
                # Passed strict requirements for 2-3 copies
                consensus['quality_grade'] = 'high'
                filtered_high.append(consensus)
            else:
                rejected.append(consensus)
                reasons = []
                if not has_tsd:
                    reasons.append('no_TSD')
                if final_score < 0.80:
                    reasons.append(f'low_quality({final_score:.2f})')
                if validation < 0.80:
                    reasons.append(f'low_validation({validation:.2f})')
                rejection_reasons[f'2-3_copies_insufficient_evidence:' + ','.join(reasons)] += 1

        # <2 copies: definitely reject
        else:
            rejected.append(consensus)
            rejection_reasons[f'too_few_copies({copies}<2)'] += 1

    # Statistics
    total = len(consensus_list)
    stats = {
        'total_input': total,
        'high_quality': len(filtered_high),
        'medium_quality': len(filtered_medium),
        'low_quality': len(filtered_low),
        'rejected': len(rejected),
        'rejection_reasons': dict(rejection_reasons)
    }

    logger.info(f"Stratified filtering complete:")
    logger.info(f"  Total input: {total}")
    logger.info(f"  High quality: {len(filtered_high)} ({len(filtered_high)/total*100:.1f}%)")
    logger.info(f"  Medium quality: {len(filtered_medium)} ({len(filtered_medium)/total*100:.1f}%)")
    logger.info(f"  Low quality: {len(filtered_low)} ({len(filtered_low)/total*100:.1f}%)")
    logger.info(f"  Rejected: {len(rejected)} ({len(rejected)/total*100:.1f}%)")

    return filtered_high, filtered_medium, filtered_low, stats
