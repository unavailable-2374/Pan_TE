#!/usr/bin/env python3
"""
Stratified Biological Filtering Module
Evidence-based tiered filtering with copy number-adjusted thresholds

v2.0: Adaptive thresholds based on guide.md (Standardized Abundance)
"""

import re
import math
from typing import List, Dict, Tuple, Optional
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
        except (OSError, AttributeError):
            genome_stats['genome_size'] = 100000000  # Default 100Mb

    # Calculate derived statistics
    if genome_stats['te_elements'] > 0:
        if genome_stats['te_coverage'] > 0:
            genome_stats['avg_te_length'] = genome_stats['te_coverage'] / genome_stats['te_elements']
        if genome_stats['genome_size'] > 0:
            genome_stats['te_percentage'] = (genome_stats['te_coverage'] / genome_stats['genome_size']) * 100

    return genome_stats


def calculate_standardized_abundance(copy_number: int, length_kb: float,
                                     genome_size_gb: float, L0_kb: float = 3.0) -> float:
    """
    Calculate standardized abundance I (guide.md v1.0)

    I = (C × (L/L₀)) / G_Gb

    Args:
        copy_number: Number of copies in genome
        length_kb: Family typical length in kb
        genome_size_gb: Genome size in Gb
        L0_kb: Reference length (default 3.0 kb)

    Returns:
        Standardized abundance I
    """
    if genome_size_gb <= 0:
        return 0.0

    I = (copy_number * (length_kb / L0_kb)) / genome_size_gb
    return I


def calculate_adaptive_copy_thresholds(genome_size_gb: float, length_kb: float,
                                       L0_kb: float = 3.0, b1: float = 5.0,
                                       b2: float = 25.0) -> Dict[str, int]:
    """
    Calculate adaptive copy number thresholds based on genome size and family length
    (guide.md v1.0 formula)

    C₁ = ⌈b₁ × G_Gb × (L₀/L)⌉
    C₂ = ⌈b₂ × G_Gb × (L₀/L)⌉

    Args:
        genome_size_gb: Genome size in Gb
        length_kb: Family typical length in kb
        L0_kb: Reference length (default 3.0 kb)
        b1: Tier B/C boundary in I space (default 5.0)
        b2: Tier A/B boundary in I space (default 25.0)

    Returns:
        Dict with 'tier_c_boundary' (C₁) and 'tier_a_boundary' (C₂)
    """
    if genome_size_gb <= 0 or length_kb <= 0:
        # Fallback to minimal thresholds
        return {
            'tier_c_boundary': 2,
            'tier_a_boundary': 6
        }

    # Calculate adaptive thresholds
    C1 = max(2, math.ceil(b1 * genome_size_gb * (L0_kb / length_kb)))
    C2 = max(4, math.ceil(b2 * genome_size_gb * (L0_kb / length_kb)))

    # Ensure C2 > C1
    if C2 <= C1:
        C2 = C1 + 2

    return {
        'tier_c_boundary': C1,  # MEDIUM/LOW boundary
        'tier_a_boundary': C2   # HIGH/MEDIUM boundary
    }


def calculate_structure_score(consensus: Dict, has_tsd: bool = False,
                              use_stct_v1: bool = True) -> float:
    """
    Calculate structure score S ∈ [0,1]

    Args:
        consensus: Consensus sequence dict with quality metrics
        has_tsd: Whether TSD is detected
        use_stct_v1: If True, use stct.md v1.0 comprehensive scoring
                     If False, use simplified scoring (backward compatible)

    Returns:
        Structure score S in [0, 1]
    """
    if use_stct_v1:
        # Use comprehensive structure scoring (stct.md v1.0)
        try:
            from structure_scoring_v1 import calculate_structure_score_v1

            # Selective HMM: only run for long sequences (>1000bp) to balance accuracy and speed
            sequence_length = len(consensus.get('sequence', ''))
            run_hmm = (sequence_length > 1000)

            # Calculate comprehensive structure score
            result = calculate_structure_score_v1(
                consensus,
                genome_hits=consensus.get('rm_hits'),
                options={'run_hmm': run_hmm}
            )

            # Store additional metadata
            consensus['structure_score_details'] = result
            consensus['structure_tier'] = result['structure_tier']

            return result['structure_score']

        except ImportError:
            logger.warning("structure_scoring_v1 not available, falling back to simple scoring")
            # Fall through to simple scoring

    # Simplified scoring (backward compatible with guide.md)
    base_quality = consensus.get('quality_score', 0.0)
    validation = consensus.get('validation_score', 0.0)
    length = len(consensus.get('sequence', ''))

    # Component scores
    scores = []
    weights = []

    # 1. Base quality (alignment quality, consensus quality)
    scores.append(base_quality)
    weights.append(0.4)

    # 2. Validation score (biological validation)
    scores.append(validation)
    weights.append(0.3)

    # 3. TSD evidence (strong structural signal)
    tsd_score = 1.0 if has_tsd else 0.3
    scores.append(tsd_score)
    weights.append(0.2)

    # 4. Length completeness (typical TE: 100-10000bp)
    if length < 100:
        length_score = length / 100.0
    elif length > 10000:
        length_score = max(0.5, 1.0 - (length - 10000) / 10000.0)
    else:
        length_score = 1.0
    scores.append(length_score)
    weights.append(0.1)

    # Weighted average
    S = sum(s * w for s, w in zip(scores, weights)) / sum(weights)

    return max(0.0, min(1.0, S))


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
    Apply stratified biological filtering with DISTRIBUTION-BASED DYNAMIC TRUST

    VERSION v3.0: Distribution-aware trust levels instead of fixed thresholds

    Key insight: RepeatScout already filtered sequences (≥10 k-mer occurrences),
    so we trust the distribution and apply validation proportionally.

    Trust levels based on copy number/I-value percentiles:
    - ≥P75: High trust - accept with minimal requirements
    - P50-P75: Medium-high trust - accept with relaxed requirements
    - P25-P50: Medium trust - light validation
    - <P25: Lower trust - need validation or biological features

    This approach adapts to each dataset's characteristics automatically.

    Returns:
        (filtered_high, filtered_medium, filtered_low, rejection_stats)
    """
    import numpy as np

    if not consensus_list:
        return [], [], [], {}

    # ========================================================================
    # STEP 1: Calculate genome size and sequence characteristics
    # ========================================================================
    genome_size_bytes = genome_stats.get('genome_size', 0)
    genome_size_gb = genome_size_bytes / (1024**3) if genome_size_bytes > 0 else 0.1

    family_lengths_kb = [len(c.get('sequence', '')) / 1000.0 for c in consensus_list]
    if family_lengths_kb:
        median_length_kb = np.median(family_lengths_kb)
        typical_length_kb = max(0.5, min(10.0, median_length_kb))
    else:
        typical_length_kb = 3.0

    L0_kb = 3.0

    # ========================================================================
    # STEP 2: Calculate distribution-based trust thresholds (NEW)
    # ========================================================================
    # Extract copy numbers and calculate I values for all sequences
    copy_numbers = []
    I_values = []

    for c in consensus_list:
        copies = c.get('copy_number', 0)
        length_kb = len(c.get('sequence', '')) / 1000.0
        I = calculate_standardized_abundance(copies, length_kb, genome_size_gb, L0_kb)
        copy_numbers.append(copies)
        I_values.append(I)

    # Calculate percentiles for distribution-based trust
    if copy_numbers:
        copy_p25 = np.percentile(copy_numbers, 25)
        copy_p50 = np.percentile(copy_numbers, 50)
        copy_p75 = np.percentile(copy_numbers, 75)
        I_p25 = np.percentile(I_values, 25)
        I_p50 = np.percentile(I_values, 50)
        I_p75 = np.percentile(I_values, 75)
    else:
        copy_p25, copy_p50, copy_p75 = 2, 5, 15
        I_p25, I_p50, I_p75 = 1, 5, 25

    # Dynamic trust thresholds based on distribution
    # Use max(percentile, minimum) to handle skewed distributions
    HIGH_TRUST_COPY = max(copy_p75, 3)
    MEDIUM_HIGH_TRUST_COPY = max(copy_p50, 2)
    MEDIUM_TRUST_COPY = max(copy_p25, 1)

    HIGH_TRUST_I = max(I_p75, 5)
    MEDIUM_HIGH_TRUST_I = max(I_p50, 2)
    MEDIUM_TRUST_I = max(I_p25, 0.5)

    # Basic quality thresholds (relaxed)
    MIN_LENGTH = 50
    MAX_LENGTH = 20000
    MAX_N_CONTENT = 15.0
    MIN_COMPLEXITY = 0.30  # Further relaxed from 0.35

    # Structure score thresholds by trust level
    # Slightly tightened to improve precision by ~0.5%
    S_HIGH_TRUST = 0.18       # High trust: minimal requirement
    S_MEDIUM_HIGH_TRUST = 0.22
    S_MEDIUM_TRUST = 0.27
    S_LOW_TRUST = 0.33        # Low trust: require slightly more evidence

    # Filter results
    filtered_high = []
    filtered_medium = []
    filtered_low = []
    rejected = []
    passed_basic_filtering = []
    rejection_reasons = Counter()
    trust_stats = {'high': 0, 'medium_high': 0, 'medium': 0, 'low': 0}

    logger.info(f"========== DISTRIBUTION-BASED TRUST FILTERING (v3.0) ==========")
    logger.info(f"Genome: {genome_size_gb:.3f} Gb, Median family length: {typical_length_kb:.2f} kb")
    logger.info(f"Copy number distribution: P25={copy_p25:.1f}, P50={copy_p50:.1f}, P75={copy_p75:.1f}")
    logger.info(f"I-value distribution: P25={I_p25:.2f}, P50={I_p50:.2f}, P75={I_p75:.2f}")
    logger.info(f"Dynamic trust thresholds (copy): High≥{HIGH_TRUST_COPY:.1f}, MedHigh≥{MEDIUM_HIGH_TRUST_COPY:.1f}, Med≥{MEDIUM_TRUST_COPY:.1f}")
    logger.info(f"Structure score thresholds: High≥{S_HIGH_TRUST}, MedHigh≥{S_MEDIUM_HIGH_TRUST}, Med≥{S_MEDIUM_TRUST}, Low≥{S_LOW_TRUST}")
    logger.info(f"Basic quality: length=[{MIN_LENGTH},{MAX_LENGTH}], N<{MAX_N_CONTENT}%, complexity>{MIN_COMPLEXITY}")
    logger.info(f"================================================================")

    for consensus in consensus_list:
        # ========================================================================
        # STEP 3: Extract metadata and calculate metrics
        # ========================================================================
        sequence = consensus.get('sequence', '')
        copies = consensus.get('copy_number', 0)
        length = len(sequence)
        length_kb = length / 1000.0
        tsd = consensus.get('tsd', 'N/A')
        quality_class = consensus.get('quality_class', 'B')

        # Calculate basic quality metrics
        n_content = calculate_n_content(sequence)
        complexity = calculate_sequence_complexity(sequence)
        has_tsd = (tsd not in ['N/A', 'NNNN', None] and not str(tsd).startswith('NNNN'))

        # Calculate standardized abundance I
        I = calculate_standardized_abundance(copies, length_kb, genome_size_gb, L0_kb)
        consensus['standardized_abundance'] = I

        # Calculate structure score S
        S = calculate_structure_score(consensus, has_tsd)
        consensus['structure_score'] = S

        # ========================================================================
        # STEP 4: Basic quality checks (relaxed)
        # ========================================================================
        reject = False
        reasons = []

        if length < MIN_LENGTH:
            reject = True
            reasons.append(f"too_short({length}<{MIN_LENGTH})")
        elif length > MAX_LENGTH:
            reject = True
            reasons.append(f"too_long({length}>{MAX_LENGTH})")

        if n_content > MAX_N_CONTENT:
            reject = True
            reasons.append(f"high_N({n_content:.1f}%)")

        # Relaxed complexity check - allow override for long sequences
        if complexity < MIN_COMPLEXITY:
            if length < 500:  # Only reject short low-complexity sequences
                reject = True
                reasons.append(f"low_complexity({complexity:.2f})")

        if reject:
            rejected.append(consensus)
            for reason in reasons:
                rejection_reasons[reason] += 1
            continue

        passed_basic_filtering.append(consensus)

        # ========================================================================
        # STEP 5: Distribution-based trust level assignment (NEW)
        # ========================================================================
        # Determine trust level based on copy number percentile
        if copies >= HIGH_TRUST_COPY or I >= HIGH_TRUST_I:
            trust_level = 'high'
            required_S = S_HIGH_TRUST
            trust_stats['high'] += 1
        elif copies >= MEDIUM_HIGH_TRUST_COPY or I >= MEDIUM_HIGH_TRUST_I:
            trust_level = 'medium_high'
            required_S = S_MEDIUM_HIGH_TRUST
            trust_stats['medium_high'] += 1
        elif copies >= MEDIUM_TRUST_COPY or I >= MEDIUM_TRUST_I:
            trust_level = 'medium'
            required_S = S_MEDIUM_TRUST
            trust_stats['medium'] += 1
        else:
            trust_level = 'low'
            required_S = S_LOW_TRUST
            trust_stats['low'] += 1

        consensus['trust_level'] = trust_level

        # ========================================================================
        # STEP 6: Apply trust-aware filtering
        # ========================================================================
        # Check for biological evidence that can override low scores
        has_structural_features = _check_structural_features(sequence)
        is_long_sequence = length >= 1000
        is_rescued = quality_class == 'C_rescued'
        is_high_quality_class = quality_class == 'A'

        # Determine if sequence passes
        pass_te = False
        pass_reason = None

        if trust_level in ['high', 'medium_high']:
            # High/medium-high trust: almost always accept
            if S >= required_S:
                pass_te = True
                pass_reason = f"{trust_level}_trust_pass"
            elif has_structural_features or is_long_sequence or has_tsd:
                # Override with biological evidence
                pass_te = True
                pass_reason = f"{trust_level}_trust_bio_override"
            elif S >= 0.10:  # Very minimal threshold for high trust
                pass_te = True
                pass_reason = f"{trust_level}_trust_minimal"

        elif trust_level == 'medium':
            # Medium trust: accept if meets threshold or has evidence
            if S >= required_S:
                pass_te = True
                pass_reason = "medium_trust_pass"
            elif has_structural_features or has_tsd:
                pass_te = True
                pass_reason = "medium_trust_bio_override"
            elif is_long_sequence and S >= 0.15:
                pass_te = True
                pass_reason = "medium_trust_long_seq"
            elif is_rescued or is_high_quality_class:
                pass_te = True
                pass_reason = f"medium_trust_{quality_class}_override"

        else:  # low trust
            # Low trust: need threshold OR strong evidence
            if S >= required_S:
                pass_te = True
                pass_reason = "low_trust_pass"
            elif has_structural_features and S >= 0.15:
                pass_te = True
                pass_reason = "low_trust_structural"
            elif has_tsd and S >= 0.20:
                pass_te = True
                pass_reason = "low_trust_tsd"
            elif is_long_sequence and S >= 0.20:
                pass_te = True
                pass_reason = "low_trust_long_seq"
            elif is_high_quality_class:
                pass_te = True
                pass_reason = "low_trust_A_class"
            elif is_rescued and S >= 0.15:
                pass_te = True
                pass_reason = "low_trust_rescued"

        # ========================================================================
        # STEP 7: Assign quality grade and categorize
        # ========================================================================
        if pass_te:
            # Assign grade based on trust level and S score
            if trust_level == 'high' or (trust_level == 'medium_high' and S >= S_HIGH_TRUST):
                consensus['quality_grade'] = 'high'
                filtered_high.append(consensus)
            elif trust_level in ['medium_high', 'medium'] or S >= S_MEDIUM_TRUST:
                consensus['quality_grade'] = 'medium'
                filtered_medium.append(consensus)
            else:
                consensus['quality_grade'] = 'low'
                filtered_low.append(consensus)

            consensus['pass_reason'] = pass_reason
            logger.debug(f"{consensus.get('id', 'unknown')}: {trust_level} → {consensus['quality_grade'].upper()} "
                        f"(C={copies}, I={I:.2f}, S={S:.2f}, reason={pass_reason})")
        else:
            rejected.append(consensus)
            rejection_reasons[f'{trust_level}_trust_rejected(S={S:.2f}<{required_S})'] += 1
            logger.debug(f"{consensus.get('id', 'unknown')}: Rejected {trust_level} (C={copies}, I={I:.2f}, S={S:.2f})")

    # Statistics
    total = len(consensus_list)
    total_accepted = len(filtered_high) + len(filtered_medium) + len(filtered_low)
    stats = {
        'total_input': total,
        'high_quality': len(filtered_high),
        'medium_quality': len(filtered_medium),
        'low_quality': len(filtered_low),
        'rejected': len(rejected),
        'passed_basic_filtering': len(passed_basic_filtering),
        'rejected_by_validation': len(passed_basic_filtering) - total_accepted,
        'rejection_reasons': dict(rejection_reasons),
        'passed_basic_filtering_sequences': passed_basic_filtering,
        'trust_distribution': trust_stats
    }

    logger.info(f"=======================================================")
    logger.info(f"Distribution-based trust filtering complete:")
    logger.info(f"  Total input: {total}")
    logger.info(f"  Trust level distribution:")
    logger.info(f"    High trust: {trust_stats['high']} ({trust_stats['high']/total*100:.1f}%)")
    logger.info(f"    Medium-high trust: {trust_stats['medium_high']} ({trust_stats['medium_high']/total*100:.1f}%)")
    logger.info(f"    Medium trust: {trust_stats['medium']} ({trust_stats['medium']/total*100:.1f}%)")
    logger.info(f"    Low trust: {trust_stats['low']} ({trust_stats['low']/total*100:.1f}%)")
    logger.info(f"  Results:")
    logger.info(f"    Passed basic quality: {len(passed_basic_filtering)} ({len(passed_basic_filtering)/total*100:.1f}%)")
    logger.info(f"    High quality: {len(filtered_high)} ({len(filtered_high)/total*100:.1f}%)")
    logger.info(f"    Medium quality: {len(filtered_medium)} ({len(filtered_medium)/total*100:.1f}%)")
    logger.info(f"    Low quality: {len(filtered_low)} ({len(filtered_low)/total*100:.1f}%)")
    logger.info(f"    Total accepted: {total_accepted} ({total_accepted/total*100:.1f}%)")
    logger.info(f"    Rejected: {len(rejected)} ({len(rejected)/total*100:.1f}%)")
    logger.info(f"=======================================================")

    return filtered_high, filtered_medium, filtered_low, stats


def _check_structural_features(sequence: str) -> bool:
    """
    Check if sequence has TE structural features.
    Used as biological evidence to override low scores.
    """
    if not sequence or len(sequence) < 50:
        return False

    sequence = sequence.upper()
    seq_len = len(sequence)

    # 1. Check for Terminal Inverted Repeats (TIR)
    for tir_len in range(10, min(41, seq_len // 4)):
        left_tir = sequence[:tir_len]
        right_tir = sequence[-tir_len:]
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        right_tir_rc = ''.join(complement.get(b, 'N') for b in reversed(right_tir))
        mismatches = sum(1 for a, b in zip(left_tir, right_tir_rc) if a != b)
        if mismatches <= 2:
            return True

    # 2. Check for Long Terminal Repeats (LTR)
    if seq_len >= 260:
        for ltr_len in range(80, min(501, seq_len // 3)):
            left_ltr = sequence[:ltr_len]
            right_ltr = sequence[-ltr_len:]
            matches = sum(1 for a, b in zip(left_ltr, right_ltr) if a == b)
            if matches / ltr_len >= 0.75:
                return True

    # 3. Check for poly-A tail (LINE/SINE signature)
    tail_region = sequence[-50:] if seq_len >= 50 else sequence
    if 'AAAAAA' in tail_region:
        return True

    # 4. Check for TG...CA pattern (LTR signature)
    if sequence[:2] == 'TG' and sequence[-2:] == 'CA':
        return True

    return False
