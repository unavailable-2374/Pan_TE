#!/usr/bin/env python3
"""
TE Structure Scoring Module v1.0 (based on stct.md)

Comprehensive structure evidence detection and scoring for TE families
from repeat finder libraries.

Scoring components:
  A. Terminal structure & boundary coherence (max +0.40)
  B. TSD evidence (max +0.20)
  C. Internal organization (max +0.25)
  D. Replication features (max +0.15)
  P. Penalties (max -0.20)

Final score: S = max(0, S_pos - 0.20 × S_neg) ∈ [0, 1]
"""

import re
import logging
from typing import Dict, List, Tuple, Optional
from collections import Counter
import numpy as np

logger = logging.getLogger(__name__)


# ============================================================================
# A. Terminal Structure & Boundary Coherence (max +0.40)
# ============================================================================

def score_terminal_pairs(consensus: Dict) -> float:
    """
    A1: Terminal paired structures (LTR/TIR) (0-0.28)

    LTR: Same-direction repeats at both ends (≥80 bp)
    TIR: Inverted repeats at both ends (≥10 bp)

    Score = min(1, L_repeat/threshold) × (identity/100) × 0.28

    Uses hybrid strategy: quick k-mer filter + Smith-Waterman precise alignment
    """
    sequence = consensus.get('sequence', '')
    if len(sequence) < 200:
        return 0.0

    # Try importing alignment utilities
    try:
        from alignment_utils import align_terminal_repeats

        # Check for LTR (direct repeats)
        ltr_result = align_terminal_repeats(
            sequence,
            terminal_window=500,
            min_repeat_length=80,
            min_identity=70.0,
            repeat_type='ltr'
        )

        # Check for TIR (inverted repeats)
        tir_result = align_terminal_repeats(
            sequence,
            terminal_window=100,
            min_repeat_length=10,
            min_identity=75.0,
            repeat_type='tir'
        )

        # Calculate scores
        ltr_score = 0.0
        if ltr_result:
            length_factor = min(1.0, ltr_result['length'] / 300.0)
            identity_factor = ltr_result['identity'] / 100.0
            ltr_score = length_factor * identity_factor
            logger.debug(f"LTR detected: length={ltr_result['length']}, "
                        f"identity={ltr_result['identity']:.1f}%, score={ltr_score:.3f}")

        tir_score = 0.0
        if tir_result:
            length_factor = min(1.0, tir_result['length'] / 40.0)
            identity_factor = tir_result['identity'] / 100.0
            tir_score = length_factor * identity_factor
            logger.debug(f"TIR detected: length={tir_result['length']}, "
                        f"identity={tir_result['identity']:.1f}%, score={tir_score:.3f}")

        # Take max (avoid double counting)
        return max(ltr_score, tir_score) * 0.28

    except ImportError:
        # align_terminal_repeats not implemented yet; use simplified detection
        # (logged once at debug level to avoid per-sequence spam)
        ltr_score = _detect_ltr_simple(sequence)
        tir_score = _detect_tir_simple(sequence)
        return max(ltr_score, tir_score) * 0.28


def _detect_ltr_simple(sequence: str) -> float:
    """Simplified LTR detection (checks 5' and 3' terminal similarity)"""
    # Check terminal 500bp windows
    window = min(500, len(sequence) // 4)
    if len(sequence) < 2 * window:
        return 0.0

    left = sequence[:window]
    right = sequence[-window:]

    # Simple similarity check (in real implementation, use alignment)
    matches = sum(1 for a, b in zip(left, right) if a == b)
    identity = matches / window * 100

    if identity < 60:  # Threshold for LTR-like
        return 0.0

    # Estimate LTR length (simplified)
    ltr_length = window * 0.5  # Assume half of window is LTR
    length_factor = min(1.0, ltr_length / 300.0)
    identity_factor = identity / 100.0

    return length_factor * identity_factor


def _detect_tir_simple(sequence: str) -> float:
    """Simplified TIR detection (checks inverted terminal repeats)"""
    # Check terminal 100bp windows
    window = min(100, len(sequence) // 4)
    if len(sequence) < 2 * window:
        return 0.0

    left = sequence[:window]
    right_rc = _reverse_complement(sequence[-window:])

    # Simple similarity check
    matches = sum(1 for a, b in zip(left, right_rc) if a == b)
    identity = matches / window * 100

    if identity < 70:  # Higher threshold for TIR
        return 0.0

    # Estimate TIR length
    tir_length = window * 0.3
    length_factor = min(1.0, tir_length / 40.0)
    identity_factor = identity / 100.0

    return length_factor * identity_factor


def _reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(b.upper(), 'N') for b in reversed(seq))


def score_terminal_motifs(consensus: Dict) -> float:
    """
    A2: Terminal conserved motifs (0-0.06)

    LTR: TG...CA motifs
    TIR: Common terminal motifs
    Helitron: 5'TC, 3'CTRR
    """
    sequence = consensus.get('sequence', '').upper()
    if len(sequence) < 50:
        return 0.0

    score = 0.0

    # LTR motifs: TG at 5', CA at 3'
    if sequence.startswith('TG'):
        score += 0.03
    if sequence.endswith('CA'):
        score += 0.03

    # TIR common motifs: TA, TTAA
    if sequence[:2] == 'TA' and sequence[-2:] == 'TA':
        score += 0.04
    elif sequence[:4] == 'TTAA' and sequence[-4:] == 'TTAA':
        score += 0.06

    # Helitron motifs
    if sequence[:2] == 'TC':
        score += 0.03
    if re.search(r'CT[AG][AG]$', sequence):  # CTRR at 3'
        score += 0.03

    return min(0.06, score)


def score_boundary_coherence(consensus: Dict, genome_hits: Optional[List] = None) -> float:
    """
    A3: Boundary coherence (0-0.06)

    Check if RepeatMasker hits have consistent start/end positions
    (≥70% within ±10bp → +0.06; 50-70% → +0.03)
    """
    if not genome_hits:
        # Use rm_hits from consensus if available
        genome_hits = consensus.get('rm_hits', [])

    if not genome_hits or len(genome_hits) < 10:
        return 0.0  # Not enough data

    # Extract start and end positions from hits
    starts = []
    ends = []

    for hit in genome_hits[:50]:  # Use top 50 non-overlapping hits
        if isinstance(hit, dict):
            starts.append(hit.get('query_start', 0))
            ends.append(hit.get('query_end', 0))
        # Handle other hit formats as needed

    if not starts or not ends:
        return 0.0

    # Calculate coherence (% of hits within ±10bp of mode)
    start_coherence = _calculate_position_coherence(starts, window=10)
    end_coherence = _calculate_position_coherence(ends, window=10)

    avg_coherence = (start_coherence + end_coherence) / 2.0

    if avg_coherence >= 0.70:
        return 0.06
    elif avg_coherence >= 0.50:
        return 0.03
    else:
        return 0.0


def _calculate_position_coherence(positions: List[int], window: int = 10) -> float:
    """Calculate what fraction of positions cluster within ±window of mode"""
    if not positions:
        return 0.0

    # Find mode (most common position)
    counter = Counter(positions)
    if not counter:
        return 0.0

    mode_pos = counter.most_common(1)[0][0]

    # Count positions within window
    in_window = sum(1 for p in positions if abs(p - mode_pos) <= window)

    return in_window / len(positions)


# ============================================================================
# B. TSD (Target Site Duplication) Evidence (max +0.20)
# ============================================================================

def score_tsd_rate(consensus: Dict, genome_hits: Optional[List] = None) -> Tuple[float, int]:
    """
    B1: TSD mode length hit rate (0-0.12)

    Returns: (score, mode_length)

    Rate ≥0.6 → +0.12
    Rate 0.4-0.6 → +0.08
    Rate 0.25-0.4 → +0.04
    """
    if not genome_hits:
        genome_hits = consensus.get('rm_hits', [])

    if not genome_hits:
        # Fallback: check if consensus has TSD annotation
        tsd = consensus.get('tsd', 'N/A')
        if tsd not in ['N/A', 'NNNN', None] and not str(tsd).startswith('NNNN'):
            # Has TSD, assume moderate rate
            return 0.08, len(tsd)
        return 0.0, 0

    # Infer TSD from genome hits (simplified)
    # In full implementation: extract flanking sequences and find repeats
    tsd_rate, mode_length = _infer_tsd_from_hits(genome_hits)

    if tsd_rate >= 0.6:
        score = 0.12
    elif tsd_rate >= 0.4:
        score = 0.08
    elif tsd_rate >= 0.25:
        score = 0.04
    else:
        score = 0.0

    return score, mode_length


def _infer_tsd_from_hits(genome_hits: List) -> Tuple[float, int]:
    """
    Infer TSD from RepeatMasker hits (simplified)

    In full implementation: extract flanking 30bp, find exact duplications
    """
    # Simplified: use TSD info if available in hits
    tsd_lengths = []

    for hit in genome_hits[:100]:
        if isinstance(hit, dict):
            tsd = hit.get('tsd')
            if tsd and tsd not in ['N/A', 'NNNN', None]:
                tsd_lengths.append(len(tsd))

    if not tsd_lengths:
        return 0.0, 0

    # Find mode length
    counter = Counter(tsd_lengths)
    mode_length, mode_count = counter.most_common(1)[0]

    # Calculate rate
    rate = mode_count / len(genome_hits[:100])

    return rate, mode_length


def score_tsd_expected(mode_length: int, consensus: Dict) -> float:
    """
    B2: TSD matches expected length for TE class (0-0.08)

    LTR: 4-6 bp → +0.06
    DNA/TIR: 2 bp (TA/TTAA +0.02 bonus) → +0.06-0.08
    LINE: 7-20 bp → +0.06
    """
    if mode_length == 0:
        return 0.0

    score = 0.0

    # Try to infer TE class from consensus features
    sequence = consensus.get('sequence', '').upper()

    # LTR-like: 4-6 bp TSD
    if 4 <= mode_length <= 6:
        score += 0.06

    # DNA/TIR: 2 bp or 4 bp (TTAA)
    elif mode_length == 2:
        score += 0.06
        # Check for TA
        tsd_seq = consensus.get('tsd', '')
        if tsd_seq == 'TA':
            score += 0.02
    elif mode_length == 4:
        tsd_seq = consensus.get('tsd', '')
        if tsd_seq == 'TTAA':
            score += 0.08

    # LINE: 7-20 bp
    elif 7 <= mode_length <= 20:
        score += 0.06

    return min(0.08, score)


# ============================================================================
# C. Internal Organization (max +0.25)
# ============================================================================

def score_hallmark_domains(consensus: Dict) -> float:
    """
    C1: Hallmark protein domains via HMM (0-0.17)

    LTR-RT: RT+0.05, RNaseH+0.05, Integrase+0.05, all three+0.02
    LINE: EN+0.05, RT+0.07, ORF1+0.05
    DNA/TIR: DDE transposase+0.10, auxiliary+0.05
    Helitron: RepHel+0.12

    Note: This is a simplified version. Full implementation requires
    PFAM/TREP/TEfam HMM searches.
    """
    # Check if HMM domain info is already annotated
    domains = consensus.get('protein_domains', [])

    if not domains:
        # Fallback: infer from validation score or quality
        # In real implementation, run hmmscan against PFAM/TREP
        return 0.0

    score = 0.0
    domain_types = set(d.lower() for d in domains)

    # LTR-RT domains
    rt_found = any('rt' in d or 'reverse_transcriptase' in d for d in domain_types)
    rnase_found = any('rnase' in d for d in domain_types)
    int_found = any('integrase' in d or 'int' in d for d in domain_types)

    if rt_found and rnase_found and int_found:
        score += 0.17  # All three + bonus
    elif rt_found:
        score += 0.05
    if rnase_found:
        score += 0.05
    if int_found:
        score += 0.05

    # LINE domains
    en_found = any('endonuclease' in d or 'ape' in d for d in domain_types)
    orf1_found = any('orf1' in d or 'rrm' in d or 'ccch' in d for d in domain_types)

    if en_found:
        score += 0.05
    if orf1_found:
        score += 0.05

    # DNA/TIR
    tpase_found = any('transposase' in d or 'dde' in d for d in domain_types)
    if tpase_found:
        score += 0.10

    # Helitron
    rephel_found = any('rephel' in d or 'helicase' in d for d in domain_types)
    if rephel_found:
        score += 0.12

    return min(0.17, score)


def score_orf_continuity(consensus: Dict) -> float:
    """
    C2: ORF continuity (0-0.08)

    Main ORF length ≥70% expected and ≤1 internal stop → +0.08
    ≥50% expected → +0.04
    """
    sequence = consensus.get('sequence', '')
    if len(sequence) < 300:
        return 0.0

    # Find longest ORF in all six frames
    longest_orf_len = _find_longest_orf(sequence)

    # Expected ORF lengths by TE class (simplified)
    # LINE ORF2: ~1200-1500 bp
    # LTR pol: ~1500-2000 bp
    # DNA transposase: ~900-1200 bp
    expected_length = 1200  # Default expectation

    orf_ratio = longest_orf_len / expected_length

    if orf_ratio >= 0.70:
        # Check internal stops (simplified: if ORF is long, assume few stops)
        return 0.08
    elif orf_ratio >= 0.50:
        return 0.04
    else:
        return 0.0


def _find_longest_orf(sequence: str) -> int:
    """Find longest ORF in all six reading frames"""
    stop_codons = {'TAA', 'TAG', 'TGA'}
    max_length = 0

    # Check all 6 frames (3 forward, 3 reverse)
    for seq in [sequence, _reverse_complement(sequence)]:
        for frame in range(3):
            orf_lengths = []
            current_length = 0

            for i in range(frame, len(seq) - 2, 3):
                codon = seq[i:i+3].upper()
                if len(codon) < 3:
                    break

                if codon in stop_codons:
                    if current_length > 0:
                        orf_lengths.append(current_length)
                    current_length = 0
                else:
                    current_length += 3

            if current_length > 0:
                orf_lengths.append(current_length)

            if orf_lengths:
                max_length = max(max_length, max(orf_lengths))

    return max_length


# ============================================================================
# D. Replication Features (max +0.15)
# ============================================================================

def score_polyA_tail(consensus: Dict) -> float:
    """
    D1: LINE 3' Poly(A) tail (0-0.07)

    3' end 30nt A-content ≥70% → +0.05
    Plus ≥10 consecutive A's → +0.02
    """
    sequence = consensus.get('sequence', '').upper()
    if len(sequence) < 30:
        return 0.0

    tail = sequence[-30:]
    a_content = tail.count('A') / 30.0

    score = 0.0
    if a_content >= 0.70:
        score += 0.05

    # Check for consecutive A's
    max_consecutive_a = max(len(m.group()) for m in re.finditer(r'A+', tail)) if 'A' in tail else 0
    if max_consecutive_a >= 10:
        score += 0.02

    return score


def score_pbs_ppt(consensus: Dict) -> float:
    """
    D2: LTR PBS/PPT signals (0-0.08)

    PBS: tRNA match in 100nt downstream of 5' LTR → +0.04
    PPT: U-rich (≥70%) in 10-30nt upstream of 3' LTR → +0.04

    Simplified implementation (full version requires tRNA database)
    """
    sequence = consensus.get('sequence', '').upper()
    if len(sequence) < 500:
        return 0.0

    score = 0.0

    # PBS check (simplified: look for tRNA-like patterns near 5' end)
    # In full implementation: match against tRNA database
    pbs_region = sequence[100:200]  # After 5' LTR region
    if _has_trna_like_pattern(pbs_region):
        score += 0.04

    # PPT check: U-rich region before 3' end
    ppt_region = sequence[-150:-120]  # Before 3' LTR region
    t_content = ppt_region.count('T') / len(ppt_region)  # T = U in DNA
    if t_content >= 0.70:
        score += 0.04

    return score


def _has_trna_like_pattern(sequence: str) -> bool:
    """Check for tRNA-like patterns (simplified)"""
    # tRNAs typically have specific motifs and stem-loop structures
    # Simplified: check for common tRNA primer sequences
    trna_primers = ['TGGCGCCCGAACAG', 'TGGTATCAGAGCGCC', 'GCCCGGATAGCTCA']

    for primer in trna_primers:
        if primer in sequence or _reverse_complement(primer) in sequence:
            return True

    return False


# ============================================================================
# P. Penalties (max -0.20)
# ============================================================================

def penalty_low_complexity(consensus: Dict) -> float:
    """
    P1: Low complexity / tandem repeat contamination (0-0.10)

    >40% low complexity → -0.10
    20-40% → -0.05
    """
    sequence = consensus.get('sequence', '')
    if not sequence:
        return 0.0

    # Calculate low complexity (simplified: use base composition entropy)
    from stratified_biological_filter import calculate_sequence_complexity

    complexity = calculate_sequence_complexity(sequence)

    # Low complexity penalty (inverse of complexity)
    if complexity < 0.30:  # Very low complexity
        return 0.10
    elif complexity < 0.50:  # Moderate low complexity
        return 0.05
    else:
        return 0.0


def penalty_host_contamination(consensus: Dict) -> float:
    """
    P2: Host component contamination (0-0.06)

    rRNA/tRNA/protein-coding exon significant hits → -0.06

    Simplified: check sequence characteristics
    """
    sequence = consensus.get('sequence', '').upper()
    if not sequence:
        return 0.0

    # Simplified heuristics (full implementation needs BLAST against host)

    # rRNA-like: very high GC content and specific motifs
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
    if gc_content > 0.65 and ('GGAAGGCCUG' in sequence or 'UGCCUGGCGG' in sequence):
        return 0.06

    # tRNA-like: specific length and structure
    if 70 <= len(sequence) <= 90 and gc_content > 0.55:
        if sequence.startswith('GC') and sequence.endswith('CCA'):
            return 0.06

    return 0.0


def penalty_incompatible_domains(consensus: Dict) -> float:
    """
    P3: Incompatible domain combinations (0-0.04)

    Mutually exclusive domains (e.g., RT + Transposase) → -0.04
    """
    domains = consensus.get('protein_domains', [])
    if not domains:
        return 0.0

    domain_types = set(d.lower() for d in domains)

    # Check for incompatible combinations
    has_rt = any('rt' in d or 'reverse_transcriptase' in d for d in domain_types)
    has_tpase = any('transposase' in d for d in domain_types)
    has_helicase = any('helicase' in d or 'rephel' in d for d in domain_types)

    # RT + Transposase is incompatible
    if has_rt and has_tpase:
        return 0.04

    # Helicase + RT is incompatible (Helitron vs LTR/LINE)
    if has_helicase and has_rt:
        return 0.04

    return 0.0


# ============================================================================
# Main Scoring Function with Normalization
# ============================================================================

def calculate_structure_score_v1(consensus: Dict, genome_hits: Optional[List] = None,
                                  options: Optional[Dict] = None) -> Dict:
    """
    Calculate comprehensive structure score S ∈ [0, 1] (stct.md v1.0)

    Args:
        consensus: Consensus sequence dict
        genome_hits: RepeatMasker hits for this family
        options: Optional configuration (which modules to run)

    Returns:
        Dict with:
            - structure_score: Final S value
            - positive_score: S_pos before penalty
            - negative_score: S_neg
            - components: Breakdown of individual scores
            - structure_tier: Strong/Medium/Weak/Poor
    """
    if options is None:
        options = {'run_hmm': False}  # Default: skip HMM (expensive)

    # =================================================================
    # Positive evidence
    # =================================================================
    positive_scores = {}

    # A. Terminal structure & boundary
    positive_scores['terminal_pairs'] = score_terminal_pairs(consensus)
    positive_scores['terminal_motifs'] = score_terminal_motifs(consensus)
    positive_scores['boundary_coherence'] = score_boundary_coherence(consensus, genome_hits)

    # B. TSD
    tsd_score, tsd_mode_len = score_tsd_rate(consensus, genome_hits)
    positive_scores['tsd_rate'] = tsd_score
    positive_scores['tsd_expected'] = score_tsd_expected(tsd_mode_len, consensus)

    # C. Internal organization
    if options.get('run_hmm', False):
        positive_scores['hmm_domains'] = score_hallmark_domains(consensus)
    else:
        positive_scores['hmm_domains'] = None  # Not run

    positive_scores['orf_continuity'] = score_orf_continuity(consensus)

    # D. Replication features
    positive_scores['polyA_tail'] = score_polyA_tail(consensus)
    positive_scores['pbs_ppt'] = score_pbs_ppt(consensus)

    # =================================================================
    # Negative evidence (penalties)
    # =================================================================
    negative_scores = {}

    negative_scores['low_complexity'] = penalty_low_complexity(consensus)
    negative_scores['host_contamination'] = penalty_host_contamination(consensus)
    negative_scores['incompatible_domains'] = penalty_incompatible_domains(consensus)

    # =================================================================
    # Normalization (stct.md formula)
    # =================================================================

    # Maximum possible positive scores (from stct.md)
    MAX_POSITIVE = {
        'terminal_pairs': 0.28,
        'terminal_motifs': 0.06,
        'boundary_coherence': 0.06,
        'tsd_rate': 0.12,
        'tsd_expected': 0.08,
        'hmm_domains': 0.17,
        'orf_continuity': 0.08,
        'polyA_tail': 0.07,
        'pbs_ppt': 0.08
    }

    # Maximum possible negative scores
    MAX_NEGATIVE = {
        'low_complexity': 0.10,
        'host_contamination': 0.06,
        'incompatible_domains': 0.04
    }

    # Sum actual scores and max scores (only for features that were run)
    P = sum(v for v in positive_scores.values() if v is not None)
    M_pos = sum(MAX_POSITIVE[k] for k, v in positive_scores.items() if v is not None)

    Q = sum(v for v in negative_scores.values() if v is not None)
    M_neg = sum(MAX_NEGATIVE[k] for k, v in negative_scores.items() if v is not None)

    # Normalized scores
    S_pos = 0.0 if M_pos == 0 else min(1.0, P / M_pos)
    S_neg = 0.0 if M_neg == 0 else min(1.0, Q / M_neg)

    # Final structure score (stct.md formula)
    S = max(0.0, S_pos - 0.20 * S_neg)

    # ===================================================================
    # Intelligent fallback: integrate base quality when structure
    # evidence is limited (e.g., no genome hits, no HMM)
    # ===================================================================
    base_quality = consensus.get('quality_score', 0.0)
    validation_score = consensus.get('validation_score', 0.0)

    # Check if we have limited structural evidence
    has_limited_evidence = (
        P < 0.15  # Very low raw positive score
        or M_pos < 0.5  # Few features actually ran
        or not genome_hits  # No genome hits for TSD/boundary
    )

    if has_limited_evidence and (base_quality > 0 or validation_score > 0):
        # Enhanced hybrid score: rely more on quality when structure is weak
        # Weight: 40% structure + 60% quality (increased from 40%)
        # Quality component: 80% base quality + 20% validation (increased from 70/30)
        quality_component = base_quality * 0.8 + validation_score * 0.2
        S_hybrid = S * 0.4 + quality_component * 0.6

        logger.debug(f"Limited structural evidence detected, using enhanced hybrid score: "
                    f"S_struct={S:.3f}, S_quality={quality_component:.3f}, "
                    f"S_hybrid={S_hybrid:.3f}")

        S = S_hybrid

    # Structure tier classification (updated to match relaxed filtering thresholds)
    if S >= 0.50:
        tier = 'Strong'
    elif S >= 0.40:
        tier = 'Medium'
    elif S >= 0.30:
        tier = 'Weak'
    else:
        tier = 'Poor'

    return {
        'structure_score': S,
        'positive_score': S_pos,
        'negative_score': S_neg,
        'structure_tier': tier,
        'components': {
            'positive': positive_scores,
            'negative': negative_scores
        },
        'raw_sums': {
            'P': P,
            'M_pos': M_pos,
            'Q': Q,
            'M_neg': M_neg
        },
        'used_hybrid': has_limited_evidence and (base_quality > 0 or validation_score > 0)
    }
