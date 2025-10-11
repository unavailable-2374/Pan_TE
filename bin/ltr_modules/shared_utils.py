#!/usr/bin/env python3
"""
Shared Utilities for LTR Modules

Contains common functions, constants, and utilities used across multiple
LTR processing modules to avoid code duplication.
"""

import re
from collections import Counter
from Bio.Seq import Seq


class LTRConstants:
    """Common constants for LTR retrotransposon analysis"""

    # Terminal motifs
    LTR5_MOTIFS = ["TG", "TGT", "TGTA", "TGTG", "TGCA"]
    LTR3_MOTIFS = ["CA", "ACA", "TACA", "CACA", "TGCA"]

    # TSD parameters
    TSD_MIN_LEN = 4
    TSD_MAX_LEN = 6

    # LTR length constraints
    MIN_LTR_LEN = 100
    MAX_LTR_LEN = 1500

    # PBS patterns - expanded for different organisms
    PBS_PATTERNS = {
        'tRNAPro': 'TGGCGCCCAACGTGGGGC',
        'tRNATrp': 'TGGCGCCGTAACAGGGAC',
        'tRNAGln': 'TGGCGCCCGAACAGGGAC',
        'tRNALys': 'TGGCGCCCAACCTGGGA',
        'tRNAIle': 'TGGTAGCAGAGCTGGGAA',
        'tRNAMet': 'TGGCAGCAGGTCAGGGC',
        'tRNAAla': 'TGGCGCAGTGGCAGCGC',
        'tRNAArg': 'TGGACCGCTAGCTCAGTGGTA',
        'tRNAAsn': 'TGGCTCCGTAGCTCAATGG',
        'tRNAAsp': 'TGGGTCCGTAGTGTAGCGGT',
        'tRNACys': 'TGGCGCAGTGGAAGCGC',
        'tRNAGlu': 'TGGTTCCATGGTGAGGCC',
        'tRNAGly': 'TGGCGCGGTGGCGCAG',
        'tRNAHis': 'TGGCCGTGATCGTATAGTG',
        'tRNAPhe': 'TGGTGCGTTTAACCACTA',
        'tRNASer': 'TGGACGAGTGGCCCGAG',
        'tRNAThr': 'TGGCCGCGTGGCCCAAT',
        'tRNATyr': 'TGGGTGACCTCCCGGGC',
        'tRNAVal': 'TGGGTGATTAGCTCAGC'
    }

    # TSD nucleotide bias (empirical data)
    TSD_NUCLEOTIDE_BIAS = {
        'A': 0.30,
        'C': 0.20,
        'G': 0.20,
        'T': 0.30
    }

    # Evidence weights for boundary detection
    EVIDENCE_WEIGHTS = {
        'terminal_motifs': 0.30,
        'alignment_boundaries': 0.25,
        'kmer_transitions': 0.20,
        'tsd_evidence': 0.15,
        'internal_features': 0.10
    }


def calculate_similarity(seq1, seq2):
    """
    Calculate similarity between two sequences.

    Args:
        seq1: First sequence string
        seq2: Second sequence string

    Returns:
        float: Similarity score (0-1)
    """
    if not seq1 or not seq2:
        return 0.0

    # Handle sequences of different lengths
    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0

    # Convert to uppercase for comparison
    seq1 = seq1.upper()[:min_len]
    seq2 = seq2.upper()[:min_len]

    # Count matches
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != 'N' and b != 'N')

    return matches / min_len


def reverse_complement(seq: str) -> str:
    """
    Calculate reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        str: Reverse complement sequence
    """
    try:
        return str(Seq(seq).reverse_complement())
    except Exception:
        # Fallback to manual method
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                     'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                     'N': 'N', 'n': 'n'}
        return ''.join(complement.get(base, base) for base in reversed(seq))


def gc_content(seq: str) -> float:
    """
    Calculate GC content of a sequence.

    Args:
        seq: DNA sequence string

    Returns:
        float: GC content (0-1)
    """
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    total = len(seq)

    if total == 0:
        return 0.0

    return gc_count / total


def check_low_complexity(seq: str, window_size: int = 50, threshold: float = 0.3) -> bool:
    """
    Check if sequence has low complexity.

    Args:
        seq: DNA sequence string
        window_size: Window size for complexity calculation
        threshold: Complexity threshold (lower = simpler)

    Returns:
        bool: True if low complexity detected
    """
    if len(seq) < window_size:
        window_size = len(seq)

    if window_size == 0:
        return True

    # Calculate entropy-based complexity
    seq = seq.upper()

    # Count nucleotide frequencies
    counts = Counter(seq)
    total = len(seq)

    # Calculate Shannon entropy
    import math
    entropy = 0
    for count in counts.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)

    # Normalize to 0-1 (max entropy for DNA is log2(4) = 2)
    max_entropy = 2.0
    complexity = entropy / max_entropy if max_entropy > 0 else 0

    return complexity < threshold


def find_tandem_repeats(seq: str, min_unit: int = 2, max_unit: int = 50, min_copies: int = 3) -> float:
    """
    Find tandem repeats in sequence.

    Args:
        seq: DNA sequence string
        min_unit: Minimum repeat unit size
        max_unit: Maximum repeat unit size
        min_copies: Minimum number of copies to consider

    Returns:
        float: Proportion of sequence covered by tandem repeats (0-1)
    """
    if len(seq) < min_unit * min_copies:
        return 0.0

    seq = seq.upper()
    max_coverage = 0

    # Try different unit sizes
    for unit_size in range(min_unit, min(max_unit + 1, len(seq) // min_copies + 1)):
        for start in range(len(seq) - unit_size * min_copies + 1):
            unit = seq[start:start + unit_size]

            # Count consecutive copies
            copies = 1
            pos = start + unit_size

            while pos + unit_size <= len(seq):
                if seq[pos:pos + unit_size] == unit:
                    copies += 1
                    pos += unit_size
                else:
                    break

            if copies >= min_copies:
                coverage = (copies * unit_size) / len(seq)
                max_coverage = max(max_coverage, coverage)

    return max_coverage


def validate_sequence(seq: str, min_length: int = 10, max_n_content: float = 0.95) -> tuple:
    """
    Validate sequence for extreme cases.

    Args:
        seq: DNA sequence string
        min_length: Minimum sequence length
        max_n_content: Maximum N content (0-1)

    Returns:
        tuple: (is_valid: bool, reason: str)
    """
    if not seq:
        return False, "empty sequence"

    seq_len = len(seq)

    # Check minimum length
    if seq_len < min_length:
        return False, f"too short ({seq_len} bp < {min_length} bp)"

    # Check N content
    n_ratio = seq.upper().count('N') / seq_len
    if n_ratio > max_n_content:
        return False, f"excessive N content ({n_ratio:.1%} > {max_n_content:.1%})"

    return True, "pass"


def compute_pairwise_similarity(seq1: str, seq2: str, sample_size: int = None) -> float:
    """
    Compute pairwise similarity between two sequences.

    Args:
        seq1: First sequence
        seq2: Second sequence
        sample_size: If specified, sample this many positions for faster computation

    Returns:
        float: Similarity score (0-1)
    """
    if not seq1 or not seq2:
        return 0.0

    # Use shorter length
    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0

    seq1 = seq1.upper()[:min_len]
    seq2 = seq2.upper()[:min_len]

    # Sample if requested
    if sample_size and sample_size < min_len:
        import random
        indices = random.sample(range(min_len), sample_size)
        matches = sum(1 for i in indices if seq1[i] == seq2[i] and seq1[i] != 'N')
        return matches / sample_size
    else:
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != 'N')
        return matches / min_len


def compute_avg_similarity(consensus: str, sequences: list, sample_size: int = 10) -> float:
    """
    Compute average similarity between consensus and member sequences.

    Args:
        consensus: Consensus sequence
        sequences: List of member sequences
        sample_size: Number of sequences to sample (if > len(sequences), use all)

    Returns:
        float: Average similarity (0-1)
    """
    if not sequences:
        return 0.0

    import random

    # Sample sequences if there are many
    if len(sequences) > sample_size:
        sampled = random.sample(sequences, sample_size)
    else:
        sampled = sequences

    similarities = [compute_pairwise_similarity(consensus, seq) for seq in sampled]

    return sum(similarities) / len(similarities) if similarities else 0.0
