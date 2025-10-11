#!/usr/bin/env python3
"""
LTR Identifier Module

Identifies LTR retrotransposon sequences using 3 core features:
1. LTR similarity (two terminal repeats)
2. TSD (Target Site Duplication)
3. Boundary motifs (TG...CA pattern)

This module focuses on IDENTIFICATION (is it an LTR?) rather than
quality control (how good is the sequence?).

Date: 2025-10-08
"""

import os
import sys
import logging
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class LTRQualityValidator:
    """
    LTR Retrotransposon Identifier

    Uses minimal sufficient criteria for LTR identification:
    - Core Feature 1: LTR similarity (>= 80% for recent, >= 50% for ancient)
    - Core Feature 2: TSD presence (3-9bp)
    - Core Feature 3: Boundary motifs (TG...CA, optional but recommended)

    Classification:
    - 'confirmed': All 3 features present
    - 'probable': Features 1+2 present
    - 'possible': Feature 1 present (high LTR similarity)
    - 'not_ltr': Does not meet criteria
    """

    def __init__(self, min_ltr_len=100, max_ltr_len=10000,
                 tsd_min=3, tsd_max=9,
                 min_ltr_similarity=0.50,
                 logger=None):
        """
        Initialize LTR identifier.

        Args:
            min_ltr_len: Minimum LTR length (default: 100)
            max_ltr_len: Maximum LTR length (default: 10000)
            tsd_min: Minimum TSD length (default: 3, expanded from 4)
            tsd_max: Maximum TSD length (default: 9, expanded from 6)
            min_ltr_similarity: Minimum LTR similarity (default: 0.50 for ancient elements)
            logger: Logger instance (optional)
        """
        self.min_ltr_len = min_ltr_len
        self.max_ltr_len = max_ltr_len
        self.tsd_min = tsd_min
        self.tsd_max = tsd_max
        self.min_ltr_similarity = min_ltr_similarity
        
        # Setup logger
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger("LTRQualityValidator")
            if not self.logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter(
                    '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
                )
                handler.setFormatter(formatter)
                self.logger.addHandler(handler)
                self.logger.setLevel(logging.INFO)
        
        # LTR motif patterns
        self.ltr5_motifs = ["TG", "TGT", "TGTA", "TGTG", "TGCA"]
        self.ltr3_motifs = ["CA", "ACA", "TACA", "CACA", "TGCA"]
        
        # PBS patterns (for complete element scoring)
        self.pbs_patterns = {
            'tRNAPro': 'TGGCGCCCAACGTGGGGC',
            'tRNATrp': 'TGGCGCCGTAACAGGGAC',
            'tRNAGln': 'TGGCGCCCGAACAGGGAC',
            'tRNALys': 'TGGCGCCCAACCTGGGA',
        }
    
    def _calculate_similarity(self, seq1, seq2):
        """Calculate sequence similarity (0-1)."""
        if not seq1 or not seq2:
            return 0.0
        
        # Align sequences using simple character-by-character comparison
        matches = sum(a == b for a, b in zip(seq1, seq2))
        max_len = max(len(seq1), len(seq2))
        
        if max_len == 0:
            return 0.0
        
        return matches / max_len


    def validate_ltr_boundaries(self, seq, strict=False):
        """
        Validate LTR boundary motifs (5' and 3' ends).

        This method checks if the sequence has canonical LTR boundary patterns.
        Moved from build_ltr_consensus.py as boundary validation is a core
        responsibility of the boundary optimizer.

        Args:
            seq: Sequence string to validate (str or Bio.Seq object)
            strict: If True, only accept canonical TG...CA pattern
                   If False, also accept common variants (CA...CA, TG...TA)

        Returns:
            Dictionary with:
            - 'valid': bool - whether boundaries are valid
            - '5p_motif': str - detected 5' motif (or None)
            - '3p_motif': str - detected 3' motif (or None)
            - 'pattern': str - detected pattern type (or None)
            - 'confidence': float - confidence score (0-1)
        """
        # Handle Bio.Seq objects
        seq_str = str(seq).upper() if seq else ''

        # Minimum length check
        if len(seq_str) < 4:
            self.logger.debug(f"Boundary validation failed: sequence too short ({len(seq_str)}bp)")
            return {
                'valid': False,
                '5p_motif': None,
                '3p_motif': None,
                'pattern': None,
                'confidence': 0.0,
                'reason': 'sequence too short'
            }

        # Extract terminal dinucleotides
        start_di = seq_str[:2]
        end_di = seq_str[-2:]

        # Define boundary patterns with confidence scores
        # Format: (5' motif, 3' motif, pattern name, confidence)
        boundary_patterns = [
            ('TG', 'CA', 'canonical', 1.0),      # Canonical LTR pattern
            ('CA', 'CA', 'variant_CA', 0.8),     # CA...CA variant
            ('TG', 'TA', 'variant_TA', 0.8),     # TG...TA variant
            ('TA', 'CA', 'variant_5TA', 0.6),    # TA...CA variant (less common)
            ('TG', 'TG', 'variant_TG', 0.5),     # TG...TG (rare)
        ]

        # Check against patterns
        for motif_5p, motif_3p, pattern_name, confidence in boundary_patterns:
            if start_di == motif_5p and end_di == motif_3p:
                # In strict mode, only accept canonical pattern
                if strict and pattern_name != 'canonical':
                    continue

                self.logger.debug(f"Boundary validation passed: {start_di}...{end_di} ({pattern_name}, confidence={confidence:.2f})")
                return {
                    'valid': True,
                    '5p_motif': motif_5p,
                    '3p_motif': motif_3p,
                    'pattern': pattern_name,
                    'confidence': confidence,
                    'reason': f'matched {pattern_name} pattern'
                }

        # Check for partial matches (useful for diagnostic purposes)
        partial_match_5p = any(start_di == p[0] for p in boundary_patterns)
        partial_match_3p = any(end_di == p[1] for p in boundary_patterns)

        reason = 'no valid boundary motifs'
        if partial_match_5p and not partial_match_3p:
            reason = f"valid 5' motif ({start_di}) but invalid 3' motif ({end_di})"
        elif partial_match_3p and not partial_match_5p:
            reason = f"valid 3' motif ({end_di}) but invalid 5' motif ({start_di})"
        elif partial_match_5p and partial_match_3p:
            reason = f"partial match: {start_di}...{end_di} (mismatched combination)"
        else:
            reason = f"invalid motifs: {start_di}...{end_di}"

        self.logger.debug(f"Boundary validation failed: {reason}")
        return {
            'valid': False,
            '5p_motif': start_di if partial_match_5p else None,
            '3p_motif': end_di if partial_match_3p else None,
            'pattern': None,
            'confidence': 0.0,
            'reason': reason
        }

    def comprehensive_quality_validation(self, seq, check_boundary=True, strict_boundary=False):
        """
        Comprehensive quality validation for LTR sequences AFTER boundary correction.

        This method should be called AFTER boundary optimization to perform complete
        quality checks. It includes all validations that were intentionally skipped
        in build_ltr_consensus.py's minimal pre-filtering.

        Args:
            seq: Sequence string to validate (str or Bio.Seq object)
            check_boundary: Whether to validate boundary motifs
            strict_boundary: If True, only accept canonical TG...CA pattern

        Returns:
            Dictionary with:
            - 'valid': bool - overall validation result
            - 'checks': dict - individual check results
            - 'score': float - quality score (0-100)
            - 'warnings': list - non-fatal quality issues
            - 'errors': list - fatal quality issues
        """
        seq_str = str(seq).upper() if seq else ''
        seq_len = len(seq_str)

        checks = {}
        warnings = []
        errors = []
        scores = []

        # 1. Boundary motif validation
        if check_boundary:
            boundary_result = self.validate_ltr_boundaries(seq_str, strict=strict_boundary)
            checks['boundary'] = boundary_result

            if boundary_result['valid']:
                scores.append(boundary_result['confidence'] * 100)
            else:
                errors.append(f"Invalid boundary motifs: {boundary_result['reason']}")
                scores.append(0)

        # 2. Length validation
        min_len = self.min_ltr_len if hasattr(self, 'min_ltr_len') else 100
        max_len = self.max_ltr_len if hasattr(self, 'max_ltr_len') else 10000

        if seq_len < min_len:
            checks['length'] = {'valid': False, 'value': seq_len, 'min': min_len}
            errors.append(f"Too short: {seq_len} bp < {min_len} bp")
            scores.append(0)
        elif seq_len > max_len:
            checks['length'] = {'valid': False, 'value': seq_len, 'max': max_len}
            errors.append(f"Too long: {seq_len} bp > {max_len} bp")
            scores.append(0)
        else:
            checks['length'] = {'valid': True, 'value': seq_len}
            # Score based on typical LTR length (100-1500 bp)
            if 100 <= seq_len <= 1500:
                scores.append(100)
            elif seq_len < 100:
                scores.append(max(50, seq_len / 2))
            else:
                scores.append(max(50, 100 - (seq_len - 1500) / 100))

        # 3. N content validation
        n_count = seq_str.count('N')
        n_ratio = n_count / seq_len if seq_len > 0 else 1.0
        max_n_ratio = 0.20  # 20% threshold

        checks['n_content'] = {'valid': n_ratio <= max_n_ratio, 'ratio': n_ratio, 'count': n_count}

        if n_ratio > max_n_ratio:
            errors.append(f"Excessive N content: {n_ratio:.1%} > {max_n_ratio:.1%}")
            scores.append(0)
        elif n_ratio > 0.10:
            warnings.append(f"Moderate N content: {n_ratio:.1%}")
            scores.append(50)
        else:
            scores.append(100)

        # 4. Low complexity check
        low_complexity_score = self._check_sequence_complexity(seq_str)
        checks['complexity'] = {'score': low_complexity_score}

        if low_complexity_score < 30:
            errors.append(f"Low sequence complexity: score {low_complexity_score:.1f}")
            scores.append(0)
        elif low_complexity_score < 50:
            warnings.append(f"Moderate sequence complexity: score {low_complexity_score:.1f}")
            scores.append(low_complexity_score)
        else:
            scores.append(low_complexity_score)

        # 5. Tandem repeat check
        # Note: Real LTR sequences often have some repetitive elements
        # Use relaxed thresholds: 70% for error, 50% for warning
        tandem_ratio = self._check_tandem_repeats(seq_str)
        checks['tandem_repeats'] = {'ratio': tandem_ratio}

        if tandem_ratio > 0.70:
            errors.append(f"Excessive tandem repeats: {tandem_ratio:.1%} of sequence")
            scores.append(0)
        elif tandem_ratio > 0.50:
            warnings.append(f"Moderate tandem repeats: {tandem_ratio:.1%}")
            scores.append(60)
        else:
            scores.append(100)

        # 6. GC content check (optional, for info)
        gc_count = seq_str.count('G') + seq_str.count('C')
        valid_bases = sum(1 for b in seq_str if b in 'ACGT')
        gc_ratio = gc_count / valid_bases if valid_bases > 0 else 0.5

        checks['gc_content'] = {'ratio': gc_ratio}

        # Extreme GC bias might indicate issues
        if gc_ratio < 0.10 or gc_ratio > 0.90:
            warnings.append(f"Extreme GC content: {gc_ratio:.1%}")

        # Calculate overall score
        overall_score = sum(scores) / len(scores) if scores else 0

        # Determine if valid (no fatal errors)
        is_valid = len(errors) == 0

        # Log validation results
        if is_valid:
            self.logger.debug(f"Quality validation passed: score={overall_score:.1f}, warnings={len(warnings)}")
        else:
            self.logger.info(f"Quality validation failed: score={overall_score:.1f}, errors={errors}")

        if warnings:
            self.logger.debug(f"Quality warnings: {warnings}")

        return {
            'valid': is_valid,
            'checks': checks,
            'score': overall_score,
            'warnings': warnings,
            'errors': errors
        }

    def _check_sequence_complexity(self, seq_str, window=50):
        """
        Check sequence complexity using sliding window.

        Returns:
            float: Complexity score (0-100), higher is more complex
        """
        if len(seq_str) < window:
            window = len(seq_str)

        if window < 10:
            return 50  # Cannot assess reliably

        complexity_scores = []

        for i in range(len(seq_str) - window + 1):
            win_seq = seq_str[i:i+window]

            # Calculate entropy-based complexity
            base_counts = {
                'A': win_seq.count('A'),
                'C': win_seq.count('C'),
                'G': win_seq.count('G'),
                'T': win_seq.count('T')
            }

            total = sum(base_counts.values())
            if total == 0:
                continue

            # Shannon entropy
            entropy = 0
            for count in base_counts.values():
                if count > 0:
                    p = count / total
                    entropy -= p * np.log2(p)

            # Normalize to 0-100 scale (max entropy for 4 bases is 2.0)
            complexity = (entropy / 2.0) * 100
            complexity_scores.append(complexity)

        if not complexity_scores:
            return 50

        # Return average complexity
        return sum(complexity_scores) / len(complexity_scores)

    def _check_tandem_repeats(self, seq_str, min_unit=2, max_unit=50, min_copies=3):
        """
        Check for tandem repeats.

        Returns:
            float: Ratio of sequence covered by tandem repeats (0-1)
        """
        seq_len = len(seq_str)
        if seq_len < min_unit * min_copies:
            return 0.0

        repeat_positions = set()

        # Check various repeat unit sizes
        for unit_size in range(min_unit, min(max_unit + 1, seq_len // min_copies + 1)):
            for start in range(seq_len - unit_size * min_copies + 1):
                unit = seq_str[start:start + unit_size]

                # Skip units with N
                if 'N' in unit:
                    continue

                # Count consecutive repeats
                copies = 1
                pos = start + unit_size

                while pos + unit_size <= seq_len and seq_str[pos:pos + unit_size] == unit:
                    copies += 1
                    pos += unit_size

                # If found significant tandem repeat
                if copies >= min_copies:
                    # Mark positions as part of repeat
                    for i in range(start, start + copies * unit_size):
                        repeat_positions.add(i)

        # Calculate ratio
        return len(repeat_positions) / seq_len if seq_len > 0 else 0.0

    def _score_ltr_structure(self, seq_record, alignments=None, genome_seqs=None):
        """
        Score LTR structural completeness using multi-dimensional criteria.

        IMPROVED: Automatically detects sequence type based on terminal LTR similarity
        rather than arbitrary length cutoffs.

        Detection strategy:
        1. Fast terminal similarity check (multi-scale)
        2. Classify as:
           - Complete element: Terminal LTR similarity >= 40%
           - Single LTR: Terminal LTR similarity < 40%

        This approach correctly handles:
        - Short complete elements (~1800bp with high terminal similarity)
        - Long single LTRs (~1800bp with low terminal similarity)
        - Prevents misclassification of edge cases

        Args:
            seq_record: SeqRecord object
            alignments: List of alignments (optional, for TSD detection)
            genome_seqs: Genome access object (optional, for TSD extraction)

        Returns:
            Dictionary with:
            - 'total_score': Overall score (0-100)
            - 'category': 'confirmed', 'probable', 'possible', 'not_ltr'
            - 'components': Individual component scores
            - 'issues': List of detected structural issues
            - 'sequence_type': 'single_ltr' or 'complete_element'
            - 'detection_method': How sequence type was determined
        """
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)

        # =================================================================
        # STEP 1: MULTI-SCENARIO SEQUENCE TYPE DETECTION
        # =================================================================

        detection_result = self._detect_sequence_type(seq_str, seq_len)

        seq_type = detection_result['type']
        detection_method = detection_result['method']
        max_ltr_similarity = detection_result['similarity']
        best_ltr_len = detection_result['ltr_length']

        # Log detection decision
        self.logger.debug(
            f"Sequence type detection for {seq_record.id}: "
            f"type={seq_type}, method={detection_method}, "
            f"similarity={max_ltr_similarity:.2f}, ltr_len={best_ltr_len}bp, "
            f"total_len={seq_len}bp"
        )

        # =================================================================
        # STEP 2: APPLY APPROPRIATE SCORING STRATEGY
        # =================================================================

        if seq_type == 'complete_element':
            # Pass detected LTR info to avoid redundant calculation
            result = self._score_complete_element(
                seq_record, alignments, genome_seqs,
                detected_ltr_similarity=max_ltr_similarity,
                detected_ltr_length=best_ltr_len
            )
        else:  # 'single_ltr'
            result = self._score_single_ltr(seq_record)

        # Add detection metadata to result
        result['detection_method'] = detection_method
        result['terminal_similarity'] = max_ltr_similarity

        return result


    def _detect_sequence_type(self, seq_str, seq_len):
        """
        Multi-scenario strategy for detecting whether sequence is a single LTR
        or a complete LTR element based on terminal repeat similarity.

        STRATEGY OVERVIEW:

        Scenario 1: Very short sequences (<300bp)
        - Decision: Single LTR (too short for complete element)
        - Method: Length-based

        Scenario 2: Short sequences (300-800bp)
        - Check terminal similarity with restricted LTR length range
        - Decision: Based on similarity (threshold: 35%)

        Scenario 3: Medium sequences (800-3000bp)
        - Check terminal similarity with standard LTR length range
        - Decision: Based on similarity (threshold: 40%)

        Scenario 4: Long sequences (>3000bp)
        - Check terminal similarity with extended LTR length range
        - Decision: Based on similarity (threshold: 40%)

        Args:
            seq_str: Uppercase sequence string
            seq_len: Length of sequence

        Returns:
            Dictionary with:
            - 'type': 'single_ltr' or 'complete_element'
            - 'method': Detection method used
            - 'similarity': Maximum LTR similarity found
            - 'ltr_length': Length of best matching LTR pair
        """

        # =================================================================
        # SCENARIO 1: Very Short Sequences (<300bp)
        # =================================================================
        if seq_len < 300:
            # Too short to be a complete element (would need LTR + internal region)
            # Minimum complete element: 100bp LTR + 100bp internal + 100bp LTR = 300bp
            return {
                'type': 'single_ltr',
                'method': 'length_based_too_short',
                'similarity': 0.0,
                'ltr_length': 0
            }

        # =================================================================
        # SCENARIO 2: Short Sequences (300-800bp)
        # =================================================================
        elif seq_len < 800:
            # Restricted LTR length range for short sequences
            # LTR can be 20-40% of total length
            max_ltr_similarity = 0.0
            best_ltr_len = 0

            for ltr_ratio in [0.20, 0.25, 0.30, 0.35, 0.40]:
                test_ltr_len = int(seq_len * ltr_ratio)

                # Ensure LTR length is reasonable
                if test_ltr_len < 50 or test_ltr_len > min(400, seq_len // 2):
                    continue

                left_ltr = seq_str[:test_ltr_len]
                right_ltr = seq_str[-test_ltr_len:]
                similarity = self._calculate_similarity(left_ltr, right_ltr)

                if similarity > max_ltr_similarity:
                    max_ltr_similarity = similarity
                    best_ltr_len = test_ltr_len

            # Decision threshold: 35% for short sequences (more lenient)
            if max_ltr_similarity >= 0.35:
                return {
                    'type': 'complete_element',
                    'method': 'similarity_based_short',
                    'similarity': max_ltr_similarity,
                    'ltr_length': best_ltr_len
                }
            else:
                return {
                    'type': 'single_ltr',
                    'method': 'similarity_based_short',
                    'similarity': max_ltr_similarity,
                    'ltr_length': best_ltr_len
                }

        # =================================================================
        # SCENARIO 3: Medium Sequences (800-3000bp)
        # =================================================================
        elif seq_len < 3000:
            # Standard LTR length range
            # LTR can be 10-40% of total length, constrained by min/max_ltr_len
            max_ltr_similarity = 0.0
            best_ltr_len = 0

            for ltr_ratio in [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
                test_ltr_len = int(seq_len * ltr_ratio)

                # Apply global LTR length constraints
                if test_ltr_len < self.min_ltr_len or test_ltr_len > self.max_ltr_len:
                    continue

                left_ltr = seq_str[:test_ltr_len]
                right_ltr = seq_str[-test_ltr_len:]
                similarity = self._calculate_similarity(left_ltr, right_ltr)

                if similarity > max_ltr_similarity:
                    max_ltr_similarity = similarity
                    best_ltr_len = test_ltr_len

            # Decision threshold: 40% for medium sequences (standard)
            if max_ltr_similarity >= 0.40:
                return {
                    'type': 'complete_element',
                    'method': 'similarity_based_medium',
                    'similarity': max_ltr_similarity,
                    'ltr_length': best_ltr_len
                }
            else:
                return {
                    'type': 'single_ltr',
                    'method': 'similarity_based_medium',
                    'similarity': max_ltr_similarity,
                    'ltr_length': best_ltr_len
                }

        # =================================================================
        # SCENARIO 4: Long Sequences (>3000bp)
        # =================================================================
        else:
            # Extended LTR length range for long sequences
            # Try broader range of LTR ratios
            max_ltr_similarity = 0.0
            best_ltr_len = 0

            for ltr_ratio in [0.05, 0.08, 0.10, 0.15, 0.20, 0.25, 0.30]:
                test_ltr_len = int(seq_len * ltr_ratio)

                # Apply global LTR length constraints
                if test_ltr_len < self.min_ltr_len or test_ltr_len > self.max_ltr_len:
                    continue

                left_ltr = seq_str[:test_ltr_len]
                right_ltr = seq_str[-test_ltr_len:]
                similarity = self._calculate_similarity(left_ltr, right_ltr)

                if similarity > max_ltr_similarity:
                    max_ltr_similarity = similarity
                    best_ltr_len = test_ltr_len

            # Decision threshold: 40% for long sequences (standard)
            if max_ltr_similarity >= 0.40:
                return {
                    'type': 'complete_element',
                    'method': 'similarity_based_long',
                    'similarity': max_ltr_similarity,
                    'ltr_length': best_ltr_len
                }
            else:
                return {
                    'type': 'single_ltr',
                    'method': 'similarity_based_long',
                    'similarity': max_ltr_similarity,
                    'ltr_length': best_ltr_len
                }

    def _score_single_ltr(self, seq_record):
        """
        Identify single LTR consensus sequence (100-1500 bp).

        Uses 2 core features (cannot assess LTR similarity/TSD for single LTR):
        - Boundary motifs (70%): TG...CA pattern - strongest LTR signature
        - Length (30%): Typical LTR range 100-1500 bp

        Returns:
            Dictionary with identification information
        """
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)

        scores = {}
        issues = []

        # === Core Feature: Boundary Motifs (Weight: 0.70) ===
        boundary_result = self.validate_ltr_boundaries(seq_str, strict=False)

        if boundary_result['valid']:
            motif_score = boundary_result['confidence'] * 100
        else:
            motif_score = 0.0
            issues.append(f"Missing LTR boundary motifs: {boundary_result['reason']}")

        scores['boundary_motifs'] = motif_score

        # === Core Feature: Length (Weight: 0.30) ===
        optimal_min = 100
        optimal_max = 1500

        if optimal_min <= seq_len <= optimal_max:
            length_score = 100.0
        elif 50 <= seq_len < optimal_min:
            length_score = 70.0  # Short but possible
        elif optimal_max < seq_len < 2000:
            length_score = 85.0  # Slightly long but acceptable
        else:
            length_score = 0.0
            if seq_len < 50:
                issues.append(f"Too short: {seq_len}bp")
            else:
                issues.append(f"Too long: {seq_len}bp")

        scores['length'] = length_score

        # === Calculate Total Score ===
        weights = {
            'boundary_motifs': 0.70,  # Primary identifier
            'length': 0.30
        }

        total_score = sum(scores[key] * weights[key] for key in scores)

        # === Classify Sequence (Identification, not quality) ===
        if total_score >= 70 and motif_score >= 80:
            category = 'confirmed'  # Strong LTR signature
        elif total_score >= 50:
            category = 'probable'   # Likely LTR
        elif total_score >= 30:
            category = 'possible'   # Weak evidence
        else:
            category = 'not_ltr'    # Does not meet criteria

        # Log identification result
        self.logger.info(f"LTR identification: {seq_record.id} = {category} (score={total_score:.1f})")
        self.logger.debug(f"  Boundary: {motif_score:.1f}, Length: {length_score:.1f}")

        return {
            'total_score': total_score,
            'category': category,
            'components': scores,
            'issues': issues,
            'sequence_type': 'single_ltr',
            'ltr_similarity': None,
            'ltr_length': seq_len,
            'tsd': None,
            'pbs': None,
            'ppt': None
        }

    def _score_complete_element(self, seq_record, alignments=None, genome_seqs=None,
                                 detected_ltr_similarity=None, detected_ltr_length=None):
        """
        Identify complete LTR retrotransposon element (>2000 bp).

        Uses 3 core identification features:
        - LTR similarity (50%): Two terminal repeats
        - TSD (35%): Target site duplication
        - Boundary motifs (15%): TG...CA pattern

        Args:
            seq_record: SeqRecord object
            alignments: Optional alignment data for TSD detection
            genome_seqs: Optional genome sequences for TSD extraction
            detected_ltr_similarity: Pre-calculated LTR similarity from _detect_sequence_type()
            detected_ltr_length: Pre-calculated LTR length from _detect_sequence_type()

        Returns:
            Dictionary with identification information
        """
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)

        scores = {}
        issues = []

        # === Core Feature 1: LTR Similarity (Weight: 0.50) ===
        # Use pre-calculated values if available to avoid redundant computation
        if detected_ltr_similarity is not None and detected_ltr_length is not None:
            best_ltr_sim = detected_ltr_similarity
            best_ltr_len = detected_ltr_length
            self.logger.debug(f"Using pre-detected LTR similarity: {best_ltr_sim:.3f} ({best_ltr_len}bp)")
        else:
            # Fallback: calculate LTR similarity if not provided
            best_ltr_sim = 0.0
            best_ltr_len = 0

            for ltr_ratio in [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]:
                test_ltr_len = int(seq_len * ltr_ratio)
                if test_ltr_len < self.min_ltr_len or test_ltr_len > self.max_ltr_len:
                    continue

                left_ltr = seq_str[:test_ltr_len]
                right_ltr = seq_str[-test_ltr_len:]
                similarity = self._calculate_similarity(left_ltr, right_ltr)

                if similarity > best_ltr_sim:
                    best_ltr_sim = similarity
                    best_ltr_len = test_ltr_len

        # Score LTR similarity with gradual scaling
        if best_ltr_sim >= 0.90:
            ltr_sim_score = 100.0
        elif best_ltr_sim >= 0.80:
            ltr_sim_score = 95.0
        elif best_ltr_sim >= 0.70:
            ltr_sim_score = 85.0
        elif best_ltr_sim >= 0.50:
            ltr_sim_score = 70.0  # Ancient but valid
        elif best_ltr_sim >= 0.40:
            ltr_sim_score = 50.0
        else:
            ltr_sim_score = 0.0
            issues.append(f"Low LTR similarity: {best_ltr_sim:.2f}")

        scores['ltr_similarity'] = ltr_sim_score

        # === Core Feature 2: TSD (Weight: 0.35) ===
        tsd_score = 0.0
        tsd_found = None

        if alignments and genome_seqs:
            tsd_evidence = self._find_tsd_evidence(alignments, genome_seqs)
            if tsd_evidence and 'tsd' in tsd_evidence:
                tsd_found = tsd_evidence['tsd']
                tsd_len = len(tsd_found)

                if self.tsd_min <= tsd_len <= self.tsd_max:
                    tsd_score = 100.0
                elif tsd_len > 0:
                    tsd_score = 70.0
        else:
            if seq_len > 20:
                left_flank = seq_str[:10]
                right_flank = seq_str[-10:]

                for tsd_len in range(self.tsd_min, min(self.tsd_max + 1, 10)):
                    if left_flank[:tsd_len] == right_flank[-tsd_len:]:
                        tsd_score = 80.0
                        tsd_found = left_flank[:tsd_len]
                        break

        scores['tsd'] = tsd_score

        if tsd_score == 0:
            issues.append("No TSD detected")

        # === Core Feature 3: Boundary Motifs (Weight: 0.15) ===
        has_5_motif = any(seq_str.startswith(motif) for motif in self.ltr5_motifs)
        has_3_motif = any(seq_str.endswith(motif) for motif in self.ltr3_motifs)

        if has_5_motif and has_3_motif:
            motif_score = 100.0
        elif has_5_motif or has_3_motif:
            motif_score = 50.0
        else:
            motif_score = 0.0
            issues.append("Missing terminal motifs")

        scores['boundary_motifs'] = motif_score

        # === Calculate Total Score ===
        weights = {
            'ltr_similarity': 0.50,
            'tsd': 0.35,
            'boundary_motifs': 0.15
        }

        total_score = sum(scores[key] * weights[key] for key in scores)

        # === Classify (Identification) ===
        if total_score >= 70 and best_ltr_sim >= 0.80:
            category = 'confirmed'
        elif total_score >= 60 and best_ltr_sim >= 0.70:
            category = 'confirmed'
        elif total_score >= 50 and best_ltr_sim >= 0.50:
            category = 'probable'  # Ancient LTR
        elif total_score >= 35:
            category = 'possible'
        else:
            category = 'not_ltr'

        # Log result
        self.logger.info(f"LTR identification: {seq_record.id} = {category} (score={total_score:.1f})")
        self.logger.debug(f"  LTR sim={best_ltr_sim:.2f} ({best_ltr_len}bp), TSD={tsd_found}, Motifs={has_5_motif}/{has_3_motif}")

        return {
            'total_score': total_score,
            'category': category,
            'components': scores,
            'issues': issues,
            'sequence_type': 'complete_element',
            'ltr_similarity': best_ltr_sim,
            'ltr_length': best_ltr_len,
            'tsd': tsd_found,
            'pbs': None,
            'ppt': None
        }

    def _classify_sequences_by_quality(self, fasta_file, alignments_by_query, genome_seqs):
        """
        Identify LTR sequences and classify by confidence.

        Args:
            fasta_file: Path to input FASTA file
            alignments_by_query: Dictionary mapping query IDs to alignments
            genome_seqs: Genome access object

        Returns:
            Dictionary with:
            - 'complete': List of confirmed+probable LTR sequences
            - 'needs_trimming': Empty (not used in identification mode)
            - 'incomplete': List of possible LTR sequences
            - 'quality_info': Dictionary mapping seq_id to identification info
        """
        sequences = list(SeqIO.parse(fasta_file, "fasta"))

        classification = {
            'complete': [],          # confirmed + probable LTRs
            'needs_trimming': [],    # Not used in identification mode
            'incomplete': [],        # possible LTRs
            'quality_info': {}
        }

        self.logger.info(f"Identifying LTR sequences: {len(sequences)} candidates")

        # Track statistics
        category_stats = {
            'confirmed': 0,
            'probable': 0,
            'possible': 0,
            'not_ltr': 0
        }
        score_ranges = {'0-35': 0, '35-50': 0, '50-70': 0, '70-100': 0}

        for seq_record in sequences:
            seq_id = seq_record.id
            alignments = alignments_by_query.get(seq_id, [])

            # Identify the sequence
            id_info = self._score_ltr_structure(seq_record, alignments, genome_seqs)
            classification['quality_info'][seq_id] = id_info

            # Track statistics
            category = id_info['category']
            total_score = id_info['total_score']
            category_stats[category] = category_stats.get(category, 0) + 1

            # Track score distribution
            if total_score < 35:
                score_ranges['0-35'] += 1
            elif total_score < 50:
                score_ranges['35-50'] += 1
            elif total_score < 70:
                score_ranges['50-70'] += 1
            else:
                score_ranges['70-100'] += 1

            # Classify based on identification confidence
            if category in ['confirmed', 'probable']:
                # High-confidence LTRs
                classification['complete'].append(seq_record)
            elif category == 'possible':
                # Low-confidence LTRs (borderline)
                classification['incomplete'].append(seq_record)
            else:  # not_ltr
                # Not an LTR - filter out
                self.logger.debug(f"Not an LTR: {seq_id} (score={total_score:.1f})")
                pass

        # Log detailed statistics
        total_retained = len(classification['complete']) + len(classification['incomplete'])
        total_filtered = category_stats['not_ltr']

        self.logger.info(f"Identification results: {total_retained} LTRs identified, {total_filtered} filtered")
        self.logger.info(f"  Confirmed: {category_stats['confirmed']} ({category_stats['confirmed']*100/len(sequences):.1f}%)")
        self.logger.info(f"  Probable: {category_stats['probable']} ({category_stats['probable']*100/len(sequences):.1f}%)")
        self.logger.info(f"  Possible: {category_stats['possible']} ({category_stats['possible']*100/len(sequences):.1f}%)")
        self.logger.info(f"  Not LTR (filtered): {category_stats['not_ltr']} ({category_stats['not_ltr']*100/len(sequences):.1f}%)")

        self.logger.info(f"Score distribution: "
                        f"[0-35]={score_ranges['0-35']}, "
                        f"[35-50]={score_ranges['35-50']}, "
                        f"[50-70]={score_ranges['50-70']}, "
                        f"[70-100]={score_ranges['70-100']}")

        return classification

    def _is_high_quality_sequence(self, seq_str, min_length=100, max_n_content=0.05):
        """
        Check if sequence meets quality criteria.

        Args:
            seq_str: Sequence string
            min_length: Minimum acceptable length
            max_n_content: Maximum proportion of N's allowed

        Returns:
            Boolean indicating if sequence is high quality
        """
        if len(seq_str) < min_length:
            return False

        n_count = seq_str.upper().count('N')
        n_proportion = n_count / len(seq_str) if len(seq_str) > 0 else 1.0

        return n_proportion <= max_n_content
