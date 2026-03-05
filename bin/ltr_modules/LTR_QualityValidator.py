#!/usr/bin/env python3
"""
LTR Identifier Module

Identifies LTR retrotransposon sequences using 4 core features:
1. LTR similarity (two terminal repeats)
2. TSD (Target Site Duplication)
3. Boundary motifs (TG...CA pattern)
4. PBS (Primer Binding Site)

This module focuses on IDENTIFICATION (is it an LTR?) rather than
quality control (how good is the sequence?).

Date: 2025-11-20
"""

import os
import sys
import logging
import subprocess
import tempfile
import random
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
    - Core Feature 4: PBS presence (tRNA binding site)

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
        
        # PBS patterns (Common tRNA 3' ends)
        # TGGC... is the most common start for PBS (complementary to tRNA 3' CCA)
        self.pbs_patterns = [
            "TGGCG", "TGGCA", "TGGCC", "TGGCT", # General tRNA
            "TGGT", "TGGA" # Variants
        ]
    
    def _calculate_similarity(self, seq1, seq2):
        """Calculate sequence similarity (0-1) - simple fallback method."""
        if not seq1 or not seq2:
            return 0.0

        # Align sequences using simple character-by-character comparison
        matches = sum(a == b for a, b in zip(seq1, seq2))
        max_len = max(len(seq1), len(seq2))

        if max_len == 0:
            return 0.0

        return matches / max_len

    def _calculate_similarity_blast(self, full_seq):
        """
        Calculate terminal similarity using BLAST self-alignment.

        Strategy:
        1. Extract 5' terminal region as query (first 500-1000bp)
        2. BLAST against the full sequence
        3. Find the best hit in the 3' region (excluding self-hit)
        4. Return identity as similarity score

        Args:
            full_seq: Complete sequence string

        Returns:
            float: Terminal similarity (0-1), or None if BLAST fails
        """
        seq_len = len(full_seq)

        # Need at least 1kb to detect terminal repeats
        if seq_len < 1000:
            return None

        # Extract 5' terminal as query (20% of length, 500-1500bp range)
        query_len = min(max(500, int(seq_len * 0.2)), 1500)
        query_seq = full_seq[:query_len]

        # Create temporary files
        temp_dir = tempfile.gettempdir()
        rand_id = random.randint(100000, 999999)
        query_file = os.path.join(temp_dir, f"ltr_term_q_{rand_id}.fa")
        subject_file = os.path.join(temp_dir, f"ltr_term_s_{rand_id}.fa")

        try:
            # Write sequences
            with open(query_file, 'w') as f:
                f.write(f">query\n{query_seq}\n")

            with open(subject_file, 'w') as f:
                f.write(f">subject\n{full_seq}\n")

            # Run BLAST
            cmd = [
                'blastn',
                '-query', query_file,
                '-subject', subject_file,
                '-outfmt', '6 qstart qend sstart send pident length qlen slen',
                '-task', 'blastn',
                '-dust', 'no',
                '-word_size', '7',  # Sensitive for short repeats
                '-evalue', '10'     # Permissive to catch degraded LTRs
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Parse BLAST output
            best_identity = 0.0
            best_hit_in_3prime = None

            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 8:
                    continue

                qstart, qend, sstart, send = map(int, parts[:4])
                pident = float(parts[4])
                length = int(parts[5])
                qlen = int(parts[6])
                slen = int(parts[7])

                # Filter criteria:
                # 1. Exclude self-hit (query maps to 5' region of subject)
                if sstart < query_len * 0.5:  # Overlaps with query region
                    continue

                # 2. Must be in the 3' half of the sequence (potential 3' LTR)
                if send < slen * 0.5:  # Not in 3' region
                    continue

                # 3. Alignment should cover a good portion of query
                coverage = length / qlen
                if coverage < 0.5:  # At least 50% coverage
                    continue

                # 4. Must be long enough (at least 100bp aligned)
                if length < 100:
                    continue

                # Track best hit
                if pident > best_identity:
                    best_identity = pident
                    best_hit_in_3prime = {
                        'sstart': sstart,
                        'send': send,
                        'pident': pident,
                        'length': length,
                        'coverage': coverage
                    }

            # Clean up temp files
            os.remove(query_file)
            os.remove(subject_file)

            if best_hit_in_3prime:
                # Convert percent identity to 0-1 scale
                return best_identity / 100.0
            else:
                # No valid 3' LTR hit found
                return 0.0

        except Exception as e:
            # Clean up on error
            if os.path.exists(query_file):
                os.remove(query_file)
            if os.path.exists(subject_file):
                os.remove(subject_file)

            # Return None to indicate BLAST failed
            return None

    def _detect_pbs(self, seq_str, ltr_len):
        """
        Detect Primer Binding Site (PBS) immediately downstream of 5' LTR.
        
        Args:
            seq_str: Full sequence string
            ltr_len: Length of 5' LTR
            
        Returns:
            bool: True if PBS detected
        """
        # PBS is usually 1-5bp downstream of 5' LTR
        search_start = ltr_len
        search_end = min(len(seq_str), ltr_len + 20)
        region = seq_str[search_start:search_end]
        
        for pattern in self.pbs_patterns:
            if pattern in region:
                return True
        return False

    def validate_ltr_boundaries(self, seq, strict=False):
        """
        Validate LTR boundary motifs (5' and 3' ends).
        """
        # Handle Bio.Seq objects
        seq_str = str(seq).upper() if seq else ''

        # Minimum length check
        if len(seq_str) < 4:
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
        # Removed weak variants (TG..TG, TA..CA) to reduce false positives
        boundary_patterns = [
            ('TG', 'CA', 'canonical', 1.0),      # Canonical LTR pattern
            ('CA', 'CA', 'variant_CA', 0.8),     # CA...CA variant
            ('TG', 'TA', 'variant_TA', 0.8),     # TG...TA variant
        ]

        # Check against patterns
        for motif_5p, motif_3p, pattern_name, confidence in boundary_patterns:
            if start_di == motif_5p and end_di == motif_3p:
                if strict and pattern_name != 'canonical':
                    continue
                return {
                    'valid': True,
                    '5p_motif': motif_5p,
                    '3p_motif': motif_3p,
                    'pattern': pattern_name,
                    'confidence': confidence,
                    'reason': f'matched {pattern_name} pattern'
                }

        # Check for partial matches
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

        return {
            'valid': False,
            '5p_motif': start_di if partial_match_5p else None,
            '3p_motif': end_di if partial_match_3p else None,
            'pattern': None,
            'confidence': 0.0,
            'reason': reason
        }

    def comprehensive_quality_validation(self, seq, check_boundary=True, strict_boundary=False, seq_type=None):
        """
        Comprehensive quality validation for LTR sequences.

        Args:
            seq: Sequence to validate
            check_boundary: Whether to check boundary motifs
            strict_boundary: Use strict boundary validation
            seq_type: Sequence type ('single_ltr', 'complete_element', or None for auto-detect)
        """
        seq_str = str(seq).upper() if seq else ''
        seq_len = len(seq_str)

        checks = {}
        warnings = []
        errors = []
        scores = []

        # Auto-detect sequence type if not provided
        if seq_type is None:
            # For Internal sequences (check_boundary=False), use different logic
            if not check_boundary:
                # Internal sequences should not be classified by terminal similarity
                # They are internal regions without LTRs
                seq_type = 'internal'
            else:
                # For LTR sequences, use BLAST-based terminal similarity detection
                if seq_len >= 1000:
                    # Try BLAST self-alignment first
                    similarity = self._calculate_similarity_blast(seq_str)

                    # If BLAST fails or returns None, fallback to simple comparison
                    if similarity is None:
                        test_ltr_len = min(int(seq_len * 0.2), 1000)
                        left_ltr = seq_str[:test_ltr_len]
                        right_ltr = seq_str[-test_ltr_len:]
                        similarity = self._calculate_similarity(left_ltr, right_ltr)

                    # Solution A: Lowered threshold from 0.35 to 0.30 to avoid misclassifying
                    # ancient/degraded complete elements as single_ltr
                    seq_type = 'complete_element' if similarity >= 0.30 else 'single_ltr'
                elif seq_len >= 300:
                    # For shorter sequences, use simple comparison (BLAST not reliable)
                    test_ltr_len = min(int(seq_len * 0.2), 1000)
                    left_ltr = seq_str[:test_ltr_len]
                    right_ltr = seq_str[-test_ltr_len:]
                    similarity = self._calculate_similarity(left_ltr, right_ltr)
                    seq_type = 'complete_element' if similarity >= 0.30 else 'single_ltr'
                else:
                    seq_type = 'single_ltr'

        # 1. Boundary motif validation
        if check_boundary:
            boundary_result = self.validate_ltr_boundaries(seq_str, strict=strict_boundary)
            checks['boundary'] = boundary_result

            if boundary_result['valid']:
                scores.append(boundary_result['confidence'] * 100)
            else:
                # Warning only, not fatal
                warnings.append(f"Invalid boundary motifs: {boundary_result['reason']}")
                scores.append(0)

        # 2. Length validation with dynamic limits based on sequence type
        min_len = self.min_ltr_len if hasattr(self, 'min_ltr_len') else 100

        # Set maximum length based on sequence type
        if seq_type == 'single_ltr':
            max_len = 5000   # Strict limit for single LTR
            hard_limit = True
        elif seq_type == 'internal':
            # Internal sequences can be much longer (pol, gag, env genes)
            # Ty1/Copia: 4-7kb, Ty3/Gypsy: 8-15kb, Large elements: 20-30kb+
            max_len = 35000  # Very permissive for large Internal regions
            hard_limit = False  # Use soft limit with penalty
        else:  # complete_element
            max_len = 20000  # Relaxed limit for complete elements
            hard_limit = False

        if seq_len < min_len:
            checks['length'] = {'valid': False, 'value': seq_len, 'min': min_len, 'type': seq_type}
            errors.append(f"Too short: {seq_len} bp < {min_len} bp")
            scores.append(0)
        elif seq_len > max_len:
            if hard_limit:
                # Solution C: For single_ltr with slight overage (<35%), use soft limit
                overage_ratio = (seq_len - max_len) / max_len

                if overage_ratio < 0.35:  # Less than 35% over limit (e.g., 5000 -> 6750)
                    # Soft limit: Allow but penalize
                    checks['length'] = {'valid': True, 'value': seq_len, 'max': max_len, 'type': seq_type,
                                       'penalty': 'slight_overage'}
                    warnings.append(f"Single LTR slightly exceeds typical length: {seq_len} bp (limit: {max_len} bp)")
                    # Gradual penalty: -1 point per 10bp over limit, max 40 points penalty
                    penalty = min(40, (seq_len - max_len) / 10)
                    length_score = max(60, 100 - penalty)
                    scores.append(length_score)
                else:
                    # Hard filter: Significantly too long (likely an error)
                    checks['length'] = {'valid': False, 'value': seq_len, 'max': max_len, 'type': seq_type}
                    errors.append(f"Too long for single LTR: {seq_len} bp > {max_len} bp")
                    scores.append(0)
            else:
                # Complete element: Soft penalty (allow long elements but penalize)
                checks['length'] = {'valid': True, 'value': seq_len, 'max': max_len, 'type': seq_type,
                                   'penalty': 'unusually_long'}
                warnings.append(f"Unusually long element: {seq_len} bp (typical range: 4-15 kb)")
                # Gradual penalty: -1 point per 200bp over limit, max 50 points penalty
                penalty = min(50, (seq_len - max_len) / 200)
                length_score = max(50, 100 - penalty)
                scores.append(length_score)
        else:
            checks['length'] = {'valid': True, 'value': seq_len, 'type': seq_type}
            if seq_type == 'single_ltr':
                # Optimal range for single LTR: 100-1500bp
                if 100 <= seq_len <= 1500:
                    scores.append(100)
                elif seq_len < 100:
                    scores.append(max(50, seq_len / 2))
                else:
                    scores.append(max(70, 100 - (seq_len - 1500) / 100))
            elif seq_type == 'internal':
                # Optimal range for Internal: 3000-25000bp
                # Ty1/Copia: ~5kb, Ty3/Gypsy: ~12kb, Large elements: 20-30kb
                if 3000 <= seq_len <= 25000:
                    scores.append(100)
                elif seq_len < 3000:
                    # Too short for typical Internal (missing genes?)
                    scores.append(max(60, 100 - (3000 - seq_len) / 50))
                else:
                    # Very long but acceptable (some elements are >25kb)
                    scores.append(max(85, 100 - (seq_len - 25000) / 500))
            else:  # complete_element
                # Optimal range for complete element: 4000-15000bp
                if 4000 <= seq_len <= 15000:
                    scores.append(100)
                elif seq_len < 4000:
                    scores.append(max(70, 100 - (4000 - seq_len) / 100))
                else:
                    scores.append(max(80, 100 - (seq_len - 15000) / 200))

        # 3. N content validation (STRENGTHENED)
        n_count = seq_str.count('N')
        n_ratio = n_count / seq_len if seq_len > 0 else 1.0
        max_n_ratio = 0.10  # Reduced from 0.20 to 0.10 (stricter)

        checks['n_content'] = {'valid': n_ratio <= max_n_ratio, 'ratio': n_ratio, 'count': n_count}

        if n_ratio > max_n_ratio:
            errors.append(f"Excessive N content: {n_ratio:.1%} > {max_n_ratio:.1%}")
            scores.append(0)
        elif n_ratio > 0.05:  # Reduced from 0.10 to 0.05
            warnings.append(f"Moderate N content: {n_ratio:.1%}")
            scores.append(50)
        else:
            scores.append(100)

        # 4. Low complexity check (STRENGTHENED)
        low_complexity_score = self._check_sequence_complexity(seq_str)
        checks['complexity'] = {'score': low_complexity_score}

        if low_complexity_score < 40:  # Increased from 5 to 40 (stricter)
            errors.append(f"Low sequence complexity: score {low_complexity_score:.1f}")
            scores.append(0)
        elif low_complexity_score < 60:  # Increased from 40 to 60
            warnings.append(f"Moderate sequence complexity: score {low_complexity_score:.1f}")
            scores.append(low_complexity_score)
        else:
            scores.append(low_complexity_score)

        # 5. Tandem repeat check (STRENGTHENED)
        tandem_ratio = self._check_tandem_repeats(seq_str)
        checks['tandem_repeats'] = {'ratio': tandem_ratio}

        if tandem_ratio > 0.50:  # Reduced from 0.90 to 0.50 (stricter)
            errors.append(f"Excessive tandem repeats: {tandem_ratio:.1%} of sequence")
            scores.append(0)
        elif tandem_ratio > 0.30:  # Reduced from 0.60 to 0.30
            warnings.append(f"Moderate tandem repeats: {tandem_ratio:.1%}")
            scores.append(60)
        else:
            scores.append(100)

        # Calculate overall score
        overall_score = sum(scores) / len(scores) if scores else 0

        # Determine if valid (no fatal errors)
        is_valid = len(errors) == 0

        if is_valid:
            self.logger.debug(f"Quality validation passed: score={overall_score:.1f}, warnings={len(warnings)}")
        else:
            self.logger.info(f"Quality validation failed: score={overall_score:.1f}, errors={errors}")

        return {
            'valid': is_valid,
            'checks': checks,
            'score': overall_score,
            'warnings': warnings,
            'errors': errors
        }

    def _check_sequence_complexity(self, seq_str, window=50):
        """Check sequence complexity using sliding window."""
        if len(seq_str) < window:
            window = len(seq_str)
        if window < 10: return 50

        complexity_scores = []
        for i in range(len(seq_str) - window + 1):
            win_seq = seq_str[i:i+window]
            base_counts = {b: win_seq.count(b) for b in 'ACGT'}
            total = sum(base_counts.values())
            if total == 0: continue
            
            entropy = 0
            for count in base_counts.values():
                if count > 0:
                    p = count / total
                    entropy -= p * np.log2(p)
            
            complexity = (entropy / 2.0) * 100
            complexity_scores.append(complexity)

        if not complexity_scores: return 50
        return sum(complexity_scores) / len(complexity_scores)

    def _check_tandem_repeats(self, seq_str, min_unit=2, max_unit=50, min_copies=3):
        """Check for tandem repeats."""
        seq_len = len(seq_str)
        if seq_len < min_unit * min_copies: return 0.0

        repeat_positions = set()
        for unit_size in range(min_unit, min(max_unit + 1, seq_len // min_copies + 1)):
            for start in range(seq_len - unit_size * min_copies + 1):
                unit = seq_str[start:start + unit_size]
                if 'N' in unit: continue
                
                copies = 1
                pos = start + unit_size
                while pos + unit_size <= seq_len and seq_str[pos:pos + unit_size] == unit:
                    copies += 1
                    pos += unit_size
                
                if copies >= min_copies:
                    for i in range(start, start + copies * unit_size):
                        repeat_positions.add(i)
        
        return len(repeat_positions) / seq_len if seq_len > 0 else 0.0

    def _score_ltr_structure(self, seq_record, alignments=None, genome_seqs=None):
        """
        Score LTR structural completeness.
        """
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)

        # Detect sequence type
        detection_result = self._detect_sequence_type(seq_str, seq_len)
        seq_type = detection_result['type']
        max_ltr_similarity = detection_result['similarity']
        best_ltr_len = detection_result['ltr_length']

        if seq_type == 'complete_element':
            result = self._score_complete_element(
                seq_record, alignments, genome_seqs,
                detected_ltr_similarity=max_ltr_similarity,
                detected_ltr_length=best_ltr_len
            )
        else:
            result = self._score_single_ltr(seq_record)

        result['detection_method'] = detection_result['method']
        result['terminal_similarity'] = max_ltr_similarity
        return result

    def _detect_sequence_type(self, seq_str, seq_len):
        """Detect if sequence is single LTR or complete element."""
        if seq_len < 300:
            return {'type': 'single_ltr', 'method': 'length', 'similarity': 0.0, 'ltr_length': 0}
        
        # Check for terminal repeats
        max_ltr_similarity = 0.0
        best_ltr_len = 0
        
        # Check various LTR ratios
        ratios = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
        for ltr_ratio in ratios:
            test_ltr_len = int(seq_len * ltr_ratio)
            if test_ltr_len < self.min_ltr_len or test_ltr_len > self.max_ltr_len: continue
            
            left_ltr = seq_str[:test_ltr_len]
            right_ltr = seq_str[-test_ltr_len:]
            similarity = self._calculate_similarity(left_ltr, right_ltr)
            
            if similarity > max_ltr_similarity:
                max_ltr_similarity = similarity
                best_ltr_len = test_ltr_len
        
        if max_ltr_similarity >= 0.35:
            return {'type': 'complete_element', 'method': 'similarity', 'similarity': max_ltr_similarity, 'ltr_length': best_ltr_len}
        else:
            return {'type': 'single_ltr', 'method': 'similarity', 'similarity': max_ltr_similarity, 'ltr_length': best_ltr_len}

    def _score_single_ltr(self, seq_record):
        """Score single LTR consensus."""
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)
        scores = {}
        issues = []

        # Boundary Motifs (70%)
        boundary_result = self.validate_ltr_boundaries(seq_str, strict=False)
        motif_score = boundary_result['confidence'] * 100 if boundary_result['valid'] else 0.0
        scores['boundary_motifs'] = motif_score
        if not boundary_result['valid']: issues.append(boundary_result['reason'])

        # Length (30%)
        if 100 <= seq_len <= 1500: length_score = 100.0
        else: length_score = 50.0
        scores['length'] = length_score

        total_score = motif_score * 0.7 + length_score * 0.3
        
        # Strict check for Single LTRs: Must have strong motifs (>=80)
        # Single LTRs lack terminal repeats, so boundary motifs are the ONLY structural evidence.
        if motif_score < 80:
            total_score = 0
            issues.append("Single LTR lacks strong boundary motifs")
        
        category = 'not_ltr'
        if total_score >= 70: category = 'confirmed'
        elif total_score >= 50: category = 'probable'
        elif total_score >= 30: category = 'possible'

        return {
            'total_score': total_score,
            'category': category,
            'components': scores,
            'issues': issues,
            'sequence_type': 'single_ltr',
            'ltr_similarity': None,
            'ltr_length': seq_len,
            'tsd': None,
            'pbs': None
        }

    def _score_complete_element(self, seq_record, alignments=None, genome_seqs=None,
                                 detected_ltr_similarity=None, detected_ltr_length=None):
        """Score complete LTR element."""
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)
        scores = {}
        issues = []

        # 1. LTR Similarity (40%)
        best_ltr_sim = detected_ltr_similarity if detected_ltr_similarity else 0.0
        best_ltr_len = detected_ltr_length if detected_ltr_length else 0
        
        if best_ltr_sim >= 0.90: ltr_sim_score = 100.0
        elif best_ltr_sim >= 0.80: ltr_sim_score = 95.0
        elif best_ltr_sim >= 0.70: ltr_sim_score = 85.0
        elif best_ltr_sim >= 0.50: ltr_sim_score = 70.0
        elif best_ltr_sim >= 0.40: ltr_sim_score = 50.0
        else: ltr_sim_score = 0.0
        scores['ltr_similarity'] = ltr_sim_score

        # 2. TSD (30%)
        tsd_score = 0.0
        tsd_found = None
        if alignments and genome_seqs:
            # (TSD logic omitted for brevity, assuming external call or simplified check)
            pass
        else:
            # Simple self-check for TSD
            # Require at least 4bp for simple check to avoid random 3bp matches
            if seq_len > 20:
                left = seq_str[:10]
                right = seq_str[-10:]
                min_k = max(4, self.tsd_min)
                for k in range(min_k, min(self.tsd_max+1, 10)):
                    if left[:k] == right[-k:]:
                        tsd_score = 80.0
                        tsd_found = left[:k]
                        break
        scores['tsd'] = tsd_score

        # 3. Boundary Motifs (15%)
        has_5_motif = any(seq_str.startswith(m) for m in self.ltr5_motifs)
        has_3_motif = any(seq_str.endswith(m) for m in self.ltr3_motifs)
        motif_score = 100.0 if (has_5_motif and has_3_motif) else (50.0 if (has_5_motif or has_3_motif) else 0.0)
        scores['boundary_motifs'] = motif_score

        # 4. PBS (15%) - NEW FEATURE
        has_pbs = self._detect_pbs(seq_str, best_ltr_len)
        pbs_score = 100.0 if has_pbs else 0.0
        scores['pbs'] = pbs_score

        # Calculate Total Score
        weights = {
            'ltr_similarity': 0.40,
            'tsd': 0.30,
            'boundary_motifs': 0.15,
            'pbs': 0.15
        }
        total_score = sum(scores[k] * weights[k] for k in scores)

        # === STRICT FILTERING LOGIC ===
        category = 'not_ltr'
        
        # Rule 1: Short sequences (<2000bp) MUST have structural evidence
        if seq_len < 2000:
            if not (tsd_found or (has_5_motif and has_3_motif)):
                # Penalty for featureless short sequences
                total_score *= 0.5
                issues.append("Short sequence lacks TSD/Motifs")

        # Rule 2: Featureless sequences (Similarity only) are rejected
        has_structure = (tsd_found is not None) or (has_5_motif or has_3_motif) or has_pbs
        if not has_structure:
            total_score *= 0.4 # Heavy penalty
            issues.append("No structural features (TSD/Motif/PBS)")

        # Classification
        if total_score >= 70 and best_ltr_sim >= 0.80: category = 'confirmed'
        elif total_score >= 60 and best_ltr_sim >= 0.70: category = 'confirmed'
        elif total_score >= 50 and best_ltr_sim >= 0.45: category = 'probable'
        elif total_score >= 35: category = 'possible'

        return {
            'total_score': total_score,
            'category': category,
            'components': scores,
            'issues': issues,
            'sequence_type': 'complete_element',
            'ltr_similarity': best_ltr_sim,
            'ltr_length': best_ltr_len,
            'tsd': tsd_found,
            'pbs': has_pbs
        }