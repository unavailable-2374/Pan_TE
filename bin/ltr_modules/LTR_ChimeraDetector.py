#!/usr/bin/env python3
"""
LTR Chimera Detector Module

Handles detection and splitting of chimeric LTR sequences.

This module contains all chimera-related functionality extracted from
LTR_Boundary_Optimizer.py for better code organization and reusability.
"""

import os
import re
import logging
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np


class LTRChimeraDetector:
    """
    LTR Chimera Detection and Splitting

    Handles all chimera-related tasks including:
    - Chimeric sequence detection
    - Pattern analysis
    - Breakpoint identification and splitting
    - Secondary chimera detection
    """

    def __init__(self, chimera_threshold=0.3, min_segment_length=500,
                 blank_region_threshold=100, orientation_aware=True,
                 logger=None):
        """
        Initialize the chimera detector.

        Args:
            chimera_threshold: Threshold for chimera detection (higher = more sensitive)
            min_segment_length: Minimum length for segments when splitting
            blank_region_threshold: Minimum length of blank regions for detection
            orientation_aware: Whether to use orientation-aware detection
            logger: Logger instance (optional)
        """
        self.chimera_threshold = chimera_threshold
        self.min_segment_length = min_segment_length
        self.blank_region_threshold = blank_region_threshold
        self.orientation_aware = orientation_aware

        # Setup logger
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger('LTRChimeraDetector')
            if not self.logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                self.logger.addHandler(handler)
                self.logger.setLevel(logging.INFO)

        # Chimera tracking
        self.chimera_info = {}

    # ========================================================================
    # PUBLIC API
    # ========================================================================

    def detect_chimeras(self, fasta_file, alignments_by_query):
        """
        PUBLIC API - Main entry point for chimera detection and splitting.

        Detects chimeric sequences and splits them at breakpoints.

        Args:
            fasta_file: Path to input FASTA file
            alignments_by_query: Dictionary mapping query IDs to their alignments

        Returns:
            dict: Results with keys:
                - 'split_records': List of split SeqRecord objects
                - 'chimera_info': Dictionary of chimera detection information
                - 'num_chimeras': Number of chimeric sequences detected
                - 'num_splits': Total number of splits performed

        Example:
            >>> detector = LTRChimeraDetector()
            >>> result = detector.detect_chimeras('consensus.fa', alignments)
            >>> split_seqs = result['split_records']
            >>> print(f"Found {result['num_chimeras']} chimeras")
        """
        self.logger.info(f"Starting chimera detection on {fasta_file}")
        self.logger.debug(f"  Chimera threshold: {self.chimera_threshold}")
        self.logger.debug(f"  Min segment length: {self.min_segment_length}bp")

        split_records, chimera_info = self._detect_and_split_chimeras(fasta_file, alignments_by_query)

        num_chimeras = sum(1 for info in chimera_info.values() if info.get('is_chimeric', False))
        num_splits = sum(info.get('num_segments', 1) - 1 for info in chimera_info.values())

        self.logger.info(f"Chimera detection complete: {num_chimeras} chimeras detected, {num_splits} splits performed")
        if num_chimeras > 0:
            self.logger.info(f"  Total sequences: {len(chimera_info)}, chimeric: {num_chimeras} ({num_chimeras*100/len(chimera_info):.1f}%)")

        return {
            'split_records': split_records,
            'chimera_info': chimera_info,
            'num_chimeras': num_chimeras,
            'num_splits': num_splits
        }

    def analyze_chimera_pattern(self, seq_record, alignments):
        """
        PUBLIC API - Analyze chimeric pattern of a sequence.

        Args:
            seq_record: BioPython SeqRecord object
            alignments: List of alignments for this sequence

        Returns:
            dict: Pattern analysis with keys:
                - 'is_chimeric': bool
                - 'breakpoints': List of breakpoint positions
                - 'confidence': Confidence score (0-1)
                - 'blank_regions': List of blank regions

        Example:
            >>> detector = LTRChimeraDetector()
            >>> pattern = detector.analyze_chimera_pattern(seq_record, alignments)
            >>> if pattern['is_chimeric']:
            >>>     print(f"Breakpoints: {pattern['breakpoints']}")
        """
        # Call private method and convert tuple to dict
        result = self._analyze_chimeric_pattern(seq_record, alignments)

        # The private method returns (is_chimeric, breakpoints, confidence, blank_regions)
        if isinstance(result, tuple) and len(result) == 4:
            is_chimeric, breakpoints, confidence, blank_regions = result
            return {
                'is_chimeric': is_chimeric,
                'breakpoints': breakpoints if breakpoints else [],
                'confidence': confidence,
                'blank_regions': blank_regions if blank_regions else []
            }
        else:
            # Fallback for unexpected return format
            return {
                'is_chimeric': False,
                'breakpoints': [],
                'confidence': 0.0,
                'blank_regions': []
            }

    def split_sequence(self, seq_record, breakpoints):
        """
        PUBLIC API - Split sequence at specified breakpoints.

        Args:
            seq_record: BioPython SeqRecord object
            breakpoints: List of breakpoint positions

        Returns:
            list: List of split SeqRecord objects

        Example:
            >>> detector = LTRChimeraDetector()
            >>> segments = detector.split_sequence(seq_record, [500, 1000])
        """
        return self._split_at_breakpoints(seq_record, breakpoints)

    # ========================================================================
    # PRIVATE METHODS
    # ========================================================================

    def _detect_and_split_chimeras(self, fasta_file, alignments_by_query):
        """
        Enhanced method to detect and split chimeric consensus sequences.
        Uses multiple evidence types:
        1. Coverage discontinuities
        2. Target chromosome changes
        3. Alignment orientation changes
        4. Blank regions (no alignment)
        5. Orientation-aware pattern analysis
        
        Args:
            fasta_file: Path to input FASTA file
            alignments_by_query: Dictionary mapping query IDs to their alignments
            
        Returns:
            (split_records, chimera_info) tuple
        """
        # Load sequences
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        
        # Track chimeric sequences and their segments
        split_records = []
        chimera_info = {}
        
        # Process each sequence
        for seq_record in sequences:
            seq_id = seq_record.id
            
            # Get alignments for this sequence
            seq_alignments = alignments_by_query.get(seq_id, [])
            
            # Analyze alignment pattern to detect chimeras
            is_chimeric, breakpoints, breakpoint_confidence, blank_regions = self._analyze_chimeric_pattern(seq_record, seq_alignments)
            
            # Record chimera information
            chimera_info[seq_id] = {
                'is_chimeric': is_chimeric,
                'breakpoints': breakpoints,
                'breakpoint_confidence': breakpoint_confidence,
                'blank_regions': blank_regions,
                'num_segments': len(breakpoints) + 1 if breakpoints else 1
            }
            
            if is_chimeric and breakpoints:
                # Split sequence at breakpoints
                segments = self._split_at_breakpoints(seq_record, breakpoints)
                split_records.extend(segments)
            else:
                # Not chimeric, keep original
                split_records.append(seq_record)
        
        return split_records, chimera_info
    

    def _analyze_chimeric_pattern(self, seq_record, alignments):
        """
        Enhanced method to analyze alignment patterns to detect chimeric sequences.
        Uses multiple evidence types:
        1. Coverage discontinuities
        2. Target chromosome changes
        3. Alignment orientation changes
        4. Blank regions (no alignment)
        5. Orientation-aware pattern analysis
        
        Args:
            seq_record: SeqRecord object
            alignments: List of alignments for this sequence
            
        Returns:
            (is_chimeric, breakpoints, breakpoint_confidence, blank_regions) tuple
        """
        seq_len = len(seq_record.seq)
        
        # Not enough alignments to analyze
        if len(alignments) < 3:
            return False, [], [], []
            
        # Create coverage profile by position
        coverage = np.zeros(seq_len)
        position_targets = [set() for _ in range(seq_len)]
        position_strands = [set() for _ in range(seq_len)]
        
        # Fill coverage profile
        for aln in alignments:
            start = aln['query_start']
            end = aln['query_end']
            
            for pos in range(start, end):
                if pos < seq_len:
                    coverage[pos] += 1
                    position_targets[pos].add(aln['target_name'])
                    position_strands[pos].add(aln['strand'])
        
        # Detect coverage discontinuities and target changes
        breakpoints = []
        breakpoint_confidence = []
        window_size = min(100, seq_len // 10)  # Adaptive window size
        if window_size < 20:  # Minimum window size
            window_size = 20
            
        threshold = self.chimera_threshold  # Threshold for significant change
        
        # Analyze for breakpoints
        for i in range(window_size, seq_len - window_size):
            # Check coverage difference
            left_avg = np.mean(coverage[i-window_size:i])
            right_avg = np.mean(coverage[i:i+window_size])
            
            # Check target composition difference
            left_targets = set().union(*position_targets[i-window_size:i])
            right_targets = set().union(*position_targets[i:i+window_size])
            
            # Check strand composition difference
            left_strands = set().union(*position_strands[i-window_size:i])
            right_strands = set().union(*position_strands[i:i+window_size])
            
            # Calculate target similarity
            target_similarity = 0
            if left_targets and right_targets:
                common_targets = left_targets.intersection(right_targets)
                target_similarity = len(common_targets) / max(len(left_targets), len(right_targets))
                
            # Check for strand changes
            strand_change = 0
            if left_strands and right_strands:
                if len(left_strands.symmetric_difference(right_strands)) > 0:
                    strand_change = 1
            
            # Combined metric
            coverage_diff = abs(left_avg - right_avg) / max(left_avg, right_avg) if max(left_avg, right_avg) > 0 else 0
            change_score = coverage_diff + (1 - target_similarity) + strand_change
            
            # Enhanced: Orientation-aware analysis if enabled
            if self.orientation_aware and strand_change > 0:
                # Check if this is a systematic inversion rather than a chimera
                # For systematic inversions, we'd expect consistent strand within regions
                # and similar target distributions
                if target_similarity > 0.7:  # Same targets but different orientation
                    # This might be an inversion rather than a chimera
                    # Reduce the change score
                    change_score *= 0.5
            
            if change_score > threshold:
                # Significant change detected
                breakpoints.append(i)
                breakpoint_confidence.append(change_score)
                
        # Detect blank regions (contiguous zero-coverage)
        blank_regions = []
        if self.blank_region_threshold > 0:
            in_blank = False
            blank_start = 0
            
            for i in range(seq_len):
                if coverage[i] == 0:
                    if not in_blank:
                        in_blank = True
                        blank_start = i
                else:
                    if in_blank:
                        in_blank = False
                        blank_end = i
                        blank_length = blank_end - blank_start
                        
                        if blank_length >= self.blank_region_threshold:
                            blank_regions.append((blank_start, blank_end))
                            
                            # Add as potential breakpoint if not already close to one
                            midpoint = (blank_start + blank_end) // 2
                            if not any(abs(bp - midpoint) < window_size for bp in breakpoints):
                                breakpoints.append(midpoint)
                                breakpoint_confidence.append(1.0)  # High confidence for blank regions
            
            # Check if sequence ends with a blank region
            if in_blank and seq_len - blank_start >= self.blank_region_threshold:
                blank_regions.append((blank_start, seq_len))
                
        # Merge nearby breakpoints
        if breakpoints:
            merged_breakpoints = []
            merged_confidence = []
            
            # Sort breakpoints and confidence together
            sorted_pairs = sorted(zip(breakpoints, breakpoint_confidence))
            breakpoints = [p[0] for p in sorted_pairs]
            breakpoint_confidence = [p[1] for p in sorted_pairs]
            
            current_group = [breakpoints[0]]
            current_conf = [breakpoint_confidence[0]]
            
            for i in range(1, len(breakpoints)):
                bp = breakpoints[i]
                conf = breakpoint_confidence[i]
                
                if bp - current_group[-1] < window_size:
                    # Add to current group
                    current_group.append(bp)
                    current_conf.append(conf)
                else:
                    # Start new group
                    # Use weighted average by confidence
                    weighted_bp = sum(b * c for b, c in zip(current_group, current_conf)) / sum(current_conf)
                    merged_breakpoints.append(int(weighted_bp))
                    merged_confidence.append(max(current_conf))  # Use max confidence
                    
                    current_group = [bp]
                    current_conf = [conf]
                    
            # Add last group
            if current_group:
                weighted_bp = sum(b * c for b, c in zip(current_group, current_conf)) / sum(current_conf)
                merged_breakpoints.append(int(weighted_bp))
                merged_confidence.append(max(current_conf))
                
            # Filter breakpoints by minimum segment length
            valid_breakpoints = []
            valid_confidence = []
            last_pos = 0
            
            breakpoint_pairs = sorted(zip(merged_breakpoints, merged_confidence))
            
            for bp, conf in breakpoint_pairs:
                if bp - last_pos >= self.min_segment_length:
                    valid_breakpoints.append(bp)
                    valid_confidence.append(conf)
                    last_pos = bp
                    
            # Check if last segment is long enough
            if valid_breakpoints and seq_len - valid_breakpoints[-1] < self.min_segment_length:
                valid_breakpoints.pop()
                valid_confidence.pop()
                
            # If we still have valid breakpoints, it's chimeric
            return len(valid_breakpoints) > 0, valid_breakpoints, valid_confidence, blank_regions
        
        return False, [], [], blank_regions
    

    def _split_at_breakpoints(self, seq_record, breakpoints):
        """
        Split a sequence at the given breakpoints.
        
        Args:
            seq_record: SeqRecord object
            breakpoints: List of positions to split at
            
        Returns:
            List of split SeqRecord objects
        """
        seq_len = len(seq_record.seq)
        segments = []
        
        # Sort breakpoints
        breakpoints = sorted(breakpoints)
        
        # Add start and end positions
        positions = [0] + breakpoints + [seq_len]
        
        # Create segments
        for i in range(len(positions) - 1):
            start = positions[i]
            end = positions[i + 1]
            
            # Create new record
            segment_id = f"{seq_record.id}_segment_{i+1}"
            segment_desc = f"Segment {i+1} from {seq_record.id} position {start+1}-{end}"
            
            segment = SeqRecord(
                seq_record.seq[start:end],
                id=segment_id,
                name=segment_id,
                description=segment_desc
            )
            
            segments.append(segment)

        return segments


    def _check_sequence_in_alignment(self, seq_str, alignment):
        """
        Check if a sequence appears in an alignment.

        Args:
            seq_str: Sequence string to search for
            alignment: Alignment dictionary

        Returns:
            Boolean indicating if sequence is covered by alignment
        """
        # Simple check: if sequence is within alignment boundaries
        return True  # Placeholder - could be more sophisticated


    def _secondary_chimera_detection(self, sequences, genome_seqs):
        """
        Perform secondary chimera detection after extension.

        This catches chimeras that may have been created or expanded during
        the aggressive extension phase.

        Args:
            sequences: List of sequence records (potentially extended)
            genome_seqs: Genome access object

        Returns:
            Tuple of (refined_sequences, chimera_info)
        """
        self.logger.info(f"Performing secondary chimera detection on {len(sequences)} sequences")

        # Write sequences to temporary file
        temp_file = os.path.join(self.temp_dir, f"secondary_chimera_check_{int(time.time())}.fa")
        SeqIO.write(sequences, temp_file, "fasta")

        # Align to genome
        alignment_file = os.path.join(self.temp_dir, f"secondary_chimera_{int(time.time())}.paf")
        self._align_to_genome(temp_file, alignment_file)

        # Filter alignments
        filtered_file = alignment_file.replace('.paf', '.filtered.paf')
        alignments_by_query = self._filter_alignments(alignment_file, filtered_file)

        # Detect and split chimeras with stricter threshold
        # Use lower threshold for secondary detection (more sensitive)
        original_threshold = self.chimera_threshold
        self.chimera_threshold = max(0.2, self.chimera_threshold - 0.1)  # More sensitive

        split_records, chimera_info = self._detect_and_split_chimeras(temp_file, alignments_by_query)

        # Restore original threshold
        self.chimera_threshold = original_threshold

        # Count how many were split
        num_chimeras = sum(1 for info in chimera_info.values() if info['is_chimeric'])

        self.logger.info(f"Secondary chimera detection: {num_chimeras} chimeras detected and split")

        return split_records, chimera_info


