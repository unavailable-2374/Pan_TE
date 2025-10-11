#!/usr/bin/env python3
"""
LTR Boundary Refiner Module

Handles LTR boundary detection, TSD identification, and boundary extension.

This module contains all boundary-related functionality extracted from
LTR_Boundary_Optimizer.py for better code organization and reusability.
"""

import os
import re
import math
import logging
import subprocess
import tempfile
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from itertools import groupby


class LTRBoundaryRefiner:
    """
    LTR Boundary Detection and Refinement

    Handles all boundary-related tasks including:
    - Boundary detection using multiple evidence sources
    - TSD (Target Site Duplication) identification
    - Boundary extension and optimization
    - Boundary validation and trimming
    """

    def __init__(self, min_ltr_len=100, max_ltr_len=1500,
                 tsd_min=4, tsd_max=6, flanking_seq=20,
                 kmer_boundary=True, advanced_tsd=True,
                 custom_motifs=None, weighted_evidence=True,
                 genome_file=None, temp_dir=None,
                 logger=None):
        """
        Initialize the boundary refiner.

        Args:
            min_ltr_len: Minimum LTR length
            max_ltr_len: Maximum LTR length
            tsd_min: Minimum TSD length
            tsd_max: Maximum TSD length
            flanking_seq: Length of flanking sequence for TSD analysis
            kmer_boundary: Whether to use k-mer based boundary detection
            advanced_tsd: Whether to use advanced TSD detection
            custom_motifs: Custom motifs for boundary detection
            weighted_evidence: Whether to use weighted evidence integration
            genome_file: Reference genome file for realignment (optional)
            temp_dir: Temporary directory for alignment files (optional)
            logger: Logger instance (optional)
        """
        self.min_ltr_len = min_ltr_len
        self.max_ltr_len = max_ltr_len
        self.tsd_min = tsd_min
        self.tsd_max = tsd_max
        self.flanking_seq = flanking_seq
        self.kmer_boundary = kmer_boundary
        self.advanced_tsd = advanced_tsd
        self.custom_motifs = custom_motifs
        self.weighted_evidence = weighted_evidence
        self.genome_file = genome_file
        self.temp_dir = temp_dir if temp_dir else tempfile.gettempdir()

        # Setup logger
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger('LTRBoundaryRefiner')
            if not self.logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                self.logger.addHandler(handler)
                self.logger.setLevel(logging.INFO)

        # LTR-specific patterns
        self.ltr5_motifs = ["TG", "TGT", "TGTA", "TGTG", "TGCA"]
        self.ltr3_motifs = ["CA", "ACA", "TACA", "CACA", "TGCA"]

        # Add custom motifs if provided
        if self.custom_motifs:
            if 'ltr5' in self.custom_motifs and isinstance(self.custom_motifs['ltr5'], list):
                self.ltr5_motifs.extend(self.custom_motifs['ltr5'])
            if 'ltr3' in self.custom_motifs and isinstance(self.custom_motifs['ltr3'], list):
                self.ltr3_motifs.extend(self.custom_motifs['ltr3'])

        # TSD patterns storage
        self.tsd_patterns = {}
        self.tsd_frequency = Counter()

        # TSD nucleotide bias from empirical data
        self.tsd_nucleotide_bias = {
            'A': 0.30,
            'C': 0.20,
            'G': 0.20,
            'T': 0.30
        }

        # PBS patterns - expanded for different organisms
        self.pbs_patterns = {
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

        # Define DEFAULT evidence weights for boundary detection
        # These will be adjusted dynamically based on sequence type
        self.default_evidence_weights = {
            'terminal_motifs': 0.30,
            'alignment_boundaries': 0.25,
            'kmer_transitions': 0.20,
            'tsd_evidence': 0.15,
            'internal_features': 0.10
        }

        # Evidence weights for SINGLE LTR sequences
        # Higher weight on terminal motifs (only direct evidence available)
        self.single_ltr_weights = {
            'terminal_motifs': 0.70,      # PRIMARY: Only clear structural signal
            'alignment_boundaries': 0.15,  # Secondary: Alignment consistency
            'kmer_transitions': 0.10,      # Tertiary: Composition changes
            'tsd_evidence': 0.03,          # Minimal: Usually not present
            'internal_features': 0.02      # Minimal: No internal structure expected
        }

        # Evidence weights for COMPLETE LTR elements
        # Lower weight on terminal motifs (can rely on LTR similarity)
        self.complete_element_weights = {
            'terminal_motifs': 0.15,       # Tertiary: Nice to have but not critical
            'alignment_boundaries': 0.35,  # PRIMARY: LTR similarity via alignments
            'kmer_transitions': 0.20,      # Secondary: LTR/internal transitions
            'tsd_evidence': 0.20,          # Secondary: Insertion site evidence
            'internal_features': 0.10      # Tertiary: PBS/PPT support
        }

        # Current active weights (will be set based on sequence type)
        self.evidence_weights = self.default_evidence_weights.copy()

    # ========================================================================
    # PUBLIC API
    # ========================================================================

    def refine_boundaries(self, seq_record, alignments=None, genome_seqs=None, iteration=1,
                         sequence_type=None):
        """
        PUBLIC API - Main entry point for boundary refinement.

        Detects and refines LTR boundaries using multiple evidence sources.
        Uses ADAPTIVE strategy based on sequence type.

        Args:
            seq_record: BioPython SeqRecord object
            alignments: List of alignments for this sequence (optional)
            genome_seqs: Genome access object (optional)
            iteration: Current iteration number (default: 1)
            sequence_type: Sequence type ('single_ltr' or 'complete_element', optional)
                          If provided, adapts boundary detection strategy accordingly

        Returns:
            dict: Boundary information with keys:
                - 'refined_record': Refined SeqRecord
                - 'boundaries': Boundary coordinates
                - 'tsd': TSD information
                - 'quality_score': Quality score (0-100)
                - 'method': Detection method used

        Example:
            >>> refiner = LTRBoundaryRefiner()
            >>> # For single LTR (relies heavily on terminal motifs)
            >>> result = refiner.refine_boundaries(seq_record, alignments, genome_seqs,
            ...                                    sequence_type='single_ltr')
            >>> # For complete element (relies on LTR similarity)
            >>> result = refiner.refine_boundaries(seq_record, alignments, genome_seqs,
            ...                                    sequence_type='complete_element')
        """
        # Handle None alignments
        if alignments is None:
            alignments = []

        seq_len = len(seq_record.seq)

        # === ADAPTIVE STRATEGY: Adjust evidence weights based on sequence type ===
        if sequence_type == 'single_ltr':
            self.evidence_weights = self.single_ltr_weights.copy()
            strategy_name = "Single LTR (motif-focused)"
            self.logger.debug(f"{seq_record.id}: Using SINGLE_LTR strategy (motif weight: 0.70)")
        elif sequence_type == 'complete_element':
            self.evidence_weights = self.complete_element_weights.copy()
            strategy_name = "Complete element (similarity-focused)"
            self.logger.debug(f"{seq_record.id}: Using COMPLETE_ELEMENT strategy (alignment weight: 0.35)")
        else:
            self.evidence_weights = self.default_evidence_weights.copy()
            strategy_name = "Default (balanced)"
            self.logger.debug(f"{seq_record.id}: Using DEFAULT strategy (balanced weights)")

        self.logger.debug(f"Refining boundaries for {seq_record.id} ({seq_len}bp), "
                         f"iteration {iteration}, strategy: {strategy_name}")

        # Detect boundaries
        boundary_info = self._detect_ltr_boundaries(seq_record, alignments, genome_seqs, iteration)

        if not boundary_info:
            self.logger.debug(f"  No boundaries detected for {seq_record.id}")
            return {
                'refined_record': seq_record,
                'boundaries': None,
                'tsd': None,
                'quality_score': 0,
                'method': 'no_boundaries_detected'
            }

        # Refine sequence
        refined_record = self._refine_sequence(seq_record, boundary_info)

        quality_score = boundary_info.get('confidence_score', 0)
        method = boundary_info.get('method', 'integrated')
        self.logger.debug(f"  Boundary refinement complete: quality={quality_score:.1f}, method={method}")

        # Log TSD information if available
        tsd_info = boundary_info.get('tsd')
        if tsd_info and tsd_info.get('found'):
            self.logger.debug(f"  TSD detected: {tsd_info.get('sequence', 'N/A')} "
                            f"({tsd_info.get('length', 0)}bp, confidence={tsd_info.get('confidence', 0):.2f})")

        return {
            'refined_record': refined_record,
            'boundaries': boundary_info,
            'tsd': tsd_info,
            'quality_score': quality_score,
            'method': method
        }

    def detect_tsd(self, left_flank, right_flank, strand='+'):
        """
        PUBLIC API - Detect Target Site Duplication (TSD).

        Args:
            left_flank: Left flanking sequence
            right_flank: Right flanking sequence
            strand: Strand orientation ('+' or '-')

        Returns:
            dict: TSD information with keys:
                - 'found': bool
                - 'sequence': TSD sequence
                - 'length': TSD length
                - 'similarity': Similarity score (0-1)
                - 'confidence': Confidence score (0-1)

        Example:
            >>> refiner = LTRBoundaryRefiner()
            >>> tsd_info = refiner.detect_tsd(left_seq, right_seq)
            >>> if tsd_info['found']:
            >>>     print(f"TSD: {tsd_info['sequence']}")
        """
        if self.advanced_tsd:
            tsd_result = self._find_enhanced_tsd(left_flank, right_flank, strand)
            # _find_enhanced_tsd returns a dict or None
            if tsd_result is None or not isinstance(tsd_result, dict):
                tsd_result = {
                    'found': False,
                    'sequence': '',
                    'length': 0,
                    'similarity': 0.0,
                    'confidence': 0.0
                }
            else:
                # Ensure 'found' key exists
                if 'found' not in tsd_result:
                    tsd_result['found'] = tsd_result.get('sequence', '') != ''
        else:
            # Use simple TSD detection
            tsd_result = {
                'found': False,
                'sequence': '',
                'length': 0,
                'similarity': 0.0,
                'confidence': 0.0
            }

        return tsd_result

    def extend_boundaries(self, seq_record, boundaries, max_extension=500):
        """
        PUBLIC API - Extend boundaries to capture complete LTR elements.

        Args:
            seq_record: BioPython SeqRecord object
            boundaries: Current boundary coordinates
            max_extension: Maximum extension in bp (default: 500)

        Returns:
            dict: Extended boundary information

        Example:
            >>> refiner = LTRBoundaryRefiner()
            >>> extended = refiner.extend_boundaries(seq_record, boundaries)
        """
        return self._aggressive_extension(seq_record, boundaries, max_extension)

    # ========================================================================
    # PRIVATE METHODS
    # ========================================================================

    def _detect_ltr_boundaries(self, seq_record, alignments, genome_seqs, iteration=1):
        """
        Enhanced method to detect LTR boundaries using multiple evidence types
        with an integrative approach that doesn't rely on PSSM.

        TWO-STAGE STRATEGY:
        - Iteration 1: Aggressive extension (500bp each side) to capture true boundaries
        - Iteration 2+: Precise refinement using multiple evidence types

        Args:
            seq_record: SeqRecord object
            alignments: List of alignments for this sequence
            genome_seqs: Genome access object
            iteration: Current iteration number (default: 1)

        Returns:
            Dictionary with boundary information or None
        """
        seq_len = len(seq_record.seq)
        seq_id = seq_record.id

        # Not enough alignments
        if len(alignments) < 2:
            return None

        # Collect evidence from multiple sources
        evidence = {
            'terminal_motifs': self._find_terminal_motifs(seq_record),
            'kmer_boundaries': self._detect_kmer_boundaries(seq_record) if self.kmer_boundary else (None, None),
            'alignment_boundaries': self._find_alignment_boundaries(alignments),
            'tsd_evidence': self._find_tsd_evidence(alignments, genome_seqs),
            'internal_features': self._find_internal_features_enhanced(seq_record)
        }

        # Integrate evidence using weighted approach
        if self.weighted_evidence:
            boundaries = self._integrate_boundary_evidence(seq_record, evidence)
        else:
            # Fallback to simpler approach
            boundaries = self._simple_boundary_detection(seq_record, evidence)

        if not boundaries:
            return None

        # === TWO-STAGE BOUNDARY EXTENSION ===
        original_boundaries = boundaries.copy()

        if iteration == 1:
            # STAGE 1 (Iteration 1): AGGRESSIVE EXTENSION
            # Goal: Ensure we capture the true boundaries even if input is incomplete
            self.logger.debug(f"{seq_id}: Iteration 1 - Applying aggressive 500bp extension")
            boundaries = self._aggressive_extension(seq_record, boundaries, max_extension=500)
        else:
            # STAGE 2 (Iteration 2+): PRECISE REFINEMENT
            # Goal: Fine-tune boundaries on the extended sequence from iteration 1
            self.logger.debug(f"{seq_id}: Iteration {iteration} - Applying precise refinement")

            # STAGE 1.5 (NEW - Only on iteration 2): REALIGNMENT AND DEPTH-BASED DETECTION
            # After aggressive extension in iteration 1, realign to genome and use depth profile
            # to get better boundary estimates before applying precise refinement
            if iteration == 2 and self.genome_file:
                self.logger.info(f"{seq_id}: Stage 1.5 - Realigning extended sequence to genome")

                # Realign the extended sequence to genome
                realignments = self._realign_extended_sequence(seq_record)

                if realignments and len(realignments) > 0:
                    # Create depth profile from realignments
                    depth_profile = self._create_depth_profile(realignments, seq_len)

                    # Detect boundaries from depth profile
                    depth_boundaries = self._detect_boundaries_from_depth(
                        depth_profile,
                        current_boundaries=boundaries,
                        min_depth_ratio=0.3,
                        smoothing_window=51
                    )

                    # If depth-based detection succeeded, use it as starting point
                    if depth_boundaries and depth_boundaries != boundaries:
                        self.logger.info(f"{seq_id}: Using depth-based boundaries as starting point for refinement")
                        self.logger.debug(f"  Original: {boundaries['start']}-{boundaries['end']}")
                        self.logger.debug(f"  Depth-based: {depth_boundaries['start']}-{depth_boundaries['end']}")
                        boundaries = depth_boundaries
                    else:
                        self.logger.debug(f"{seq_id}: Depth-based detection did not improve boundaries, using original")
                else:
                    self.logger.debug(f"{seq_id}: Realignment did not produce valid alignments, skipping depth-based detection")

            # Step 1: Terminal motif-guided precise extension (highest priority, relaxed mode)
            boundaries = self._extend_boundaries_by_motifs(seq_record, boundaries, max_extension=100)

            # Step 2: Conservative extension based on alignment distribution
            boundaries = self._extend_boundaries_by_alignment_distribution(boundaries, alignments)

            # Step 3: Safe extension based on complexity (fill small gaps)
            boundaries = self._safe_extend_by_complexity(seq_record, boundaries, max_extension=30)

            # Step 4: For complete LTR elements, optimize using symmetry
            if self._is_potential_complete_ltr(boundaries):
                boundaries = self._extend_by_ltr_symmetry(seq_record, boundaries)

            # Step 5: Validate extension is reasonable
            boundaries = self._validate_extension(seq_record, original_boundaries,
                                                 boundaries, alignments, genome_seqs)

        # Validate final boundaries
        if boundaries['start'] >= boundaries['end']:
            return None

        # Check if boundaries define a reasonable LTR length
        ltr_length = boundaries['end'] - boundaries['start']
        if ltr_length < self.min_ltr_len or ltr_length > self.max_ltr_len:
            return None

        # Add additional information
        boundaries['seq_id'] = seq_id
        boundaries['original_length'] = seq_len

        return boundaries


    def _find_terminal_motifs(self, seq_record):
        """
        Find LTR terminal motifs (TG...CA patterns) in the sequence.
        
        Args:
            seq_record: SeqRecord object
            
        Returns:
            Dictionary with terminal motif information
        """
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)
        
        # Search for 5' motifs in the first 20 bp
        start_motifs = []
        for motif in self.ltr5_motifs:
            for i in range(min(20, seq_len - len(motif))):
                if seq_str[i:i+len(motif)] == motif:
                    start_motifs.append({
                        'motif': motif,
                        'position': i,
                        'confidence': 1.0 - (i / 20.0)  # Higher confidence for motifs closer to start
                    })
        
        # Search for 3' motifs in the last 20 bp
        end_motifs = []
        for motif in self.ltr3_motifs:
            for i in range(max(0, seq_len - 20 - len(motif)), seq_len - len(motif)):
                if seq_str[i:i+len(motif)] == motif:
                    end_motifs.append({
                        'motif': motif,
                        'position': i + len(motif),  # End position after motif
                        'confidence': 1.0 - ((seq_len - (i + len(motif))) / 20.0)  # Higher confidence for motifs closer to end
                    })
        
        return {
            'start_motifs': start_motifs,
            'end_motifs': end_motifs
        }


    def _find_alignment_boundaries(self, alignments):
        """
        Find potential boundaries based on alignment start/end positions.
        
        Args:
            alignments: List of alignments for this sequence
            
        Returns:
            Dictionary with alignment boundary information
        """
        if not alignments:
            return {'starts': [], 'ends': []}
            
        starts = [aln['query_start'] for aln in alignments]
        ends = [aln['query_end'] for aln in alignments]
        
        # Calculate frequency of each start/end position
        start_counter = Counter(starts)
        end_counter = Counter(ends)
        
        # Convert to format with position and confidence
        start_positions = []
        for pos, count in start_counter.items():
            confidence = count / len(alignments)
            start_positions.append({
                'position': pos,
                'count': count,
                'confidence': confidence
            })
            
        end_positions = []
        for pos, count in end_counter.items():
            confidence = count / len(alignments)
            end_positions.append({
                'position': pos,
                'count': count,
                'confidence': confidence
            })
            
        # Sort by confidence (highest first)
        start_positions.sort(key=lambda x: x['confidence'], reverse=True)
        end_positions.sort(key=lambda x: x['confidence'], reverse=True)
        
        return {
            'starts': start_positions,
            'ends': end_positions
        }


    def _detect_kmer_boundaries(self, seq_record, k=5, window_size=50):
        """
        Detect potential LTR boundaries based on k-mer frequency transitions.
        
        Args:
            seq_record: SeqRecord object
            k: k-mer size
            window_size: Size of sliding window for k-mer analysis
            
        Returns:
            (5'_boundary, 3'_boundary) tuple or (None, None) if not found
        """
        seq_str = str(seq_record.seq)
        seq_len = len(seq_str)
        
        # Too short for meaningful analysis
        if seq_len < window_size * 2 + k:
            return None, None
            
        # Compute k-mer frequencies in sliding windows
        window_kmers = []
        
        for i in range(0, seq_len - window_size + 1):
            window = seq_str[i:i+window_size]
            kmers = {}
            
            for j in range(len(window) - k + 1):
                kmer = window[j:j+k]
                kmers[kmer] = kmers.get(kmer, 0) + 1
                
            window_kmers.append(kmers)
            
        # Compute Jaccard distance between adjacent windows
        distances = []
        
        for i in range(1, len(window_kmers)):
            kmers1 = set(window_kmers[i-1].keys())
            kmers2 = set(window_kmers[i].keys())
            
            union = len(kmers1.union(kmers2))
            intersection = len(kmers1.intersection(kmers2))
            
            if union > 0:
                distance = 1 - (intersection / union)
            else:
                distance = 1.0
                
            distances.append(distance)
            
        # Find peaks in the distance profile
        # These represent potential boundaries with significant k-mer composition changes
        peaks = []
        min_peak_height = 0.3  # Minimum peak height
        
        for i in range(1, len(distances) - 1):
            if (distances[i] > distances[i-1] and 
                distances[i] > distances[i+1] and 
                distances[i] >= min_peak_height):
                
                # Adjust index to center of window
                peak_pos = i + window_size // 2
                peaks.append((peak_pos, distances[i]))
                
        # No significant peaks found
        if not peaks:
            return None, None
            
        # Sort peaks by height (highest first)
        peaks.sort(key=lambda x: x[1], reverse=True)
        
        # Look for peaks that could correspond to 5' and 3' boundaries
        # Typically, we need two peaks separated by a reasonable distance
        boundary5 = None
        boundary3 = None
        
        # Try different peak pairs
        for i in range(len(peaks)):
            pos1, _ = peaks[i]
            
            for j in range(i+1, len(peaks)):
                pos2, _ = peaks[j]
                
                # Make sure pos1 < pos2
                if pos1 > pos2:
                    pos1, pos2 = pos2, pos1
                    
                # Check if distance between peaks is reasonable for an LTR
                if self.min_ltr_len <= pos2 - pos1 <= self.max_ltr_len:
                    boundary5 = pos1
                    boundary3 = pos2
                    break
                    
            if boundary5 is not None:
                break
                
        return boundary5, boundary3


    def _integrate_boundary_evidence(self, seq_record, evidence):
        """
        Integrate multiple evidence types to determine LTR boundaries.
        Uses weighted evidence approach without relying on PSSM.
        
        Args:
            seq_record: SeqRecord object
            evidence: Dictionary with evidence from multiple sources
            
        Returns:
            Dictionary with boundary information or None
        """
        seq_len = len(seq_record.seq)
        seq_str = str(seq_record.seq)
        
        # Check for minimum evidence
        if not evidence['alignment_boundaries']['starts'] and not evidence['terminal_motifs']['start_motifs']:
            return None
            
        # 1. Process START boundary evidence
        start_candidates = []
        
        # Add terminal motif evidence
        for motif_info in evidence['terminal_motifs']['start_motifs']:
            start_candidates.append({
                'position': motif_info['position'],
                'source': 'terminal_motif',
                'weight': self.evidence_weights['terminal_motifs'],
                'confidence': motif_info['confidence'],
                'motif': motif_info['motif']
            })
            
        # Add alignment boundary evidence
        for aln_info in evidence['alignment_boundaries']['starts']:
            start_candidates.append({
                'position': aln_info['position'],
                'source': 'alignment',
                'weight': self.evidence_weights['alignment_boundaries'],
                'confidence': aln_info['confidence']
            })
            
        # Add k-mer transition evidence
        if evidence['kmer_boundaries'][0] is not None:
            start_candidates.append({
                'position': evidence['kmer_boundaries'][0],
                'source': 'kmer',
                'weight': self.evidence_weights['kmer_transitions'],
                'confidence': 0.8  # Reasonable confidence for k-mer boundaries
            })
            
        # 2. Process END boundary evidence
        end_candidates = []
        
        # Add terminal motif evidence
        for motif_info in evidence['terminal_motifs']['end_motifs']:
            end_candidates.append({
                'position': motif_info['position'],
                'source': 'terminal_motif',
                'weight': self.evidence_weights['terminal_motifs'],
                'confidence': motif_info['confidence'],
                'motif': motif_info['motif']
            })
            
        # Add alignment boundary evidence
        for aln_info in evidence['alignment_boundaries']['ends']:
            end_candidates.append({
                'position': aln_info['position'],
                'source': 'alignment',
                'weight': self.evidence_weights['alignment_boundaries'],
                'confidence': aln_info['confidence']
            })
            
        # Add k-mer transition evidence
        if evidence['kmer_boundaries'][1] is not None:
            end_candidates.append({
                'position': evidence['kmer_boundaries'][1],
                'source': 'kmer',
                'weight': self.evidence_weights['kmer_transitions'],
                'confidence': 0.8  # Reasonable confidence for k-mer boundaries
            })
            
        # 3. Handle TSD evidence (can affect both boundaries)
        if evidence['tsd_evidence']:
            # Find best TSD evidence
            best_tsd = max(evidence['tsd_evidence'], key=lambda x: x['score'])
            
            # TSD can provide indirect evidence for boundaries
            # in typical insertions, TSDs flank the element
            tsd_weight = self.evidence_weights['tsd_evidence']
            
            # Since we don't know exact mapping of TSD to the element boundaries,
            # we just increase confidence for boundaries near the ends
            for candidate in start_candidates:
                if candidate['position'] < 20:  # Near start
                    candidate['confidence'] = min(1.0, candidate['confidence'] + 0.1)
                    candidate['tsd_support'] = True
                    
            for candidate in end_candidates:
                if candidate['position'] > seq_len - 20:  # Near end
                    candidate['confidence'] = min(1.0, candidate['confidence'] + 0.1)
                    candidate['tsd_support'] = True
        
        # 4. Use internal features to adjust boundary confidence
        if 'pbs' in evidence['internal_features']:
            pbs_pos = evidence['internal_features']['pbs']['position']
            # PBS should be just after 5' LTR
            for candidate in start_candidates:
                # If PBS is close to candidate start position
                if 5 <= pbs_pos - candidate['position'] <= 20:
                    candidate['confidence'] = min(1.0, candidate['confidence'] + 0.1)
                    candidate['pbs_support'] = True
                    
        if 'ppt' in evidence['internal_features']:
            ppt_pos = evidence['internal_features']['ppt']['position']
            # PPT should be just before 3' LTR
            for candidate in end_candidates:
                # If PPT is close to candidate end position
                if 5 <= candidate['position'] - ppt_pos <= 20:
                    candidate['confidence'] = min(1.0, candidate['confidence'] + 0.1)
                    candidate['ppt_support'] = True
        
        # 5. Process candidate boundaries to find best start/end combination
        if not start_candidates or not end_candidates:
            return None
            
        # Calculate weighted scores for candidates
        for candidate in start_candidates:
            candidate['score'] = candidate['weight'] * candidate['confidence']
            
        for candidate in end_candidates:
            candidate['score'] = candidate['weight'] * candidate['confidence']
            
        # Sort by score (highest first)
        start_candidates.sort(key=lambda x: x['score'], reverse=True)
        end_candidates.sort(key=lambda x: x['score'], reverse=True)
        
        # Find the best combination that gives a valid LTR length
        best_start = None
        best_end = None
        best_combined_score = 0
        
        for start_candidate in start_candidates:
            start_pos = start_candidate['position']
            
            for end_candidate in end_candidates:
                end_pos = end_candidate['position']
                
                # Check if this gives a valid LTR length
                if end_pos <= start_pos or end_pos - start_pos < self.min_ltr_len or end_pos - start_pos > self.max_ltr_len:
                    continue
                    
                # Calculate combined score
                combined_score = start_candidate['score'] + end_candidate['score']
                
                # Apply bonus for consistent evidence types
                if start_candidate['source'] == end_candidate['source']:
                    combined_score *= 1.1  # 10% bonus
                    
                # Check terminal motifs
                start_motif = start_candidate.get('motif')
                end_motif = end_candidate.get('motif')
                
                # Apply bonus for canonical TG...CA pattern
                if start_motif and end_motif:
                    if start_motif.startswith('TG') and end_motif.endswith('CA'):
                        combined_score *= 1.2  # 20% bonus
                
                # Check if better than current best
                if combined_score > best_combined_score:
                    best_combined_score = combined_score
                    best_start = start_candidate
                    best_end = end_candidate
        
        # No valid combination found
        if not best_start or not best_end:
            return None
            
        # 6. Create final boundary information
        start_feature = best_start.get('motif')
        end_feature = best_end.get('motif')
        
        # Find TSD from evidence
        tsd = None
        if evidence['tsd_evidence']:
            # Get best TSD
            best_tsd_evidence = max(evidence['tsd_evidence'], key=lambda x: x['score'])
            tsd = {
                'found': True,
                'sequence': best_tsd_evidence.get('sequence', ''),
                'length': len(best_tsd_evidence.get('sequence', '')),
                'confidence': best_tsd_evidence.get('score', 0.0)
            }
            
        # Calculate confidence values
        start_confidence = best_start['confidence']
        end_confidence = best_end['confidence']
        
        # Check for k-mer boundary support
        kmer_support = (best_start['source'] == 'kmer' or best_end['source'] == 'kmer')
        
        return {
            'start': best_start['position'],
            'end': best_end['position'],
            'start_confidence': start_confidence,
            'end_confidence': end_confidence,
            'start_feature': start_feature,
            'end_feature': end_feature,
            'tsd': tsd,
            'kmer_boundary_support': kmer_support,
            'internal_features': evidence['internal_features']
        }


    def _simple_boundary_detection(self, seq_record, evidence):
        """
        Simplified boundary detection when weighted evidence integration is disabled.
        
        Args:
            seq_record: SeqRecord object
            evidence: Dictionary with evidence from multiple sources
            
        Returns:
            Dictionary with boundary information or None
        """
        seq_len = len(seq_record.seq)
        
        # Get alignment boundaries (most reliable source)
        start_positions = evidence['alignment_boundaries']['starts']
        end_positions = evidence['alignment_boundaries']['ends']
        
        # Not enough evidence
        if not start_positions or not end_positions:
            return None
            
        # Use most common alignment boundaries
        best_start = start_positions[0]['position']
        best_end = end_positions[0]['position']
        
        # Try to refine using terminal motifs
        if evidence['terminal_motifs']['start_motifs']:
            start_motif = evidence['terminal_motifs']['start_motifs'][0]
            # If motif is close to alignment boundary, use it
            if abs(start_motif['position'] - best_start) <= 20:
                best_start = start_motif['position']
                start_feature = start_motif['motif']
            else:
                start_feature = None
        else:
            start_feature = None
            
        if evidence['terminal_motifs']['end_motifs']:
            end_motif = evidence['terminal_motifs']['end_motifs'][0]
            # If motif is close to alignment boundary, use it
            if abs(end_motif['position'] - best_end) <= 20:
                best_end = end_motif['position']
                end_feature = end_motif['motif']
            else:
                end_feature = None
        else:
            end_feature = None
        
        # Get TSD if available
        tsd = None
        if evidence['tsd_evidence']:
            best_tsd = max(evidence['tsd_evidence'], key=lambda x: x['score'])
            tsd = {
                'found': True,
                'sequence': best_tsd.get('sequence', ''),
                'length': len(best_tsd.get('sequence', '')),
                'confidence': best_tsd.get('score', 0.0)
            }
        
        # Calculate confidence values
        start_confidence = start_positions[0]['confidence']
        end_confidence = end_positions[0]['confidence']
        
        # Check for k-mer boundary support
        kmer5, kmer3 = evidence['kmer_boundaries']
        kmer_support = False
        
        if kmer5 is not None and abs(kmer5 - best_start) <= 30:
            kmer_support = True
            
        if kmer3 is not None and abs(kmer3 - best_end) <= 30:
            kmer_support = True
            
        return {
            'start': best_start,
            'end': best_end,
            'start_confidence': start_confidence,
            'end_confidence': end_confidence,
            'start_feature': start_feature,
            'end_feature': end_feature,
            'tsd': tsd,
            'kmer_boundary_support': kmer_support,
            'internal_features': evidence['internal_features']
        }


    def _find_tsd_evidence(self, alignments, genome_seqs):
        """
        Enhanced TSD detection using local sequence context and alignment patterns
        instead of PSSM.

        Args:
            alignments: List of alignments for this sequence
            genome_seqs: Genome access object (optional)

        Returns:
            List of TSD evidence entries
        """
        tsd_evidence = []

        # If no genome sequences available, return empty list
        if not genome_seqs:
            return tsd_evidence

        # Process each alignment
        for aln in alignments:
            target_name = aln['target_name']
            target_start = aln['target_start']
            target_end = aln['target_end']
            strand = aln['strand']

            # Check if we can access the target sequence
            if not genome_seqs.has_sequence(target_name):
                continue

            # Extract flanking sequences
            left_flank_start = max(0, target_start - self.flanking_seq)
            left_flank = genome_seqs.get_sequence_region(target_name, left_flank_start, target_start)
            
            right_flank_end = target_end + self.flanking_seq
            right_flank = genome_seqs.get_sequence_region(target_name, target_end, right_flank_end)
            
            if not left_flank or not right_flank:
                continue
            
            # Find potential TSDs using advanced methods
            tsd_result = self._find_enhanced_tsd(left_flank, right_flank, strand)
            
            if tsd_result:
                tsd_evidence.append(tsd_result)
                
                # Learn from this TSD by tracking its pattern
                self.tsd_frequency[tsd_result['sequence']] += 1
        
        return tsd_evidence


    def _find_enhanced_tsd(self, left_flank, right_flank, strand):
        """
        Enhanced method to find Target Site Duplication (TSD) without using PSSM.
        Uses multiple techniques:
        1. Exact match search
        2. Approximate match with mismatches
        3. Common TSD pattern matching
        4. Nucleotide bias incorporation
        
        Args:
            left_flank: Left flanking sequence
            right_flank: Right flanking sequence
            strand: Strand of the alignment
            
        Returns:
            Dictionary with TSD information or None
        """
        # Check if flanking sequences are long enough
        if len(left_flank) < self.tsd_min or len(right_flank) < self.tsd_min:
            return None
            
        # Convert to uppercase
        left_flank = left_flank.upper()
        right_flank = right_flank.upper()
        
        # 1. Look for exact match TSDs
        best_exact_tsd = None
        best_exact_score = 0
        
        for tsd_len in range(self.tsd_min, min(self.tsd_max + 1, min(len(left_flank), len(right_flank)) + 1)):
            left_tsd = left_flank[-tsd_len:]
            right_tsd = right_flank[:tsd_len]
            
            if left_tsd == right_tsd:
                # Calculate base composition score
                # TSDs with balanced base composition are more reliable
                bases = Counter(left_tsd)
                composition_entropy = 0
                for base, count in bases.items():
                    p = count / tsd_len
                    composition_entropy -= p * math.log2(p) if p > 0 else 0
                
                # Normalize entropy to 0-1 range (max entropy for 4 bases is 2)
                composition_score = min(1.0, composition_entropy / 2.0)
                
                # Score based on TSD length (longer exact matches are better)
                length_factor = min(1.0, tsd_len / self.tsd_max)
                
                # Calculate bias-adjusted score using nucleotide preferences
                bias_score = 0
                for base in left_tsd:
                    if base in self.tsd_nucleotide_bias:
                        bias_score += self.tsd_nucleotide_bias[base]
                bias_score /= tsd_len
                
                # Combined score
                exact_score = 0.7 + (0.1 * composition_score) + (0.1 * length_factor) + (0.1 * bias_score)
                
                if exact_score > best_exact_score:
                    best_exact_score = exact_score
                    best_exact_tsd = {
                        'sequence': left_tsd,
                        'length': tsd_len,
                        'type': 'exact',
                        'score': exact_score,
                        'position': (len(left_flank) - tsd_len, 0)  # Positions in left and right flanks
                    }
        
        # 2. Look for approximate match TSDs (allow mismatches)
        best_approx_tsd = None
        best_approx_score = 0
        
        for tsd_len in range(self.tsd_min, min(self.tsd_max + 1, min(len(left_flank), len(right_flank)) + 1)):
            left_tsd = left_flank[-tsd_len:]
            right_tsd = right_flank[:tsd_len]
            
            # Count matches
            matches = sum(1 for a, b in zip(left_tsd, right_tsd) if a == b)
            match_ratio = matches / tsd_len
            
            # Only consider if match ratio is high enough
            if match_ratio >= 0.7:
                # Calculate mismatch-aware consensus TSD
                consensus_tsd = ""
                for a, b in zip(left_tsd, right_tsd):
                    if a == b:
                        consensus_tsd += a
                    else:
                        # Choose based on nucleotide bias
                        a_bias = self.tsd_nucleotide_bias.get(a, 0.25)
                        b_bias = self.tsd_nucleotide_bias.get(b, 0.25)
                        consensus_tsd += a if a_bias > b_bias else b
                
                # Adjust score based on match quality and length
                approx_score = match_ratio * 0.6  # Base score from match ratio
                
                # Bonus for longer TSD
                length_factor = min(1.0, tsd_len / self.tsd_max)
                approx_score += length_factor * 0.2
                
                # Bonus for common TSD patterns
                if self.tsd_frequency and consensus_tsd in self.tsd_frequency:
                    frequency_bonus = min(0.2, self.tsd_frequency[consensus_tsd] * 0.05)
                    approx_score += frequency_bonus
                
                if approx_score > best_approx_score:
                    best_approx_score = approx_score
                    best_approx_tsd = {
                        'sequence': consensus_tsd,
                        'length': tsd_len,
                        'type': 'approximate',
                        'score': approx_score,
                        'match_ratio': match_ratio,
                        'position': (len(left_flank) - tsd_len, 0)
                    }
        
        # 3. Prefer exact match over approximate if it's good
        if best_exact_tsd and best_exact_score >= 0.8:
            return best_exact_tsd
        
        # 4. Use approximate match if it's reasonable
        if best_approx_tsd and best_approx_score >= 0.7:
            return best_approx_tsd
        
        # 5. Check if there are common known TSDs from the learned patterns
        if self.tsd_frequency:
            for tsd, freq in self.tsd_frequency.most_common(5):  # Check top 5 common TSDs
                tsd_len = len(tsd)
                
                # Skip if too long for flanking sequences
                if tsd_len > min(len(left_flank), len(right_flank)):
                    continue
                    
                left_match = left_flank[-tsd_len:].count(tsd)
                right_match = right_flank[:tsd_len].count(tsd)
                
                if left_match > 0 and right_match > 0:
                    return {
                        'sequence': tsd,
                        'length': tsd_len,
                        'type': 'common_pattern',
                        'score': 0.7,  # Lower confidence for pattern match
                        'position': (len(left_flank) - tsd_len, 0)
                    }
        
        # No valid TSD found
        return None


    def _aggressive_extension(self, seq_record, boundaries, max_extension=500):
        """
        Aggressive boundary extension for iteration 1.
        Unconditionally extends both boundaries by max_extension to ensure true boundaries are captured.

        Strategy:
        1. Extend 5' boundary by max_extension (or to sequence start)
        2. Extend 3' boundary by max_extension (or to sequence end)
        3. NO motif-based early stopping - full extension for recall
        4. Iteration 2 will refine based on motifs and other evidence

        This ensures incomplete input sequences are extended enough to capture true boundaries.

        Args:
            seq_record: SeqRecord object
            boundaries: Dictionary with current boundary information
            max_extension: Maximum extension in bp (default: 500)

        Returns:
            Updated boundaries dictionary with aggressive extension applied
        """
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)
        current_start = boundaries['start']
        current_end = boundaries['end']

        # Calculate extension limits (full extension, no early stopping)
        new_start = max(0, current_start - max_extension)
        new_end = min(seq_len, current_end + max_extension)

        # Calculate actual extension amounts
        extension_5prime = current_start - new_start
        extension_3prime = new_end - current_end

        # Apply full extensions
        boundaries['start'] = new_start
        boundaries['end'] = new_end
        boundaries['aggressive_5prime_extension'] = extension_5prime
        boundaries['aggressive_3prime_extension'] = extension_3prime

        # Log extensions
        if extension_5prime > 0:
            self.logger.info(f"Aggressively extended 5' boundary by {extension_5prime} bp (full extension)")
        if extension_3prime > 0:
            self.logger.info(f"Aggressively extended 3' boundary by {extension_3prime} bp (full extension)")

        # Optional: Log if we found motifs in extended regions (for debugging)
        if extension_5prime > 0:
            search_region_5 = seq_str[new_start:current_start]
            found_5_motifs = [m for m in self.ltr5_motifs if m in search_region_5]
            if found_5_motifs:
                self.logger.debug(f"Found 5' motifs in extended region: {found_5_motifs}")

        if extension_3prime > 0:
            search_region_3 = seq_str[current_end:new_end]
            found_3_motifs = [m for m in self.ltr3_motifs if m in search_region_3]
            if found_3_motifs:
                self.logger.debug(f"Found 3' motifs in extended region: {found_3_motifs}")

        return boundaries


    def _extend_boundaries_by_motifs(self, seq_record, boundaries, max_extension=100):
        """
        Extend boundaries based on terminal motif signals outside current boundaries.

        Args:
            seq_record: SeqRecord object
            boundaries: Dictionary with current boundary information
            max_extension: Maximum number of bases to extend (default: 100, relaxed mode)

        Returns:
            Updated boundaries dictionary
        """
        seq_str = str(seq_record.seq).upper()
        start = boundaries['start']
        end = boundaries['end']

        # 5' end extension: search for stronger start motifs
        best_5prime_extension = 0
        best_5prime_score = 0

        for ext in range(0, min(max_extension, start)):
            search_start = max(0, start - ext - 20)
            search_end = start - ext
            if search_start >= search_end:
                continue
            search_region = seq_str[search_start:search_end]

            # Search for 5' motifs
            for motif in self.ltr5_motifs:
                if motif in search_region:
                    pos = search_region.rfind(motif)  # Rightmost motif
                    # Score: consider motif length and distance
                    # Relaxed mode: decay over 120bp (max_extension + search_window)
                    score = len(motif) * (1.0 - (ext + (20 - pos)) / 120.0)
                    if score > best_5prime_score:
                        best_5prime_score = score
                        best_5prime_extension = ext + (20 - pos) - len(motif)

        # 3' end extension: search for stronger end motifs
        best_3prime_extension = 0
        best_3prime_score = 0

        for ext in range(0, min(max_extension, len(seq_str) - end)):
            search_start = end + ext
            search_end = min(len(seq_str), end + ext + 20)
            if search_start >= search_end:
                continue
            search_region = seq_str[search_start:search_end]

            for motif in self.ltr3_motifs:
                if motif in search_region:
                    pos = search_region.find(motif)
                    # Relaxed mode: decay over 120bp (max_extension + search_window)
                    score = len(motif) * (1.0 - (ext + pos) / 120.0)
                    if score > best_3prime_score:
                        best_3prime_score = score
                        best_3prime_extension = ext + pos + len(motif)

        # Only extend if significantly better motif found
        if best_5prime_score > 0.5:  # Threshold can be adjusted
            boundaries['start'] = max(0, start - best_5prime_extension)
            boundaries['5prime_extended'] = best_5prime_extension
            self.logger.debug(f"Extended 5' boundary by {best_5prime_extension} bp (score: {best_5prime_score:.2f})")

        if best_3prime_score > 0.5:
            boundaries['end'] = min(len(seq_str), end + best_3prime_extension)
            boundaries['3prime_extended'] = best_3prime_extension
            self.logger.debug(f"Extended 3' boundary by {best_3prime_extension} bp (score: {best_3prime_score:.2f})")

        return boundaries


    def _extend_boundaries_by_alignment_distribution(self, boundaries, alignments, percentile=10):
        """
        Extend boundaries using outer percentiles of alignment boundary distribution.

        Args:
            boundaries: Dictionary with current boundary information
            alignments: List of alignments for this sequence
            percentile: Percentile threshold for extension (default: 10)

        Returns:
            Updated boundaries dictionary
        """
        if len(alignments) < 3:
            return boundaries

        starts = [aln['query_start'] for aln in alignments]
        ends = [aln['query_end'] for aln in alignments]

        # Calculate boundary distribution
        start_p10 = int(np.percentile(starts, percentile))  # 10% of alignments start earlier
        end_p90 = int(np.percentile(ends, 100 - percentile))      # 90% of alignments end later

        current_start = boundaries['start']
        current_end = boundaries['end']

        # Conservative extension: only extend to percentiles with sufficient support
        if start_p10 < current_start:
            # Check support count
            support_count = sum(1 for s in starts if s <= start_p10)
            if support_count >= len(starts) * (percentile / 100.0):  # At least percentile% support
                boundaries['start'] = start_p10
                boundaries['start_alignment_extended'] = True
                self.logger.debug(f"Extended start boundary from {current_start} to {start_p10} based on alignment distribution")

        if end_p90 > current_end:
            support_count = sum(1 for e in ends if e >= end_p90)
            if support_count >= len(ends) * (percentile / 100.0):
                boundaries['end'] = end_p90
                boundaries['end_alignment_extended'] = True
                self.logger.debug(f"Extended end boundary from {current_end} to {end_p90} based on alignment distribution")

        return boundaries


    def _safe_extend_by_complexity(self, seq_record, boundaries, max_extension=30):
        """
        Safely extend boundaries but stop at low complexity regions.

        Args:
            seq_record: SeqRecord object
            boundaries: Dictionary with current boundary information
            max_extension: Maximum number of bases to extend (default: 30)

        Returns:
            Updated boundaries dictionary
        """
        seq_str = str(seq_record.seq).upper()

        def is_high_complexity(subseq, min_entropy=1.2):
            """Check if subsequence is high complexity"""
            if len(subseq) < 10:
                return False
            counter = Counter(subseq)
            if sum(counter.values()) == 0:
                return False
            entropy = -sum((c/len(subseq)) * math.log2(c/len(subseq))
                          for c in counter.values() if c > 0)
            return entropy >= min_entropy and 'NNNN' not in subseq

        # 5' end extension
        start = boundaries['start']
        extended_5prime = 0
        for ext in range(5, max_extension + 1, 5):  # Extend in 5bp steps
            if start - ext < 0:
                break
            extension_region = seq_str[start-ext:start]
            if not is_high_complexity(extension_region):
                break
            extended_5prime = ext

        if extended_5prime > 0:
            boundaries['start'] = start - extended_5prime
            self.logger.debug(f"Extended 5' boundary by {extended_5prime} bp based on complexity")

        # 3' end extension
        end = boundaries['end']
        seq_len = len(seq_str)
        extended_3prime = 0
        for ext in range(5, max_extension + 1, 5):
            if end + ext > seq_len:
                break
            extension_region = seq_str[end:end+ext]
            if not is_high_complexity(extension_region):
                break
            extended_3prime = ext

        if extended_3prime > 0:
            boundaries['end'] = end + extended_3prime
            self.logger.debug(f"Extended 3' boundary by {extended_3prime} bp based on complexity")

        return boundaries


    def _extend_by_ltr_symmetry(self, seq_record, boundaries, ltr_similarity_threshold=0.7):
        """
        Fine-tune boundaries for complete LTR elements using LTR similarity.

        Args:
            seq_record: SeqRecord object
            boundaries: Dictionary with current boundary information
            ltr_similarity_threshold: Minimum similarity required (default: 0.7)

        Returns:
            Updated boundaries dictionary
        """
        seq_str = str(seq_record.seq)
        start = boundaries['start']
        end = boundaries['end']
        element_length = end - start

        # Estimate LTR length (typically 200-600bp)
        estimated_ltr_len = min(600, element_length // 4)

        if estimated_ltr_len < 100:
            return boundaries  # Too short for meaningful comparison

        ltr5_seq = seq_str[start:start+estimated_ltr_len]
        ltr3_seq = seq_str[end-estimated_ltr_len:end]

        # Calculate current similarity
        current_similarity = self._calculate_similarity(ltr5_seq, ltr3_seq)

        # Try micro-adjustments to improve LTR similarity
        best_similarity = current_similarity
        best_adjustment = (0, 0)

        for start_adj in range(-20, 21, 5):
            for end_adj in range(-20, 21, 5):
                new_start = start + start_adj
                new_end = end + end_adj

                if new_start < 0 or new_end > len(seq_str) or new_end - new_start < 100:
                    continue

                test_ltr5 = seq_str[new_start:new_start+estimated_ltr_len]
                test_ltr3 = seq_str[new_end-estimated_ltr_len:new_end]

                if len(test_ltr5) < estimated_ltr_len or len(test_ltr3) < estimated_ltr_len:
                    continue

                similarity = self._calculate_similarity(test_ltr5, test_ltr3)

                if similarity > best_similarity:
                    best_similarity = similarity
                    best_adjustment = (start_adj, end_adj)

        # Apply adjustment if significant improvement found
        if best_similarity > current_similarity + 0.05:  # At least 5% improvement
            boundaries['start'] += best_adjustment[0]
            boundaries['end'] += best_adjustment[1]
            boundaries['ltr_similarity'] = best_similarity
            boundaries['symmetry_optimized'] = True
            self.logger.debug(f"Optimized boundaries for LTR symmetry: similarity improved from {current_similarity:.2f} to {best_similarity:.2f}")

        return boundaries


    def _aggressive_extension_enhanced(self, seq_record, boundaries, alignments,
                                       genome_seqs, max_extension=1000):
        """
        Enhanced aggressive extension for incomplete sequences.

        Uses top-K high-quality alignments to guide extension and attempts
        to complete structural features (LTRs, TSDs, PBS/PPT).

        Strategy:
        1. Select top-K alignments based on quality
        2. For each alignment, extract extended genomic context
        3. Attempt to identify complete LTR boundaries
        4. Validate extension using structural features
        5. Return best extension candidate

        Args:
            seq_record: SeqRecord object
            boundaries: Initial boundary estimates
            alignments: List of alignments for this sequence
            genome_seqs: Genome access object
            max_extension: Maximum extension in bp (default: 1000)

        Returns:
            Updated boundaries with extension information
        """
        if not alignments or not genome_seqs:
            # Fallback to original aggressive extension
            return self._aggressive_extension(seq_record, boundaries, max_extension=max_extension)

        # Select top alignments
        top_alignments = self._select_top_alignments(alignments, k=5)

        if not top_alignments:
            return boundaries

        best_extension = None
        best_score = -1

        for aln in top_alignments:
            try:
                # Extract alignment coordinates
                target_id = aln['target_name']
                target_start = aln['target_start']
                target_end = aln['target_end']
                strand = aln.get('strand', '+')

                # Calculate extension range on genome
                genome_ext_start = max(0, target_start - max_extension)
                genome_ext_end = target_end + max_extension

                # Extract extended genomic region
                if not genome_seqs.has_sequence(target_id):
                    continue

                extended_region = genome_seqs.get_sequence_region(
                    target_id, genome_ext_start, genome_ext_end
                )

                if not extended_region:
                    continue

                # Reverse complement if needed
                if strand == '-':
                    extended_region = str(Seq(extended_region).reverse_complement())

                # Create extended sequence record
                extended_record = SeqRecord(
                    Seq(extended_region),
                    id=f"{seq_record.id}_extended",
                    description="Extended from genome"
                )

                # Score this extension using structural completeness
                extension_quality = self._score_ltr_structure(
                    extended_record,
                    alignments=None,  # Don't use alignments for extension scoring
                    genome_seqs=None
                )

                extension_score = extension_quality['total_score']

                # Track best extension
                if extension_score > best_score:
                    best_score = extension_score
                    best_extension = {
                        'sequence': extended_region,
                        'quality': extension_quality,
                        'source_alignment': aln,
                        'boundaries': {
                            'start': 0,
                            'end': len(extended_region),
                            'original_length': len(seq_record.seq),
                            'extension_left': target_start - genome_ext_start,
                            'extension_right': genome_ext_end - target_end,
                            'start_confidence': 0.9,  # High confidence from genome
                            'end_confidence': 0.9,
                            'start_feature': 'genomic_extension',
                            'end_feature': 'genomic_extension'
                        }
                    }

                self.logger.debug(f"{seq_record.id}: Extension candidate score={extension_score:.1f} "
                                f"from {target_id}:{target_start}-{target_end}")

            except Exception as e:
                self.logger.warning(f"Failed to extract extension from alignment: {e}")
                continue

        # If we found a good extension, use it
        if best_extension and best_score >= 60:
            self.logger.info(f"{seq_record.id}: Applied genomic extension "
                           f"(score: {best_score:.1f}, "
                           f"length: {len(seq_record.seq)} -> {len(best_extension['sequence'])}bp)")

            # Update the sequence record
            seq_record.seq = Seq(best_extension['sequence'])
            seq_record.description += f" | Genomically extended (score={best_score:.1f})"

            return best_extension['boundaries']
        else:
            # No good extension found, fallback to original method
            self.logger.debug(f"{seq_record.id}: No suitable genomic extension found, "
                            f"using fallback extension")
            return self._aggressive_extension(seq_record, boundaries, max_extension=max_extension)


    def _validate_extension(self, seq_record, old_boundaries, new_boundaries,
                           alignments, genome_seqs):
        """
        Validate that boundary extension is reasonable.

        Args:
            seq_record: SeqRecord object
            old_boundaries: Original boundaries dictionary
            new_boundaries: Extended boundaries dictionary
            alignments: List of alignments
            genome_seqs: Genome access object

        Returns:
            Validated boundaries (either new or old if validation fails)
        """
        # Extract extension regions
        seq_str = str(seq_record.seq)

        validation_score = 0.0

        # 1. Check if extension improved terminal motifs
        new_seq_str = seq_str[new_boundaries['start']:new_boundaries['end']]
        has_5prime_motif = any(new_seq_str[:20].find(m) != -1 for m in self.ltr5_motifs)
        has_3prime_motif = any(new_seq_str[-20:].rfind(m) != -1 for m in self.ltr3_motifs)

        if has_5prime_motif:
            validation_score += 0.3
        if has_3prime_motif:
            validation_score += 0.3

        # 2. Check sequence quality
        if self._is_high_quality_sequence(new_seq_str):
            validation_score += 0.2

        # 3. Check length is reasonable
        new_length = new_boundaries['end'] - new_boundaries['start']
        if self.min_ltr_len <= new_length <= self.max_ltr_len:
            validation_score += 0.2

        # Accept extension if validation score is sufficient
        if validation_score >= 0.5:
            self.logger.debug(f"Extension validated with score {validation_score:.2f}")
            return new_boundaries
        else:
            self.logger.debug(f"Extension rejected (score {validation_score:.2f}), reverting to original boundaries")
            return old_boundaries


    def _trim_excess_flanking(self, seq_record, quality_info):
        """
        Trim excess flanking sequences based on LTR boundaries.

        Args:
            seq_record: SeqRecord object
            quality_info: Quality information from _score_ltr_structure

        Returns:
            Trimmed SeqRecord object
        """
        seq_str = str(seq_record.seq)
        seq_len = len(seq_str)
        ltr_len = quality_info.get('ltr_length', 0)

        if ltr_len == 0:
            # Can't trim without LTR info
            return seq_record

        # Estimate reasonable total length
        # LTR retrotransposons: 2 * LTR + internal region (typically 2-8kb)
        expected_internal_min = 2000
        expected_internal_max = 8000

        expected_total_min = 2 * ltr_len + expected_internal_min
        expected_total_max = 2 * ltr_len + expected_internal_max

        # If sequence is within expected range, don't trim
        if expected_total_min <= seq_len <= expected_total_max:
            return seq_record

        # If too long, trim flanking regions
        if seq_len > expected_total_max:
            excess = seq_len - expected_total_max

            # Distribute trimming equally to both sides
            trim_left = excess // 2
            trim_right = excess - trim_left

            # Trim
            trimmed_seq = seq_str[trim_left:seq_len - trim_right]

            trimmed_record = SeqRecord(
                Seq(trimmed_seq),
                id=seq_record.id,
                name=seq_record.name,
                description=f"{seq_record.description} | Trimmed: {trim_left}bp left, {trim_right}bp right"
            )

            self.logger.debug(f"{seq_record.id}: Trimmed {excess}bp excess flanking sequence")

            return trimmed_record

        # If too short, return as-is (will be handled by extension)
        return seq_record


    def _refine_sequence(self, seq_record, boundaries):
        """
        Refine sequence based on detected boundaries.
        
        Args:
            seq_record: SeqRecord object
            boundaries: Dictionary with boundary information
            
        Returns:
            Refined SeqRecord object
        """
        # Extract boundaries
        start = boundaries['start']
        end = boundaries['end']
        
        # Safety check
        if start < 0:
            start = 0
        if end > len(seq_record.seq):
            end = len(seq_record.seq)
            
        # Extract refined sequence
        refined_seq = seq_record.seq[start:end]
        
        # Create new record
        refined_record = SeqRecord(
            refined_seq,
            id=seq_record.id,
            name=seq_record.name,
            description=f"Refined: {start+1}-{end}, " + 
                      f"confidence: {boundaries['start_confidence']:.2f}/{boundaries['end_confidence']:.2f}, " +
                      f"features: {boundaries['start_feature'] or 'NA'}/{boundaries['end_feature'] or 'NA'}"
        )
        
        return refined_record
    

    def _calculate_similarity(self, seq1, seq2):
        """
        Calculate similarity between two sequences using simple identity metric.

        Args:
            seq1: First sequence string
            seq2: Second sequence string

        Returns:
            Similarity score between 0 and 1
        """
        if len(seq1) == 0 or len(seq2) == 0:
            return 0.0

        # Use the shorter length for comparison
        min_len = min(len(seq1), len(seq2))
        matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])

        return matches / min_len


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


    def _is_potential_complete_ltr(self, boundaries):
        """
        Check if boundaries suggest a complete LTR element.

        Args:
            boundaries: Dictionary with boundary information

        Returns:
            Boolean indicating if this appears to be a complete LTR element
        """
        element_length = boundaries['end'] - boundaries['start']

        # Complete LTR elements typically have:
        # - Reasonable length (at least 1kb)
        # - Both terminal motifs
        # - Internal features (PBS/PPT)

        if element_length < 1000:
            return False

        has_terminal_motifs = boundaries.get('start_feature') and boundaries.get('end_feature')
        has_internal_features = bool(boundaries.get('internal_features'))

        return has_terminal_motifs and has_internal_features


    def _find_internal_features_enhanced(self, seq_record):
        """
        Enhanced method to find LTR internal features like PBS and PPT.
        
        Args:
            seq_record: SeqRecord object
            
        Returns:
            Dictionary with internal feature information
        """
        seq_str = str(seq_record.seq)
        seq_len = len(seq_str)
        
        results = {}
        
        # Search for PBS patterns
        pbs_results = []
        pbs_region = seq_str[:min(200, seq_len)]  # Search in first 200 bp
        
        for pattern_name, pattern in self.pbs_patterns.items():
            pattern_len = len(pattern)
            
            # Slide through the region looking for matches
            for i in range(len(pbs_region) - pattern_len + 1):
                window = pbs_region[i:i+pattern_len].upper()
                
                # Calculate similarity
                matches = sum(1 for a, b in zip(window, pattern) if a == b)
                similarity = matches / pattern_len
                
                if similarity >= 0.8:  # At least 80% match
                    pbs_results.append({
                        'position': i,
                        'sequence': window,
                        'pattern': pattern_name,
                        'similarity': similarity
                    })
        
        # Take best PBS hit based on similarity
        if pbs_results:
            pbs_results.sort(key=lambda x: x['similarity'], reverse=True)
            results['pbs'] = pbs_results[0]
        
        # Search for PPT (polypurine tract)
        # PPT is typically a stretch of purines (A/G) near 3' LTR
        ppt_results = []
        
        # Use regex for PPT detection
        import re
        ppt_pattern = r'[AG]{8,15}'
        
        # Search in the second half of the sequence
        search_start = max(0, seq_len // 2)
        search_region = seq_str[search_start:]
        
        for match in re.finditer(ppt_pattern, search_region.upper()):
            ppt_results.append({
                'position': search_start + match.start(),
                'sequence': match.group(),
                'length': len(match.group())
            })
            
        # Take PPT hit closest to the end
        if ppt_results:
            # Sort by position (descending) to get the one closest to 3' end
            ppt_results.sort(key=lambda x: x['position'], reverse=True)
            results['ppt'] = ppt_results[0]
            
        # Enhanced: Look for inverted repeats that might indicate terminal repeats
        terminal_repeats = self._find_inverted_repeats(seq_record)
        if terminal_repeats:
            results['terminal_repeats'] = terminal_repeats
            
        return results


    def _find_inverted_repeats(self, seq_record, min_len=5, max_len=15):
        """
        Find inverted repeats that might indicate terminal structures.
        
        Args:
            seq_record: SeqRecord object
            min_len: Minimum repeat length
            max_len: Maximum repeat length
            
        Returns:
            List of inverted repeat information or None
        """
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)
        
        # Too short
        if seq_len < min_len * 2:
            return None
            
        repeats = []
        
        # Look for inverted repeats at the ends
        for repeat_len in range(min_len, min(max_len + 1, seq_len // 2)):
            # Get sequences from the ends
            left_seq = seq_str[:repeat_len]
            right_seq = seq_str[-repeat_len:]
            
            # Get reverse complement of right sequence
            right_rc = str(Seq(right_seq).reverse_complement())
            
            # Calculate similarity
            matches = sum(1 for a, b in zip(left_seq, right_rc) if a == b)
            similarity = matches / repeat_len
            
            if similarity >= 0.8:  # At least 80% match
                repeats.append({
                    'left_position': 0,
                    'right_position': seq_len - repeat_len,
                    'length': repeat_len,
                    'similarity': similarity,
                    'left_seq': left_seq,
                    'right_seq': right_seq
                })
                
        if not repeats:
            return None
            
        # Return the best match (by similarity and length)
        repeats.sort(key=lambda x: (x['similarity'], x['length']), reverse=True)
        return repeats[0]


    # ========================================================================
    # REALIGNMENT AND DEPTH-BASED BOUNDARY DETECTION
    # ========================================================================

    def _realign_extended_sequence(self, seq_record):
        """
        Realign extended sequence to genome and generate new alignments.
        
        This method is called after aggressive extension (iteration 1) to get
        fresh alignment data that includes the extended regions.
        
        Args:
            seq_record: Extended SeqRecord object
            
        Returns:
            List of alignment dictionaries or None if realignment fails
        """
        if not self.genome_file or not os.path.exists(self.genome_file):
            self.logger.debug("Genome file not available for realignment")
            return None
            
        seq_id = seq_record.id
        seq_len = len(seq_record.seq)
        
        self.logger.info(f"{seq_id}: Realigning extended sequence ({seq_len}bp) to genome")
        
        try:
            # Write extended sequence to temp file
            import time
            import random
            temp_seq_file = os.path.join(self.temp_dir, f"{seq_id}_extended_{int(time.time())}_{random.randint(1000, 9999)}.fa")
            SeqIO.write([seq_record], temp_seq_file, "fasta")
            
            # Run minimap2 alignment
            temp_paf_file = temp_seq_file.replace('.fa', '.paf')
            
            cmd = [
                'minimap2',
                '-x', 'map-ont',  # Long-read preset
                '-c',  # Generate CIGAR
                '-t', '4',  # Use 4 threads
                self.genome_file,
                temp_seq_file
            ]
            
            self.logger.debug(f"Running minimap2: {' '.join(cmd)}")
            
            with open(temp_paf_file, 'w') as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True
                )
            
            if result.returncode != 0:
                self.logger.warning(f"minimap2 failed with return code {result.returncode}")
                return None
            
            # Parse PAF alignments
            alignments = self._parse_paf_alignments(temp_paf_file)
            
            # Clean up temp files
            try:
                os.remove(temp_seq_file)
                os.remove(temp_paf_file)
            except:
                pass
            
            if alignments:
                self.logger.info(f"{seq_id}: Realignment complete, found {len(alignments)} alignments")
            else:
                self.logger.warning(f"{seq_id}: Realignment produced no valid alignments")
            
            return alignments

        except Exception as e:
            self.logger.error(f"{seq_id}: Realignment failed: {e}")
            return None
    
    
    def _parse_paf_alignments(self, paf_file):
        """
        Parse PAF format alignment file.
        
        PAF format (12 columns minimum):
        0: Query sequence name
        1: Query sequence length
        2: Query start (0-based)
        3: Query end (0-based)
        4: Strand ('+' or '-')
        5: Target sequence name
        6: Target sequence length
        7: Target start (0-based)
        8: Target end (0-based)
        9: Number of matching bases
        10: Alignment block length
        11: Mapping quality (0-255)
        
        Args:
            paf_file: Path to PAF file
            
        Returns:
            List of alignment dictionaries
        """
        alignments = []
        
        if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
            return alignments
        
        try:
            with open(paf_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 12:
                        continue
                    
                    # Parse basic alignment fields
                    aln = {
                        'query_name': fields[0],
                        'query_len': int(fields[1]),
                        'query_start': int(fields[2]),
                        'query_end': int(fields[3]),
                        'strand': fields[4],
                        'target_name': fields[5],
                        'target_len': int(fields[6]),
                        'target_start': int(fields[7]),
                        'target_end': int(fields[8]),
                        'num_matches': int(fields[9]),
                        'block_len': int(fields[10]),
                        'mapq': int(fields[11])
                    }
                    
                    # Calculate additional metrics
                    aln['query_cov'] = (aln['query_end'] - aln['query_start']) / aln['query_len']
                    aln['identity'] = aln['num_matches'] / aln['block_len'] if aln['block_len'] > 0 else 0
                    
                    # Filter low-quality alignments
                    if aln['query_cov'] >= 0.5 and aln['identity'] >= 0.7:
                        alignments.append(aln)
        
        except Exception as e:
            self.logger.warning(f"Error parsing PAF file: {e}")
        
        return alignments
    
    
    def _create_depth_profile(self, alignments, seq_len):
        """
        Create coverage depth profile from alignments.
        
        Args:
            alignments: List of alignment dictionaries
            seq_len: Length of query sequence
            
        Returns:
            numpy array of depth values (one per position)
        """
        depth = np.zeros(seq_len, dtype=np.float32)
        
        for aln in alignments:
            start = aln['query_start']
            end = aln['query_end']
            
            # Weight by alignment quality
            weight = aln['identity'] * min(1.0, aln['mapq'] / 60.0)
            
            # Increment depth
            depth[start:end] += weight
        
        return depth
    
    
    def _detect_boundaries_from_depth(self, depth_profile, current_boundaries=None, 
                                     min_depth_ratio=0.3, smoothing_window=51):
        """
        Detect precise boundaries from depth profile using signal processing.
        
        Strategy:
        1. Smooth depth profile to reduce noise
        2. Detect significant coverage transitions
        3. Identify stable high-depth regions (core LTR)
        4. Refine boundaries to transition points
        
        Args:
            depth_profile: numpy array of depth values
            current_boundaries: Current boundary estimates (optional)
            min_depth_ratio: Minimum depth ratio to consider valid (default: 0.3)
            smoothing_window: Window size for smoothing (default: 51, must be odd)
            
        Returns:
            Dictionary with refined boundary information or None
        """
        seq_len = len(depth_profile)
        
        # Check if we have coverage
        max_depth = np.max(depth_profile)
        if max_depth == 0:
            self.logger.debug("No coverage found in depth profile")
            return current_boundaries
        
        # Normalize depth
        normalized_depth = depth_profile / max_depth
        
        # Smooth depth profile using Savitzky-Golay filter
        try:
            from scipy.signal import savgol_filter
            
            # Ensure window size is appropriate
            window = min(smoothing_window, seq_len)
            if window % 2 == 0:
                window -= 1  # Must be odd
            if window < 5:
                window = 5
            
            if seq_len >= window:
                smoothed = savgol_filter(normalized_depth, window_length=window, polyorder=3)
            else:
                smoothed = normalized_depth
        except ImportError:
            self.logger.warning("scipy not available, using simple moving average")
            # Fallback to simple moving average
            window = min(51, seq_len)
            smoothed = np.convolve(normalized_depth, np.ones(window)/window, mode='same')
        
        # Find regions with sufficient coverage
        high_cov_mask = smoothed >= min_depth_ratio
        
        if not np.any(high_cov_mask):
            self.logger.debug(f"No regions with sufficient coverage (>={min_depth_ratio})")
            return current_boundaries
        
        # Find contiguous high-coverage regions
        high_cov_regions = []
        in_region = False
        region_start = 0
        
        for i in range(seq_len):
            if high_cov_mask[i] and not in_region:
                # Start of new region
                in_region = True
                region_start = i
            elif not high_cov_mask[i] and in_region:
                # End of region
                in_region = False
                high_cov_regions.append((region_start, i))
        
        # Handle case where region extends to end
        if in_region:
            high_cov_regions.append((region_start, seq_len))
        
        if not high_cov_regions:
            self.logger.debug("No contiguous high-coverage regions found")
            return current_boundaries
        
        # Find the longest high-coverage region (most likely the core LTR)
        longest_region = max(high_cov_regions, key=lambda x: x[1] - x[0])
        core_start, core_end = longest_region
        
        self.logger.debug(f"Core high-coverage region: {core_start}-{core_end} ({core_end - core_start}bp)")
        
        # Refine boundaries by finding sharp transitions
        # Look for the steepest drops on each side
        
        # 5' boundary refinement
        refined_start = core_start
        if core_start > 20:
            # Search window before core
            search_start = max(0, core_start - 100)
            search_region = smoothed[search_start:core_start]
            
            # Find position with steepest increase (derivative peak)
            if len(search_region) > 5:
                gradient = np.gradient(search_region)
                # Find position with maximum positive gradient
                peak_idx = np.argmax(gradient)
                refined_start = search_start + peak_idx
        
        # 3' boundary refinement
        refined_end = core_end
        if core_end < seq_len - 20:
            # Search window after core
            search_end = min(seq_len, core_end + 100)
            search_region = smoothed[core_end:search_end]
            
            # Find position with steepest decrease (derivative trough)
            if len(search_region) > 5:
                gradient = np.gradient(search_region)
                # Find position with maximum negative gradient
                trough_idx = np.argmin(gradient)
                refined_end = core_end + trough_idx
        
        # Calculate confidence based on depth drop magnitude
        if refined_start > 0 and refined_end < seq_len:
            left_drop = smoothed[refined_start] - np.mean(smoothed[max(0, refined_start-10):refined_start])
            right_drop = smoothed[refined_end] - np.mean(smoothed[refined_end:min(seq_len, refined_end+10)])
            
            # Normalize to 0-1
            depth_confidence = min(1.0, (abs(left_drop) + abs(right_drop)) / 2.0)
        else:
            depth_confidence = 0.5
        
        # Create boundary info
        boundaries = {
            'start': int(refined_start),
            'end': int(refined_end),
            'start_confidence': min(1.0, depth_confidence + 0.2),  # Boost confidence slightly
            'end_confidence': min(1.0, depth_confidence + 0.2),
            'start_feature': 'depth_transition',
            'end_feature': 'depth_transition',
            'method': 'depth_profile',
            'confidence_score': depth_confidence * 100,
            'max_depth': float(max_depth),
            'core_region': (int(core_start), int(core_end))
        }
        
        # Validate boundaries
        boundary_len = boundaries['end'] - boundaries['start']
        if boundary_len < self.min_ltr_len or boundary_len > self.max_ltr_len:
            self.logger.debug(f"Depth-based boundaries outside valid range: {boundary_len}bp")
            return current_boundaries
        
        self.logger.info(f"Depth-based boundaries: {refined_start}-{refined_end} "
                        f"({boundary_len}bp, confidence={depth_confidence:.2f})")
        
        return boundaries


