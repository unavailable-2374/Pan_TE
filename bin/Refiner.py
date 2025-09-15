#!/usr/bin/env python3
"""
Enhanced Refiner.py based on bin/Refiner/phase2 logic
Replaces the original Refiner.py with advanced consensus building capabilities
"""

import os
import sys
import logging
import argparse
import tempfile
import subprocess
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict, Counter
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import combinations

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)

class EnhancedTERefiner:
    """Enhanced TE Consensus Builder with 5 key improvements:
    1. Enhanced sequence similarity scoring with k-mer analysis
    2. Adaptive threshold determination based on sequence divergence
    3. Multi-round consensus building with iterative refinement
    4. TSD-aware boundary refinement for accurate element detection
    5. Quality-based subfamily detection and classification
    """
    
    def __init__(self, min_score=150, gap_init=20, gap_ext=5, threads=4, kmer_size=21):
        self.min_score = min_score
        self.gap_init = gap_init
        self.gap_ext = gap_ext
        self.threads = threads or 4
        self.kmer_size = kmer_size
        
        # Adaptive thresholds - will be determined based on sequence characteristics
        self.adaptive_identity_threshold = 85.0
        self.adaptive_length_threshold = 0.2
        
        # Multi-round parameters
        self.max_refinement_rounds = 3
        self.convergence_threshold = 0.95  # Stop if 95% sequences unchanged between rounds
        
        # TSD detection parameters
        self.tsd_min_length = 2
        self.tsd_max_length = 20
        self.tsd_similarity_threshold = 0.8
        
        # Try to find required tools
        self.mafft_path = self._find_tool('mafft')
        self.blast_path = self._find_tool('blastn')
        if not self.mafft_path:
            logger.warning("MAFFT not found in PATH, will use basic consensus")
        
        logger.info(f"Enhanced TE Refiner initialized: min_score={min_score}, threads={threads}, kmer_size={kmer_size}")
    
    def _find_tool(self, tool_name: str) -> Optional[str]:
        """Find tool in PATH"""
        import shutil
        return shutil.which(tool_name)
    
    def build_consensus(self, input_file: str, output_file: str) -> bool:
        """Main consensus building function"""
        try:
            # Read input sequences
            sequences = list(SeqIO.parse(input_file, "fasta"))
            if not sequences:
                logger.error(f"No sequences found in {input_file}")
                return False
            
            logger.info(f"Processing {len(sequences)} sequences from {input_file}")
            
            # Convert to internal format
            copies = []
            for i, seq in enumerate(sequences):
                copies.append({
                    'id': seq.id,
                    'sequence': str(seq.seq),
                    'length': len(seq.seq),
                    'description': seq.description
                })
            
            # Build consensus using phase2 logic
            consensus_result = self._build_enhanced_consensus(copies)
            
            if not consensus_result:
                logger.error("Failed to build consensus")
                return False
            
            # Write output
            self._write_consensus(consensus_result, output_file)
            logger.info(f"Consensus written to {output_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error building consensus: {e}")
            return False
    
    def _build_enhanced_consensus(self, copies: List[Dict]) -> Optional[Dict]:
        """Multi-round consensus building with iterative refinement"""
        
        if not copies:
            return None
        
        if len(copies) == 1:
            # Single sequence - apply TSD detection and boundary refinement
            refined_copy = self._refine_sequence_boundaries([copies[0]])[0]
            return self._create_single_consensus(refined_copy)
        
        # Initial analysis
        characteristics = self._analyze_sequence_characteristics(copies)
        logger.info(f"Initial sequence characteristics: {characteristics}")
        
        # Apply TSD-aware boundary refinement to all sequences
        refined_copies = self._refine_sequence_boundaries(copies)
        
        # Re-analyze after boundary refinement
        refined_characteristics = self._analyze_sequence_characteristics(refined_copies)
        logger.info(f"Post-refinement characteristics: {refined_characteristics}")
        
        # Quality-based subfamily detection
        if self._should_split_subfamilies(refined_characteristics, len(refined_copies)):
            logger.info("Applying quality-based subfamily detection")
            subfamilies = self._identify_subfamilies_enhanced(refined_copies, refined_characteristics)
            
            # Build consensus for each subfamily and select best
            if subfamilies and len(subfamilies) > 1:
                subfamily_consensuses = []
                for i, subfamily in enumerate(subfamilies):
                    logger.info(f"Building consensus for subfamily {i+1} with {len(subfamily)} sequences")
                    consensus = self._multi_round_consensus(subfamily, refined_characteristics)
                    if consensus:
                        subfamily_consensuses.append((len(subfamily), consensus))
                
                # Return consensus from largest high-quality subfamily
                if subfamily_consensuses:
                    subfamily_consensuses.sort(key=lambda x: x[0], reverse=True)
                    return subfamily_consensuses[0][1]
        
        # Multi-round consensus building for all sequences
        return self._multi_round_consensus(refined_copies, refined_characteristics)
    
    def _multi_round_consensus(self, copies: List[Dict], characteristics: Dict) -> Optional[Dict]:
        """Multi-round iterative consensus refinement"""
        
        current_copies = copies.copy()
        previous_consensus = None
        
        for round_num in range(1, self.max_refinement_rounds + 1):
            logger.info(f"Consensus refinement round {round_num}/{self.max_refinement_rounds}")
            
            # Build consensus for current round
            current_consensus = self._build_msa_consensus(current_copies, characteristics)
            
            if not current_consensus:
                logger.warning(f"Failed to build consensus in round {round_num}")
                return previous_consensus if previous_consensus else None
            
            # Check convergence
            if previous_consensus and self._check_convergence(previous_consensus, current_consensus):
                logger.info(f"Consensus converged after round {round_num}")
                return current_consensus
            
            # Filter sequences based on similarity to current consensus for next round
            if round_num < self.max_refinement_rounds:
                filtered_copies = self._filter_sequences_by_consensus(
                    current_copies, current_consensus, characteristics
                )
                
                # If we lose too many sequences, stop refinement
                retention_rate = len(filtered_copies) / len(current_copies)
                if retention_rate < 0.5:
                    logger.info(f"Low retention rate ({retention_rate:.2f}), stopping refinement")
                    return current_consensus
                
                current_copies = filtered_copies
                logger.info(f"Round {round_num}: retained {len(current_copies)}/{len(copies)} sequences")
            
            previous_consensus = current_consensus
        
        return previous_consensus
    
    def _check_convergence(self, prev_consensus: Dict, curr_consensus: Dict) -> bool:
        """Check if consensus has converged between rounds"""
        if not prev_consensus or not curr_consensus:
            return False
        
        prev_seq = prev_consensus['sequence']
        curr_seq = curr_consensus['sequence']
        
        # Calculate similarity between consecutive consensus sequences
        similarity = self._calculate_identity(prev_seq, curr_seq)
        
        return similarity >= self.convergence_threshold * 100
    
    def _filter_sequences_by_consensus(self, copies: List[Dict], consensus: Dict, 
                                     characteristics: Dict) -> List[Dict]:
        """Filter sequences based on similarity to consensus"""
        consensus_seq = consensus['sequence']
        filtered_copies = []
        
        # Dynamic threshold based on sequence characteristics
        min_similarity = max(
            self.adaptive_identity_threshold - 10,  # Allow some flexibility
            50.0  # Minimum threshold
        )
        
        for copy in copies:
            similarity = self._calculate_identity(copy['sequence'], consensus_seq)
            if similarity >= min_similarity:
                filtered_copies.append(copy)
            else:
                logger.debug(f"Filtered out sequence {copy['id']} (similarity: {similarity:.1f}%)")
        
        return filtered_copies if len(filtered_copies) >= 2 else copies
    
    def _analyze_sequence_characteristics(self, copies: List[Dict]) -> Dict:
        """Analyze sequence characteristics with adaptive threshold determination"""
        lengths = [copy['length'] for copy in copies]
        sequences = [copy['sequence'] for copy in copies]
        
        # Calculate pairwise identities (enhanced with k-mer analysis)
        identities = []
        gc_contents = []
        complexity_scores = []
        
        # Analyze individual sequences
        for seq in sequences:
            gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100 if len(seq) > 0 else 0
            gc_contents.append(gc_content)
            
            # Calculate sequence complexity (Shannon entropy)
            complexity = self._calculate_sequence_complexity(seq)
            complexity_scores.append(complexity)
        
        # Sample pairwise comparisons efficiently
        if len(sequences) > 1:
            max_pairs = min(100, len(sequences) * (len(sequences) - 1) // 2)
            pairs_sampled = 0
            
            for i in range(len(sequences)):
                for j in range(i + 1, len(sequences)):
                    if pairs_sampled >= max_pairs:
                        break
                    identity = self._calculate_identity(sequences[i], sequences[j])
                    identities.append(identity)
                    pairs_sampled += 1
                if pairs_sampled >= max_pairs:
                    break
        
        characteristics = {
            'count': len(copies),
            'avg_length': np.mean(lengths),
            'length_std': np.std(lengths),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'length_cv': np.std(lengths) / np.mean(lengths) if np.mean(lengths) > 0 else 0,
            'avg_identity': np.mean(identities) if identities else 100.0,
            'identity_std': np.std(identities) if identities else 0.0,
            'avg_gc_content': np.mean(gc_contents),
            'gc_std': np.std(gc_contents),
            'avg_complexity': np.mean(complexity_scores),
            'complexity_std': np.std(complexity_scores)
        }
        
        # Adaptive threshold determination
        self._determine_adaptive_thresholds(characteristics)
        
        return characteristics
    
    def _calculate_sequence_complexity(self, sequence: str) -> float:
        """Calculate sequence complexity using Shannon entropy"""
        if len(sequence) == 0:
            return 0.0
        
        # Count nucleotide frequencies
        counts = Counter(sequence.upper())
        total = sum(counts.values())
        
        if total == 0:
            return 0.0
        
        # Calculate Shannon entropy
        entropy = 0.0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)
        
        # Normalize to 0-1 range (max entropy for DNA is 2 bits)
        return entropy / 2.0
    
    def _determine_adaptive_thresholds(self, characteristics: Dict):
        """Determine adaptive thresholds based on sequence divergence patterns"""
        avg_identity = characteristics['avg_identity']
        identity_std = characteristics['identity_std']
        length_cv = characteristics['length_cv']
        count = characteristics['count']
        
        # Adaptive identity threshold based on sequence divergence
        if avg_identity > 95:
            # Very similar sequences - use stricter threshold
            self.adaptive_identity_threshold = 92.0
        elif avg_identity > 85:
            # Moderately similar - standard threshold
            self.adaptive_identity_threshold = 80.0
        elif avg_identity > 70:
            # Divergent sequences - relaxed threshold
            self.adaptive_identity_threshold = 65.0
        else:
            # Very divergent - very relaxed threshold
            self.adaptive_identity_threshold = 50.0
        
        # Adjust based on variability
        if identity_std > 15:
            # High variability - relax threshold
            self.adaptive_identity_threshold -= 5.0
        elif identity_std < 5:
            # Low variability - can be more stringent
            self.adaptive_identity_threshold += 5.0
        
        # Adaptive length threshold
        if count < 5:
            # Few sequences - be more permissive with length variation
            self.adaptive_length_threshold = 0.3
        elif length_cv > 0.4:
            # High length variation
            self.adaptive_length_threshold = 0.5
        else:
            # Standard threshold
            self.adaptive_length_threshold = 0.2
        
        logger.info(f"Adaptive thresholds: identity={self.adaptive_identity_threshold:.1f}%, length_cv={self.adaptive_length_threshold:.2f}")
    
    def _calculate_identity(self, seq1: str, seq2: str) -> float:
        """Enhanced sequence similarity scoring with k-mer analysis"""
        if len(seq1) == 0 or len(seq2) == 0:
            return 0.0
        
        # Use k-mer based similarity for better accuracy
        kmer_similarity = self._calculate_kmer_similarity(seq1, seq2)
        
        # Combine with direct alignment for short sequences
        if min(len(seq1), len(seq2)) < 100:
            direct_similarity = self._calculate_direct_similarity(seq1, seq2)
            # Weighted average favoring k-mer for longer sequences
            weight = min(len(seq1), len(seq2)) / 100.0
            return kmer_similarity * weight + direct_similarity * (1 - weight)
        else:
            return kmer_similarity
    
    def _calculate_kmer_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate k-mer based similarity (Jaccard index)"""
        # Generate k-mers
        kmers1 = self._get_kmers(seq1, self.kmer_size)
        kmers2 = self._get_kmers(seq2, self.kmer_size)
        
        if not kmers1 or not kmers2:
            return 0.0
        
        # Calculate Jaccard similarity
        intersection = len(kmers1 & kmers2)
        union = len(kmers1 | kmers2)
        
        if union == 0:
            return 0.0
        
        jaccard = intersection / union
        return jaccard * 100.0
    
    def _get_kmers(self, sequence: str, k: int) -> Set[str]:
        """Extract k-mers from sequence"""
        kmers = set()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k].upper()
            if 'N' not in kmer:  # Skip k-mers with ambiguous nucleotides
                kmers.add(kmer)
        return kmers
    
    def _calculate_direct_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate direct alignment similarity for short sequences"""
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))
        
        # Simple alignment with sliding window
        best_matches = 0
        best_overlap = 0
        
        for offset in range(-min_len//2, min_len//2 + 1):
            matches = 0
            overlap = 0
            
            start1 = max(0, offset)
            start2 = max(0, -offset)
            
            end1 = min(len(seq1), len(seq2) + offset)
            end2 = min(len(seq2), len(seq1) - offset)
            
            for i in range(max(start1, start2), min(end1, end2)):
                pos1 = i - offset if offset >= 0 else i
                pos2 = i + offset if offset < 0 else i
                
                if 0 <= pos1 < len(seq1) and 0 <= pos2 < len(seq2):
                    overlap += 1
                    if seq1[pos1] == seq2[pos2]:
                        matches += 1
            
            if overlap > best_overlap or (overlap == best_overlap and matches > best_matches):
                best_matches = matches
                best_overlap = overlap
        
        return (best_matches / best_overlap * 100.0) if best_overlap > 0 else 0.0
    
    def _should_split_subfamilies(self, characteristics: Dict, num_copies: int) -> bool:
        """Quality-based subfamily detection decision using adaptive thresholds"""
        if num_copies < 4:
            return False
        
        # Use adaptive thresholds instead of fixed ones
        identity_threshold = self.adaptive_identity_threshold
        length_threshold = self.adaptive_length_threshold
        
        # Multiple criteria for subfamily splitting
        conditions = [
            # High sequence count with diversity
            num_copies >= 10 and characteristics['avg_identity'] < identity_threshold,
            # Length variation exceeds adaptive threshold
            characteristics['length_cv'] > length_threshold,
            # High identity variance suggests subfamilies
            characteristics['identity_std'] > 20.0,
            # GC content variation suggests different origins
            characteristics.get('gc_std', 0) > 10.0,
            # Complexity variation suggests structural differences
            characteristics.get('complexity_std', 0) > 0.3
        ]
        
        # Require at least 2 conditions for splitting (more conservative)
        split_score = sum(conditions)
        
        # Additional quality checks
        if split_score >= 2:
            # Ensure we have enough sequences for meaningful subfamilies
            if num_copies >= 8:
                logger.info(f"Subfamily splitting recommended: {split_score}/5 conditions met")
                return True
        
        return False
    
    def _identify_subfamilies_enhanced(self, copies: List[Dict], characteristics: Dict) -> List[List[Dict]]:
        """Enhanced subfamily identification with quality-based clustering"""
        
        if len(copies) < 4:
            return [copies]
        
        try:
            # Multi-dimensional clustering using multiple features
            features = []
            
            for copy in copies:
                seq = copy['sequence']
                feature_vector = [
                    copy['length'],
                    (seq.count('G') + seq.count('C')) / len(seq) * 100,  # GC content
                    self._calculate_sequence_complexity(seq),  # Complexity
                    self._count_tandem_repeats(seq),  # Tandem repeat content
                ]
                features.append(feature_vector)
            
            # Normalize features
            features = np.array(features)
            features = (features - features.mean(axis=0)) / (features.std(axis=0) + 1e-8)
            
            # Build combined distance matrix
            n = len(copies)
            distances = np.zeros((n, n))
            
            for i in range(n):
                for j in range(i + 1, n):
                    # Sequence similarity component
                    seq_identity = self._calculate_identity(copies[i]['sequence'], copies[j]['sequence'])
                    seq_distance = 100.0 - seq_identity
                    
                    # Feature similarity component
                    feature_distance = np.linalg.norm(features[i] - features[j])
                    
                    # Combined distance (weighted)
                    combined_distance = 0.7 * seq_distance + 0.3 * feature_distance * 20
                    distances[i, j] = combined_distance
                    distances[j, i] = combined_distance
            
            # Adaptive threshold based on characteristics
            base_threshold = 100.0 - self.adaptive_identity_threshold
            
            # Dynamic clustering
            subfamilies = self._hierarchical_clustering(copies, distances, base_threshold)
            
            # Quality filtering of subfamilies
            quality_subfamilies = []
            for subfamily in subfamilies:
                if len(subfamily) >= 2:
                    # Calculate subfamily cohesion
                    cohesion = self._calculate_subfamily_cohesion(subfamily)
                    if cohesion > 0.7:  # High cohesion threshold
                        quality_subfamilies.append(subfamily)
                    else:
                        logger.debug(f"Filtered out subfamily with low cohesion: {cohesion:.2f}")
            
            # If no quality subfamilies, return largest clusters
            if not quality_subfamilies:
                subfamilies.sort(key=len, reverse=True)
                quality_subfamilies = subfamilies[:3]  # Top 3 largest
            
            logger.info(f"Identified {len(quality_subfamilies)} quality subfamilies")
            return quality_subfamilies if quality_subfamilies else [copies]
            
        except Exception as e:
            logger.warning(f"Enhanced subfamily identification failed: {e}, using simple clustering")
            return self._identify_subfamilies_simple(copies, characteristics)
    
    def _identify_subfamilies_simple(self, copies: List[Dict], characteristics: Dict) -> List[List[Dict]]:
        """Simple fallback subfamily identification"""
        
        if len(copies) < 4:
            return [copies]
        
        try:
            n = len(copies)
            distances = np.zeros((n, n))
            
            for i in range(n):
                for j in range(i + 1, n):
                    identity = self._calculate_identity(copies[i]['sequence'], copies[j]['sequence'])
                    distance = 100.0 - identity
                    distances[i, j] = distance
                    distances[j, i] = distance
            
            threshold = 100.0 - self.adaptive_identity_threshold
            return self._hierarchical_clustering(copies, distances, threshold)
            
        except Exception as e:
            logger.warning(f"Simple subfamily identification failed: {e}, using all sequences")
            return [copies]
    
    def _hierarchical_clustering(self, copies: List[Dict], distances: np.ndarray, threshold: float) -> List[List[Dict]]:
        """Simple hierarchical clustering"""
        
        n = len(copies)
        visited = [False] * n
        subfamilies = []
        
        for i in range(n):
            if visited[i]:
                continue
            
            subfamily = [copies[i]]
            visited[i] = True
            
            # Find all sequences within threshold
            for j in range(i + 1, n):
                if not visited[j] and distances[i, j] <= threshold:
                    subfamily.append(copies[j])
                    visited[j] = True
            
            if len(subfamily) >= 1:  # Keep all clusters
                subfamilies.append(subfamily)
        
        return subfamilies if subfamilies else [copies]
    
    def _calculate_subfamily_cohesion(self, subfamily: List[Dict]) -> float:
        """Calculate cohesion score for a subfamily"""
        
        if len(subfamily) < 2:
            return 1.0
        
        similarities = []
        for i in range(len(subfamily)):
            for j in range(i + 1, len(subfamily)):
                sim = self._calculate_identity(subfamily[i]['sequence'], subfamily[j]['sequence'])
                similarities.append(sim)
        
        return np.mean(similarities) / 100.0 if similarities else 0.0
    
    def _count_tandem_repeats(self, sequence: str) -> float:
        """Count tandem repeat content in sequence"""
        
        if len(sequence) < 20:
            return 0.0
        
        repeat_count = 0
        window_size = min(10, len(sequence) // 4)
        
        for i in range(len(sequence) - window_size * 2):
            motif = sequence[i:i + window_size]
            next_motif = sequence[i + window_size:i + window_size * 2]
            
            if motif == next_motif:
                repeat_count += 1
        
        return repeat_count / (len(sequence) - window_size * 2) if len(sequence) > window_size * 2 else 0.0
    
    def _refine_sequence_boundaries(self, copies: List[Dict]) -> List[Dict]:
        """TSD-aware boundary refinement for accurate element detection"""
        
        refined_copies = []
        
        for copy in copies:
            refined_copy = copy.copy()
            
            # Detect and trim terminal repeats/TSDs
            refined_seq = self._detect_and_trim_tsds(copy['sequence'])
            
            if refined_seq and len(refined_seq) >= 50:  # Minimum length threshold
                refined_copy['sequence'] = refined_seq
                refined_copy['length'] = len(refined_seq)
                refined_copy['boundary_refined'] = True
                
                # Calculate improvement metrics
                original_len = len(copy['sequence'])
                refined_len = len(refined_seq)
                refined_copy['length_change'] = refined_len - original_len
                refined_copy['length_ratio'] = refined_len / original_len
            else:
                refined_copy['boundary_refined'] = False
            
            refined_copies.append(refined_copy)
        
        refined_count = sum(1 for c in refined_copies if c.get('boundary_refined', False))
        logger.info(f"Boundary refinement: {refined_count}/{len(copies)} sequences refined")
        
        return refined_copies
    
    def _detect_and_trim_tsds(self, sequence: str) -> Optional[str]:
        """Detect Target Site Duplications and trim boundaries"""
        
        if len(sequence) < self.tsd_max_length * 4:  # Too short for TSD analysis
            return sequence
        
        # Check for TSDs at both ends
        best_tsd_length = 0
        best_start = 0
        best_end = len(sequence)
        
        for tsd_len in range(self.tsd_min_length, min(self.tsd_max_length + 1, len(sequence) // 4)):
            left_tsd = sequence[:tsd_len]
            right_tsd = sequence[-tsd_len:]
            
            # Calculate similarity between terminal sequences
            similarity = self._calculate_tsd_similarity(left_tsd, right_tsd)
            
            if similarity >= self.tsd_similarity_threshold:
                # Found potential TSD
                if tsd_len > best_tsd_length:
                    best_tsd_length = tsd_len
                    best_start = tsd_len
                    best_end = len(sequence) - tsd_len
        
        # Additional boundary refinement - remove low complexity regions
        if best_tsd_length > 0:
            trimmed_seq = sequence[best_start:best_end]
            
            # Further refine by removing terminal low-complexity regions
            refined_seq = self._trim_low_complexity_ends(trimmed_seq)
            
            if len(refined_seq) >= 50:
                logger.debug(f"TSD detected: length={best_tsd_length}, refined length: {len(sequence)} -> {len(refined_seq)}")
                return refined_seq
        
        # No significant TSDs found, apply light trimming
        return self._trim_low_complexity_ends(sequence)
    
    def _calculate_tsd_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate similarity between potential TSDs"""
        
        if len(seq1) != len(seq2) or len(seq1) == 0:
            return 0.0
        
        matches = sum(1 for i in range(len(seq1)) if seq1[i] == seq2[i])
        return matches / len(seq1)
    
    def _trim_low_complexity_ends(self, sequence: str, window_size: int = 20) -> str:
        """Trim low-complexity regions from sequence ends"""
        
        if len(sequence) <= window_size * 2:
            return sequence
        
        # Trim from start
        start_pos = 0
        for i in range(0, len(sequence) - window_size, window_size // 2):
            window = sequence[i:i + window_size]
            complexity = self._calculate_sequence_complexity(window)
            
            if complexity > 0.3:  # Found reasonable complexity
                start_pos = max(0, i - window_size // 2)
                break
        
        # Trim from end
        end_pos = len(sequence)
        for i in range(len(sequence) - window_size, window_size, -window_size // 2):
            window = sequence[i - window_size:i]
            complexity = self._calculate_sequence_complexity(window)
            
            if complexity > 0.3:  # Found reasonable complexity
                end_pos = min(len(sequence), i + window_size // 2)
                break
        
        return sequence[start_pos:end_pos] if end_pos > start_pos else sequence
    
    def _build_msa_consensus(self, copies: List[Dict], characteristics: Dict) -> Optional[Dict]:
        """Build consensus using multiple sequence alignment"""
        
        if not copies:
            return None
        
        if len(copies) == 1:
            return self._create_single_consensus(copies[0])
        
        # Try MAFFT if available
        if self.mafft_path:
            consensus_seq = self._run_mafft_consensus(copies, characteristics)
            if consensus_seq:
                return self._create_consensus_record(consensus_seq, copies, "mafft")
        
        # Fallback to simple consensus
        logger.info("Using simple consensus method")
        consensus_seq = self._build_simple_consensus(copies)
        if consensus_seq:
            return self._create_consensus_record(consensus_seq, copies, "simple")
        
        # Last resort - use longest sequence
        longest = max(copies, key=lambda x: x['length'])
        return self._create_single_consensus(longest)
    
    def _run_mafft_consensus(self, copies: List[Dict], characteristics: Dict) -> Optional[str]:
        """Run MAFFT and build consensus from MSA"""
        
        try:
            # Create temporary input file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_input:
                for i, copy in enumerate(copies):
                    tmp_input.write(f">seq_{i}\n{copy['sequence']}\n")
                tmp_input_path = tmp_input.name
            
            # Choose MAFFT algorithm based on characteristics
            if len(copies) > 100:
                algorithm = ['--auto']
            elif characteristics.get('avg_identity', 70) > 80:
                algorithm = ['--localpair', '--maxiterate', '1000']
            else:
                algorithm = ['--genafpair', '--maxiterate', '1000']
            
            # Run MAFFT
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_output:
                tmp_output_path = tmp_output.name
            
            cmd = [self.mafft_path] + algorithm + ['--quiet', '--thread', str(self.threads), tmp_input_path]
            
            logger.info(f"Running MAFFT: {' '.join(cmd)}")
            
            with open(tmp_output_path, 'w') as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
            
            if result.returncode != 0:
                logger.warning(f"MAFFT failed: {result.stderr}")
                return None
            
            # Build consensus from MSA
            consensus_seq = self._consensus_from_msa(tmp_output_path)
            
            # Cleanup
            os.unlink(tmp_input_path)
            os.unlink(tmp_output_path)
            
            return consensus_seq
            
        except Exception as e:
            logger.warning(f"MAFFT consensus failed: {e}")
            return None
    
    def _consensus_from_msa(self, msa_file: str) -> Optional[str]:
        """Build consensus from MSA file"""
        
        try:
            aligned_seqs = []
            for record in SeqIO.parse(msa_file, "fasta"):
                aligned_seqs.append(str(record.seq))
            
            if not aligned_seqs:
                return None
            
            # Build consensus position by position
            consensus = []
            msa_length = len(aligned_seqs[0])
            
            for pos in range(msa_length):
                column = [seq[pos] for seq in aligned_seqs if pos < len(seq)]
                column_clean = [c for c in column if c not in '-']
                
                if not column_clean:
                    continue  # Skip gap-only columns
                
                # Get most common nucleotide
                counter = Counter(column_clean)
                most_common = counter.most_common(1)[0]
                
                # Use majority rule with minimum coverage
                if most_common[1] >= len(column_clean) * 0.3:  # At least 30% support
                    consensus.append(most_common[0])
            
            consensus_seq = ''.join(consensus)
            
            # Remove excessive gaps and clean up
            consensus_seq = consensus_seq.replace('-', '')
            
            return consensus_seq if len(consensus_seq) >= 50 else None
            
        except Exception as e:
            logger.warning(f"Error building consensus from MSA: {e}")
            return None
    
    def _build_simple_consensus(self, copies: List[Dict]) -> Optional[str]:
        """Build simple consensus without MSA"""
        
        if not copies:
            return None
        
        # For simple consensus, just use the longest sequence
        # In a more sophisticated implementation, we could do position-wise consensus
        longest = max(copies, key=lambda x: x['length'])
        longest = max(copies, key=lambda x: x['length'])
        return longest['sequence']
    
    def _create_single_consensus(self, copy: Dict) -> Dict:
        """Create consensus record from single sequence"""
        return {
            'id': copy['id'] + '_consensus',
            'sequence': copy['sequence'],
            'source_id': copy['id'],
            'copy_number': 1,
            'length': copy['length'],
            'method': 'single'
        }
    
    def _create_consensus_record(self, consensus_seq: str, copies: List[Dict], method: str) -> Dict:
        """Create consensus record"""
        
        # Use first sequence ID as base
        base_id = copies[0]['id']
        if '_' in base_id:
            base_id = base_id.split('_')[0]
        
        return {
            'id': f"{base_id}_consensus",
            'sequence': consensus_seq,
            'source_id': base_id,
            'copy_number': len(copies),
            'length': len(consensus_seq),
            'method': method,
            'original_length': np.mean([c['length'] for c in copies]),
            'improvement_ratio': len(consensus_seq) / np.mean([c['length'] for c in copies])
        }
    
    def _write_consensus(self, consensus_result: Dict, output_file: str):
        """Write consensus to output file"""
        
        record = SeqRecord(
            Seq(consensus_result['sequence']),
            id=consensus_result['id'],
            description=f"consensus from {consensus_result['copy_number']} copies, method={consensus_result['method']}"
        )
        
        with open(output_file, 'w') as handle:
            SeqIO.write([record], handle, "fasta")


def main():
    """Command line interface - maintains compatibility with original Refiner.py"""
    
    parser = argparse.ArgumentParser(
        description='Enhanced TE Consensus Builder based on advanced phase2 logic'
    )
    
    parser.add_argument('input', help='Input FASTA file with family sequences')
    parser.add_argument('output', help='Output FASTA file for consensus')
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='Number of threads (default: 4)')
    parser.add_argument('--min-score', type=float, default=150,
                       help='Minimum alignment score (default: 150)')
    parser.add_argument('--gap-init', type=int, default=20,
                       help='Gap initiation penalty (default: 20)')
    parser.add_argument('--gap-ext', type=int, default=5,
                       help='Gap extension penalty (default: 5)')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate input file
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        # Create refiner and build consensus
        refiner = EnhancedTERefiner(
            min_score=args.min_score,
            gap_init=args.gap_init,
            gap_ext=args.gap_ext,
            threads=args.threads
        )
        
        success = refiner.build_consensus(args.input, args.output)
        
        if success:
            logger.info("Consensus building completed successfully")
            sys.exit(0)
        else:
            logger.error("Consensus building failed")
            sys.exit(1)
            
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()