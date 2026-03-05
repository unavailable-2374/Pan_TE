#!/usr/bin/env python3
"""
LTR Boundary Refiner Module - Simplified Strategy

New Strategy:
1. Physical extraction: ±500bp
2. Self-alignment using blastn/rmblastn to detect LTR pairs
3. Branch strategy:
   - If both LTRs present → use LTR coordinates to update boundaries → search for motifs
   - If LTR missing → align to genome → use depth to update boundaries → search for motifs
"""

import os
import logging
import subprocess
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np


class LTRBoundaryRefiner:
    """
    LTR Boundary Detection and Refinement - Simplified Strategy

    Key improvements:
    - Self-alignment based LTR detection
    - Depth-based boundary refinement when LTRs are incomplete
    - Motif-based validation
    """

    def __init__(self, min_ltr_len=100, max_ltr_len=5000,
                 min_ltr_similarity=0.70, tsd_min=3, tsd_max=6, temp_dir=None, logger=None):
        """
        Initialize the boundary refiner.

        Args:
            min_ltr_len: Minimum LTR length (default: 100bp)
            max_ltr_len: Maximum LTR length (default: 5000bp)
            min_ltr_similarity: Minimum LTR similarity (default: 0.70)
            tsd_min: Minimum TSD length (default: 3)
            tsd_max: Maximum TSD length (default: 6)
            temp_dir: Temporary directory for alignment files
            logger: Logger instance
        """
        self.min_ltr_len = min_ltr_len
        self.max_ltr_len = max_ltr_len
        self.min_ltr_similarity = min_ltr_similarity
        self.tsd_min = tsd_min
        self.tsd_max = tsd_max
        self.temp_dir = temp_dir if temp_dir else tempfile.gettempdir()

        # Setup logger
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger('LTRBoundaryRefiner')
            if not self.logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                self.logger.addHandler(handler)
                self.logger.setLevel(logging.INFO)

        # LTR terminal motifs
        self.ltr5_motifs = ["TG", "TGT", "TGTA", "TGTG", "TGCA"]
        self.ltr3_motifs = ["CA", "ACA", "TACA", "CACA", "TGCA"]

    # ========================================================================
    # PUBLIC API
    # ========================================================================

    def refine_boundaries(self, seq_record, alignments=None, genome_seqs=None, iteration=1):
        """
        Main entry point for boundary refinement using simplified strategy.

        Args:
            seq_record: BioPython SeqRecord object
            alignments: List of genome alignments (optional)
            genome_seqs: Genome access object (optional)
            iteration: Iteration number (not used in new strategy)

        Returns:
            dict with keys:
                - 'refined_record': Refined SeqRecord
                - 'boundaries': Boundary information
                - 'ltr_pair': LTR pair information if detected
                - 'quality_score': Quality score (0-100)
                - 'method': Detection method used
        """
        seq_id = seq_record.id
        seq_len = len(seq_record.seq)

        self.logger.info(f"Refining boundaries for {seq_id} ({seq_len}bp)")

        # Step 1: Self-alignment to detect LTR pairs
        ltr_pair = self._detect_ltr_pair_by_self_alignment(seq_record)

        if ltr_pair and ltr_pair['both_ltrs_found']:
            # Branch A: Both LTRs detected
            self.logger.info(f"{seq_id}: Both LTRs detected via self-alignment")
            boundaries = self._refine_with_ltr_pair(seq_record, ltr_pair)
            method = 'ltr_pair_based'
        else:
            # Branch B: Missing LTR(s), use genome alignment + depth
            self.logger.info(f"{seq_id}: LTRs incomplete, using genome alignment + depth")
            boundaries = self._refine_with_genome_depth(seq_record, alignments, genome_seqs)
            method = 'depth_based'

        if not boundaries:
            self.logger.warning(f"{seq_id}: Failed to detect boundaries")
            return {
                'refined_record': seq_record,
                'boundaries': None,
                'ltr_pair': ltr_pair,
                'quality_score': 0,
                'method': 'failed'
            }

        # Step 2: Search for terminal motifs around boundaries
        boundaries = self._search_motifs_around_boundaries(seq_record, boundaries)

        # Step 3: Extract refined sequence
        refined_record = self._extract_refined_sequence(seq_record, boundaries)

        # Calculate quality score
        quality_score = self._calculate_quality_score(boundaries, ltr_pair)

        self.logger.info(f"{seq_id}: Refinement complete (method={method}, quality={quality_score:.1f})")

        return {
            'refined_record': refined_record,
            'boundaries': boundaries,
            'ltr_pair': ltr_pair,
            'quality_score': quality_score,
            'method': method
        }

    # ========================================================================
    # SELF-ALIGNMENT LTR DETECTION
    # ========================================================================

    def _detect_ltr_pair_by_self_alignment(self, seq_record):
        """
        Detect LTR pairs using self-alignment with rmblastn/blastn.

        Strategy:
        1. Write sequence to temp file
        2. Run rmblastn (or blastn fallback) self-alignment
        3. Parse alignments to identify LTR pairs
        4. Validate LTR similarity and structure

        Args:
            seq_record: SeqRecord object

        Returns:
            dict with LTR pair information or None:
                - 'both_ltrs_found': bool
                - 'ltr5': {start, end, seq}
                - 'ltr3': {start, end, seq}
                - 'similarity': float (0-1)
                - 'internal_start': int
                - 'internal_end': int
        """
        seq_id = seq_record.id
        seq_len = len(seq_record.seq)

        # Write sequence to temp file
        import time
        import random
        temp_fasta = os.path.join(self.temp_dir, f"{seq_id}_selfalign_{int(time.time())}_{random.randint(1000, 9999)}.fa")
        SeqIO.write([seq_record], temp_fasta, "fasta")

        # Try rmblastn first (optimized for repeat detection)
        alignments = self._run_self_alignment_rmblastn(temp_fasta, seq_id)

        if not alignments:
            # Fallback to regular blastn
            self.logger.debug(f"{seq_id}: rmblastn failed, trying blastn")
            alignments = self._run_self_alignment_blastn(temp_fasta, seq_id)

        # Clean up temp file
        try:
            os.remove(temp_fasta)
        except:
            pass

        if not alignments:
            self.logger.debug(f"{seq_id}: Self-alignment produced no valid alignments")
            return {'both_ltrs_found': False}

        self.logger.debug(f"{seq_id}: Found {len(alignments)} self-alignments")

        # Parse alignments to find LTR pairs
        ltr_pair = self._identify_ltr_pair_from_alignments(alignments, seq_len)

        return ltr_pair

    def _run_self_alignment_rmblastn(self, fasta_file, seq_id):
        """
        Run rmblastn self-alignment.

        Args:
            fasta_file: Path to FASTA file
            seq_id: Sequence ID for filtering

        Returns:
            List of alignment dicts or None
        """
        try:
            # Create blast database
            makeblastdb_cmd = [
                'makeblastdb',
                '-in', fasta_file,
                '-dbtype', 'nucl',
                '-out', fasta_file + '.db'
            ]

            subprocess.run(makeblastdb_cmd, check=True,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Run rmblastn
            output_file = fasta_file + '.blast'
            rmblastn_cmd = [
                'rmblastn',
                '-query', fasta_file,
                '-subject', fasta_file,
                '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen',
                '-out', output_file,
                '-gapopen', '12',  # Higher gap penalty for LTR detection
                '-gapextend', '2'
            ]

            subprocess.run(rmblastn_cmd, check=True,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Parse output
            alignments = self._parse_blast_output(output_file, seq_id)

            # Clean up
            try:
                os.remove(output_file)
                for ext in ['.db.nhr', '.db.nin', '.db.nsq']:
                    try:
                        os.remove(fasta_file + ext)
                    except:
                        pass
            except:
                pass

            return alignments

        except (subprocess.CalledProcessError, FileNotFoundError):
            return None

    def _run_self_alignment_blastn(self, fasta_file, seq_id):
        """
        Run blastn self-alignment (fallback).

        Args:
            fasta_file: Path to FASTA file
            seq_id: Sequence ID for filtering

        Returns:
            List of alignment dicts or None
        """
        try:
            # Create blast database
            makeblastdb_cmd = [
                'makeblastdb',
                '-in', fasta_file,
                '-dbtype', 'nucl',
                '-out', fasta_file + '.db'
            ]

            subprocess.run(makeblastdb_cmd, check=True,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Run blastn
            output_file = fasta_file + '.blast'
            blastn_cmd = [
                'blastn',
                '-query', fasta_file,
                '-subject', fasta_file,
                '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen',
                '-out', output_file,
                '-word_size', '11',  # Sensitive for repeats
                '-reward', '2',
                '-penalty', '-3',
                '-gapopen', '5',
                '-gapextend', '2'
            ]

            subprocess.run(blastn_cmd, check=True,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Parse output
            alignments = self._parse_blast_output(output_file, seq_id)

            # Clean up
            try:
                os.remove(output_file)
                for ext in ['.db.nhr', '.db.nin', '.db.nsq']:
                    try:
                        os.remove(fasta_file + ext)
                    except:
                        pass
            except:
                pass

            return alignments

        except (subprocess.CalledProcessError, FileNotFoundError):
            return None

    def _parse_blast_output(self, blast_file, seq_id):
        """
        Parse BLAST tabular output.

        Format: qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen

        Args:
            blast_file: Path to BLAST output file
            seq_id: Sequence ID to filter

        Returns:
            List of alignment dicts
        """
        alignments = []

        if not os.path.exists(blast_file):
            return alignments

        try:
            with open(blast_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 12:
                        continue

                    qseqid = fields[0]
                    sseqid = fields[1]
                    qstart = int(fields[4]) - 1  # Convert to 0-based
                    qend = int(fields[5])
                    sstart = int(fields[6]) - 1
                    send = int(fields[7])

                    # Skip perfect self-hits (same position, same length)
                    if qseqid == sseqid:
                        if abs(qstart - sstart) < 10 and abs(qend - send) < 10:
                            continue

                    aln = {
                        'qseqid': qseqid,
                        'sseqid': sseqid,
                        'pident': float(fields[2]),
                        'length': int(fields[3]),
                        'qstart': qstart,
                        'qend': qend,
                        'sstart': sstart,
                        'send': send,
                        'evalue': float(fields[8]),
                        'bitscore': float(fields[9]),
                        'qlen': int(fields[10]),
                        'slen': int(fields[11])
                    }

                    # Filter by quality
                    if aln['pident'] >= 70 and aln['length'] >= self.min_ltr_len:
                        alignments.append(aln)

        except Exception as e:
            self.logger.warning(f"Error parsing BLAST output: {e}")

        return alignments

    def _identify_ltr_pair_from_alignments(self, alignments, seq_len):
        """
        Identify LTR pairs from self-alignment results.

        Strategy:
        - Look for two alignments that suggest LTR structure:
          1. One near the start (5' LTR)
          2. One near the end (3' LTR)
          3. Similar lengths and high identity

        Args:
            alignments: List of alignment dicts from BLAST
            seq_len: Total sequence length

        Returns:
            dict with LTR pair information
        """
        if len(alignments) < 2:
            return {'both_ltrs_found': False}

        # Group alignments by query position (5' vs 3')
        left_alignments = []  # Alignments starting in first third
        right_alignments = []  # Alignments starting in last third

        for aln in alignments:
            qstart = aln['qstart']
            qend = aln['qend']
            sstart = aln['sstart']
            send = aln['send']

            # Ensure qstart < qend and sstart < send
            if qstart > qend:
                qstart, qend = qend, qstart
            if sstart > send:
                sstart, send = send, sstart

            # Check if this is a potential LTR pair
            # LTR pair: one alignment in 5' region matches 3' region
            if qstart < seq_len / 3 and sstart > seq_len * 2 / 3:
                left_alignments.append({
                    'ltr5_start': qstart,
                    'ltr5_end': qend,
                    'ltr3_start': sstart,
                    'ltr3_end': send,
                    'similarity': aln['pident'] / 100.0,
                    'length': aln['length']
                })
            elif sstart < seq_len / 3 and qstart > seq_len * 2 / 3:
                # Reverse case
                left_alignments.append({
                    'ltr5_start': sstart,
                    'ltr5_end': send,
                    'ltr3_start': qstart,
                    'ltr3_end': qend,
                    'similarity': aln['pident'] / 100.0,
                    'length': aln['length']
                })

        if not left_alignments:
            return {'both_ltrs_found': False}

        # Find best LTR pair (highest similarity, longest length)
        best_pair = max(left_alignments,
                       key=lambda x: (x['similarity'], x['length']))

        # Validate LTR pair
        if best_pair['similarity'] < self.min_ltr_similarity:
            return {'both_ltrs_found': False}

        ltr5_len = best_pair['ltr5_end'] - best_pair['ltr5_start']
        ltr3_len = best_pair['ltr3_end'] - best_pair['ltr3_start']

        if ltr5_len < self.min_ltr_len or ltr3_len < self.min_ltr_len:
            return {'both_ltrs_found': False}

        if ltr5_len > self.max_ltr_len or ltr3_len > self.max_ltr_len:
            return {'both_ltrs_found': False}

        # Build LTR pair info
        return {
            'both_ltrs_found': True,
            'ltr5': {
                'start': best_pair['ltr5_start'],
                'end': best_pair['ltr5_end'],
                'length': ltr5_len
            },
            'ltr3': {
                'start': best_pair['ltr3_start'],
                'end': best_pair['ltr3_end'],
                'length': ltr3_len
            },
            'similarity': best_pair['similarity'],
            'internal_start': best_pair['ltr5_end'],
            'internal_end': best_pair['ltr3_start']
        }

    # ========================================================================
    # LTR-BASED BOUNDARY REFINEMENT
    # ========================================================================

    def _refine_with_ltr_pair(self, seq_record, ltr_pair):
        """
        Refine boundaries using detected LTR pair coordinates.

        Strategy:
        - Use LTR coordinates as initial boundaries
        - Validate with motif search
        - Adjust boundaries if motifs found nearby

        Args:
            seq_record: SeqRecord object
            ltr_pair: LTR pair information from self-alignment

        Returns:
            dict with boundary information
        """
        ltr5 = ltr_pair['ltr5']
        ltr3 = ltr_pair['ltr3']

        # Initial boundaries from LTR coordinates
        boundaries = {
            'start': ltr5['start'],
            'end': ltr3['end'],
            'method': 'ltr_pair',
            'ltr5_end': ltr5['end'],
            'ltr3_start': ltr3['start'],
            'ltr_similarity': ltr_pair['similarity'],
            'confidence_score': ltr_pair['similarity'] * 100
        }

        self.logger.debug(f"LTR-based boundaries: {boundaries['start']}-{boundaries['end']}")

        return boundaries

    # ========================================================================
    # GENOME DEPTH-BASED BOUNDARY REFINEMENT
    # ========================================================================

    def _refine_with_genome_depth(self, seq_record, alignments, genome_seqs):
        """
        Refine boundaries using genome alignment depth profile.

        Strategy:
        1. Re-align sequence to genome (if not already aligned)
        2. Build depth profile from alignments
        3. Detect boundaries from depth transitions

        Args:
            seq_record: SeqRecord object
            alignments: List of genome alignments
            genome_seqs: Genome access object

        Returns:
            dict with boundary information
        """
        seq_len = len(seq_record.seq)

        if not alignments or len(alignments) < 2:
            self.logger.warning(f"{seq_record.id}: Not enough alignments for depth-based refinement")
            return None

        # Build depth profile
        depth_profile = np.zeros(seq_len, dtype=np.float32)

        for aln in alignments:
            qstart = aln.get('query_start', 0)
            qend = aln.get('query_end', seq_len)
            identity = aln.get('identity', 1.0)

            # Weight by identity
            depth_profile[qstart:qend] += identity

        # Detect boundaries from depth
        boundaries = self._detect_boundaries_from_depth_simple(depth_profile)

        if boundaries:
            boundaries['method'] = 'depth_profile'
            self.logger.debug(f"Depth-based boundaries: {boundaries['start']}-{boundaries['end']}")

        return boundaries

    def _detect_boundaries_from_depth_simple(self, depth_profile):
        """
        Simple depth-based boundary detection.

        Strategy:
        1. Find core high-depth region
        2. Extend to depth drop-off points

        Args:
            depth_profile: numpy array of depth values

        Returns:
            dict with boundary information or None
        """
        seq_len = len(depth_profile)
        max_depth = np.max(depth_profile)

        if max_depth == 0:
            return None

        # Normalize
        norm_depth = depth_profile / max_depth

        # Find region with depth > 50% of max
        threshold = 0.5
        high_depth_mask = norm_depth >= threshold

        if not np.any(high_depth_mask):
            return None

        # Find first and last high-depth positions
        high_indices = np.where(high_depth_mask)[0]
        start = int(high_indices[0])
        end = int(high_indices[-1])

        # Validate length
        length = end - start
        if length < self.min_ltr_len or length > self.max_ltr_len:
            return None

        return {
            'start': start,
            'end': end,
            'confidence_score': 70.0,  # Moderate confidence for depth-based
            'max_depth': float(max_depth)
        }

    # ========================================================================
    # MOTIF SEARCH AND VALIDATION
    # ========================================================================

    def _search_motifs_around_boundaries(self, seq_record, boundaries):
        """
        Search for terminal motifs AND Target Site Duplications (TSDs) around boundaries.

        Strategy:
        1. Search for TSDs (4-6bp direct repeats flanking the element).
           This is the "Gold Standard" for LTR boundaries.
        2. If TSD found, snap boundaries to TSD.
        3. If no TSD, fallback to Motif search (TG..CA).

        Args:
            seq_record: SeqRecord object
            boundaries: Current boundary information

        Returns:
            Updated boundaries with motif/TSD information
        """
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)
        start = boundaries['start']
        end = boundaries['end']
        
        # --- Strategy 1: TSD Detection (The Gold Standard) ---
        # Look for 4-6bp identical sequences flanking the current boundaries
        # Search window: ±20bp around current start/end
        
        tsd_found = False
        best_tsd_score = 0
        best_tsd_coords = None
        
        # Define search windows
        # Expanded window size to catch boundaries that are further off
        window_size = 50
        start_window_min = max(0, start - window_size)
        start_window_max = min(seq_len, start + window_size)
        
        end_window_min = max(0, end - window_size)
        end_window_max = min(seq_len, end + window_size)
        
        # Extract regions to search for TSDs
        # We look for: Sequence[i : i+k] == Sequence[j : j+k]
        # where i is near start, j is near end
        
        # Optimization: Only check k=4,5,6
        for k in [4, 5, 6]:
            # Get candidate TSDs from the start region
            # We assume the TSD is *immediately* before the LTR start
            # So we iterate possible "real" starts
            for s_pos in range(start_window_min, start_window_max - k):
                candidate_tsd = seq_str[s_pos : s_pos + k]
                
                # Check if this TSD exists near the end region
                # The TSD should be *immediately* after the LTR end
                # So we look for it in the end window
                e_pos = seq_str.find(candidate_tsd, end_window_min, end_window_max)
                
                if e_pos != -1:
                    # Found a match!
                    # Calculate "LTR length" implied by this TSD
                    # LTR would be from s_pos + k to e_pos
                    implied_len = e_pos - (s_pos + k)
                    
                    if self.min_ltr_len <= implied_len <= self.max_ltr_len:
                        # Validate with Motifs inside the TSDs
                        # LTR starts at s_pos + k, ends at e_pos
                        ltr_seq = seq_str[s_pos + k : e_pos]
                        
                        has_5_motif = any(ltr_seq.startswith(m) for m in self.ltr5_motifs)
                        has_3_motif = any(ltr_seq.endswith(m) for m in self.ltr3_motifs)
                        
                        score = k * 10
                        if has_5_motif: score += 5
                        if has_3_motif: score += 5
                        
                        # Prefer TSDs closer to original boundaries
                        dist_penalty = abs(s_pos + k - start) + abs(e_pos - end)
                        score -= dist_penalty * 0.5
                        
                        if score > best_tsd_score:
                            best_tsd_score = score
                            best_tsd_coords = (s_pos + k, e_pos, candidate_tsd)
                            tsd_found = True

        if tsd_found:
            new_start, new_end, tsd_seq = best_tsd_coords
            boundaries['start'] = new_start
            boundaries['end'] = new_end
            boundaries['tsd'] = tsd_seq
            boundaries['method'] = 'TSD_refined'
            self.logger.debug(f"TSD found: {tsd_seq} at {new_start}-{new_end}")
            return boundaries

        # --- Strategy 2: Motif Search (Fallback) ---
        # Search for 5' motifs near start
        # Reverted search window to 20bp to avoid snapping to random noise
        search_window = 20
        search_start = max(0, start - search_window)
        search_end = min(seq_len, start + search_window)
        search_region_5 = seq_str[search_start:search_end]

        best_5_motif = None
        best_5_pos = start

        for motif in self.ltr5_motifs:
            pos = search_region_5.find(motif)
            if pos != -1:
                actual_pos = search_start + pos
                if best_5_motif is None or len(motif) > len(best_5_motif):
                    best_5_motif = motif
                    best_5_pos = actual_pos

        # Search for 3' motifs near end
        search_start = max(0, end - search_window)
        search_end = min(seq_len, end + search_window)
        search_region_3 = seq_str[search_start:search_end]

        best_3_motif = None
        best_3_pos = end

        for motif in self.ltr3_motifs:
            pos = search_region_3.rfind(motif)
            if pos != -1:
                actual_pos = search_start + pos + len(motif)
                if best_3_motif is None or len(motif) > len(best_3_motif):
                    best_3_motif = motif
                    best_3_pos = actual_pos

        # Update boundaries if motifs found close to original
        # Relaxed proximity check
        if best_5_motif and abs(best_5_pos - start) <= search_window:
            boundaries['start'] = best_5_pos
            boundaries['start_motif'] = best_5_motif
            self.logger.debug(f"Adjusted start boundary to motif {best_5_motif} at {best_5_pos}")

        if best_3_motif and abs(best_3_pos - end) <= search_window:
            boundaries['end'] = best_3_pos
            boundaries['end_motif'] = best_3_motif
            self.logger.debug(f"Adjusted end boundary to motif {best_3_motif} at {best_3_pos}")

        return boundaries

    # ========================================================================
    # UTILITY METHODS
    # ========================================================================

    def _extract_refined_sequence(self, seq_record, boundaries):
        """
        Extract refined sequence based on boundaries.

        Args:
            seq_record: Original SeqRecord
            boundaries: Boundary information

        Returns:
            Refined SeqRecord
        """
        start = max(0, boundaries['start'])
        end = min(len(seq_record.seq), boundaries['end'])

        refined_seq = seq_record.seq[start:end]

        refined_record = SeqRecord(
            refined_seq,
            id=seq_record.id,
            name=seq_record.name,
            description=f"Refined: {start+1}-{end} | Method: {boundaries.get('method', 'unknown')}"
        )

        return refined_record

    def _calculate_quality_score(self, boundaries, ltr_pair):
        """
        Calculate quality score for refined boundaries.

        Args:
            boundaries: Boundary information
            ltr_pair: LTR pair information (may be None)

        Returns:
            Quality score (0-100)
        """
        score = boundaries.get('confidence_score', 50.0)

        # Bonus for LTR pair detection
        if ltr_pair and ltr_pair.get('both_ltrs_found'):
            score += 20

        # Bonus for motif detection
        if boundaries.get('start_motif'):
            score += 10
        if boundaries.get('end_motif'):
            score += 10

        return min(100.0, score)
