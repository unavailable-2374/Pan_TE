#!/usr/bin/env python3
"""
LTR Chimera Detector Module

Handles detection and splitting of chimeric LTR sequences using structural analysis.
Detects internal LTR copies that indicate chimeric fusion of multiple elements.
"""

import os
import logging
import subprocess
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class LTRChimeraDetector:
    """
    LTR Chimera Detection and Splitting based on Structural Analysis.
    
    Strategy:
    1. Extract the 5' LTR sequence (putative).
    2. Blast this LTR against the rest of the sequence.
    3. If a significant hit is found in the middle, it indicates a chimera 
       (e.g., LTR1-Internal-LTR1-LTR2-Internal-LTR2).
    4. Split the sequence at the internal LTR boundary.
    """

    def __init__(self, chimera_threshold=0.8, min_segment_length=500,
                 orientation_aware=True, logger=None):
        """
        Initialize the chimera detector.

        Args:
            chimera_threshold: Minimum identity to consider an internal match as LTR (default: 0.8)
            min_segment_length: Minimum length for segments when splitting
            orientation_aware: (Not used in structural mode but kept for API compatibility)
            logger: Logger instance
        """
        self.chimera_threshold = chimera_threshold
        self.min_segment_length = min_segment_length
        self.temp_dir = tempfile.gettempdir()

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

    def detect_chimeras(self, fasta_file, alignments_by_query=None):
        """
        Main entry point for chimera detection.
        
        Args:
            fasta_file: Path to input FASTA file
            alignments_by_query: (Ignored in structural mode)
            
        Returns:
            dict with split records and statistics
        """
        self.logger.info(f"Starting structural chimera detection on {fasta_file}")
        
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        split_records = []
        chimera_info = {}
        
        for seq_record in sequences:
            is_chimeric, segments = self._analyze_and_split_structural(seq_record)
            
            chimera_info[seq_record.id] = {
                'is_chimeric': is_chimeric,
                'num_segments': len(segments)
            }
            
            split_records.extend(segments)
            
        num_chimeras = sum(1 for info in chimera_info.values() if info['is_chimeric'])
        num_splits = sum(info['num_segments'] - 1 for info in chimera_info.values())
        
        self.logger.info(f"Chimera detection complete: {num_chimeras} chimeras detected, {num_splits} splits performed")
        
        return {
            'split_records': split_records,
            'chimera_info': chimera_info,
            'num_chimeras': num_chimeras,
            'num_splits': num_splits
        }

    def _analyze_and_split_structural(self, seq_record):
        """
        Analyze sequence structure to detect internal LTRs and split if necessary.
        """
        seq_len = len(seq_record.seq)
        if seq_len < 2000: # Too short to be a chimera of two complete elements
            return False, [seq_record]
            
        # Heuristic: Assume 5' LTR is within the first 20% or max 1000bp
        ltr_guess_len = min(1000, int(seq_len * 0.2))
        ltr_query = seq_record.seq[:ltr_guess_len]
        
        # Write query and subject to temp files
        import random
        rand_id = random.randint(10000, 99999)
        query_file = os.path.join(self.temp_dir, f"chimera_q_{rand_id}.fa")
        subject_file = os.path.join(self.temp_dir, f"chimera_s_{rand_id}.fa")
        
        SeqIO.write(SeqRecord(ltr_query, id="query"), query_file, "fasta")
        SeqIO.write(seq_record, subject_file, "fasta")
        
        # Run blastn
        # Look for hits starting AFTER the query region
        try:
            cmd = [
                'blastn',
                '-query', query_file,
                '-subject', subject_file,
                '-outfmt', '6 qstart qend sstart send pident length',
                '-task', 'blastn', # Standard blastn for high similarity
                '-dust', 'no'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            breakpoints = []
            
            for line in result.stdout.strip().split('\n'):
                if not line: continue
                parts = line.split('\t')
                qstart, qend, sstart, send = map(int, parts[:4])
                pident = float(parts[4])
                length = int(parts[5])
                
                # Filter hits
                # 1. Must be high identity (it's a copy of the LTR)
                if pident < self.chimera_threshold * 100: continue
                
                # 2. Must be long enough (at least 100bp)
                if length < 100: continue
                
                # 3. Must be "internal" - significantly after the query region
                # The query is at 0-ltr_guess_len. 
                # A hit representing the 3' LTR of the FIRST element should be around:
                # Start of 2nd element = End of 1st element
                
                # We are looking for the START of the 2nd element's 5' LTR.
                # Or the 3' LTR of the 1st element.
                
                # If we find a hit starting at 'sstart', it means there is an LTR-like sequence there.
                # If this is a chimera: LTR1 ... LTR1 [Break] LTR2 ... LTR2
                # The hit at [Break] corresponds to the 3' LTR of element 1.
                # Or LTR1 ... [Break] LTR2 ... LTR2
                
                # Let's simplify: Any significant LTR hit in the middle (excluding the very end)
                # is a potential breakpoint.
                
                # Exclude self-hit (start ~ 0)
                if sstart < ltr_guess_len: continue
                
                # Exclude the true 3' LTR of the full sequence (end ~ seq_len)
                if send > seq_len - 100: continue
                
                # Found an internal LTR!
                # The breakpoint is likely around 'send' (end of 1st element) 
                # or 'sstart' (start of 2nd element if they share LTRs? No, usually tandem)
                
                # If it's LTR1 ... LTR1-LTR2 ... LTR2
                # The hit is the 2nd LTR1. So the split should be after it.
                breakpoints.append(send)
                
            # Clean up
            os.remove(query_file)
            os.remove(subject_file)
            
            if not breakpoints:
                return False, [seq_record]
                
            # Sort and filter breakpoints
            breakpoints.sort()
            final_breakpoints = []
            last_pos = 0
            
            for bp in breakpoints:
                # Ensure segments are long enough
                if bp - last_pos > self.min_segment_length and seq_len - bp > self.min_segment_length:
                    final_breakpoints.append(bp)
                    last_pos = bp
            
            if not final_breakpoints:
                return False, [seq_record]
                
            # Split
            segments = []
            positions = [0] + final_breakpoints + [seq_len]
            
            for i in range(len(positions) - 1):
                start = positions[i]
                end = positions[i+1]
                seg_seq = seq_record.seq[start:end]

                # Preserve original description (including copy number info)
                # Format: "Original_description | Chimera_segment from position X-Y"
                original_desc = seq_record.description if seq_record.description else ""
                seg_rec = SeqRecord(
                    seg_seq,
                    id=f"{seq_record.id}_seg{i+1}",
                    description=f"{original_desc} | Chimera_segment_{i+1} ({start}-{end})"
                )
                segments.append(seg_rec)
                
            self.logger.info(f"Detected chimera: {seq_record.id} -> {len(segments)} segments")
            return True, segments
            
        except Exception as e:
            self.logger.warning(f"Chimera detection failed for {seq_record.id}: {e}")
            if os.path.exists(query_file): os.remove(query_file)
            if os.path.exists(subject_file): os.remove(subject_file)
            return False, [seq_record]