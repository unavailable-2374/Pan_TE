#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
import tempfile
import time
import re
import shutil
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from itertools import groupby

class LTREnhancedOptimizer:
    """
    Comprehensive tool for refining LTR consensus sequences with:
    1. Chimeric sequence detection and resolution
    2. Precise boundary determination using structural features
    3. Non-transposon sequence filtering
    4. Efficient processing of GB-scale genomes
    
    Can be used as a standalone script or imported as a module in workflows.
    """
    
    def __init__(self, consensus_file=None, genome_file=None, output_dir=None, 
                 threads=1, min_identity=80, min_coverage=80, 
                 min_ltr_len=100, max_ltr_len=1500, flanking_seq=20,
                 chimera_threshold=0.3, min_segment_length=500,
                 tsd_min=4, tsd_max=6, iterations=2, 
                 chunk_size=10, batch_size=50, 
                 keep_temp=False, log_level="INFO"):
        """
        Initialize the LTR enhanced optimizer.
        
        Args:
            consensus_file: Path to consensus sequences from Refiner_for_LTR
            genome_file: Path to reference genome
            output_dir: Directory for output files
            threads: Number of CPU threads to use
            min_identity: Minimum identity percentage for valid alignments
            min_coverage: Minimum coverage percentage for valid alignments
            min_ltr_len: Minimum length of an LTR
            max_ltr_len: Maximum length of an LTR
            flanking_seq: Length of flanking sequence to analyze for TSD
            chimera_threshold: Threshold for chimera detection (higher = more sensitive)
            min_segment_length: Minimum length for segments when splitting chimeras
            tsd_min: Minimum TSD length
            tsd_max: Maximum TSD length
            iterations: Number of optimization iterations
            chunk_size: Number of genome sequences to process at once
            batch_size: Number of consensus sequences to process in each batch
            keep_temp: Whether to keep temporary files
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        """
        # Core parameters
        self.consensus_file = consensus_file
        self.genome_file = genome_file
        self.output_dir = output_dir
        
        # Processing parameters
        self.threads = threads
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.min_ltr_len = min_ltr_len
        self.max_ltr_len = max_ltr_len
        self.flanking_seq = flanking_seq
        self.chimera_threshold = chimera_threshold
        self.min_segment_length = min_segment_length
        self.tsd_min = tsd_min
        self.tsd_max = tsd_max
        self.iterations = iterations
        self.chunk_size = chunk_size
        self.batch_size = batch_size
        self.keep_temp = keep_temp
        
        # Initialize logger
        self._setup_logging(log_level)
        
        # Status tracking
        self.initialized = False
        self.genome_indexed = False
        self.temp_dir = None
        
        # LTR-specific patterns
        self.ltr5_motifs = ["TG", "TGT", "TGTA", "TGTG", "TGCA"]
        self.ltr3_motifs = ["CA", "ACA", "TACA", "CACA", "TGCA"]
        self.tsd_patterns = {}  # Will store found TSDs
        
    def _setup_logging(self, log_level):
        """Setup logging configuration"""
        numeric_level = getattr(logging, log_level.upper(), None)
        if not isinstance(numeric_level, int):
            numeric_level = logging.INFO
            
        logging.basicConfig(
            level=numeric_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler("ltr_optimizer.log"),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger('LTROptimizer')
        
    def initialize(self, consensus_file=None, genome_file=None, output_dir=None):
        """
        Initialize the optimizer with files and directories.
        Call this if you didn't provide files in the constructor.
        """
        # Update file paths if provided
        if consensus_file:
            self.consensus_file = consensus_file
        if genome_file:
            self.genome_file = genome_file
        if output_dir:
            self.output_dir = output_dir
            
        # Validate required files
        if not self.consensus_file or not os.path.exists(self.consensus_file):
            raise ValueError(f"Consensus file not found: {self.consensus_file}")
        if not self.genome_file or not os.path.exists(self.genome_file):
            raise ValueError(f"Genome file not found: {self.genome_file}")
            
        # Create output directory
        if not self.output_dir:
            self.output_dir = os.path.join(os.path.dirname(self.consensus_file), "optimized")
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Create temporary directory
        self.temp_dir = os.path.join(self.output_dir, f"temp_{int(time.time())}")
        os.makedirs(self.temp_dir, exist_ok=True)
        
        # Create result directories
        self.alignment_dir = os.path.join(self.output_dir, "alignments")
        self.chimera_dir = os.path.join(self.output_dir, "chimera_detection")
        self.boundary_dir = os.path.join(self.output_dir, "boundaries")
        self.stats_dir = os.path.join(self.output_dir, "statistics")
        
        for directory in [self.alignment_dir, self.chimera_dir, self.boundary_dir, self.stats_dir]:
            os.makedirs(directory, exist_ok=True)
            
        # Check for required external tools
        if not self._check_dependencies():
            raise EnvironmentError("Missing required external tools")
            
        self.initialized = True
        self.logger.info("LTR Enhanced Optimizer initialized successfully")
        return self
    
    def _check_dependencies(self):
        """Check if required external tools are available"""
        required_tools = ["minimap2", "bedtools"]
        optional_tools = ["seqkit", "blast+"]
        
        missing_tools = []
        missing_optional = []
        
        for tool in required_tools:
            if not shutil.which(tool):
                missing_tools.append(tool)
                
        for tool in optional_tools:
            if not shutil.which(tool):
                missing_optional.append(tool)
                
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            return False
            
        if missing_optional:
            self.logger.warning(f"Missing optional tools (some features may be limited): {', '.join(missing_optional)}")
            
        return True
    
    def prepare_genome_index(self):
        """Prepare genome index for efficient alignment"""
        if self.genome_indexed:
            return True
            
        if not self.initialized:
            raise RuntimeError("Optimizer not initialized. Call initialize() first.")
            
        index_file = os.path.join(self.temp_dir, "genome.mmi")
        
        # Create genome index with minimap2
        self.logger.info(f"Creating genome index for {self.genome_file}")
        cmd = ["minimap2", "-d", index_file, self.genome_file]
        
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.genome_index = index_file
            self.genome_indexed = True
            self.logger.info(f"Genome index created: {index_file}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to create genome index: {e.stderr.decode()}")
            return False
    
    def optimize(self, consensus_file=None, genome_file=None, output_dir=None, **kwargs):
        """
        Main method to optimize LTR consensus sequences.
        
        Args:
            consensus_file: Path to consensus sequences (if not provided during initialization)
            genome_file: Path to reference genome (if not provided during initialization)
            output_dir: Directory for output files (if not provided during initialization)
            **kwargs: Additional parameters to override defaults
            
        Returns:
            Path to the optimized consensus file
        """
        # Update parameters if provided
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
                
        # Initialize if not already done
        if not self.initialized:
            self.initialize(consensus_file, genome_file, output_dir)
            
        # Prepare genome index
        if not self.prepare_genome_index():
            raise RuntimeError("Failed to prepare genome index")
            
        self.logger.info("Starting LTR consensus optimization")
        
        # Load genome sequences - efficiently chunked for large genomes
        self.logger.info("Preparing genome sequence access")
        genome_seqs = self._prepare_genome_access()
        
        # Process in iterations for incremental improvement
        current_input = self.consensus_file
        for iteration in range(1, self.iterations + 1):
            self.logger.info(f"Starting iteration {iteration}/{self.iterations}")
            
            # Step 1: Process consensus sequences in batches
            refined_sequences = []
            
            # Process batches in parallel
            with ProcessPoolExecutor(max_workers=max(1, self.threads // 2)) as executor:
                batch_futures = []
                
                for batch_id, batch_file in self._create_sequence_batches(current_input):
                    # Submit batch for processing
                    future = executor.submit(
                        self._process_batch,
                        batch_id, 
                        batch_file,
                        genome_seqs,
                        iteration
                    )
                    batch_futures.append((batch_id, future))
                
                # Collect results as they complete
                for batch_id, future in sorted(batch_futures, key=lambda x: x[0]):
                    try:
                        batch_results = future.result()
                        refined_sequences.extend(batch_results)
                        self.logger.info(f"Completed batch {batch_id} with {len(batch_results)} sequences")
                    except Exception as e:
                        self.logger.error(f"Error processing batch {batch_id}: {str(e)}")
            
            # Write sequences from this iteration
            iter_output = os.path.join(self.output_dir, f"iteration_{iteration}_refined.fa")
            SeqIO.write(refined_sequences, iter_output, "fasta")
            
            self.logger.info(f"Iteration {iteration} completed with {len(refined_sequences)} sequences")
            
            # Prepare for next iteration
            current_input = iter_output
            
        # Create final output
        final_output = os.path.join(self.output_dir, "optimized_consensus.fa")
        shutil.copy(current_input, final_output)
        
        # Generate summary statistics
        self._generate_statistics(final_output)
        
        # Clean up temporary files
        if not self.keep_temp:
            self.logger.info(f"Cleaning up temporary files in {self.temp_dir}")
            shutil.rmtree(self.temp_dir)
            
        self.logger.info(f"LTR consensus optimization completed. Results saved to {final_output}")
        return final_output
    
    def _prepare_genome_access(self):
        """
        Prepare efficient access to genome sequences.
        For large genomes, this uses a dictionary of sequence IDs to file positions,
        loading sequences only when needed.
        
        Returns:
            Object that can be used to access genome sequences
        """
        # Check if seqkit is available for efficient handling
        if shutil.which("seqkit"):
            self.logger.info("Using seqkit for efficient genome handling")
            return self._prepare_genome_index_seqkit()
        else:
            self.logger.warning("seqkit not available. Using simple dictionary for genome access (may use more memory)")
            return self._prepare_genome_simple()
    
    def _prepare_genome_index_seqkit(self):
        """
        Create an efficient genome access index using seqkit.
        This allows us to extract sequences on demand without loading
        the entire genome into memory.
        
        Returns:
            GenomeAccess object
        """
        # Create index file
        index_file = os.path.join(self.temp_dir, "genome.seqkit.idx")
        
        # Run seqkit index
        cmd = ["seqkit", "index", "-j", str(self.threads), self.genome_file]
        
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Get sequence names and lengths
            cmd = ["seqkit", "seq", "--name", "--only-id", self.genome_file]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            seq_names = result.stdout.strip().split('\n')
            
            return GenomeSeqkitAccess(self.genome_file, seq_names, self.chunk_size, self.temp_dir)
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to create seqkit index: {e.stderr.decode()}")
            # Fall back to simple access
            return self._prepare_genome_simple()
    
    def _prepare_genome_simple(self):
        """
        Create a simple genome access object.
        This loads sequences in chunks when needed.
        
        Returns:
            GenomeAccess object
        """
        return GenomeSimpleAccess(self.genome_file, self.chunk_size)
    
    def _create_sequence_batches(self, fasta_file):
        """
        Split sequences into batches for parallel processing.
        
        Args:
            fasta_file: Path to FASTA file with sequences
            
        Yields:
            (batch_id, batch_file) tuples
        """
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        self.logger.info(f"Creating batches from {len(sequences)} sequences")
        
        for i in range(0, len(sequences), self.batch_size):
            batch = sequences[i:i+self.batch_size]
            batch_id = i // self.batch_size + 1
            batch_file = os.path.join(self.temp_dir, f"batch_{batch_id}.fa")
            
            # Write batch to file
            SeqIO.write(batch, batch_file, "fasta")
            
            yield batch_id, batch_file
    
    def _process_batch(self, batch_id, batch_file, genome_seqs, iteration):
        """
        Process a batch of sequences.
        This is designed to be run in parallel.
        
        Args:
            batch_id: Batch identifier
            batch_file: Path to batch FASTA file
            genome_seqs: Genome access object
            iteration: Current iteration number
            
        Returns:
            List of optimized sequence records
        """
        # Configure batch-specific output files
        alignment_file = os.path.join(self.alignment_dir, f"batch_{batch_id}_iter_{iteration}.paf")
        filtered_file = os.path.join(self.alignment_dir, f"batch_{batch_id}_iter_{iteration}.filtered.paf")
        
        # Step 1: Align sequences to genome
        self._align_to_genome(batch_file, alignment_file)
        
        # Step 2: Filter alignments
        alignments_by_query = self._filter_alignments(alignment_file, filtered_file)
        
        # Step 3: Detect and split chimeric sequences
        split_records, chimera_info = self._detect_and_split_chimeras(batch_file, alignments_by_query)
        
        # Write chimera info
        chimera_file = os.path.join(self.chimera_dir, f"batch_{batch_id}_iter_{iteration}_chimeras.tsv")
        with open(chimera_file, "w") as f:
            f.write("seq_id\tis_chimeric\tnum_segments\tbreakpoints\n")
            for seq_id, info in chimera_info.items():
                breakpoints_str = ",".join(map(str, info["breakpoints"])) if info["breakpoints"] else "NA"
                f.write(f"{seq_id}\t{info['is_chimeric']}\t{info['num_segments']}\t{breakpoints_str}\n")
        
        # Step 4: Re-align split sequences if needed
        if any(info["is_chimeric"] for info in chimera_info.values()):
            # Write split sequences to temporary file
            split_file = os.path.join(self.temp_dir, f"batch_{batch_id}_iter_{iteration}_split.fa")
            SeqIO.write(split_records, split_file, "fasta")
            
            # Re-align split sequences
            split_aln_file = os.path.join(self.alignment_dir, f"batch_{batch_id}_iter_{iteration}_split.paf")
            split_filtered_file = os.path.join(self.alignment_dir, f"batch_{batch_id}_iter_{iteration}_split.filtered.paf")
            
            self._align_to_genome(split_file, split_aln_file)
            alignments_by_query = self._filter_alignments(split_aln_file, split_filtered_file)
        
        # Step 5: Detect boundaries for each sequence
        boundary_info = {}
        for seq in split_records:
            seq_alignments = alignments_by_query.get(seq.id, [])
            if seq_alignments:
                boundaries = self._detect_ltr_boundaries(seq, seq_alignments, genome_seqs)
                if boundaries:
                    boundary_info[seq.id] = boundaries
        
        # Write boundary info
        boundary_file = os.path.join(self.boundary_dir, f"batch_{batch_id}_iter_{iteration}_boundaries.tsv")
        with open(boundary_file, "w") as f:
            f.write("seq_id\toriginal_length\tstart\tend\tstart_confidence\tend_confidence\t" +
                   "start_feature\tend_feature\ttsd\n")
            
            for seq_id, boundaries in boundary_info.items():
                f.write(f"{seq_id}\t{boundaries['original_length']}\t" +
                       f"{boundaries['start']}\t{boundaries['end']}\t" +
                       f"{boundaries['start_confidence']:.2f}\t{boundaries['end_confidence']:.2f}\t" +
                       f"{boundaries['start_feature'] or 'NA'}\t{boundaries['end_feature'] or 'NA'}\t" +
                       f"{boundaries['tsd'] or 'NA'}\n")
        
        # Step 6: Refine sequences based on boundaries
        refined_records = []
        for seq in split_records:
            if seq.id in boundary_info:
                # Refine sequence based on detected boundaries
                refined_seq = self._refine_sequence(seq, boundary_info[seq.id])
                # Filter out non-TE regions
                filtered_seq = self._filter_non_te_regions(refined_seq, alignments_by_query.get(seq.id, []))
                refined_records.append(filtered_seq)
            else:
                # Keep original if no boundaries detected
                refined_records.append(seq)
        
        return refined_records
    
    def _align_to_genome(self, fasta_file, output_file):
        """
        Align sequences to the genome using minimap2.
        
        Args:
            fasta_file: Path to input FASTA file
            output_file: Path to output PAF file
            
        Returns:
            True if successful, False otherwise
        """
        # Ensure genome is indexed
        if not self.genome_indexed:
            if not self.prepare_genome_index():
                return False
        
        # Run minimap2 alignment
        cmd = [
            "minimap2",
            "-cx", "map-ont",  # Use more sensitive mode for consensus sequences
            "--secondary=yes",  # Allow secondary alignments
            "-p", "0.5",  # Lower secondary alignment threshold
            "-N", "50",   # Report up to 50 alignments per sequence
            "--cs",       # Output CIGAR string for detailed analysis
            "-t", str(self.threads),
            self.genome_index,
            fasta_file
        ]
        
        try:
            with open(output_file, "w") as f:
                subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE)
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Alignment failed: {e.stderr.decode() if e.stderr else str(e)}")
            return False
    
    def _filter_alignments(self, paf_file, filtered_file):
        """
        Filter alignments based on identity and coverage thresholds.
        
        Args:
            paf_file: Path to PAF alignment file
            filtered_file: Path to output filtered PAF file
            
        Returns:
            Dictionary mapping query IDs to their alignments
        """
        # Parse PAF file and filter alignments
        alignments_by_query = defaultdict(list)
        filtered_alignments = []
        total_count = 0
        passed_count = 0
        
        with open(paf_file, "r") as f:
            for line in f:
                total_count += 1
                fields = line.strip().split("\t")
                
                # Basic PAF fields
                query_name = fields[0]
                query_len = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                strand = fields[4]
                target_name = fields[5]
                target_len = int(fields[6])
                target_start = int(fields[7])
                target_end = int(fields[8])
                matches = int(fields[9])
                aln_len = int(fields[10])
                mapq = int(fields[11])
                
                # Calculate identity and coverage
                identity = (matches / aln_len) * 100 if aln_len > 0 else 0
                coverage = ((query_end - query_start) / query_len) * 100
                
                # Filter based on thresholds
                if identity >= self.min_identity and coverage >= self.min_coverage:
                    passed_count += 1
                    filtered_alignments.append(line)
                    
                    # Create alignment record
                    alignment = {
                        'query_name': query_name,
                        'query_len': query_len,
                        'query_start': query_start,
                        'query_end': query_end,
                        'strand': strand,
                        'target_name': target_name,
                        'target_len': target_len,
                        'target_start': target_start,
                        'target_end': target_end,
                        'matches': matches,
                        'aln_len': aln_len,
                        'mapq': mapq,
                        'identity': identity,
                        'coverage': coverage
                    }
                    
                    # Additional fields (CIGAR, etc.) if available
                    for i in range(12, len(fields)):
                        if fields[i].startswith("cg:Z:"):
                            alignment['cigar'] = fields[i][5:]
                        elif fields[i].startswith("cs:Z:"):
                            alignment['cs'] = fields[i][5:]
                    
                    alignments_by_query[query_name].append(alignment)
        
        # Write filtered alignments to file
        with open(filtered_file, "w") as f:
            f.writelines(filtered_alignments)
            
        self.logger.info(f"Filtered alignments: {passed_count}/{total_count} passed thresholds")
        return alignments_by_query
    
    def _detect_and_split_chimeras(self, fasta_file, alignments_by_query):
        """
        Detect and split chimeric consensus sequences.
        
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
            is_chimeric, breakpoints = self._analyze_chimeric_pattern(seq_record, seq_alignments)
            
            # Record chimera information
            chimera_info[seq_id] = {
                'is_chimeric': is_chimeric,
                'breakpoints': breakpoints,
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
        Analyze alignment patterns to detect chimeric sequences.
        This enhanced version uses multiple evidence types:
        1. Coverage discontinuities
        2. Target chromosome changes
        3. Alignment orientation changes
        
        Args:
            seq_record: SeqRecord object
            alignments: List of alignments for this sequence
            
        Returns:
            (is_chimeric, breakpoints) tuple
        """
        seq_len = len(seq_record.seq)
        
        # Not enough alignments to analyze
        if len(alignments) < 3:
            return False, []
            
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
        window_size = min(100, seq_len // 10)  # Adaptive window size
        if window_size < 20:  # Minimum window size
            window_size = 20
            
        threshold = self.chimera_threshold  # Threshold for significant change
        
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
            
            if change_score > threshold:
                # Significant change detected
                breakpoints.append(i)
                
        # Merge nearby breakpoints
        if breakpoints:
            merged_breakpoints = []
            current_group = [breakpoints[0]]
            
            for bp in breakpoints[1:]:
                if bp - current_group[-1] < window_size:
                    # Add to current group
                    current_group.append(bp)
                else:
                    # Start new group
                    merged_breakpoints.append(sum(current_group) // len(current_group))
                    current_group = [bp]
                    
            # Add last group
            if current_group:
                merged_breakpoints.append(sum(current_group) // len(current_group))
                
            # Filter breakpoints by minimum segment length
            valid_breakpoints = []
            last_pos = 0
            
            for bp in sorted(merged_breakpoints):
                if bp - last_pos >= self.min_segment_length:
                    valid_breakpoints.append(bp)
                    last_pos = bp
                    
            # Check if last segment is long enough
            if seq_len - last_pos < self.min_segment_length:
                if valid_breakpoints:  # Only pop if there are breakpoints
                    valid_breakpoints.pop()
                
            # If we still have valid breakpoints, it's chimeric
            return len(valid_breakpoints) > 0, valid_breakpoints
        
        return False, []
    
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
    
    def _detect_ltr_boundaries(self, seq_record, alignments, genome_seqs):
        """
        Enhanced method to detect LTR boundaries using multiple evidence types:
        1. Alignment boundaries
        2. Structure motifs (TG...CA)
        3. Target Site Duplications (TSDs)
        4. PBS/PPT features
        
        Args:
            seq_record: SeqRecord object
            alignments: List of alignments for this sequence
            genome_seqs: Genome access object
            
        Returns:
            Dictionary with boundary information or None
        """
        seq_len = len(seq_record.seq)
        seq_id = seq_record.id
        
        # Not enough alignments
        if len(alignments) < 2:
            return None
            
        # Analyze alignment boundaries
        starts = []
        ends = []
        
        # Track motif evidence
        motif5_positions = []
        motif3_positions = []
        tsd_evidence = []
        
        # Examine each alignment
        for aln in alignments:
            q_start = aln['query_start']
            q_end = aln['query_end']
            
            starts.append(q_start)
            ends.append(q_end)
            
            # Check for 5' and 3' motifs at boundaries
            seq_str = str(seq_record.seq)
            
            # Check 5' end
            for motif in self.ltr5_motifs:
                motif_pos = seq_str.find(motif, q_start, q_start + 10)
                if motif_pos >= 0:
                    motif5_positions.append(motif_pos)
                    
            # Check 3' end
            for motif in self.ltr3_motifs:
                motif_pos = seq_str.rfind(motif, q_end - 10, q_end)
                if motif_pos >= 0:
                    motif3_positions.append(motif_pos + len(motif))  # End position after motif
                    
            # Try to find TSDs from genome context
            tsd = self._find_tsd_from_alignment(aln, genome_seqs)
            if tsd:
                tsd_evidence.append(tsd)
        
        # Calculate consensus boundaries
        start_counter = Counter(starts)
        end_counter = Counter(ends)
        motif5_counter = Counter(motif5_positions)
        motif3_counter = Counter(motif3_positions)
        
        # Find most common boundaries
        best_start = start_counter.most_common(1)[0][0] if start_counter else 0
        best_end = end_counter.most_common(1)[0][0] if end_counter else seq_len
        
        # Adjust boundaries based on motif evidence if available
        if motif5_counter:
            best_motif5 = motif5_counter.most_common(1)[0][0]
            # If motif is close to alignment boundary, use it
            if abs(best_motif5 - best_start) <= 20:
                best_start = best_motif5
                
        if motif3_counter:
            best_motif3 = motif3_counter.most_common(1)[0][0]
            # If motif is close to alignment boundary, use it
            if abs(best_motif3 - best_end) <= 20:
                best_end = best_motif3
                
        # Calculate confidence scores
        total_alns = len(alignments)
        start_confidence = start_counter[best_start] / total_alns if total_alns > 0 else 0
        end_confidence = end_counter[best_end] / total_alns if total_alns > 0 else 0
        
        # Search for internal features (PBS/PPT)
        internal_features = self._find_internal_features(seq_record, best_start, best_end)
        
        # Use internal features to further refine boundaries if possible
        if internal_features['pbs'] and best_start < internal_features['pbs']['position']:
            # PBS should be just inside the 5' LTR, so adjust if needed
            if internal_features['pbs']['position'] - best_start > 20:
                # PBS is too far from current start, might need adjustment
                pbs_pos = internal_features['pbs']['position']
                if pbs_pos > 20:
                    # Look for LTR start motifs near PBS
                    for motif in self.ltr5_motifs:
                        motif_pos = str(seq_record.seq).rfind(motif, max(0, pbs_pos - 30), pbs_pos - 5)
                        if motif_pos >= 0:
                            best_start = motif_pos
                            break
        
        if internal_features['ppt'] and best_end > internal_features['ppt']['position']:
            # PPT should be just before the 3' LTR, so adjust if needed
            if best_end - internal_features['ppt']['position'] > 20:
                # PPT is too far from current end, might need adjustment
                ppt_pos = internal_features['ppt']['position'] + len(internal_features['ppt']['sequence'])
                if ppt_pos < seq_len - 20:
                    # Look for LTR end motifs near PPT
                    for motif in self.ltr3_motifs:
                        motif_pos = str(seq_record.seq).find(motif, ppt_pos, ppt_pos + 30)
                        if motif_pos >= 0:
                            best_end = motif_pos + len(motif)
                            break
        
        # Determine LTR structural features
        start_feature = None
        end_feature = None
        
        for motif in self.ltr5_motifs:
            if str(seq_record.seq[best_start:best_start+len(motif)]).upper() == motif:
                start_feature = motif
                break
                
        for motif in self.ltr3_motifs:
            end_pos = best_end - len(motif)
            if end_pos >= 0 and str(seq_record.seq[end_pos:best_end]).upper() == motif:
                end_feature = motif
                break
                
        # Summarize TSD evidence
        tsd = None
        if tsd_evidence:
            tsd_counter = Counter(tsd_evidence)
            best_tsd, count = tsd_counter.most_common(1)[0]
            if count >= 2:  # Require at least 2 supporting alignments
                tsd = best_tsd
                
        # Validate boundaries
        if best_start >= best_end:
            return None
            
        # Check if boundaries define a reasonable LTR length
        if best_end - best_start < self.min_ltr_len:
            return None
            
        if best_end - best_start > self.max_ltr_len:
            # Try to find the true LTR sequences within this range
            # Look for similar terminal repeat patterns
            ltr_start, ltr_end = self._find_terminal_repeats(seq_record, best_start, best_end)
            if ltr_start is not None and ltr_end is not None:
                best_start = ltr_start
                best_end = ltr_end
        
        # Final boundaries check
        if best_end - best_start < self.min_ltr_len or best_end - best_start > self.max_ltr_len:
            return None
                
        return {
            'seq_id': seq_id,
            'original_length': seq_len,
            'start': best_start,
            'end': best_end,
            'start_confidence': start_confidence,
            'end_confidence': end_confidence,
            'start_feature': start_feature,
            'end_feature': end_feature,
            'tsd': tsd,
            'internal_features': internal_features
        }

    def _find_terminal_repeats(self, seq_record, start, end):
        """
        Find terminal repeats within a sequence range to identify true LTR boundaries.
        LTRs typically have similar sequences at both ends.
        
        Args:
            seq_record: SeqRecord object
            start: Start position
            end: End position
            
        Returns:
            (ltr_start, ltr_end) tuple or (None, None) if not found
        """
        seq = seq_record.seq[start:end]
        seq_str = str(seq)
        
        # Try to identify similar patterns at the beginning and end
        for terminal_len in range(100, 500, 50):  # Try different lengths
            if len(seq) < terminal_len * 2:
                continue
                
            # Get terminal sequences
            left_terminal = seq_str[:terminal_len]
            right_terminal = seq_str[-terminal_len:]
            
            # Calculate similarity
            matches = sum(a == b for a, b in zip(left_terminal, right_terminal))
            similarity = matches / terminal_len
            
            if similarity > 0.7:  # At least 70% similarity
                return start, end
                
        # If no repeats found, try to look for LTR motifs
        left_boundary = None
        right_boundary = None
        
        # Look for 5' motifs near the start
        for i in range(start, min(start + 200, end)):
            for motif in self.ltr5_motifs:
                if str(seq_record.seq[i:i+len(motif)]).upper() == motif:
                    left_boundary = i
                    break
            if left_boundary is not None:
                break
                
        # Look for 3' motifs near the end
        for i in range(max(left_boundary + 100 if left_boundary is not None else start, end - 200), end):
            for motif in self.ltr3_motifs:
                end_pos = i - len(motif)
                if end_pos >= 0 and str(seq_record.seq[end_pos:i]).upper() == motif:
                    right_boundary = i
                    break
            if right_boundary is not None:
                break
                
        if left_boundary is not None and right_boundary is not None:
            return left_boundary, right_boundary
            
        return None, None
    
    def _find_tsd_from_alignment(self, alignment, genome_seqs):
        """
        Try to find Target Site Duplication (TSD) from genomic context.
        
        Args:
            alignment: Alignment information
            genome_seqs: Genome access object
            
        Returns:
            TSD sequence if found, None otherwise
        """
        target_name = alignment['target_name']
        target_start = alignment['target_start']
        target_end = alignment['target_end']
        
        # Check if we can access the target sequence
        if not genome_seqs.has_sequence(target_name):
            return None
            
        # Extract flanking sequences
        left_flank_start = max(0, target_start - self.flanking_seq)
        left_flank = genome_seqs.get_sequence_region(target_name, left_flank_start, target_start)
        
        right_flank_end = target_end + self.flanking_seq
        right_flank = genome_seqs.get_sequence_region(target_name, target_end, right_flank_end)
        
        if not left_flank or not right_flank:
            return None
            
        # Look for TSD patterns
        best_tsd = None
        best_score = 0
        
        for tsd_len in range(self.tsd_min, self.tsd_max + 1):
            # Check if flanks are long enough
            if len(left_flank) < tsd_len or len(right_flank) < tsd_len:
                continue
                
            # Compare end of left flank to start of right flank
            left_tsd = left_flank[-tsd_len:].upper()
            right_tsd = right_flank[:tsd_len].upper()
            
            # Count matches
            matches = sum(1 for a, b in zip(left_tsd, right_tsd) if a == b)
            score = matches / tsd_len
            
            if score > best_score and score >= 0.8:  # At least 80% identity
                best_score = score
                best_tsd = left_tsd
                
        return best_tsd
    
    def _find_internal_features(self, seq_record, start, end):
        """
        Find LTR internal features like PBS and PPT.
        
        Args:
            seq_record: SeqRecord object
            start: Start position
            end: End position
            
        Returns:
            Dictionary with internal feature information
        """
        # PBS (Primer Binding Site) patterns
        pbs_patterns = [
            'TGGCGCCCAACGTGGGGC',  # tRNAPro
            'TGGCGCCGTAACAGGGAC',  # tRNATrp
            'TGGCGCCCGAACAGGGAC',  # tRNAGln
            'TGGCGCCCAACCTGGGA',   # tRNALys
            'TGGTAGCAGAGCTGGGAA'   # tRNAIle
        ]
        
        # PPT (PolyPurine Tract) pattern - A/G rich sequences
        ppt_pattern = r'[AG]{8,15}'
        
        seq_str = str(seq_record.seq)
        
        # Look for PBS near the 5' LTR (5-20 bp after LTR)
        pbs_result = None
        for pbs_pos in range(start + 5, min(start + 30, end - 100, len(seq_str) - 18)):
            window = seq_str[pbs_pos:pbs_pos+18]
            
            for pattern in pbs_patterns:
                # Calculate similarity
                matches = sum(1 for a, b in zip(window.upper(), pattern) if a == b)
                similarity = matches / len(pattern)
                
                if similarity >= 0.8:  # At least 80% match
                    pbs_result = {
                        'position': pbs_pos,
                        'sequence': window,
                        'pattern': pattern,
                        'similarity': similarity
                    }
                    break
                    
            if pbs_result:
                break
                
        # Look for PPT near the 3' LTR (5-20 bp before LTR)
        ppt_result = None
        if end > 20:
            # Use regex to find PPT patterns
            import re
            ppt_matches = list(re.finditer(ppt_pattern, seq_str[max(start + 100, end - 50):end]))
            
            if ppt_matches:
                # Use the last one (closest to 3' LTR)
                match = ppt_matches[-1]
                match_start = max(start + 100, end - 50) + match.start()
                
                ppt_result = {
                    'position': match_start,
                    'sequence': match.group(),
                    'length': len(match.group())
                }
                
        return {
            'pbs': pbs_result,
            'ppt': ppt_result
        }
    
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
    
    def _filter_non_te_regions(self, seq_record, alignments):
        """
        Filter out regions that don't appear to be part of the transposon.
        This method uses alignment patterns and coding potential.
        
        Args:
            seq_record: SeqRecord object
            alignments: List of alignments for this sequence
            
        Returns:
            Filtered SeqRecord object
        """
        seq_len = len(seq_record.seq)
        
        # Not enough alignments to analyze
        if len(alignments) < 3:
            return seq_record
            
        # Create coverage profile
        coverage = np.zeros(seq_len)
        
        # Fill coverage profile
        for aln in alignments:
            q_start = max(0, aln['query_start'])
            q_end = min(seq_len, aln['query_end'])
            
            coverage[q_start:q_end] += 1
            
        # Find low-coverage regions
        # First, smooth the coverage profile
        window_size = min(50, seq_len // 10)
        if window_size < 10:
            window_size = 10
            
        smoothed = np.zeros(seq_len)
        for i in range(seq_len):
            start = max(0, i - window_size // 2)
            end = min(seq_len, i + window_size // 2)
            smoothed[i] = np.mean(coverage[start:end])
            
        # Calculate coverage threshold
        mean_coverage = np.mean(smoothed)
        threshold = mean_coverage * 0.3  # 30% of mean coverage
        
        # Find regions to keep (above threshold)
        keep_regions = []
        in_region = False
        region_start = 0
        
        for i in range(seq_len):
            if smoothed[i] >= threshold:
                if not in_region:
                    in_region = True
                    region_start = i
            else:
                if in_region:
                    in_region = False
                    # Only keep regions of reasonable size
                    if i - region_start >= 100:
                        keep_regions.append((region_start, i))
                        
        # Handle last region
        if in_region and seq_len - region_start >= 100:
            keep_regions.append((region_start, seq_len))
            
        # If no regions to keep, return original
        if not keep_regions:
            return seq_record
            
        # If everything is kept, return original
        if len(keep_regions) == 1 and keep_regions[0][0] == 0 and keep_regions[0][1] == seq_len:
            return seq_record
            
        # Create filtered sequence
        filtered_seq = ""
        for start, end in keep_regions:
            filtered_seq += str(seq_record.seq[start:end])
            
        # Create new record
        regions_str = ",".join([f"{start+1}-{end}" for start, end in keep_regions])
        filtered_record = SeqRecord(
            Seq(filtered_seq),
            id=seq_record.id,
            name=seq_record.name,
            description=f"{seq_record.description} | Filtered regions: {regions_str}"
        )
        
        return filtered_record
    
    def _generate_statistics(self, final_file):
        """
        Generate summary statistics for the optimization process.
        
        Args:
            final_file: Path to the final optimized consensus file
        """
        # Compare original and optimized sequences
        original_seqs = list(SeqIO.parse(self.consensus_file, "fasta"))
        optimized_seqs = list(SeqIO.parse(final_file, "fasta"))
        
        # Create sequence dictionaries for lookup
        original_dict = {record.id.split('_segment_')[0]: record for record in original_seqs}
        
        # Collect statistics
        stats = []
        for record in optimized_seqs:
            # Handle segment IDs
            base_id = record.id.split('_segment_')[0]
            original = original_dict.get(base_id)
            
            if original:
                stat = {
                    'seq_id': record.id,
                    'base_id': base_id,
                    'original_length': len(original.seq),
                    'optimized_length': len(record.seq),
                    'difference': len(original.seq) - len(record.seq),
                    'percent_change': ((len(original.seq) - len(record.seq)) / len(original.seq) * 100) 
                              if len(original.seq) > 0 else 0
                }
                stats.append(stat)
        
        # Calculate summary statistics
        total_original = sum(len(record.seq) for record in original_seqs)
        total_optimized = sum(len(record.seq) for record in optimized_seqs)
        
        avg_original = total_original / len(original_seqs) if original_seqs else 0
        avg_optimized = total_optimized / len(optimized_seqs) if optimized_seqs else 0
        
        # Save summary
        summary_file = os.path.join(self.stats_dir, "optimization_summary.txt")
        with open(summary_file, "w") as f:
            f.write(f"Original sequences: {len(original_seqs)}\n")
            f.write(f"Optimized sequences: {len(optimized_seqs)}\n")
            f.write(f"Total original length: {total_original} bp\n")
            f.write(f"Total optimized length: {total_optimized} bp\n")
            f.write(f"Average original length: {avg_original:.2f} bp\n")
            f.write(f"Average optimized length: {avg_optimized:.2f} bp\n")
            f.write(f"Total reduction: {total_original - total_optimized} bp ({(total_original - total_optimized) / total_original * 100:.2f}%)\n")
            f.write(f"Chimeric sequences split: {len(optimized_seqs) - len(original_seqs)}\n")
        
        # Save detailed statistics
        stats_file = os.path.join(self.stats_dir, "sequence_statistics.tsv")
        with open(stats_file, "w") as f:
            f.write("seq_id\tbase_id\toriginal_length\toptimized_length\tdifference\tpercent_change\n")
            for stat in stats:
                f.write(f"{stat['seq_id']}\t{stat['base_id']}\t{stat['original_length']}\t" +
                       f"{stat['optimized_length']}\t{stat['difference']}\t{stat['percent_change']:.2f}\n")
                
        self.logger.info(f"Statistics saved to {self.stats_dir}")


class GenomeAccess:
    """Base class for genome sequence access"""
    def has_sequence(self, seq_id):
        raise NotImplementedError
        
    def get_sequence_region(self, seq_id, start, end):
        raise NotImplementedError


class GenomeSeqkitAccess(GenomeAccess):
    """Efficient genome access using seqkit"""
    def __init__(self, genome_file, seq_names, chunk_size, temp_dir):
        self.genome_file = genome_file
        self.seq_names = set(seq_names)
        self.chunk_size = chunk_size
        self.temp_dir = temp_dir
        self.cache = {}
        
    def has_sequence(self, seq_id):
        return seq_id in self.seq_names
        
    def get_sequence_region(self, seq_id, start, end):
        if not self.has_sequence(seq_id):
            return None
            
        # Check if in cache
        if seq_id in self.cache:
            cached_seq = self.cache[seq_id]
            if start >= 0 and end <= len(cached_seq):
                return cached_seq[start:end]
                
        # Extract sequence region using seqkit
        region = f"{seq_id}:{start+1}-{end}"
        temp_file = os.path.join(self.temp_dir, f"temp_region_{int(time.time())}_{random.randint(1000, 9999)}.fa")
        
        try:
            cmd = ["seqkit", "subseq", "-r", region, self.genome_file]
            with open(temp_file, "w") as f:
                subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE)
                
            # Read sequence
            for record in SeqIO.parse(temp_file, "fasta"):
                seq_str = str(record.seq)
                # Cache small sequences
                if len(seq_str) <= 10000:
                    self.cache[seq_id] = seq_str
                return seq_str
                
            return None
            
        except Exception:
            return None
            
        finally:
            # Clean up
            if os.path.exists(temp_file):
                os.remove(temp_file)


class GenomeSimpleAccess(GenomeAccess):
    """Simple genome access using BioPython"""
    def __init__(self, genome_file, chunk_size):
        self.genome_file = genome_file
        self.chunk_size = chunk_size
        self.cache = {}
        self._scan_genome()
        
    def _scan_genome(self):
        """Scan genome file to build index of sequence IDs"""
        self.seq_dict = {}
        self.seq_names = set()
        
        try:
            from Bio import SeqIO
            for record in SeqIO.parse(self.genome_file, "fasta"):
                self.seq_names.add(record.id)
                # For small sequences, cache them fully
                if len(record.seq) <= 10000:
                    self.cache[record.id] = str(record.seq)
        except Exception as e:
            raise RuntimeError(f"Error scanning genome file: {str(e)}")
            
    def has_sequence(self, seq_id):
        return seq_id in self.seq_names
        
    def get_sequence_region(self, seq_id, start, end):
        if not self.has_sequence(seq_id):
            return None
            
        # Check if in cache
        if seq_id in self.cache:
            cached_seq = self.cache[seq_id]
            if start >= 0 and end <= len(cached_seq):
                return cached_seq[start:end]
                
        # Parse genome to find sequence
        for record in SeqIO.parse(self.genome_file, "fasta"):
            if record.id == seq_id:
                return str(record.seq[start:end])
                
        return None


def main():
    """Command-line entry point"""
    parser = argparse.ArgumentParser(
        description="Enhanced LTR Consensus Optimizer: Refines LTR consensus sequences by "
                   "detecting precise boundaries, resolving chimeric sequences, "
                   "and removing non-TE regions. Optimized for GB-scale genomes."
    )
    
    parser.add_argument("consensus_file", help="Path to consensus sequences from Refiner_for_LTR")
    parser.add_argument("genome_file", help="Path to reference genome")
    parser.add_argument("output_dir", help="Directory for output files")
    
    parser.add_argument("--threads", type=int, default=1,
                       help="Number of CPU threads to use (default: 1)")
    parser.add_argument("--min-identity", type=float, default=80,
                       help="Minimum identity percentage for valid alignments (default: 80)")
    parser.add_argument("--min-coverage", type=float, default=80,
                       help="Minimum coverage percentage for valid alignments (default: 80)")
    parser.add_argument("--min-ltr-len", type=int, default=100,
                       help="Minimum length of an LTR (default: 100)")
    parser.add_argument("--max-ltr-len", type=int, default=1500,
                       help="Maximum length of an LTR (default: 1500)")
    parser.add_argument("--chimera-threshold", type=float, default=0.3,
                       help="Threshold for chimera detection (default: 0.3)")
    parser.add_argument("--min-segment-length", type=int, default=500,
                       help="Minimum length for segments when splitting chimeras (default: 500)")
    parser.add_argument("--iterations", type=int, default=2,
                       help="Number of optimization iterations (default: 2)")
    parser.add_argument("--chunk-size", type=int, default=10,
                       help="Number of genome sequences to process at once (default: 10)")
    parser.add_argument("--batch-size", type=int, default=50,
                       help="Number of consensus sequences to process in each batch (default: 50)")
    parser.add_argument("--keep-temp", action="store_true",
                       help="Keep temporary files")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       default="INFO", help="Logging level")
    
    args = parser.parse_args()
    
    # Initialize optimizer
    optimizer = LTREnhancedOptimizer(
        args.consensus_file,
        args.genome_file,
        args.output_dir,
        threads=args.threads,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        min_ltr_len=args.min_ltr_len,
        max_ltr_len=args.max_ltr_len,
        chimera_threshold=args.chimera_threshold,
        min_segment_length=args.min_segment_length,
        iterations=args.iterations,
        chunk_size=args.chunk_size,
        batch_size=args.batch_size,
        keep_temp=args.keep_temp,
        log_level=args.log_level
    )
    
    # Run optimization
    optimizer.optimize()


if __name__ == "__main__":
    main()
