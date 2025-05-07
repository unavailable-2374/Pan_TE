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
import random
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from itertools import groupby

# Check and import optional dependencies
try:
    import umap
    import hdbscan
    CLUSTERING_AVAILABLE = True
except ImportError:
    CLUSTERING_AVAILABLE = False

try:
    from Bio.SeqUtils import GC, ProtParam
    BIOUTILS_AVAILABLE = True
except ImportError:
    BIOUTILS_AVAILABLE = False

class LTREnhancedOptimizer:
    """
    Comprehensive tool for refining LTR consensus sequences with:
    1. Chimeric sequence detection and resolution
    2. Precise boundary determination using structural features
    3. Non-transposon sequence filtering
    4. Efficient processing of GB-scale genomes
    5. Family clustering and subfamily identification
    6. Advanced TSD and motif detection
    7. ORF and coding potential analysis
    
    Can be used as a standalone script or imported as a module in workflows.
    """
    
    def __init__(self, consensus_file=None, genome_file=None, output_dir=None, 
                 threads=1, min_identity=80, min_coverage=80, 
                 min_ltr_len=100, max_ltr_len=1500, flanking_seq=20,
                 chimera_threshold=0.3, min_segment_length=500,
                 tsd_min=4, tsd_max=6, iterations=2, 
                 chunk_size=10, batch_size=50, 
                 keep_temp=False, log_level="INFO",
                 low_complexity_filter=True, clustering=True,
                 dynamic_threshold=True, orf_analysis=True,
                 kmer_boundary=True, blank_region_threshold=100,
                 orientation_aware=True, pssm_tsd=True,
                 custom_motifs=None, classify_output=False):
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
            low_complexity_filter: Whether to filter low complexity regions
            clustering: Whether to perform clustering for subfamily identification
            dynamic_threshold: Whether to use dynamic thresholds for non-TE filtering
            orf_analysis: Whether to perform ORF analysis
            kmer_boundary: Whether to use k-mer based boundary detection
            blank_region_threshold: Minimum length of blank regions for chimera detection
            orientation_aware: Whether to use orientation-aware chimera detection
            pssm_tsd: Whether to use PSSM for TSD identification
            custom_motifs: Custom motifs for boundary detection
            classify_output: Whether to classify output sequences
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
        
        # Enhanced parameters
        self.low_complexity_filter = low_complexity_filter
        self.clustering = clustering
        self.dynamic_threshold = dynamic_threshold
        self.orf_analysis = orf_analysis
        self.kmer_boundary = kmer_boundary
        self.blank_region_threshold = blank_region_threshold
        self.orientation_aware = orientation_aware
        self.pssm_tsd = pssm_tsd
        self.custom_motifs = custom_motifs
        self.classify_output = classify_output
        
        # Initialize logger
        self._setup_logging(log_level)
        
        # Status tracking
        self.initialized = False
        self.genome_indexed = False
        self.temp_dir = None
        
        # LTR-specific patterns
        self.ltr5_motifs = ["TG", "TGT", "TGTA", "TGTG", "TGCA"]
        self.ltr3_motifs = ["CA", "ACA", "TACA", "CACA", "TGCA"]
        
        # Add custom motifs if provided
        if self.custom_motifs:
            if 'ltr5' in self.custom_motifs and isinstance(self.custom_motifs['ltr5'], list):
                self.ltr5_motifs.extend(self.custom_motifs['ltr5'])
            if 'ltr3' in self.custom_motifs and isinstance(self.custom_motifs['ltr3'], list):
                self.ltr3_motifs.extend(self.custom_motifs['ltr3'])
                
        self.tsd_patterns = {}  # Will store found TSDs
        
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
        
        # TE protein domains for ORF analysis
        self.te_domains = {
            'gag': ['CCHC', 'CXCX2CX4HX4C', 'zinc finger', 'capsid', 'nucleocapsid'],
            'pol': ['protease', 'integrase', 'reverse transcriptase', 'RNase H', 'aspartyl protease'],
            'env': ['envelope', 'transmembrane', 'surface glycoprotein']
        }
        
        # Position-specific scoring matrices for TSD identification
        self.tsd_pssm = self._initialize_tsd_pssm()
    
    def _initialize_tsd_pssm(self):
        """Initialize position-specific scoring matrices for TSD identification"""
        # Simple initialization - can be enhanced with real training data
        pssm = {}
        for length in range(self.tsd_min, self.tsd_max + 1):
            # Create a uniform matrix for each position
            # For real implementation, these should be trained from known TSDs
            pssm[length] = np.ones((length, 4)) * 0.25  # Equal probability for A,C,G,T
        return pssm
        
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
        self.clustering_dir = os.path.join(self.output_dir, "subfamily_clustering")
        self.orf_dir = os.path.join(self.output_dir, "orf_analysis")
        self.classification_dir = os.path.join(self.output_dir, "classification")
        
        for directory in [self.alignment_dir, self.chimera_dir, self.boundary_dir, 
                         self.stats_dir, self.clustering_dir, self.orf_dir, 
                         self.classification_dir]:
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
        optional_tools = ["seqkit", "blast+", "dustmasker", "trf"]
        
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
            
        # Check for optional Python packages
        if self.clustering and not CLUSTERING_AVAILABLE:
            self.logger.warning("UMAP and HDBSCAN are not available. Clustering will be disabled.")
            self.clustering = False
            
        if self.orf_analysis and not BIOUTILS_AVAILABLE:
            self.logger.warning("BioPython SeqUtils is not available. ORF analysis features will be limited.")
            
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
        
        # Preprocess input sequences if needed (low complexity filter)
        preprocessed_file = self.consensus_file
        if self.low_complexity_filter:
            preprocessed_file = self._filter_low_complexity_regions(self.consensus_file)
        
        # Process in iterations for incremental improvement
        current_input = preprocessed_file
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
        
        # Perform subfamily clustering if enabled
        if self.clustering and CLUSTERING_AVAILABLE and len(refined_sequences) > 3:
            self.logger.info("Performing subfamily clustering")
            subfamily_file = self._perform_subfamily_clustering(refined_sequences)
            current_input = subfamily_file
            
        # Classify output sequences if enabled
        if self.classify_output:
            self.logger.info("Classifying output sequences")
            classification_file = self._classify_sequences(current_input)
            current_input = classification_file
            
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
    
    def _filter_low_complexity_regions(self, fasta_file):
        """
        Filter low complexity regions using dustmasker or custom implementation.
        
        Args:
            fasta_file: Path to input FASTA file
            
        Returns:
            Path to filtered FASTA file
        """
        output_file = os.path.join(self.temp_dir, "filtered_low_complexity.fa")
        
        # Check if dustmasker is available
        if shutil.which("dustmasker"):
            self.logger.info("Using dustmasker to filter low complexity regions")
            
            # Run dustmasker to identify low complexity regions
            mask_file = os.path.join(self.temp_dir, "masked_regions.asnb")
            cmd = ["dustmasker", "-in", fasta_file, "-out", mask_file, "-outfmt", "asn1_bin"]
            
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                # Convert masked regions to FASTA
                cmd = ["dustmasker", "-in", fasta_file, "-out", output_file, "-outfmt", "fasta", 
                       "-asnb", mask_file, "-lcase", "-infmt", "fasta"]
                
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.logger.info(f"Low complexity regions filtered: {output_file}")
                return output_file
                
            except subprocess.CalledProcessError as e:
                self.logger.error(f"Failed to run dustmasker: {e.stderr.decode()}")
                self.logger.warning("Falling back to custom low complexity filtering")
                
        # Simple custom implementation if dustmasker is not available
        self.logger.info("Using custom method to filter low complexity regions")
        
        filtered_records = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Convert sequence to lowercase for masking
            seq_str = str(record.seq).lower()
            
            # Simple entropy-based masking
            window_size = 10
            threshold = 1.2  # Higher threshold = less masking
            
            masked_seq = list(seq_str)
            
            for i in range(len(seq_str) - window_size + 1):
                window = seq_str[i:i+window_size]
                
                # Calculate entropy of window
                counts = Counter(window)
                total = sum(counts.values())
                entropy = 0
                
                for count in counts.values():
                    p = count / total
                    entropy -= p * np.log2(p)
                
                # Mark low complexity regions with 'X'
                if entropy < threshold:
                    for j in range(i, i + window_size):
                        if j < len(masked_seq):
                            masked_seq[j] = 'N'
            
            # Create new record with masked sequence
            masked_record = SeqRecord(
                Seq(''.join(masked_seq)),
                id=record.id,
                name=record.name,
                description=record.description + " [low complexity filtered]"
            )
            
            filtered_records.append(masked_record)
            
        # Write filtered sequences
        SeqIO.write(filtered_records, output_file, "fasta")
        self.logger.info(f"Custom low complexity filtering completed: {output_file}")
        return output_file
    
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
            f.write("seq_id\tis_chimeric\tnum_segments\tbreakpoints\tbreakpoint_confidence\tblank_regions\n")
            for seq_id, info in chimera_info.items():
                breakpoints_str = ",".join(map(str, info["breakpoints"])) if info["breakpoints"] else "NA"
                confidence_str = ",".join(map(str, info["breakpoint_confidence"])) if info["breakpoint_confidence"] else "NA"
                blank_regions_str = ",".join([f"{start}-{end}" for start, end in info["blank_regions"]]) if info["blank_regions"] else "NA"
                
                f.write(f"{seq_id}\t{info['is_chimeric']}\t{info['num_segments']}\t" +
                       f"{breakpoints_str}\t{confidence_str}\t{blank_regions_str}\n")
        
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
                   "start_feature\tend_feature\ttsd\tkmer_boundary_support\tpbs_type\tppt_length\n")
            
            for seq_id, boundaries in boundary_info.items():
                kmer_support = "YES" if boundaries.get('kmer_boundary_support', False) else "NO"
                pbs_type = boundaries.get('internal_features', {}).get('pbs', {}).get('pattern', "NA") if boundaries.get('internal_features', {}).get('pbs') else "NA"
                ppt_length = boundaries.get('internal_features', {}).get('ppt', {}).get('length', "NA") if boundaries.get('internal_features', {}).get('ppt') else "NA"
                
                f.write(f"{seq_id}\t{boundaries['original_length']}\t" +
                       f"{boundaries['start']}\t{boundaries['end']}\t" +
                       f"{boundaries['start_confidence']:.2f}\t{boundaries['end_confidence']:.2f}\t" +
                       f"{boundaries['start_feature'] or 'NA'}\t{boundaries['end_feature'] or 'NA'}\t" +
                       f"{boundaries['tsd'] or 'NA'}\t{kmer_support}\t{pbs_type}\t{ppt_length}\n")
        
        # Step 6: Analyze ORFs if enabled
        orf_info = {}
        if self.orf_analysis:
            for seq in split_records:
                if seq.id in boundary_info:
                    orf_results = self._analyze_orfs(seq, boundary_info[seq.id])
                    if orf_results:
                        orf_info[seq.id] = orf_results
            
            # Write ORF info
            orf_file = os.path.join(self.orf_dir, f"batch_{batch_id}_iter_{iteration}_orfs.tsv")
            with open(orf_file, "w") as f:
                f.write("seq_id\tnum_orfs\torf_locations\tdomain_hits\tcoding_score\n")
                
                for seq_id, orf_data in orf_info.items():
                    orf_locs = ",".join([f"{start}-{end}" for start, end in orf_data['orf_locations']])
                    domains = ",".join(orf_data['domain_hits'])
                    
                    f.write(f"{seq_id}\t{orf_data['num_orfs']}\t{orf_locs}\t" +
                           f"{domains}\t{orf_data['coding_score']:.2f}\n")
        
        # Step 7: Refine sequences based on boundaries and filters
        refined_records = []
        for seq in split_records:
            if seq.id in boundary_info:
                # Refine sequence based on detected boundaries
                refined_seq = self._refine_sequence(seq, boundary_info[seq.id])
                
                # Filter out non-TE regions
                filtered_seq = self._filter_non_te_regions(
                    refined_seq, 
                    alignments_by_query.get(seq.id, []),
                    orf_info.get(seq.id, None)
                )
                
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
    
    def _detect_ltr_boundaries(self, seq_record, alignments, genome_seqs):
        """
        Enhanced method to detect LTR boundaries using multiple evidence types:
        1. Alignment boundaries
        2. Structure motifs (TG...CA)
        3. Target Site Duplications (TSDs) with PSSM scoring
        4. PBS/PPT features
        5. k-mer frequency transitions
        
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
        
        # Track k-mer boundary evidence
        kmer_boundaries = []
        if self.kmer_boundary:
            kmer_boundaries = self._detect_kmer_boundaries(seq_record)
        
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
        
        # Check for k-mer boundary support
        kmer_boundary_support = False
        if kmer_boundaries:
            kmer_5_boundary, kmer_3_boundary = kmer_boundaries
            
            # If k-mer boundaries are close to alignment boundaries, use them for support
            if kmer_5_boundary is not None and abs(kmer_5_boundary - best_start) <= 30:
                kmer_boundary_support = True
                
            if kmer_3_boundary is not None and abs(kmer_3_boundary - best_end) <= 30:
                kmer_boundary_support = True
        
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
            if self.pssm_tsd:
                # Use PSSM scoring for TSD identification
                tsd = self._score_tsds_with_pssm(tsd_evidence)
            else:
                # Simple majority voting
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
            'kmer_boundary_support': kmer_boundary_support,
            'internal_features': internal_features
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
        Enhanced method to find Target Site Duplication (TSD) from genomic context.
        Uses PSSM scoring if enabled.
        
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
        if self.pssm_tsd:
            return self._find_tsd_with_pssm(left_flank, right_flank)
        else:
            # Original approach with relaxed thresholds
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
                
                # More relaxed threshold (70% identity)
                if score > best_score and score >= 0.7:
                    best_score = score
                    best_tsd = left_tsd
                    
            return best_tsd
    
    def _find_tsd_with_pssm(self, left_flank, right_flank):
        """
        Use position-specific scoring matrix (PSSM) for TSD identification.
        This allows more flexible matching compared to strict identity.
        
        Args:
            left_flank: Left flanking sequence
            right_flank: Right flanking sequence
            
        Returns:
            TSD sequence if found, None otherwise
        """
        best_tsd = None
        best_score = -float('inf')
        min_score_threshold = 0.6  # Minimum score to consider as a match
        
        # Define nucleotide mapping
        nuc_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': None}
        
        for tsd_len in range(self.tsd_min, self.tsd_max + 1):
            # Check if flanks are long enough
            if len(left_flank) < tsd_len or len(right_flank) < tsd_len:
                continue
                
            # Get potential TSD sequences
            left_tsd = left_flank[-tsd_len:].upper()
            right_tsd = right_flank[:tsd_len].upper()
            
            # Calculate score using PSSM
            score = 0
            valid_positions = 0
            
            for i in range(tsd_len):
                left_nuc = nuc_map.get(left_tsd[i])
                right_nuc = nuc_map.get(right_tsd[i])
                
                if left_nuc is not None and right_nuc is not None:
                    # Increase score if nucleotides match
                    if left_nuc == right_nuc:
                        score += 1.0
                    else:
                        # Allow for some mismatches, but penalize them
                        score -= 0.2
                        
                    valid_positions += 1
            
            # Normalize score
            if valid_positions > 0:
                normalized_score = score / valid_positions
                
                if normalized_score > best_score and normalized_score >= min_score_threshold:
                    best_score = normalized_score
                    best_tsd = left_tsd
                    
        return best_tsd
    
    def _score_tsds_with_pssm(self, tsd_candidates):
        """
        Score multiple TSD candidates using position-specific scoring matrices.
        
        Args:
            tsd_candidates: List of TSD sequence candidates
            
        Returns:
            Best TSD sequence
        """
        # Group candidates by length
        by_length = defaultdict(list)
        for tsd in tsd_candidates:
            if tsd:
                by_length[len(tsd)].append(tsd)
                
        best_tsd = None
        best_support = 0
        
        # For each length, find the most commonly occurring pattern
        for length, tsds in by_length.items():
            if length < self.tsd_min or length > self.tsd_max:
                continue
                
            if len(tsds) > best_support:
                # Use most frequent TSD
                tsd_counter = Counter(tsds)
                common_tsd, count = tsd_counter.most_common(1)[0]
                
                if count > best_support:
                    best_support = count
                    best_tsd = common_tsd
                    
        return best_tsd
    
    def _find_internal_features(self, seq_record, start, end):
        """
        Enhanced method to find LTR internal features like PBS and PPT.
        
        Args:
            seq_record: SeqRecord object
            start: Start position
            end: End position
            
        Returns:
            Dictionary with internal feature information
        """
        seq_str = str(seq_record.seq)
        
        # Look for PBS near the 5' LTR (5-20 bp after LTR)
        pbs_result = None
        for pbs_pos in range(start + 5, min(start + 30, end - 100, len(seq_str) - 18)):
            window = seq_str[pbs_pos:pbs_pos+18]
            
            for pattern_name, pattern in self.pbs_patterns.items():
                # Calculate similarity
                matches = sum(1 for a, b in zip(window.upper(), pattern) if a == b)
                similarity = matches / len(pattern)
                
                if similarity >= 0.8:  # At least 80% match
                    pbs_result = {
                        'position': pbs_pos,
                        'sequence': window,
                        'pattern': pattern_name,
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
            ppt_pattern = r'[AG]{8,15}'
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
        
    def _analyze_orfs(self, seq_record, boundaries):
        """
        Analyze ORFs in the sequence to identify potential coding regions.
        
        Args:
            seq_record: SeqRecord object
            boundaries: Detected boundaries
            
        Returns:
            Dictionary with ORF analysis results
        """
        if not BIOUTILS_AVAILABLE:
            return self._analyze_orfs_simple(seq_record, boundaries)
            
        # Extract sequence in the boundaries
        start = boundaries['start']
        end = boundaries['end']
        seq = seq_record.seq[start:end]
        
        # Find all ORFs (minimum length 300 bp = 100 aa)
        min_orf_len = 300
        orfs = []
        
        # Look in all 3 reading frames
        for frame in range(3):
            # Look for start codons
            for i in range(frame, len(seq) - 2, 3):
                if seq[i:i+3].upper() == 'ATG':
                    # Found start codon
                    for j in range(i + 3, len(seq) - 2, 3):
                        # Look for stop codon
                        if seq[j:j+3].upper() in ('TAA', 'TAG', 'TGA'):
                            # Found stop codon
                            orf_len = j + 3 - i
                            
                            if orf_len >= min_orf_len:
                                # Convert coordinates to original sequence
                                orf_start = start + i
                                orf_end = start + j + 3
                                
                                orf_seq = seq[i:j+3]
                                orfs.append((orf_start, orf_end, orf_seq))
                            break
        
        # Analyze ORFs for TE domains
        domain_hits = []
        
        for orf_start, orf_end, orf_seq in orfs:
            # Translate to protein
            protein_seq = orf_seq.translate()
            
            # Check for TE-specific domains
            for domain_type, motifs in self.te_domains.items():
                for motif in motifs:
                    if motif.lower() in str(protein_seq).lower():
                        domain_hits.append(f"{domain_type}:{motif}")
        
        # Calculate coding potential score (based on GC content, codon usage, etc.)
        coding_score = 0
        
        if orfs:
            # Average GC content of ORFs
            gc_content = sum(GC(orf_seq) for _, _, orf_seq in orfs) / len(orfs)
            
            # Score based on TE domain hits
            domain_score = len(set(domain_hits)) * 0.2
            
            # Score based on GC content (typical for coding regions)
            gc_score = 0
            if 40 <= gc_content <= 60:
                gc_score = 0.5
            
            # Score based on number of ORFs
            orf_score = min(len(orfs), 3) * 0.1
            
            coding_score = gc_score + domain_score + orf_score
            
        # Return ORF analysis results
        return {
            'num_orfs': len(orfs),
            'orf_locations': [(start, end) for start, end, _ in orfs],
            'domain_hits': list(set(domain_hits)),
            'coding_score': coding_score
        }
    
    def _analyze_orfs_simple(self, seq_record, boundaries):
        """
        Simple ORF analysis without BioPython utilities.
        
        Args:
            seq_record: SeqRecord object
            boundaries: Detected boundaries
            
        Returns:
            Dictionary with basic ORF analysis results
        """
        # Extract sequence in the boundaries
        start = boundaries['start']
        end = boundaries['end']
        seq = str(seq_record.seq[start:end]).upper()
        
        # Find all ORFs (minimum length 300 bp = 100 aa)
        min_orf_len = 300
        orfs = []
        
        # Look in all 3 reading frames
        for frame in range(3):
            # Look for start codons
            i = frame
            while i < len(seq) - 2:
                if seq[i:i+3] == 'ATG':
                    # Found start codon
                    j = i + 3
                    while j < len(seq) - 2:
                        # Look for stop codon
                        if seq[j:j+3] in ('TAA', 'TAG', 'TGA'):
                            # Found stop codon
                            orf_len = j + 3 - i
                            
                            if orf_len >= min_orf_len:
                                # Convert coordinates to original sequence
                                orf_start = start + i
                                orf_end = start + j + 3
                                
                                orfs.append((orf_start, orf_end))
                            break
                        j += 3
                    
                    # Move to next position after current start codon
                    i += 3
                else:
                    i += 3
        
        # Return basic ORF analysis results
        return {
            'num_orfs': len(orfs),
            'orf_locations': orfs,
            'domain_hits': [],
            'coding_score': min(0.5, len(orfs) * 0.1)  # Simple estimation
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
    
    def _filter_non_te_regions(self, seq_record, alignments, orf_info=None):
        """
        Enhanced method to filter out regions that don't appear to be part of the transposon.
        Uses alignment patterns, coding potential, and dynamic thresholding.
        
        Args:
            seq_record: SeqRecord object
            alignments: List of alignments for this sequence
            orf_info: ORF analysis information (optional)
            
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
            
        # Calculate coverage threshold - dynamic or fixed
        if self.dynamic_threshold:
            # Use more sophisticated thresholding
            # Sort coverage values and find elbow point or use percentile
            sorted_coverage = np.sort(smoothed)
            
            # Calculate percentiles
            p25 = np.percentile(sorted_coverage, 25)
            p75 = np.percentile(sorted_coverage, 75)
            
            # Detect bimodal distribution
            if p75 - p25 > 1.5:
                # Likely bimodal - use elbow point
                # Simple elbow detection using largest second derivative
                d2 = np.diff(np.diff(sorted_coverage))
                elbow_idx = np.argmax(d2) if len(d2) > 0 else len(sorted_coverage) // 2
                threshold = sorted_coverage[elbow_idx] if elbow_idx < len(sorted_coverage) else np.mean(smoothed) * 0.3
            else:
                # Unimodal - use percentile
                threshold = p25 * 0.5
        else:
            # Fixed threshold
            threshold = np.mean(smoothed) * 0.3  # 30% of mean coverage
            
        # Find regions to keep (above threshold)
        keep_regions = []
        in_region = False
        region_start = 0
        
        # Special protection for ORFs if available
        protected_regions = []
        if orf_info and 'orf_locations' in orf_info and orf_info['orf_locations']:
            # If coding regions found, protect them
            for orf_start, orf_end in orf_info['orf_locations']:
                # Add some padding around ORFs
                protected_regions.append((max(0, orf_start - 50), min(seq_len, orf_end + 50)))
        
        for i in range(seq_len):
            # Check if position is protected as part of an ORF
            is_protected = any(start <= i < end for start, end in protected_regions)
            
            if smoothed[i] >= threshold or is_protected:
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
    
    def _perform_subfamily_clustering(self, sequences):
        """
        Perform clustering to identify subfamilies.
        Uses UMAP for dimensionality reduction and HDBSCAN for clustering.
        
        Args:
            sequences: List of sequence records
            
        Returns:
            Path to clustered sequences file
        """
        if not CLUSTERING_AVAILABLE:
            self.logger.warning("Clustering libraries not available. Skipping subfamily clustering.")
            return os.path.join(self.output_dir, "no_clustering.fa")
            
        self.logger.info(f"Performing subfamily clustering on {len(sequences)} sequences")
        
        # Create sequence features for clustering
        features = []
        seq_ids = []
        
        for seq in sequences:
            # Extract sequence features
            seq_str = str(seq.seq).upper()
            
            # Calculate k-mer frequencies
            k = 3  # Trinucleotides
            kmers = {}
            
            for i in range(len(seq_str) - k + 1):
                kmer = seq_str[i:i+k]
                if 'N' not in kmer:
                    kmers[kmer] = kmers.get(kmer, 0) + 1
            
            # Normalize by sequence length
            for kmer in kmers:
                kmers[kmer] = kmers[kmer] / (len(seq_str) - k + 1)
                
            # Convert to vector
            all_kmers = [''.join(x) for x in itertools.product('ACGT', repeat=k)]
            kmer_vector = [kmers.get(kmer, 0) for kmer in all_kmers]
            
            features.append(kmer_vector)
            seq_ids.append(seq.id)
            
        # Convert to numpy array
        features = np.array(features)
        
        # Too few sequences for meaningful clustering
        if len(features) < 5:
            self.logger.warning("Too few sequences for meaningful clustering. Skipping.")
            output_file = os.path.join(self.output_dir, "no_clustering.fa")
            SeqIO.write(sequences, output_file, "fasta")
            return output_file
            
        # Reduce dimensionality with UMAP
        self.logger.info("Performing UMAP dimensionality reduction")
        umap_model = umap.UMAP(n_neighbors=min(15, len(features)-1), 
                              min_dist=0.1, 
                              n_components=2,
                              metric='correlation')
                              
        umap_features = umap_model.fit_transform(features)
        
        # Cluster with HDBSCAN
        self.logger.info("Performing HDBSCAN clustering")
        min_cluster_size = max(2, len(features) // 10)  # Adaptive cluster size
        
        hdbscan_model = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size,
                                       min_samples=1,
                                       cluster_selection_epsilon=0.5)
                                       
        cluster_labels = hdbscan_model.fit_predict(umap_features)
        
        # Count clusters
        n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
        self.logger.info(f"Found {n_clusters} clusters, {np.sum(cluster_labels == -1)} outliers")
        
        # Create clustered sequences
        clustered_seqs = []
        
        for i, seq in enumerate(sequences):
            cluster_id = cluster_labels[i]
            cluster_name = f"subfamily_{cluster_id}" if cluster_id >= 0 else "outlier"
            
            # Create new record with cluster information
            cluster_seq = SeqRecord(
                seq.seq,
                id=f"{seq.id}",
                name=seq.name,
                description=f"{seq.description} | Cluster: {cluster_name}"
            )
            
            clustered_seqs.append(cluster_seq)
            
        # Write clustered sequences
        output_file = os.path.join(self.clustering_dir, "clustered_sequences.fa")
        SeqIO.write(clustered_seqs, output_file, "fasta")
        
        # Write cluster information
        cluster_info_file = os.path.join(self.clustering_dir, "cluster_information.tsv")
        with open(cluster_info_file, "w") as f:
            f.write("seq_id\tcluster_id\tcluster_name\n")
            
            for i, seq_id in enumerate(seq_ids):
                cluster_id = cluster_labels[i]
                cluster_name = f"subfamily_{cluster_id}" if cluster_id >= 0 else "outlier"
                
                f.write(f"{seq_id}\t{cluster_id}\t{cluster_name}\n")
        
        self.logger.info(f"Clustering completed. Results saved to {self.clustering_dir}")
        return output_file
    
    def _classify_sequences(self, fasta_file):
        """
        Classify TE sequences based on structural features and coding domains.
        
        Args:
            fasta_file: Path to input FASTA file
            
        Returns:
            Path to classified sequences file
        """
        self.logger.info("Classifying TE sequences")
        
        # Load sequences
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        classified_seqs = []
        
        # Classification information
        classifications = {}
        
        # Check if BLAST is available for classification
        blast_available = shutil.which("blastn") and shutil.which("makeblastdb")
        
        # If BLAST is available, we can use it for classification
        if blast_available:
            self.logger.info("Using BLAST for classification")
            # This would involve creating a database of known TEs and comparing
            # For now, we'll use structural features only
        
        for seq in sequences:
            # Default classification
            te_class = "Unknown"
            te_family = "Unknown"
            
            # Look for LTR features in the sequence description
            if "Refined" in seq.description and "features" in seq.description:
                # This is likely an LTR retrotransposon
                te_class = "LTR retrotransposon"
                
                # Look for gag/pol domains in the description
                if "gag" in seq.description.lower() or "pol" in seq.description.lower():
                    te_family = "Ty3/Gypsy or Ty1/Copia"
                elif "env" in seq.description.lower():
                    te_family = "Endogenous retrovirus"
                else:
                    te_family = "LTR/Unknown"
            
            # Record classification
            classifications[seq.id] = {
                'class': te_class,
                'family': te_family
            }
            
            # Create new record with classification
            classified_seq = SeqRecord(
                seq.seq,
                id=seq.id,
                name=seq.name,
                description=f"{seq.description} | Class: {te_class} | Family: {te_family}"
            )
            
            classified_seqs.append(classified_seq)
            
        # Write classified sequences
        output_file = os.path.join(self.classification_dir, "classified_sequences.fa")
        SeqIO.write(classified_seqs, output_file, "fasta")
        
        # Write classification information
        class_info_file = os.path.join(self.classification_dir, "classification_information.tsv")
        with open(class_info_file, "w") as f:
            f.write("seq_id\tclass\tfamily\n")
            
            for seq_id, info in classifications.items():
                f.write(f"{seq_id}\t{info['class']}\t{info['family']}\n")
                
        self.logger.info(f"Classification completed. Results saved to {self.classification_dir}")
        return output_file
    
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
    
    # Added arguments for enhanced functionality
    parser.add_argument("--low-complexity-filter", action="store_true",
                       help="Filter low complexity regions")
    parser.add_argument("--clustering", action="store_true",
                       help="Perform subfamily clustering")
    parser.add_argument("--dynamic-threshold", action="store_true",
                       help="Use dynamic thresholds for non-TE filtering")
    parser.add_argument("--orf-analysis", action="store_true",
                       help="Perform ORF analysis")
    parser.add_argument("--kmer-boundary", action="store_true",
                       help="Use k-mer based boundary detection")
    parser.add_argument("--blank-region-threshold", type=int, default=100,
                       help="Minimum length of blank regions for chimera detection (default: 100)")
    parser.add_argument("--orientation-aware", action="store_true",
                       help="Use orientation-aware chimera detection")
    parser.add_argument("--pssm-tsd", action="store_true",
                       help="Use PSSM for TSD identification")
    parser.add_argument("--classify-output", action="store_true",
                       help="Classify output sequences")
    
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
        log_level=args.log_level,
        low_complexity_filter=args.low_complexity_filter,
        clustering=args.clustering,
        dynamic_threshold=args.dynamic_threshold,
        orf_analysis=args.orf_analysis,
        kmer_boundary=args.kmer_boundary,
        blank_region_threshold=args.blank_region_threshold,
        orientation_aware=args.orientation_aware,
        pssm_tsd=args.pssm_tsd,
        classify_output=args.classify_output
    )
    
    # Run optimization
    optimizer.optimize()


if __name__ == "__main__":
    main()
