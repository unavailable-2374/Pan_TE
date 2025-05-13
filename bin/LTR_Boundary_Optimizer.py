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
import glob
import itertools
import math
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
except (ImportError, TypeError):
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
                 orientation_aware=True, advanced_tsd=True,
                 custom_motifs=None, weighted_evidence=True):
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
            advanced_tsd: Whether to use advanced TSD detection instead of PSSM
            custom_motifs: Custom motifs for boundary detection
            weighted_evidence: Whether to use weighted evidence integration
        """
        # Add tracking for false positive detection
        self.chimera_info = {}  # Store chimera detection results for later use
        self.fp_score_threshold = 0.4  # Threshold for false positive classification

        # Create false positive confidence categories
        self.confidence_categories = {
            0.8: "High confidence LTR",
            0.6: "Moderate confidence LTR", 
            0.4: "Low confidence LTR",
            0.0: "Potential false positive"
        }

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
        self.advanced_tsd = advanced_tsd  # Replace PSSM with advanced methods
        self.custom_motifs = custom_motifs
        self.weighted_evidence = weighted_evidence
        
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
        
        # Define evidence weights for boundary detection
        self.evidence_weights = {
            'terminal_motifs': 0.30,    # TG...CA motifs
            'alignment_boundaries': 0.25,  # Based on alignment data
            'kmer_transitions': 0.20,    # K-mer composition changes
            'tsd_evidence': 0.15,        # Target site duplications
            'internal_features': 0.10    # PBS/PPT evidence
        }
        
        # Store previously found TSDs for adaptive learning
        self.tsd_frequency = Counter()  # Tracks TSD patterns by frequency
        
        # Initialize TSD nucleotide bias from empirical data
        # Based on common observed patterns in LTR transposons
        self.tsd_nucleotide_bias = {
            'A': 0.30,  # Slight preference for A and T in TSDs
            'C': 0.20,
            'G': 0.20,
            'T': 0.30
        }
    
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
        self.false_positive_dir = os.path.join(self.output_dir, "false_positive_analysis")
        
        for directory in [self.alignment_dir, self.chimera_dir, self.boundary_dir, 
                        self.stats_dir, self.clustering_dir, self.orf_dir, 
                        self.false_positive_dir]:
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
        optional_tools = ["seqkit", "blastn", "dustmasker", "trf"]
        
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
    
    def _evaluate_false_positive_likelihood(self, seq_record, boundary_info, alignments_by_query, orf_info=None):
        """
        Evaluate the likelihood that a sequence is a false positive LTR element.
        
        Args:
            seq_record: SeqRecord object
            boundary_info: Dictionary with boundary information
            alignments_by_query: Dictionary mapping query IDs to their alignments
            orf_info: ORF analysis information (optional)
            
        Returns:
            Dictionary with scoring information and final classification
        """
        seq_id = seq_record.id
        seq_str = str(seq_record.seq)
        seq_len = len(seq_str)
        alignments = alignments_by_query.get(seq_id, [])
        
        # Initialize scores for different categories
        # Each category will have a score between 0-1, where higher scores mean more likely to be a true LTR
        scores = {
            'structural_features': 0.0,  # Terminal motifs, TSD, etc.
            'alignment_quality': 0.0,    # Coverage, identity, copy number
            'coding_potential': 0.0,     # Presence of ORFs and TE domains
            'sequence_composition': 0.0,  # GC content, complexity
            'homogeneity': 0.0           # Non-chimeric nature, clustering
        }
        
        # 1. Evaluate structural features
        structural_score = 0.0
        feature_weights = {
            'terminal_motifs': 0.4,  # TG...CA
            'tsd': 0.3,              # Target Site Duplications
            'pbs_ppt': 0.3           # Primer Binding Site and Polypurine Tract
        }
        
        # Terminal motifs check
        motif_score = 0.0
        if boundary_info and 'start_feature' in boundary_info and boundary_info['start_feature']:
            motif_score += 0.5
        if boundary_info and 'end_feature' in boundary_info and boundary_info['end_feature']:
            motif_score += 0.5
        
        # TSD check
        tsd_score = 0.0
        if boundary_info and 'tsd' in boundary_info and boundary_info['tsd']:
            tsd_score = 1.0
        
        # PBS/PPT check
        internal_score = 0.0
        if boundary_info and 'internal_features' in boundary_info:
            if boundary_info['internal_features'].get('pbs'):
                internal_score += 0.5
            if boundary_info['internal_features'].get('ppt'):
                internal_score += 0.5
        
        # Calculate weighted structural feature score
        structural_score = (
            feature_weights['terminal_motifs'] * motif_score + 
            feature_weights['tsd'] * tsd_score + 
            feature_weights['pbs_ppt'] * internal_score
        )
        scores['structural_features'] = structural_score
        
        # 2. Evaluate alignment quality
        alignment_score = 0.0
        if alignments:
            # Calculate average identity and coverage
            avg_identity = sum(aln['identity'] for aln in alignments) / len(alignments)
            avg_coverage = sum(aln['coverage'] for aln in alignments) / len(alignments)
            
            # Normalize to 0-1 scale (80% and above is considered good)
            identity_score = min(1.0, avg_identity / 80.0)
            coverage_score = min(1.0, avg_coverage / 80.0)
            
            # Consider number of alignments (more is better, up to a point)
            copy_number = len(alignments)
            copy_score = min(1.0, copy_number / 10.0)  # 10+ copies gets full score
            
            # Weight the factors
            alignment_score = 0.4 * identity_score + 0.4 * coverage_score + 0.2 * copy_score
        
        scores['alignment_quality'] = alignment_score
        
        # 3. Evaluate coding potential
        coding_score = 0.0
        if orf_info:
            # Consider number of ORFs
            num_orfs = orf_info.get('num_orfs', 0)
            
            # Score based on ORF presence (solo LTRs won't have this, so we limit the weight)
            orf_presence_score = min(1.0, num_orfs / 2.0)  # 2+ ORFs gets full score
            
            # Score based on domain hits
            domain_hits = orf_info.get('domain_hits', [])
            domain_score = min(1.0, len(domain_hits) / 3.0)  # 3+ domain hits gets full score
            
            # Use the provided coding score if available
            provided_score = orf_info.get('coding_score', 0.0)
            
            # Weight the factors 
            coding_score = 0.3 * orf_presence_score + 0.5 * domain_score + 0.2 * provided_score
        
        scores['coding_potential'] = coding_score
        
        # 4. Evaluate sequence composition
        composition_score = 0.0
        
        # GC content check - most LTRs have GC content between 30-70%
        if 'BIOUTILS_AVAILABLE' in globals() and BIOUTILS_AVAILABLE:
            from Bio.SeqUtils import GC
            gc_content = GC(seq_str)
            
            # Score higher for GC content in the expected range
            if 30 <= gc_content <= 70:
                gc_score = 1.0
            else:
                # Lower score for extreme GC content
                gc_score = 1.0 - min(1.0, abs(gc_content - 50) / 20.0)
        else:
            # Simple GC calculation if BioPython utils not available
            gc_count = seq_str.count('G') + seq_str.count('C')
            gc_content = gc_count / seq_len if seq_len > 0 else 0
            
            if 0.3 <= gc_content <= 0.7:
                gc_score = 1.0
            else:
                gc_score = 1.0 - min(1.0, abs(gc_content - 0.5) / 0.2)
        
        # Low complexity check - count proportion of sequence that is not N or simple repeats
        low_complexity_proportion = 0.0
        n_count = seq_str.upper().count('N')
        
        # Simple repeat detection (e.g., ATATATATAT)
        simple_repeat_count = 0
        for repeat_len in range(1, 6):  # Check for repeats of length 1-5
            for i in range(len(seq_str) - repeat_len * 5):
                substr = seq_str[i:i+repeat_len]
                if substr * 5 == seq_str[i:i+repeat_len*5]:
                    simple_repeat_count += repeat_len * 5
                    i += repeat_len * 5  # Skip ahead
        
        low_complexity_proportion = (n_count + simple_repeat_count) / seq_len if seq_len > 0 else 0
        complexity_score = 1.0 - min(1.0, low_complexity_proportion / 0.3)  # Penalize if >30% low complexity
        
        # Length check - very short or very long sequences might be suspicious
        if self.min_ltr_len <= seq_len <= self.max_ltr_len:
            length_score = 1.0
        else:
            # Decrease score for lengths outside the expected range
            deviation = min(
                abs(seq_len - self.min_ltr_len) / self.min_ltr_len,
                abs(seq_len - self.max_ltr_len) / self.max_ltr_len
            )
            length_score = max(0.0, 1.0 - deviation)
        
        # Calculate weighted composition score
        composition_score = 0.4 * gc_score + 0.4 * complexity_score + 0.2 * length_score
        scores['sequence_composition'] = composition_score
        
        # 5. Evaluate homogeneity (non-chimeric nature)
        homogeneity_score = 0.0
        
        # Check if sequence was flagged as chimeric
        is_chimeric = False
        for seq_id_check, info in self.chimera_info.items():
            if seq_record.id.startswith(seq_id_check.split('_segment_')[0]):
                is_chimeric = info.get('is_chimeric', False)
                break
        
        # Higher score for non-chimeric sequences
        if not is_chimeric:
            homogeneity_score += 0.7
        else:
            # For chimeric sequences, check how many segments
            for seq_id_check, info in self.chimera_info.items():
                if seq_record.id.startswith(seq_id_check.split('_segment_')[0]):
                    # Less segments is better (less fragmented)
                    num_segments = info.get('num_segments', 1)
                    if num_segments <= 2:
                        homogeneity_score += 0.3  # Some LTRs might have 2 segments
                    else:
                        # Highly fragmented is suspicious
                        homogeneity_score += max(0.0, 0.3 - (num_segments - 2) * 0.1)
        
        # Check for blank regions
        blank_regions = []
        for seq_id_check, info in self.chimera_info.items():
            if seq_record.id.startswith(seq_id_check.split('_segment_')[0]):
                blank_regions = info.get('blank_regions', [])
                break
        
        # Penalize if there are many blank regions
        blank_penalty = min(0.3, len(blank_regions) * 0.1)
        homogeneity_score = max(0.0, homogeneity_score - blank_penalty)
        
        scores['homogeneity'] = homogeneity_score
        
        # Calculate final score with category weights
        category_weights = {
            'structural_features': 0.35,
            'alignment_quality': 0.25,
            'coding_potential': 0.15,
            'sequence_composition': 0.15,
            'homogeneity': 0.10
        }
        
        final_score = sum(scores[category] * category_weights[category] for category in scores)
        
        # Classify based on final score
        classification = "Unknown"
        for threshold, category in sorted(self.confidence_categories.items(), reverse=True):
            if final_score >= threshold:
                classification = category
                break
        
        # Store detailed evidence
        evidence = {
            'terminal_motifs': {
                'has_5prime_motif': bool(boundary_info and boundary_info.get('start_feature')),
                'has_3prime_motif': bool(boundary_info and boundary_info.get('end_feature')),
                'motif_5prime': boundary_info.get('start_feature', 'NA') if boundary_info else 'NA',
                'motif_3prime': boundary_info.get('end_feature', 'NA') if boundary_info else 'NA'
            },
            'tsd': {
                'present': bool(boundary_info and boundary_info.get('tsd')),
                'sequence': boundary_info.get('tsd', 'NA') if boundary_info else 'NA'
            },
            'internal_features': {
                'has_pbs': bool(boundary_info and boundary_info.get('internal_features', {}).get('pbs')),
                'has_ppt': bool(boundary_info and boundary_info.get('internal_features', {}).get('ppt')),
                'pbs_type': boundary_info.get('internal_features', {}).get('pbs', {}).get('pattern', 'NA') if boundary_info else 'NA',
                'ppt_length': boundary_info.get('internal_features', {}).get('ppt', {}).get('length', 'NA') if boundary_info else 'NA'
            },
            'alignment': {
                'num_alignments': len(alignments),
                'avg_identity': sum(aln['identity'] for aln in alignments) / len(alignments) if alignments else 0,
                'avg_coverage': sum(aln['coverage'] for aln in alignments) / len(alignments) if alignments else 0
            },
            'coding': {
                'num_orfs': orf_info.get('num_orfs', 0) if orf_info else 0,
                'domain_hits': orf_info.get('domain_hits', []) if orf_info else [],
                'coding_score': orf_info.get('coding_score', 0.0) if orf_info else 0.0
            },
            'composition': {
                'gc_content': gc_content,
                'length': seq_len,
                'low_complexity_proportion': low_complexity_proportion
            },
            'chimeric': {
                'is_chimeric': is_chimeric,
                'num_segments': 0,
                'num_blank_regions': len(blank_regions)
            }
        }
        
        # Return the evaluation results
        return {
            'seq_id': seq_id,
            'category_scores': scores,
            'final_score': final_score,
            'classification': classification,
            'evidence': evidence
        }

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
            all_chimera_info = {}  # Collect all chimera info here

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
                        batch_results, batch_chimera_info = future.result()
                        refined_sequences.extend(batch_results)
                        all_chimera_info.update(batch_chimera_info)

                        self.logger.info(f"Completed batch {batch_id} with {len(batch_results)} sequences")
                    except Exception as e:
                        self.logger.error(f"Error processing batch {batch_id}: {str(e)}")
            
            # Write sequences from this iteration
            self.chimera_info = all_chimera_info
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
            cmd = ["dustmasker", "-in", fasta_file, "-out", output_file, "-outfmt","fasta"]
            
            try:
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
        cmd = ["seqkit", "faidx", self.genome_file]
        
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
        
        # Store chimera info at the class level for later use
        for seq_id, info in chimera_info.items():
            self.chimera_info[seq_id] = info

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
                    
        # Step 6b: Evaluate false positive likelihood
        false_positive_evaluations = {}
        for seq in split_records:
            if seq.id in boundary_info:
                evaluation = self._evaluate_false_positive_likelihood(
                    seq, 
                    boundary_info[seq.id], 
                    alignments_by_query, 
                    orf_info.get(seq.id, None)
                )
                false_positive_evaluations[seq.id] = evaluation
            else:
                # For sequences without detected boundaries, assign low confidence
                evaluation = {
                    'seq_id': seq.id,
                    'category_scores': {k: 0.0 for k in ['structural_features', 'alignment_quality', 
                                                        'coding_potential', 'sequence_composition', 'homogeneity']},
                    'final_score': 0.2,  # Very low score
                    'classification': "Potential false positive",
                    'evidence': {
                        'terminal_motifs': {'has_5prime_motif': False, 'has_3prime_motif': False, 
                                            'motif_5prime': 'NA', 'motif_3prime': 'NA'},
                        'tsd': {'present': False, 'sequence': 'NA'},
                        'internal_features': {'has_pbs': False, 'has_ppt': False, 'pbs_type': 'NA', 'ppt_length': 'NA'},
                        'alignment': {'num_alignments': len(alignments_by_query.get(seq.id, [])), 
                                    'avg_identity': 0, 'avg_coverage': 0},
                        'coding': {'num_orfs': 0, 'domain_hits': [], 'coding_score': 0.0},
                        'composition': {'gc_content': 0, 'length': len(seq.seq), 'low_complexity_proportion': 0},
                        'chimeric': {'is_chimeric': False, 'num_segments': 0, 'num_blank_regions': 0}
                    }
                }
                false_positive_evaluations[seq.id] = evaluation
        
        # Write false positive evaluation results
        fp_file = os.path.join(self.false_positive_dir, f"batch_{batch_id}_iter_{iteration}_fp_evaluation.tsv")
        with open(fp_file, "w") as f:
            f.write("seq_id\tfinal_score\tclassification\tstructural_score\talignment_score\t" +
                "coding_score\tcomposition_score\thomogeneity_score\tevidence_summary\n")
            
            for seq_id, eval_data in false_positive_evaluations.items():
                scores = eval_data['category_scores']
                
                # Create a brief evidence summary
                evidence = eval_data['evidence']
                summary_parts = []
                
                if evidence['terminal_motifs']['has_5prime_motif'] and evidence['terminal_motifs']['has_3prime_motif']:
                    summary_parts.append("Has terminal motifs")
                if evidence['tsd']['present']:
                    summary_parts.append(f"Has TSD: {evidence['tsd']['sequence']}")
                if evidence['internal_features']['has_pbs']:
                    summary_parts.append(f"Has PBS: {evidence['internal_features']['pbs_type']}")
                if evidence['internal_features']['has_ppt']:
                    summary_parts.append("Has PPT")
                if evidence['coding']['num_orfs'] > 0:
                    summary_parts.append(f"Has {evidence['coding']['num_orfs']} ORFs")
                if evidence['coding']['domain_hits']:
                    domain_info = ','.join(evidence['coding']['domain_hits'])
                    if len(domain_info) > 30:
                        domain_info = domain_info[:27] + "..."
                    summary_parts.append(f"Has domains: {domain_info}")
                if evidence['chimeric']['is_chimeric']:
                    summary_parts.append("Is chimeric")
                
                evidence_summary = "; ".join(summary_parts) if summary_parts else "No supporting evidence"
                
                f.write(f"{seq_id}\t{eval_data['final_score']:.3f}\t{eval_data['classification']}\t" +
                    f"{scores['structural_features']:.3f}\t{scores['alignment_quality']:.3f}\t" +
                    f"{scores['coding_potential']:.3f}\t{scores['sequence_composition']:.3f}\t" +
                    f"{scores['homogeneity']:.3f}\t{evidence_summary}\n")
        
        # Step 7: Refine sequences based on boundaries and filters
        refined_records = []
        for seq in split_records:
            if seq.id in boundary_info:
                # Refine sequence based on detected boundaries
                refined_seq = self._refine_sequence(seq, boundary_info[seq.id])
                
                # Add classification to description if available
                if seq.id in false_positive_evaluations:
                    eval_data = false_positive_evaluations[seq.id]
                    refined_seq.description += f" | Confidence: {eval_data['classification']} ({eval_data['final_score']:.2f})"
                
                # Filter out non-TE regions
                filtered_seq = self._filter_non_te_regions(
                    refined_seq, 
                    alignments_by_query.get(seq.id, []),
                    orf_info.get(seq.id, None)
                )
                
                refined_records.append(filtered_seq)
            else:
                # Keep original if no boundaries detected
                # Add low confidence note
                seq.description += " | Confidence: Potential false positive (no boundaries detected)"
                refined_records.append(seq)
        
        return refined_records, chimera_info
    
    # Add a method to create an HTML report
    def _create_fp_html_report(self, evaluations):
        """
        Create an HTML report for false positive analysis.
        
        Args:
            evaluations: List of evaluation dictionaries with 'seq_id', 'final_score', and 'classification' keys
        """
        html_file = os.path.join(self.false_positive_dir, "false_positive_report.html")
        
        # Sort evaluations by score (ascending)
        sorted_evals = sorted(evaluations, key=lambda x: x['final_score'])
        
        with open(html_file, "w") as f:
            f.write("""
            <!DOCTYPE html>
            <html>
            <head>
                <title>LTR Consensus False Positive Analysis</title>
                <style>
                    body { font-family: Arial, sans-serif; margin: 20px; }
                    h1, h2 { color: #333; }
                    .container { max-width: 1200px; margin: 0 auto; }
                    table { border-collapse: collapse; width: 100%; margin-top: 20px; }
                    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
                    th { background-color: #f2f2f2; }
                    tr:nth-child(even) { background-color: #f9f9f9; }
                    .high { background-color: #d4edda; }
                    .moderate { background-color: #fff3cd; }
                    .low { background-color: #f8d7da; }
                    .false-positive { background-color: #f5c6cb; }
                    .summary { margin: 20px 0; padding: 15px; background-color: #f8f9fa; border-radius: 5px; }
                </style>
            </head>
            <body>
                <div class="container">
                    <h1>LTR Consensus False Positive Analysis</h1>
            """)
            
            # Add summary
            class_counts = {}
            for eval_data in evaluations:
                classification = eval_data['classification']
                class_counts[classification] = class_counts.get(classification, 0) + 1
            
            f.write('<div class="summary">')
            f.write(f'<h2>Summary</h2>')
            f.write(f'<p>Total sequences evaluated: {len(evaluations)}</p>')
            f.write('<ul>')
            
            for classification, count in sorted(class_counts.items(), key=lambda x: x[1], reverse=True):
                percentage = (count / len(evaluations)) * 100 if evaluations else 0
                f.write(f'<li>{classification}: {count} ({percentage:.1f}%)</li>')
            
            f.write('</ul>')
            f.write('</div>')
            
            # Add sequence table
            f.write('<h2>Sequence Confidence Scores</h2>')
            f.write('<p>Sorted by confidence score (lowest to highest)</p>')
            f.write('<table>')
            f.write('<tr><th>Sequence ID</th><th>Classification</th><th>Confidence Score</th></tr>')
            
            for eval_data in sorted_evals:
                seq_id = eval_data['seq_id']
                classification = eval_data['classification']
                score = eval_data['final_score']
                
                # Add CSS class based on classification
                css_class = ""
                if classification == "High confidence LTR":
                    css_class = "high"
                elif classification == "Moderate confidence LTR":
                    css_class = "moderate"
                elif classification == "Low confidence LTR":
                    css_class = "low"
                elif classification == "Potential false positive":
                    css_class = "false-positive"
                
                f.write(f'<tr class="{css_class}"><td>{seq_id}</td><td>{classification}</td><td>{score:.3f}</td></tr>')
            
            f.write('</table>')
            
            # Close HTML
            f.write("""
                </div>
            </body>
            </html>
            """)
        
        self.logger.info(f"Created HTML report: {html_file}")

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
        Enhanced method to detect LTR boundaries using multiple evidence types
        with an integrative approach that doesn't rely on PSSM.
        
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
            
        # Validate boundaries
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

    def _find_tsd_evidence(self, alignments, genome_seqs):
        """
        Enhanced TSD detection using local sequence context and alignment patterns
        instead of PSSM.
        
        Args:
            alignments: List of alignments for this sequence
            genome_seqs: Genome access object
            
        Returns:
            List of TSD evidence entries
        """
        tsd_evidence = []
        
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
            tsd = best_tsd_evidence['sequence']
            
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
            tsd = best_tsd['sequence']
        
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
            
            # Find original sequence if available
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
    
        # Add false positive analysis summary
        false_positive_files = glob.glob(os.path.join(self.false_positive_dir, '*_fp_evaluation.tsv'))
        
        if false_positive_files:
            # Collect all evaluations
            all_evaluations = []
            
            for fp_file in false_positive_files:
                with open(fp_file, 'r') as f:
                    # Skip header
                    next(f)
                    
                    for line in f:
                        fields = line.strip().split('\t')
                        if len(fields) >= 3:
                            seq_id = fields[0]
                            score = float(fields[1])
                            classification = fields[2]
                            
                            all_evaluations.append({
                                'seq_id': seq_id,
                                'final_score': score,
                                'classification': classification
                            })
            
            # Count by classification
            class_counts = {}
            for eval_data in all_evaluations:
                classification = eval_data['classification']
                class_counts[classification] = class_counts.get(classification, 0) + 1
            
            # Add to summary
            with open(summary_file, "a") as f:
                f.write("\n--- False Positive Analysis ---\n")
                f.write(f"Total sequences evaluated: {len(all_evaluations)}\n")
                
                for classification, count in sorted(class_counts.items(), key=lambda x: x[1], reverse=True):
                    percentage = (count / len(all_evaluations)) * 100 if all_evaluations else 0
                    f.write(f"{classification}: {count} ({percentage:.1f}%)\n")
                
                # Calculate average confidence score
                avg_score = sum(eval_data['final_score'] for eval_data in all_evaluations) / len(all_evaluations) if all_evaluations else 0
                f.write(f"Average confidence score: {avg_score:.3f}\n")
                
                # List potential false positives
                false_positives = [eval_data['seq_id'] for eval_data in all_evaluations 
                                if eval_data['classification'] == "Potential false positive"]
                
                if false_positives:
                    f.write(f"\nPotential false positives ({len(false_positives)}):\n")
                    for i, seq_id in enumerate(false_positives[:10]):  # List first 10
                        f.write(f"  {seq_id}\n")
                    if len(false_positives) > 10:
                        f.write(f"  ... and {len(false_positives) - 10} more\n")
            
            # Create a separate file listing all potential false positives
            fp_list_file = os.path.join(self.false_positive_dir, "potential_false_positives.txt")
            with open(fp_list_file, "w") as f:
                f.write("# Potential false positive sequences\n")
                f.write("# Format: sequence_id\tconfidence_score\tevidence\n\n")
                
                for eval_data in sorted(all_evaluations, key=lambda x: x['final_score']):
                    if eval_data['classification'] == "Potential false positive":
                        # Find evidence from the original TSV file for this sequence
                        evidence = "No evidence found"
                        for fp_file in false_positive_files:
                            with open(fp_file, 'r') as fp_f:
                                # Skip header
                                next(fp_f)
                                for line in fp_f:
                                    fields = line.strip().split('\t')
                                    if len(fields) >= 9 and fields[0] == eval_data['seq_id']:
                                        evidence = fields[8]
                                        break
                        
                        f.write(f"{eval_data['seq_id']}\t{eval_data['final_score']:.3f}\t{evidence}\n")
            
            # Create detailed HTML report
            self._create_fp_html_report(all_evaluations)


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
    parser.add_argument("--advanced-tsd", action="store_true",
                       help="Use advanced TSD detection instead of PSSM")
    parser.add_argument("--weighted-evidence", action="store_true",
                       help="Use weighted evidence integration for boundary detection")
    
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
        advanced_tsd=args.advanced_tsd,
        weighted_evidence=args.weighted_evidence
    )
    
    # Run optimization
    optimizer.optimize()


if __name__ == "__main__":
    main()
