#!/usr/bin/env python3
"""
LTR Clustering Tool

This script performs clustering analysis on LTR sequences using various
algorithms and refinement techniques.
"""

import os
import logging
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class ClusteringConfig:
    """Configuration for clustering parameters."""
    min_identity: float = 0.8
    min_coverage: float = 0.5
    batch_size: int = 10000
    max_family_size: int = 100
    min_family_size: int = 15
    threads: Optional[int] = None

@dataclass
class ClusteringResult:
    """Results from clustering analysis."""
    clusters: List[List[str]]
    stats: Dict[str, float]
    consensus_file: str

class LTRClusterer:
    """Main class for LTR sequence clustering."""
    
    def __init__(self, config: ClusteringConfig):
        """
        Initialize the clusterer.
        
        Args:
            config: Configuration parameters for clustering
        """
        self.config = config
        self.config.threads = config.threads or multiprocessing.cpu_count()
        self.sequences = {}
        self.distance_matrix = None
        
    def cluster_sequences(self, input_file: Path) -> ClusteringResult:
        """
        Perform clustering on input sequences.
        
        Args:
            input_file: Path to input FASTA file
            
        Returns:
            ClusteringResult object containing clusters and statistics
        """
        logger.info(f"Starting clustering analysis on {input_file}")
        
        # Load sequences
        self.sequences = self._load_sequences(input_file)
        
        # Calculate distance matrix
        self.distance_matrix = self._calculate_distance_matrix()
        
        # Perform clustering
        clusters = self._perform_clustering()
        
        # Generate consensus sequences
        consensus_file = self._generate_consensus(clusters)
        
        # Calculate statistics
        stats = self._calculate_statistics(clusters)
        
        return ClusteringResult(
            clusters=clusters,
            stats=stats,
            consensus_file=consensus_file
        )
    
    def _load_sequences(self, input_file: Path) -> Dict[str, SeqRecord]:
        """Load sequences from FASTA file."""
        logger.info("Loading sequences...")
        return SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    
    def _calculate_distance_matrix(self) -> np.ndarray:
        """Calculate distance matrix between sequences."""
        logger.info("Calculating distance matrix...")
        
        seq_count = len(self.sequences)
        distances = np.zeros((seq_count, seq_count))
        
        def process_batch(batch_idx):
            batch_distances = np.zeros((len(batch_idx), seq_count))
            for i, idx1 in enumerate(batch_idx):
                seq1 = self.sequences[idx1]
                for j, seq2 in enumerate(self.sequences.values()):
                    if idx1 == seq2.id:
                        continue
                    dist = self._calculate_sequence_distance(seq1, seq2)
                    batch_distances[i, j] = dist
            return batch_distances
        
        # Process in parallel
        with ThreadPoolExecutor(max_workers=self.config.threads) as executor:
            batch_results = list(executor.map(
                process_batch,
                np.array_split(list(self.sequences.keys()), 
                             self.config.threads)
            ))
        
        # Combine results
        for i, result in enumerate(batch_results):
            start_idx = i * (seq_count // self.config.threads)
            distances[start_idx:start_idx + result.shape[0]] = result
        
        return distances
    
    def _calculate_sequence_distance(
        self, seq1: SeqRecord, seq2: SeqRecord
    ) -> float:
        """Calculate distance between two sequences."""
        # Implement your distance calculation here
        # This is a placeholder for demonstration
        return 0.0
    
    def _perform_clustering(self) -> List[List[str]]:
        """Perform hierarchical clustering."""
        logger.info("Performing clustering...")
        
        # Convert distance matrix to condensed form
        condensed_distances = squareform(self.distance_matrix)
        
        # Perform hierarchical clustering
        linkage_matrix = linkage(condensed_distances, method='average')
        
        # Cut tree to get clusters
        labels = fcluster(
            linkage_matrix,
            t=1 - self.config.min_identity,
            criterion='distance'
        )
        
        # Group sequences into clusters
        clusters = [[] for _ in range(max(labels))]
        for seq_id, label in zip(self.sequences.keys(), labels):
            clusters[label - 1].append(seq_id)
        
        return [c for c in clusters if len(c) >= self.config.min_family_size]
    
    def _generate_consensus(self, clusters: List[List[str]]) -> str:
        """Generate consensus sequences for clusters."""
        logger.info("Generating consensus sequences...")
        
        output_file = "consensus_sequences.fasta"
        
        with open(output_file, 'w') as out_fh:
            for i, cluster in enumerate(clusters):
                consensus = self._generate_cluster_consensus(cluster)
                # Continued from previous block
                SeqIO.write(consensus, out_fh, "fasta")
                
        return output_file
    
    def _generate_cluster_consensus(self, cluster: List[str]) -> SeqRecord:
        """Generate consensus sequence for a single cluster."""
        if len(cluster) == 1:
            return self.sequences[cluster[0]]
            
        # For multiple sequences, perform multiple alignment
        sequences = [self.sequences[seq_id] for seq_id in cluster]
        alignment = self._align_sequences(sequences)
        consensus_seq = self._calculate_consensus(alignment)
        
        return SeqRecord(
            Seq(consensus_seq),
            id=f"cluster_{len(cluster)}_consensus",
            description=f"consensus from {len(cluster)} sequences"
        )
    
    def _align_sequences(self, sequences: List[SeqRecord]) -> List[SeqRecord]:
        """Perform multiple sequence alignment."""
        try:
            from Bio.Align.Applications import MafftCommandline
            
            # Create temporary file for MAFFT input
            temp_file = "temp_seqs.fasta"
            SeqIO.write(sequences, temp_file, "fasta")
            
            # Run MAFFT
            mafft_cline = MafftCommandline(
                input=temp_file,
                adjustdirection=True,
                thread=self.config.threads
            )
            stdout, stderr = mafft_cline()
            
            # Parse results
            from io import StringIO
            align = list(SeqIO.parse(StringIO(stdout), "fasta"))
            
            # Cleanup
            os.remove(temp_file)
            
            return align
            
        except Exception as e:
            logger.error(f"Alignment failed: {e}")
            return sequences
    
    def _calculate_consensus(self, alignment: List[SeqRecord]) -> str:
        """Calculate consensus sequence from alignment."""
        if not alignment:
            return ""
            
        # Convert alignment to matrix
        alignment_length = len(alignment[0].seq)
        sequences = [str(record.seq) for record in alignment]
        
        consensus = []
        for i in range(alignment_length):
            # Count bases at this position
            bases = [seq[i] for seq in sequences if i < len(seq)]
            base_counts = {}
            
            for base in bases:
                if base != '-':
                    base_counts[base] = base_counts.get(base, 0) + 1
                    
            # Select most common base
            if base_counts:
                max_base = max(base_counts.items(), key=lambda x: x[1])[0]
                consensus.append(max_base)
                
        return ''.join(consensus)
    
    def _calculate_statistics(
        self, clusters: List[List[str]]
    ) -> Dict[str, float]:
        """Calculate clustering statistics."""
        return {
            'total_clusters': len(clusters),
            'avg_cluster_size': np.mean([len(c) for c in clusters]),
            'max_cluster_size': max(len(c) for c in clusters),
            'min_cluster_size': min(len(c) for c in clusters),
            'total_sequences': sum(len(c) for c in clusters)
        }

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Cluster LTR sequences and generate consensus sequences"
    )
    
    parser.add_argument(
        '--fasta',
        required=True,
        type=Path,
        help='Input FASTA file'
    )
    
    parser.add_argument(
        '--min-identity',
        type=float,
        default=0.8,
        help='Minimum sequence identity (default: 0.8)'
    )
    
    parser.add_argument(
        '--min-coverage',
        type=float,
        default=0.5,
        help='Minimum sequence coverage (default: 0.5)'
    )
    
    parser.add_argument(
        '--batch-size',
        type=int,
        default=10000,
        help='Batch size for processing (default: 10000)'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=None,
        help='Number of threads (default: all available)'
    )
    
    parser.add_argument(
        '--max-family-size',
        type=int,
        default=100,
        help='Maximum sequences per family (default: 100)'
    )
    
    parser.add_argument(
        '--min-family-size',
        type=int,
        default=15,
        help='Minimum sequences to form family (default: 15)'
    )
    
    return parser.parse_args()

def main() -> None:
    """Main execution function."""
    args = parse_arguments()
    
    # Create configuration
    config = ClusteringConfig(
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        batch_size=args.batch_size,
        threads=args.threads,
        max_family_size=args.max_family_size,
        min_family_size=args.min_family_size
    )
    
    try:
        # Initialize clusterer
        clusterer = LTRClusterer(config)
        
        # Perform clustering
        result = clusterer.cluster_sequences(args.fasta)
        
        # Log results
        logger.info("Clustering completed successfully:")
        for stat, value in result.stats.items():
            logger.info(f"  {stat}: {value}")
        logger.info(f"Consensus sequences written to: {result.consensus_file}")
        
    except Exception as e:
        logger.error(f"Clustering failed: {e}")
        raise

if __name__ == "__main__":
    main()