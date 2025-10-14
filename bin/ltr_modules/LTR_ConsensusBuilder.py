#!/usr/bin/env python3
"""
LTR Consensus Builder Module

Handles LTR consensus sequence building using multiple strategies.

This module contains core consensus building functionality extracted from
build_ltr_consensus.py for better code organization and reusability.
"""

import os
import re
import logging
import subprocess
import shutil
from typing import List, Tuple, Optional
from .shared_utils import LTRConstants, validate_sequence


class LTRConsensusBuilder:
    """
    LTR Consensus Sequence Builder

    Provides multiple strategies for building consensus sequences:
    - Simple majority-rule consensus
    - Weighted consensus (quality-based)
    - MSA-based consensus (MAFFT)
    - CD-HIT clustering
    """

    def __init__(self, min_frac_major: float = 0.60, threads: int = 1, logger=None):
        """
        Initialize the consensus builder.

        Args:
            min_frac_major: Minimum fraction for majority base (0-1)
            threads: Number of threads for external tools
            logger: Logger instance (optional)
        """
        self.min_frac_major = min_frac_major
        self.threads = threads

        # Setup logger
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger('LTRConsensusBuilder')
            if not self.logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                self.logger.addHandler(handler)
                self.logger.setLevel(logging.INFO)

    def build_consensus(self, sequences: List[str], method: str = 'weighted') -> str:
        """
        Build consensus sequence using specified method.

        PUBLIC API - Main entry point for consensus building.

        Args:
            sequences: List of sequences to build consensus from
            method: Consensus method ('simple', 'weighted', 'msa')

        Returns:
            str: Consensus sequence
        """
        if not sequences:
            self.logger.debug("No sequences provided for consensus building")
            return ''

        if len(sequences) == 1:
            self.logger.debug(f"Single sequence, no consensus needed (length={len(sequences[0])}bp)")
            return sequences[0]

        # Log consensus building start
        seq_lens = [len(s) for s in sequences]
        avg_len = sum(seq_lens) / len(seq_lens)
        self.logger.info(f"Building consensus from {len(sequences)} sequences using '{method}' method")
        self.logger.debug(f"  Sequence lengths: avg={avg_len:.0f}bp, min={min(seq_lens)}bp, max={max(seq_lens)}bp")

        if method == 'simple':
            consensus = self.simple_consensus(sequences, self.min_frac_major)
        elif method == 'weighted':
            consensus = self.weighted_consensus(sequences, self.min_frac_major)
        elif method == 'msa':
            # MSA requires temp directory
            import tempfile
            with tempfile.TemporaryDirectory() as tmpdir:
                seq_tuples = [(f"seq_{i}", seq) for i, seq in enumerate(sequences)]
                consensus, _ = self.msa_consensus(seq_tuples, tmpdir)
        else:
            raise ValueError(f"Unknown consensus method: {method}")

        # Log consensus result
        n_ratio = consensus.count('N') / len(consensus) if consensus else 0
        self.logger.info(f"Consensus built: {len(consensus)}bp, N_ratio={n_ratio:.1%}")

        return consensus

    def simple_consensus(self, sequences: List[str], min_frac_major: float = 0.60) -> str:
        """
        Build consensus using majority rule.

        Args:
            sequences: List of sequences
            min_frac_major: Minimum fraction for majority base (0-1)

        Returns:
            str: Consensus sequence
        """
        if not sequences:
            return ''

        # Get maximum length
        max_len = max(len(s) for s in sequences)
        if max_len == 0:
            return ''

        # Pad sequences to same length
        padded = [s.ljust(max_len, 'N') for s in sequences]
        consensus = []

        for i in range(max_len):
            bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            total = 0

            for seq in padded:
                if i < len(seq):
                    base = seq[i].upper()
                    if base in bases:
                        bases[base] += 1
                        total += 1

            if total == 0:
                consensus.append('N')
                continue

            # Find most frequent base - always use the most common, no threshold
            most_freq_base = max(bases.items(), key=lambda x: x[1])
            consensus.append(most_freq_base[0])

        return ''.join(consensus).rstrip('N')

    def weighted_consensus(self, sequences: List[str], min_frac_major: float = 0.60) -> str:
        """
        Build consensus with sequence quality weighting.

        Weights are based on:
        1. Sequence completeness (fewer N's = higher weight)
        2. Sequence length relative to average

        Args:
            sequences: List of sequences
            min_frac_major: Minimum fraction for majority base (0-1)

        Returns:
            str: Consensus sequence
        """
        if not sequences:
            return ''

        if len(sequences) == 1:
            return sequences[0]

        # Calculate weights for each sequence
        weights = []
        for seq in sequences:
            weight = 1.0

            # Completeness weight: penalize sequences with many N's
            n_ratio = seq.count('N') / len(seq) if len(seq) > 0 else 1.0
            if n_ratio < 0.10:
                weight *= 1.2  # Bonus for high quality
            elif n_ratio > 0.30:
                weight *= 0.5  # Penalty for low quality

            # Length weight: penalize very short sequences
            avg_len = sum(len(s) for s in sequences) / len(sequences)
            if len(seq) < avg_len * 0.5:
                weight *= 0.6  # Penalty for truncated sequences

            weights.append(weight)

        # Normalize weights
        total_weight = sum(weights)
        if total_weight > 0:
            weights = [w / total_weight for w in weights]
        else:
            weights = [1.0 / len(sequences)] * len(sequences)

        # Get maximum length
        max_len = max(len(s) for s in sequences)
        if max_len == 0:
            return ''

        # Pad sequences to same length
        padded = [s.ljust(max_len, 'N') for s in sequences]
        consensus = []

        for i in range(max_len):
            base_weights = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}
            total_valid_weight = 0.0

            for seq, weight in zip(padded, weights):
                if i < len(seq):
                    base = seq[i].upper()
                    if base in base_weights:
                        base_weights[base] += weight
                        total_valid_weight += weight

            if total_valid_weight == 0:
                consensus.append('N')
                continue

            # Find most weighted base - always use the most weighted, no threshold
            most_weighted_base = max(base_weights.items(), key=lambda x: x[1])
            consensus.append(most_weighted_base[0])

        return ''.join(consensus).rstrip('N')

    def run_mafft(self, in_fa: str, out_fa: str) -> bool:
        """
        Run MAFFT alignment.

        Args:
            in_fa: Input FASTA file
            out_fa: Output aligned FASTA file

        Returns:
            bool: True if successful
        """
        exe = shutil.which('mafft')
        if not exe:
            self.logger.warning("MAFFT not found in PATH")
            return False

        if not os.path.exists(in_fa) or os.path.getsize(in_fa) == 0:
            return False

        # Use thread parameter for parallel execution
        cmd = [exe, '--auto', '--quiet', '--thread', str(max(1, self.threads)), in_fa]
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  text=True, check=True)
            with open(out_fa, 'w') as f:
                f.write(result.stdout)
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"MAFFT failed: {e}")
            return False

    def msa_consensus(self, sequences: List[Tuple[str, str]], tempdir: str) -> Tuple[str, int]:
        """
        MSA-based consensus with fallback to simple consensus.

        Args:
            sequences: List of (header, sequence) tuples
            tempdir: Temporary directory for intermediate files

        Returns:
            tuple: (consensus_sequence, num_sequences)
        """
        if not sequences:
            return ('', 0)

        if len(sequences) == 1:
            return (sequences[0][1], 1)

        # Write input file
        in_fa = os.path.join(tempdir, 'input.fa')
        with open(in_fa, 'w') as f:
            for header, seq in sequences:
                f.write(f'>{header}\n{seq}\n')

        # Try MAFFT alignment
        aln_fa = os.path.join(tempdir, 'aligned.fa')
        success = self.run_mafft(in_fa, aln_fa)

        if success and os.path.exists(aln_fa):
            # Parse aligned sequences
            aligned_seqs = []
            with open(aln_fa) as f:
                seq = ''
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if seq:
                            aligned_seqs.append(seq)
                            seq = ''
                    else:
                        seq += line
                if seq:
                    aligned_seqs.append(seq)

            if aligned_seqs:
                # Use weighted consensus for better representativeness
                consensus = self.weighted_consensus(aligned_seqs, self.min_frac_major)
                return (consensus, len(sequences))

        # Fallback to simple consensus without alignment
        seqs_only = [seq for _, seq in sequences]
        consensus = self.simple_consensus(seqs_only, min_frac_major=0.50)  # Lower threshold for unaligned
        return (consensus, len(sequences))

    def cluster_sequences(self, fasta_file: str, output_prefix: str,
                         identity: float = 0.85, coverage: float = 0.80) -> str:
        """
        Cluster sequences using CD-HIT-EST.

        PUBLIC API - Main entry point for clustering.

        Args:
            fasta_file: Input FASTA file
            output_prefix: Output file prefix
            identity: Sequence identity threshold (0-1)
            coverage: Coverage threshold (0-1)

        Returns:
            str: Path to cluster file (.clstr)
        """
        return self.run_cdhit_est(fasta_file, output_prefix, identity, coverage)

    def run_cdhit_est(self, fa_path: str, out_prefix: str,
                      identity: float = 0.85, coverage: float = 0.80) -> str:
        """
        Run CD-HIT-EST clustering.

        Args:
            fa_path: Input FASTA file
            out_prefix: Output prefix
            identity: Sequence identity threshold (0-1)
            coverage: Coverage threshold (0-1)

        Returns:
            str: Path to cluster file
        """
        exe = shutil.which('cd-hit-est')
        if not exe:
            raise RuntimeError("cd-hit-est not found in PATH")

        cmd = [exe, '-i', fa_path, '-o', out_prefix, '-c', str(identity),
               '-aS', str(coverage), '-T', str(max(1, self.threads)),
               '-M', '0', '-d', '0']

        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"CD-HIT-EST failed for {fa_path}")
            raise

        return out_prefix + '.clstr'

    def parse_clusters(self, clstr_file: str) -> List[List[str]]:
        """
        Parse CD-HIT cluster file.

        Args:
            clstr_file: Path to .clstr file

        Returns:
            list: List of clusters, each cluster is a list of sequence IDs
        """
        clusters = []
        current_cluster = []

        with open(clstr_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>Cluster'):
                    if current_cluster:
                        clusters.append(current_cluster)
                        current_cluster = []
                else:
                    # Extract sequence ID
                    match = re.search(r'>([^\.]+)\.\.\.', line)
                    if match:
                        current_cluster.append(match.group(1))

        if current_cluster:
            clusters.append(current_cluster)

        return clusters

    def validate_consensus(self, consensus: str) -> Tuple[bool, str]:
        """
        Validate consensus sequence for extreme cases.

        Args:
            consensus: Consensus sequence

        Returns:
            tuple: (is_valid, reason)
        """
        return validate_sequence(consensus, min_length=10, max_n_content=0.95)
