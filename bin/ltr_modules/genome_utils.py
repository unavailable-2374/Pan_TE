#!/usr/bin/env python3
"""
Genome Access Utilities Module

Provides efficient genome sequence access classes for LTR processing.

This module contains genome access classes extracted from
LTR_Boundary_Optimizer.py for better code organization and reusability.
"""

import os
import time
import random
import subprocess
from Bio import SeqIO


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
    """
    Simple genome access using BioPython with optimized initialization.

    Supports two modes:
    1. Fast mode (with seq_names provided): Skip full genome scan
    2. Fallback mode (no seq_names): Scan genome to get sequence IDs
    """
    def __init__(self, genome_file, chunk_size, seq_names=None):
        self.genome_file = genome_file
        self.chunk_size = chunk_size
        self.cache = {}

        # OPTIMIZATION: Accept pre-computed sequence names to avoid scanning
        if seq_names is not None:
            self.seq_names = set(seq_names) if not isinstance(seq_names, set) else seq_names
            self.seq_dict = {}
        else:
            # Fallback: Scan genome file (slower)
            self._scan_genome()

    def _scan_genome(self):
        """
        Scan genome file to build index of sequence IDs.

        WARNING: This is expensive for large genomes. Prefer passing seq_names to __init__.
        """
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
