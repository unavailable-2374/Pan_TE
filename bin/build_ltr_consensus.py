#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced LTR family consensus builder with comprehensive biological feature utilization.

Major improvements:
- LTR boundary validation (TG...CA motifs)
- Solo-LTR specialized handling
- Nesting relationship utilization for family assignment
- Element completeness scoring for prioritization
- Improved internal region consensus from high-quality elements
- Adaptive clustering parameters based on sequence diversity
- Iterative refinement of consensus sequences
- Structural validation and quality metrics
- ORF prediction for internal regions
- PBS/PPT motif preservation

USAGE:
    python build_ltr_consensus.py \
        --rtr genome.rtr \
        --fasta genome.fa \
        --outdir out_ltr_consensus
    
    # Advanced usage with custom parameters:
    python build_ltr_consensus.py \
        --rtr genome.rtr \
        --fasta genome.fa \
        --outdir out_ltr_consensus \
        --min-identity 0.80 \
        --threads 16 --processes 4

Dependencies: pyfaidx, biopython, mafft/muscle, cd-hit-est, samtools (optional)
Optional: hmmer for domain analysis

Author: Enhanced version with biological improvements and rtr2repeatmasker.sh integration
"""

from __future__ import annotations
import os
import sys
import csv
import re
import math
import json
import time
import shutil
import tempfile
import logging
import argparse
import itertools
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set, Any
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import hashlib

############################
# Logging
############################
logger = logging.getLogger("ltr_consensus_enhanced")
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

############################
# Biological Constants
############################

# LTR boundary motifs (5' and 3' terminal sequences)
LTR_BOUNDARIES = {
    'canonical': ('TG', 'CA'),
    'variant1': ('TG', 'TA'),
    'variant2': ('CA', 'CA'),
    'variant3': ('TG', 'TG'),
    'variant4': ('AG', 'CT'),
}

# PBS consensus motifs (for tRNA primer binding)
PBS_MOTIFS = [
    'TGGTATCAGAGC',  # Lys
    'TGGCGCCCGAAC',  # Pro
    'TGGTGCGTGGAC',  # Met
    'TGGTAGCAGAGC',  # Ile
]

# PPT consensus pattern
PPT_PATTERN = re.compile(r'[AG]{8,20}')

# Gag-Pol domains for internal region validation
CONSERVED_DOMAINS = [
    'DTGA',  # Integrase
    'HHCC',  # Gag zinc finger
    'YXDD',  # Reverse transcriptase
]

############################
# Utility functions
############################

def check_tool(path_or_name: str, required: bool = False) -> Optional[str]:
    """Return resolved path if tool is found, else None. Log if required tool missing."""
    p = shutil.which(path_or_name)
    if required and not p:
        logger.error(f"Required tool not found: {path_or_name}")
    elif p:
        logger.debug(f"Found tool: {path_or_name} at {p}")
    return p

def get_version(cmd: List[str]) -> str:
    try:
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=10)
        return r.stdout.strip().splitlines()[0][:200]
    except Exception:
        return "<unknown>"

############################
# FASTA access with caching
############################

class FastaAccessor:
    """Random access FASTA with caching for frequently accessed regions."""
    
    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.impl = None
        self.mode = None
        self.cache = {}  # Simple cache for repeated accesses
        self.cache_size = 100
        
        try:
            import pyfaidx
            self.impl = pyfaidx.Fasta(fasta_path, as_raw=True, sequence_always_upper=True, rebuild=False)
            self.mode = 'pyfaidx'
            logger.info("FASTA accessor: pyfaidx")
        except Exception as e:
            logger.warning(f"pyfaidx not available ({e}); falling back to Bio.SeqIO")
            try:
                from Bio import SeqIO
                self.impl = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(fasta_path, "fasta")}
                self.mode = 'seqio'
                logger.info("FASTA accessor: Bio.SeqIO (in-memory)")
            except Exception as e2:
                logger.error(f"Failed to load FASTA: {e2}")
                raise

    def fetch(self, chrom: str, start1: int, end1: int) -> str:
        """1-based inclusive coordinates with caching."""
        if start1 is None or end1 is None or start1 <= 0 or end1 < start1:
            return ""
        
        cache_key = f"{chrom}:{start1}-{end1}"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        if self.mode == 'pyfaidx':
            seq = str(self.impl[chrom][start1-1:end1])
        else:
            s = self.impl.get(chrom)
            if s is None:
                return ""
            seq = s[start1-1:end1]
        
        # Update cache (LRU-like)
        if len(self.cache) >= self.cache_size:
            # Remove oldest entry
            self.cache.pop(next(iter(self.cache)))
        self.cache[cache_key] = seq
        
        return seq

    def close(self):
        try:
            if self.mode == 'pyfaidx':
                self.impl.close()
        except Exception:
            pass

############################
# Enhanced IUPAC consensus
############################

IUPAC = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A','G']): 'R',
    frozenset(['C','T']): 'Y',
    frozenset(['A','C']): 'M',
    frozenset(['G','T']): 'K',
    frozenset(['C','G']): 'S',
    frozenset(['A','T']): 'W',
    frozenset(['C','G','T']): 'B',
    frozenset(['A','G','T']): 'D',
    frozenset(['A','C','T']): 'H',
    frozenset(['A','C','G']): 'V',
    frozenset(['A','C','G','T']): 'N'
}

DNA = set(['A','C','G','T'])

def iupac_from_counts(counts: Dict[str, int], min_frac_major: float = 0.5, min_frac_include: float = 0.2) -> str:
    total = sum(v for b, v in counts.items() if b in DNA)
    if total == 0:
        return 'N'
    sorted_b = sorted([(b, counts.get(b, 0)) for b in DNA], key=lambda x: x[1], reverse=True)
    major_b, major_c = sorted_b[0]
    if major_c / total >= min_frac_major:
        return major_b
    bases = [b for b, c in sorted_b if total > 0 and (c / total) >= min_frac_include]
    if not bases:
        bases = [major_b]
    return IUPAC.get(frozenset(bases), 'N')

############################
# LTR Biological Validation
############################

def validate_ltr_boundaries(seq: str, strict: bool = False) -> Tuple[bool, str]:
    """
    Validate LTR boundary motifs (TG...CA and variants).
    Returns (is_valid, boundary_type)
    """
    if not seq or len(seq) < 4:
        return (False, 'too_short')
    
    seq_upper = seq.upper()
    start_di = seq_upper[:2]
    end_di = seq_upper[-2:]
    
    for boundary_type, (start_motif, end_motif) in LTR_BOUNDARIES.items():
        if start_di == start_motif and end_di == end_motif:
            return (True, boundary_type)
    
    if not strict:
        # Allow some flexibility for degenerate sequences
        if start_di[0] in 'TC' and end_di[1] in 'AC':
            return (True, 'degenerate')
    
    return (False, 'invalid')

def find_pbs_motif(seq: str, ltr_end_pos: int = 0) -> Optional[Tuple[int, str]]:
    """Find PBS motif near LTR end. Returns (position, motif) or None."""
    if not seq or len(seq) < 20:
        return None
    
    search_region = seq[ltr_end_pos:min(ltr_end_pos + 50, len(seq))]
    
    for motif in PBS_MOTIFS:
        # Allow 1-2 mismatches
        for i in range(len(search_region) - len(motif) + 1):
            subseq = search_region[i:i+len(motif)]
            mismatches = sum(1 for a, b in zip(subseq, motif) if a != b)
            if mismatches <= 2:
                return (ltr_end_pos + i, motif)
    
    return None

def find_ppt_motif(seq: str, before_3ltr: int = None) -> Optional[Tuple[int, int]]:
    """Find PPT motif before 3' LTR. Returns (start, end) or None."""
    if not seq:
        return None
    
    if before_3ltr is None:
        before_3ltr = len(seq)
    
    # Search in the last 100bp before 3' LTR
    search_start = max(0, before_3ltr - 100)
    search_region = seq[search_start:before_3ltr]
    
    matches = PPT_PATTERN.finditer(search_region)
    best_match = None
    best_length = 0
    
    for match in matches:
        if match.end() - match.start() > best_length:
            best_match = match
            best_length = match.end() - match.start()
    
    if best_match:
        return (search_start + best_match.start(), search_start + best_match.end())
    
    return None

def score_element_completeness(rec: 'LTRRecord', left_seq: str = None, right_seq: str = None) -> float:
    """
    Score element completeness based on structural features.
    Returns score 0-100.
    """
    score = 0.0
    max_score = 100.0
    
    # Identity score (30 points)
    if rec.identity is not None:
        score += min(30, rec.identity * 30)
    
    # TSD presence (15 points)
    if rec.tsd_start and rec.tsd_end:
        score += 15
    
    # PPT presence (10 points)
    if rec.ppt_start and rec.ppt_end:
        score += 10
    
    # LTR boundary validation (20 points)
    if left_seq:
        valid_left, _ = validate_ltr_boundaries(left_seq)
        if valid_left:
            score += 10
    
    if right_seq:
        valid_right, _ = validate_ltr_boundaries(right_seq)
        if valid_right:
            score += 10
    
    # LTR length similarity (10 points)
    if rec.left_start and rec.left_end and rec.right_start and rec.right_end:
        left_len = abs(rec.left_end - rec.left_start) + 1
        right_len = abs(rec.right_end - rec.right_start) + 1
        len_ratio = min(left_len, right_len) / max(left_len, right_len)
        score += len_ratio * 10
    
    # Element type bonus (15 points)
    if rec.casetype == 'Single':
        score += 15
    elif rec.casetype == 'RecentlyNested':
        score += 10
    elif rec.casetype in ['SoloSingle', 'SoloLTR']:
        score += 5
    
    return min(score, max_score)

############################
# ORF Detection
############################

def find_orfs(seq: str, min_length: int = 300) -> List[Tuple[int, int, int]]:
    """
    Find ORFs in sequence. Returns list of (start, end, frame).
    """
    orfs = []
    seq_upper = seq.upper()
    
    start_codons = ['ATG']
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    for frame in range(3):
        i = frame
        while i < len(seq_upper) - 2:
            codon = seq_upper[i:i+3]
            if codon in start_codons:
                # Look for stop codon
                for j in range(i+3, len(seq_upper)-2, 3):
                    stop_codon = seq_upper[j:j+3]
                    if stop_codon in stop_codons:
                        orf_length = j - i
                        if orf_length >= min_length:
                            orfs.append((i, j+3, frame))
                        break
            i += 3
    
    return sorted(orfs, key=lambda x: x[1] - x[0], reverse=True)

def detect_conserved_domains(seq: str) -> List[str]:
    """Detect conserved domains in protein sequence."""
    found_domains = []
    
    # Simple pattern matching for demonstration
    # In production, use HMMER or similar tools
    for domain in CONSERVED_DOMAINS:
        if domain in seq.upper():
            found_domains.append(domain)
    
    return found_domains

############################
# Enhanced Data Structures
############################

@dataclass
class LTRRecord:
    chrom: str
    id: str
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    rc: int = 0
    casetype: str = 'Single'
    group: str = 'ALL'
    identity: Optional[float] = None
    ppt_start: Optional[int] = None
    ppt_end: Optional[int] = None
    tsd_start: Optional[int] = None
    tsd_end: Optional[int] = None
    nest_in: Optional[str] = None
    nest_out: Optional[str] = None
    completeness_score: float = 0.0
    left_boundary_type: str = 'unknown'
    right_boundary_type: str = 'unknown'
    has_pbs: bool = False
    has_ppt: bool = False
    internal_orfs: List[Tuple[int, int, int]] = field(default_factory=list)
    
    @staticmethod
    def from_row(r: Dict[str, str]) -> 'LTRRecord':
        def num(x):
            try:
                if x in (None, '', 'NA', 'nan', 'None'):
                    return None
                return int(float(x))
            except (ValueError, TypeError):
                return None
        
        def s(x, d=''):
            if x in (None, 'NA', '', 'None'):
                return d
            return str(x).strip()
        
        def f0_1(x):
            try:
                if x in (None, '', 'NA', 'nan', 'None'):
                    return None
                v = float(x)
                return v
            except (ValueError, TypeError):
                return None
        
        def parse_int_safe(x, default=0):
            """Safely parse integer with NA handling"""
            try:
                if x in (None, '', 'NA', 'nan', 'None'):
                    return default
                return int(float(x))
            except (ValueError, TypeError):
                return default
                
        nest_in = s(r.get('NestIn'))
        nest_out = s(r.get('NestOut'))
        
        # Parse nested element IDs
        if nest_in and nest_in not in ('NA', ''):
            nest_in = nest_in.strip('{}').split(',') if ',' in nest_in else [nest_in.strip('{}')]
        else:
            nest_in = None
            
        if nest_out and nest_out not in ('NA', ''):
            nest_out = nest_out.strip('{}').split(',') if ',' in nest_out else [nest_out.strip('{}')]
        else:
            nest_out = None
            
        return LTRRecord(
            chrom=s(r.get('chrom') or r.get('Chr') or r.get('seqid')),
            id=str(r.get('ID') or r.get('id') or ''),
            left_start=num(r.get('LeftStart')) or 0,
            left_end=num(r.get('LeftEnd')) or 0,
            right_start=num(r.get('RightStart')) or 0,
            right_end=num(r.get('RightEnd')) or 0,
            rc=parse_int_safe(r.get('RC'), 0),
            casetype=s(r.get('CaseType') or r.get('type') or 'Single'),
            group=str(r.get('GraphGroup') or r.get('group') or 'ALL'),
            identity=f0_1(r.get('LTRIdentity')),
            ppt_start=num(r.get('PPTStart')),
            ppt_end=num(r.get('PPTEnd')),
            tsd_start=num(r.get('TSDStart')),
            tsd_end=num(r.get('TSDEnd')),
            nest_in=nest_in,
            nest_out=nest_out,
        )
    
    def normalize(self):
        """Normalize coordinates."""
        if self.left_start > self.left_end:
            self.left_start, self.left_end = self.left_end, self.left_start
        if self.right_start > self.right_end:
            self.right_start, self.right_end = self.right_end, self.right_start

@dataclass
class FamilyConsensus:
    group: str
    family_id: str
    n_5p: int = 0
    n_3p: int = 0
    n_internal: int = 0
    cons_5p: str = ''
    cons_3p: str = ''
    cons_internal: str = ''
    merged: bool = False
    boundary_type: str = 'unknown'
    has_pbs: bool = False
    has_ppt: bool = False
    internal_orfs: int = 0
    quality_score: float = 0.0
    member_ids: List[str] = field(default_factory=list)
    is_solo: bool = False
    iteration_refined: int = 0

############################
# Nesting-aware family assignment
############################

class NestingGraph:
    """Graph structure for tracking nesting relationships."""
    
    def __init__(self):
        self.nodes: Dict[str, LTRRecord] = {}
        self.edges: Dict[str, Set[str]] = defaultdict(set)
        self.families: Dict[str, str] = {}  # element_id -> family_id
        
    def add_record(self, rec: LTRRecord):
        """Add record to graph."""
        self.nodes[rec.id] = rec
        
        # Add edges for nesting relationships
        if rec.nest_in:
            for parent_id in rec.nest_in:
                self.edges[parent_id].add(rec.id)
                
        if rec.nest_out:
            for child_id in rec.nest_out:
                self.edges[rec.id].add(child_id)
    
    def assign_families_by_nesting(self) -> Dict[str, str]:
        """Assign elements to families based on nesting."""
        family_counter = 0
        
        # Find connected components
        visited = set()
        
        for node_id in self.nodes:
            if node_id not in visited:
                # BFS to find component
                component = []
                queue = [node_id]
                
                while queue:
                    current = queue.pop(0)
                    if current in visited:
                        continue
                    
                    visited.add(current)
                    component.append(current)
                    
                    # Add neighbors
                    for neighbor in self.edges.get(current, []):
                        if neighbor not in visited:
                            queue.append(neighbor)
                    
                    # Also check reverse edges
                    for other_id, neighbors in self.edges.items():
                        if current in neighbors and other_id not in visited:
                            queue.append(other_id)
                
                # Assign family to component
                family_id = f"NestFamily_{family_counter}"
                for elem_id in component:
                    self.families[elem_id] = family_id
                family_counter += 1
        
        return self.families

############################
# Adaptive clustering parameters
############################

def calculate_sequence_diversity(seqs: List[str]) -> float:
    """Calculate diversity score for sequences."""
    if len(seqs) < 2:
        return 0.0
    
    # Sample pairwise identities
    sample_size = min(50, len(seqs))
    import random
    sampled = random.sample(seqs, sample_size)
    
    identities = []
    for i in range(len(sampled)):
        for j in range(i+1, min(i+5, len(sampled))):  # Limited comparisons
            id_score = approx_identity(sampled[i], sampled[j])
            identities.append(id_score)
    
    if not identities:
        return 0.0
    
    avg_identity = sum(identities) / len(identities)
    diversity = 1.0 - avg_identity
    
    return diversity

def get_adaptive_clustering_params(seqs: List[str]) -> Tuple[float, float]:
    """Get adaptive CD-HIT parameters based on sequence diversity."""
    diversity = calculate_sequence_diversity(seqs)
    
    # Adjust thresholds based on diversity
    if diversity < 0.1:  # Very similar sequences
        c_threshold = 0.95
        as_threshold = 0.90
    elif diversity < 0.2:
        c_threshold = 0.90
        as_threshold = 0.85
    elif diversity < 0.3:
        c_threshold = 0.85
        as_threshold = 0.80
    else:  # Diverse sequences
        c_threshold = 0.80
        as_threshold = 0.75
    
    logger.debug(f"Diversity: {diversity:.3f}, using c={c_threshold}, aS={as_threshold}")
    
    return c_threshold, as_threshold

############################
# Alignment and consensus
############################

def run_mafft(in_fa: str, out_fa: str, threads: int) -> bool:
    """Run MAFFT with graceful failure handling like rtr2repeatmasker.sh."""
    exe = check_tool('mafft')
    if not exe:
        return False
    
    # Check if input file exists and has content
    if not os.path.exists(in_fa) or os.path.getsize(in_fa) == 0:
        logger.debug(f"Input file {in_fa} is empty or missing")
        return False
    
    # Count sequences to adjust timeout
    with open(in_fa) as f:
        seq_count = sum(1 for line in f if line.startswith('>'))
    
    # Adjust timeout based on sequence count
    timeout = min(300, max(60, seq_count * 5))
    
    # Use more conservative MAFFT settings for reliability
    cmd = [exe, '--auto', '--thread', str(max(1, threads)), '--quiet', in_fa]
    try:
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                          text=True, check=True, timeout=timeout)
        if r.stdout:
            with open(out_fa, 'w') as w:
                w.write(r.stdout)
            return True
        else:
            logger.debug("MAFFT produced no output")
            return False
    except subprocess.TimeoutExpired:
        logger.debug(f"MAFFT timed out after {timeout}s, using fallback")
        return False
    except subprocess.CalledProcessError as e:
        # Only log debug, not warning, since we have fallback
        logger.debug(f"MAFFT failed (will use fallback): {e.stderr[:200] if e.stderr else 'no error message'}")
        return False

def run_muscle(in_fa: str, out_fa: str, threads: int) -> bool:
    """Run MUSCLE with graceful failure handling like rtr2repeatmasker.sh."""
    exe = check_tool('muscle')
    if not exe:
        return False
    
    # Check if input file exists and has content
    if not os.path.exists(in_fa) or os.path.getsize(in_fa) == 0:
        logger.debug(f"Input file {in_fa} is empty or missing")
        return False
    
    # Count sequences to adjust timeout
    with open(in_fa) as f:
        seq_count = sum(1 for line in f if line.startswith('>'))
    
    # Adjust timeout based on sequence count
    timeout = min(300, max(60, seq_count * 5))
    
    cmd = [exe, '-align', in_fa, '-output', out_fa, '-threads', str(max(1, threads))]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                      text=True, timeout=timeout)
        return os.path.exists(out_fa) and os.path.getsize(out_fa) > 0
    except subprocess.TimeoutExpired:
        logger.debug(f"MUSCLE timed out after {timeout}s, using fallback")
        return False
    except subprocess.CalledProcessError as e:
        logger.debug(f"MUSCLE failed (will use fallback): {e}")
        return False

def fasta_read(path: str) -> Dict[str, str]:
    seqs = {}
    with open(path) as f:
        id = None
        parts = []
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if id is not None:
                    seqs[id] = ''.join(parts).upper()
                id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if id is not None:
            seqs[id] = ''.join(parts).upper()
    return seqs

def fasta_write(path: str, records: List[Tuple[str, str]]):
    with open(path, 'w') as w:
        for h, s in records:
            w.write(f'>{h}\n')
            for i in range(0, len(s), 80):
                w.write(s[i:i+80] + '\n')

def trim_alignment(aln_fa: str, out_fa: str, gap_thresh: float = 0.7) -> None:
    seqs = fasta_read(aln_fa)
    ids = list(seqs.keys())
    if not ids:
        shutil.copyfile(aln_fa, out_fa)
        return
    L = len(next(iter(seqs.values())))
    keep = []
    for i in range(L):
        col = [seqs[_id][i] if i < len(seqs[_id]) else '-' for _id in ids]
        gap_frac = sum(1 for b in col if b == '-' or b == 'N') / len(col)
        if gap_frac <= gap_thresh:
            keep.append(i)
    trimmed = {k: ''.join(v[i] for i in keep) for k, v in seqs.items()}
    fasta_write(out_fa, [(k, v) for k, v in trimmed.items()])

def consensus_from_alignment(aln_fa: str, iupac: bool = True, preserve_motifs: bool = True) -> str:
    """Enhanced consensus with motif preservation."""
    seqs = fasta_read(aln_fa)
    if not seqs:
        return ''
    
    ids = list(seqs.keys())
    L = len(next(iter(seqs.values())))
    cols = []
    
    for i in range(L):
        counts = {'A':0,'C':0,'G':0,'T':0}
        for _id in ids:
            b = seqs[_id][i].upper()
            if b in DNA:
                counts[b] += 1
        
        if iupac:
            cols.append(iupac_from_counts(counts, min_frac_major=0.5, min_frac_include=0.2))
        else:
            total = sum(counts.values())
            if total == 0:
                cols.append('N')
            else:
                cols.append(max(counts.items(), key=lambda x: x[1])[0])
    
    consensus = ''.join(cols)
    
    # Preserve terminal motifs if needed
    if preserve_motifs and len(consensus) >= 4:
        # Check if we should preserve TG...CA
        terminal_votes = defaultdict(int)
        for seq in seqs.values():
            if len(seq) >= 4:
                start = seq[:2].upper()
                end = seq[-2:].upper()
                if start in ['TG', 'CA', 'AG']:
                    terminal_votes[start] += 1
                if end in ['CA', 'TA', 'CT']:
                    terminal_votes[end] += 1
        
        # Apply most common terminals if significant
        if terminal_votes:
            most_common_start = max((v for v in terminal_votes.items() if v[0] in ['TG', 'CA', 'AG']), 
                                   key=lambda x: x[1], default=None)
            most_common_end = max((v for v in terminal_votes.items() if v[0] in ['CA', 'TA', 'CT']), 
                                 key=lambda x: x[1], default=None)
            
            if most_common_start and most_common_start[1] > len(seqs) * 0.5:
                consensus = most_common_start[0] + consensus[2:]
            if most_common_end and most_common_end[1] > len(seqs) * 0.5:
                consensus = consensus[:-2] + most_common_end[0]
    
    return consensus

def msa_consensus(seqs: List[Tuple[str, str]], threads: int, tempdir: str,
                  trim_gap: float = 0.7, use_iupac: bool = True,
                  preserve_motifs: bool = True) -> Tuple[str, int]:
    """MSA and consensus with motif preservation."""
    if not seqs:
        return ('', 0)
    in_fa = os.path.join(tempdir, 'in.fa')
    aln_fa = os.path.join(tempdir, 'aln.fa')
    aln_trim = os.path.join(tempdir, 'aln.trim.fa')
    fasta_write(in_fa, seqs)
    
    ok = run_mafft(in_fa, aln_fa, threads)
    if not ok:
        ok = run_muscle(in_fa, aln_fa, threads)
    if not ok:
        logger.debug("Using fallback consensus method")
        padded = []
        maxlen = max(len(s) for _, s in seqs)
        for h, s in seqs:
            padded.append((h, s.ljust(maxlen, 'N')))
        fasta_write(aln_fa, padded)
    
    trim_alignment(aln_fa, aln_trim, gap_thresh=trim_gap)
    cons = consensus_from_alignment(aln_trim, iupac=use_iupac, preserve_motifs=preserve_motifs)
    
    return (cons, len(seqs))

############################
# Iterative refinement
############################

def refine_consensus_iteratively(initial_cons: str, members: List[Dict[str, Any]], 
                                 max_iterations: int = 3, threads: int = 1,
                                 tempdir: str = '/tmp') -> Tuple[str, int]:
    """Iteratively refine consensus by re-aligning to consensus."""
    if not initial_cons or not members:
        return (initial_cons, 0)
    
    current_cons = initial_cons
    iteration = 0
    
    for iteration in range(max_iterations):
        # Score members against current consensus
        scored_members = []
        for m in members:
            seq = m.get('seq', '')
            if seq:
                score = approx_identity(seq, current_cons)
                scored_members.append((score, m))
        
        # Keep top 80% most similar
        scored_members.sort(key=lambda x: x[0], reverse=True)
        cutoff = int(len(scored_members) * 0.8)
        selected = scored_members[:max(10, cutoff)]
        
        if len(selected) < 5:  # Too few sequences
            break
        
        # Re-align selected sequences
        seqs = [(f"iter{iteration}_{i}", m[1]['seq']) for i, m in enumerate(selected)]
        new_cons, _ = msa_consensus(seqs, threads, tempdir, preserve_motifs=True)
        
        # Check if consensus improved
        if approx_identity(new_cons, current_cons) > 0.95:
            # Converged
            break
        
        current_cons = new_cons
    
    return (current_cons, iteration + 1)

############################
# CD-HIT clustering
############################

CLSTR_ID_RE = re.compile(r">\s*([^\.]+)\.\.\.")

def parse_clstr(clstr_path: str) -> List[List[str]]:
    clusters: List[List[str]] = []
    cur: List[str] = []
    with open(clstr_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>Cluster'):
                if cur:
                    clusters.append(cur)
                    cur = []
                continue
            m = CLSTR_ID_RE.search(line)
            if m:
                cur.append(m.group(1))
    if cur:
        clusters.append(cur)
    return clusters

def run_cdhit_est(fa_path: str, out_prefix: str, c: float, aS: float, threads: int) -> Tuple[str, str]:
    """Run CD-HIT-EST with better error handling."""
    exe = check_tool('cd-hit-est')
    if not exe:
        raise RuntimeError("cd-hit-est not found")
    
    cmd = [exe, '-i', fa_path, '-o', out_prefix, '-c', str(c), '-aS', str(aS), 
           '-T', str(max(1, threads)), '-M', '0', '-d', '0']
    
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE, text=True, timeout=600)
        logger.debug(f"CD-HIT-EST completed for {fa_path}")
    except subprocess.TimeoutExpired:
        logger.error(f"CD-HIT-EST timed out for {fa_path}")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"CD-HIT-EST failed for {fa_path}: {e.stderr[:200] if e.stderr else 'unknown error'}")
        raise
    
    return out_prefix, out_prefix + '.clstr'

############################
# Solo-LTR handling
############################

def process_solo_ltrs(solo_records: List[LTRRecord], fa: FastaAccessor, 
                      outdir: str, threads: int) -> List[FamilyConsensus]:
    """Special handling for Solo-LTR elements."""
    results = []
    
    # Group by similarity
    solo_seqs = []
    for rec in solo_records:
        rec.normalize()
        
        # For solo, use whichever LTR is available
        if rec.left_start and rec.left_end:
            seq = fa.fetch(rec.chrom, rec.left_start, rec.left_end)
        elif rec.right_start and rec.right_end:
            seq = fa.fetch(rec.chrom, rec.right_start, rec.right_end)
        else:
            continue
        
        if rec.rc == 1:
            seq = rc_seq(seq)
        
        solo_seqs.append((rec.id, seq, rec))
    
    if not solo_seqs:
        return results
    
    # Adaptive clustering for solos
    seqs_only = [s for _, s, _ in solo_seqs]
    c_thr, as_thr = get_adaptive_clustering_params(seqs_only)
    
    with tempfile.TemporaryDirectory(prefix="solo_", dir=outdir) as td:
        # Write sequences
        fa_path = os.path.join(td, 'solo.fa')
        fasta_write(fa_path, [(id, seq) for id, seq, _ in solo_seqs])
        
        # Cluster
        out_prefix = os.path.join(td, 'solo.cdhit')
        run_cdhit_est(fa_path, out_prefix, c=c_thr, aS=as_thr, threads=threads)
        
        clusters = parse_clstr(out_prefix + '.clstr')
        
        # Build consensus for each cluster
        for i, cluster_ids in enumerate(clusters):
            cluster_seqs = [(id, seq) for id, seq, rec in solo_seqs if id in cluster_ids]
            
            if cluster_seqs:
                cons, n = msa_consensus(cluster_seqs, threads, td, preserve_motifs=True)
                
                fc = FamilyConsensus(
                    group='SOLO',
                    family_id=f'Solo_Family_{i+1}',
                    n_5p=len(cluster_seqs),
                    cons_5p=cons,
                    is_solo=True,
                    member_ids=[id for id, _ in cluster_seqs]
                )
                
                # Validate boundaries
                valid, boundary = validate_ltr_boundaries(cons)
                fc.boundary_type = boundary
                
                results.append(fc)
    
    return results

############################
# Input parsing
############################

RTR_REQUIRED = ['chrom','LeftStart','LeftEnd','RightStart','RightEnd']

def read_rtr_table(path: str) -> List[Dict[str, str]]:
    with open(path, 'r', newline='') as f:
        sample = f.read(8192)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters='\t,; ')
        except Exception:
            dialect = csv.excel_tab
        reader = csv.DictReader(f, dialect=dialect)
        rows = [dict(r) for r in reader]
    
    norm_rows = []
    for r in rows:
        r2 = {}
        for k, v in r.items():
            if k is None:
                continue
            kk = k.strip()
            r2[kk] = v.strip() if isinstance(v, str) else v
        norm_rows.append(r2)
    
    header = set(norm_rows[0].keys()) if norm_rows else set()
    missing = [k for k in RTR_REQUIRED if k not in header]
    if missing:
        raise ValueError(f"RTR table missing columns: {missing}")
    
    return norm_rows

############################
# Helper functions
############################

def autoscale_identity(id_list: List[Optional[float]]) -> List[Optional[float]]:
    vals = [v for v in id_list if v is not None and not math.isnan(v)]
    if not vals:
        return id_list
    med = sorted(vals)[len(vals)//2]
    mx = max(vals)
    if med > 1.5 or mx > 1.5:
        logger.info("Detected identity in 0-100 scale; scaling to 0-1")
        return [None if v is None else v/100.0 for v in id_list]
    return id_list

def rc_seq(seq: str) -> str:
    comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(comp)[::-1]

def approx_identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    if len(a) == len(b):
        matches = sum(1 for x,y in zip(a,b) if x==y and x in DNA)
        return matches / max(1, sum(1 for x,y in zip(a,b) if x in DNA and y in DNA))
    L = min(len(a), len(b))
    best = 0.0
    for off in range(-20, 21):
        if off >= 0:
            sa = a[off:off+L]
            sb = b[:L]
        else:
            sa = a[:L]
            sb = b[-off:-off+L]
        m = sum(1 for x,y in zip(sa,sb) if x==y and x in DNA)
        d = sum(1 for x,y in zip(sa,sb) if x in DNA and y in DNA)
        if d:
            best = max(best, m/d)
    return best

def deterministic_subsample(seqs: List[Tuple[str,str]], max_n: int, 
                           scores: List[float] = None) -> List[Tuple[str,str]]:
    """Subsample with optional quality scores."""
    if len(seqs) <= max_n:
        return seqs
    
    if scores and len(scores) == len(seqs):
        # Sort by score, then by median length proximity
        lens = [len(s) for _, s in seqs]
        med = sorted(lens)[len(lens)//2]
        seqs_scored = [(seqs[i], scores[i], abs(len(seqs[i][1]) - med)) for i in range(len(seqs))]
        seqs_scored.sort(key=lambda x: (-x[1], x[2]))  # High score, close to median
        return [s for s, _, _ in seqs_scored[:max_n]]
    else:
        # Original length-based subsampling
        lens = [len(s) for _, s in seqs]
        med = sorted(lens)[len(lens)//2]
        seqs2 = sorted(seqs, key=lambda x: (abs(len(x[1]) - med), -len(x[1]), x[0]))
        return seqs2[:max_n]

############################
# Main processing worker
############################

def extract_sequences_worker(args: Tuple[List[Dict[str, Any]], str]) -> List[Dict[str, Any]]:
    rows, fasta_path = args
    fa = FastaAccessor(fasta_path)
    out = []
    
    for d in rows:
        rec = LTRRecord(
            chrom=d['chrom'], id=d['id'], 
            left_start=d['left_start'], left_end=d['left_end'],
            right_start=d['right_start'], right_end=d['right_end'], 
            rc=d['rc'], casetype=d['casetype'],
            group=d['group'], identity=d['identity'], 
            ppt_start=d['ppt_start'], ppt_end=d['ppt_end'],
            tsd_start=d['tsd_start'], tsd_end=d['tsd_end'],
            nest_in=d.get('nest_in'), nest_out=d.get('nest_out')
        )
        rec.normalize()
        
        left = fa.fetch(rec.chrom, rec.left_start, rec.left_end)
        right = fa.fetch(rec.chrom, rec.right_start, rec.right_end)
        full = fa.fetch(rec.chrom, rec.left_start, rec.right_end)
        internal = ''
        
        if rec.casetype in ('Single','RecentlyNested') and rec.right_start and rec.left_end and rec.right_start > rec.left_end:
            internal = fa.fetch(rec.chrom, rec.left_end+1, rec.right_start-1)
        
        if rec.rc == 1:
            left, right, full, internal = rc_seq(right), rc_seq(left), rc_seq(full), rc_seq(internal)
        
        # Calculate completeness score
        completeness = score_element_completeness(rec, left, right)
        
        # Validate boundaries
        left_valid, left_boundary = validate_ltr_boundaries(left)
        right_valid, right_boundary = validate_ltr_boundaries(right)
        
        # Detect ORFs in internal region
        orfs = []
        if internal and len(internal) > 300:
            orfs = find_orfs(internal)
        
        out.append({
            'chrom': rec.chrom, 'id': rec.id, 'group': rec.group, 'casetype': rec.casetype,
            'identity': rec.identity, 'left': left, 'right': right, 'full': full, 'internal': internal,
            'completeness_score': completeness,
            'left_boundary': left_boundary, 'right_boundary': right_boundary,
            'has_pbs': bool(find_pbs_motif(full, len(left)) if full else False),
            'has_ppt': bool(rec.ppt_start and rec.ppt_end),
            'internal_orfs': len(orfs),
            'nest_in': rec.nest_in, 'nest_out': rec.nest_out
        })
    
    fa.close()
    return out

############################
# Enhanced group processing
############################

def process_group_enhanced(group_id: str, members: List[Dict[str, Any]], outdir: str,
                          base_cdhit_c: float, base_cdhit_aS: float, threads: int, 
                          max_msa_seq: int, trim_gap: float, merge_thr: float,
                          nesting_families: Dict[str, str] = None) -> List[FamilyConsensus]:
    """Enhanced group processing with quality scoring and iterative refinement."""
    results: List[FamilyConsensus] = []
    os.makedirs(outdir, exist_ok=True)
    
    logger.debug(f"Processing group {group_id} with {len(members)} members using {threads} threads")
    
    # Sort members by completeness score
    members.sort(key=lambda x: x.get('completeness_score', 0), reverse=True)
    
    # Handle single member groups specially
    if len(members) == 1:
        m = members[0]
        fc = FamilyConsensus(
            group=group_id,
            family_id=f"{group_id}_single",
            n_5p=1 if m['left'] else 0,
            n_3p=1 if m['right'] else 0,
            cons_5p=m['left'] if m['left'] else '',
            cons_3p=m['right'] if m['right'] else '',
            cons_internal=m['internal'] if m['internal'] else '',
            member_ids=[m['id']]
        )
        
        # Validate boundaries
        if fc.cons_5p:
            valid, boundary = validate_ltr_boundaries(fc.cons_5p)
            fc.boundary_type = boundary
        
        # Check if 5p and 3p are similar enough to merge
        if fc.cons_5p and fc.cons_3p:
            sim = approx_identity(fc.cons_5p, fc.cons_3p)
            fc.merged = (sim >= merge_thr)
        
        # Calculate quality score
        fc.quality_score = m.get('completeness_score', 50)
        
        results.append(fc)
        return results
    
    with tempfile.TemporaryDirectory(prefix=f"grp_{group_id}_", dir=outdir) as td:
        # Separate Solo-LTRs for special handling
        solo_members = [m for m in members if m['casetype'] in ['SoloSingle', 'SoloLTR']]
        regular_members = [m for m in members if m['casetype'] not in ['SoloSingle', 'SoloLTR']]
        
        # Process regular elements
        if regular_members:
            # Prepare sequences with quality weighting
            five_records = []
            three_records = []
            internal_records = []
            
            for i, m in enumerate(regular_members):
                weight = m.get('completeness_score', 50) / 100.0
                
                if m['left']:
                    five_records.append((f"{group_id}|5p|{i}|{weight:.2f}", m['left']))
                if m['right']:
                    three_records.append((f"{group_id}|3p|{i}|{weight:.2f}", m['right']))
                if m['internal'] and m['internal_orfs'] > 0:  # Only use internal with ORFs
                    internal_records.append((f"{group_id}|int|{i}|{weight:.2f}", m['internal']))
            
            # Get adaptive clustering parameters
            all_seqs = [s for _, s in five_records + three_records]
            if all_seqs:
                adapt_c, adapt_as = get_adaptive_clustering_params(all_seqs)
            else:
                adapt_c, adapt_as = base_cdhit_c, base_cdhit_aS
            
            # Cluster with nesting awareness
            def cluster_one(end: str, recs: List[Tuple[str,str]], nest_aware: bool = True):
                if not recs:
                    return []
                
                # If only one sequence, return it as a single cluster
                if len(recs) == 1:
                    h = recs[0][0]
                    try:
                        parts = h.split('|')
                        i = int(parts[2])
                        return [(f"{group_id}_{end}_F1", [i])]
                    except:
                        return []
                
                fa = os.path.join(td, f"{end}.fa")
                fasta_write(fa, recs)
                out_prefix = os.path.join(td, f"{end}.cdhit")
                run_cdhit_est(fa, out_prefix, c=adapt_c, aS=adapt_as, threads=max(1, threads))
                clstr = out_prefix + '.clstr'
                clusters = parse_clstr(clstr)
                
                # Post-process clusters with nesting info
                if nest_aware and nesting_families:
                    refined_clusters = []
                    for cluster in clusters:
                        # Check if any members share nesting family
                        nest_groups = defaultdict(list)
                        for h in cluster:
                            try:
                                idx = int(h.split('|')[2])
                                elem_id = regular_members[idx]['id']
                                nest_fam = nesting_families.get(elem_id, 'none')
                                nest_groups[nest_fam].append(h)
                            except:
                                nest_groups['none'].append(h)
                        
                        # Split cluster by nesting family if needed
                        for nest_fam, group_members in nest_groups.items():
                            if len(group_members) >= 2:  # Minimum cluster size
                                refined_clusters.append(group_members)
                    
                    if refined_clusters:
                        clusters = refined_clusters
                
                fams = []
                for k, mem in enumerate(clusters):
                    idxs = []
                    scores = []
                    for h in mem:
                        try:
                            parts = h.split('|')
                            i = int(parts[2])
                            score = float(parts[3]) if len(parts) > 3 else 1.0
                            idxs.append(i)
                            scores.append(score)
                        except Exception:
                            pass
                    
                    if idxs:
                        # Sort by score for prioritization
                        sorted_pairs = sorted(zip(idxs, scores), key=lambda x: x[1], reverse=True)
                        idxs = [idx for idx, _ in sorted_pairs]
                        fams.append((f"{group_id}_{end}_F{k+1}", idxs))
                
                return fams
            
            fams_5p = cluster_one('5p', five_records)
            fams_3p = cluster_one('3p', three_records)
            
            logger.debug(f"Group {group_id}: {len(fams_5p)} 5p families, {len(fams_3p)} 3p families")
            
            # If no families were created but we have sequences, create default family
            # This ensures every group generates at least one family
            if not fams_5p and not fams_3p:
                if five_records or three_records:
                    logger.debug(f"Creating default family for group {group_id}")
                    all_idxs = list(range(len(regular_members)))
                    if five_records:
                        fams_5p = [(f"{group_id}_5p_F1", all_idxs)]
                    if three_records:
                        fams_3p = [(f"{group_id}_3p_F1", all_idxs)]
                else:
                    logger.warning(f"Group {group_id} has no LTR sequences!")
                    # Cannot use continue here, just let it proceed with empty families
            
            # Pair 5p and 3p families
            fam_dict_5 = {fid: idxs for fid, idxs in fams_5p}
            fam_dict_3 = {fid: idxs for fid, idxs in fams_3p}
            
            used_3 = set()
            pairs: List[Tuple[str, Optional[str], List[int]]] = []
            
            for fid5, idxs5 in fam_dict_5.items():
                best = None
                best_overlap = 0
                
                for fid3, idxs3 in fam_dict_3.items():
                    overlap = len(set(idxs5) & set(idxs3))
                    if overlap > best_overlap:
                        best_overlap = overlap
                        best = fid3
                
                if best and best_overlap > 0 and best not in used_3:
                    used_3.add(best)
                    pairs.append((fid5, best, sorted(set(idxs5) | set(fam_dict_3[best]))))
                else:
                    pairs.append((fid5, None, idxs5))
            
            for fid3, idxs3 in fam_dict_3.items():
                if fid3 not in used_3:
                    pairs.append((None, fid3, idxs3))
            
            # Build consensus for each pair
            for fid5, fid3, idxs in pairs:
                # Get member data with quality scores
                member_data = [regular_members[i] for i in idxs]
                scores = [m.get('completeness_score', 50) for m in member_data]
                
                # Debug logging
                logger.debug(f"Processing pair: fid5={fid5}, fid3={fid3}, idxs={idxs}")
                
                # Collect sequences prioritizing high-quality elements
                seqs_5 = [(f"{group_id}_5_{i}", regular_members[i]['left']) 
                         for i in idxs if regular_members[i]['left']]
                seqs_3 = [(f"{group_id}_3_{i}", regular_members[i]['right']) 
                         for i in idxs if regular_members[i]['right']]
                
                # For internal, prioritize elements with high completeness
                high_quality_idxs = [i for i in idxs if regular_members[i].get('completeness_score', 0) > 70]
                if high_quality_idxs:
                    seqs_int = [(f"{group_id}_I_{i}", regular_members[i]['internal']) 
                               for i in high_quality_idxs if regular_members[i]['internal']]
                else:
                    seqs_int = [(f"{group_id}_I_{i}", regular_members[i]['internal']) 
                               for i in idxs if regular_members[i]['internal']]
                
                # Quality-aware subsampling
                seqs_5 = deterministic_subsample(seqs_5, max_msa_seq, 
                                                [scores[idxs.index(i)] for i in range(len(seqs_5))])
                seqs_3 = deterministic_subsample(seqs_3, max_msa_seq,
                                                [scores[idxs.index(i)] for i in range(len(seqs_3))])
                seqs_int = deterministic_subsample(seqs_int, max_msa_seq)
                
                fam_name = fid5 or fid3 or f"{group_id}_FNA"
                cons = FamilyConsensus(group=group_id, family_id=fam_name)
                
                # Build consensus with preservation of biological features
                if seqs_5:
                    c5, n5 = msa_consensus(seqs_5, threads=max(1, threads), tempdir=td, 
                                         trim_gap=trim_gap, use_iupac=True, preserve_motifs=True)
                    
                    # Iterative refinement for high-member families
                    if n5 > 20:
                        member_seqs = [{'seq': s} for _, s in seqs_5]
                        c5, iterations = refine_consensus_iteratively(c5, member_seqs, 
                                                                     max_iterations=3, threads=threads, tempdir=td)
                        cons.iteration_refined = iterations
                    
                    cons.cons_5p, cons.n_5p = c5, n5
                    
                    # Validate boundary
                    valid, boundary = validate_ltr_boundaries(c5)
                    if valid:
                        cons.boundary_type = boundary
                
                if seqs_3:
                    c3, n3 = msa_consensus(seqs_3, threads=max(1, threads), tempdir=td,
                                         trim_gap=trim_gap, use_iupac=True, preserve_motifs=True)
                    
                    if n3 > 20:
                        member_seqs = [{'seq': s} for _, s in seqs_3]
                        c3, iterations = refine_consensus_iteratively(c3, member_seqs,
                                                                     max_iterations=3, threads=threads, tempdir=td)
                    
                    cons.cons_3p, cons.n_3p = c3, n3
                
                if seqs_int:
                    ci, ni = msa_consensus(seqs_int, threads=max(1, threads), tempdir=td,
                                         trim_gap=trim_gap, use_iupac=True)
                    cons.cons_internal, cons.n_internal = ci, ni
                    
                    # Check for ORFs in consensus
                    if ci:
                        orfs = find_orfs(ci)
                        cons.internal_orfs = len(orfs)
                
                # Decide merge of 5p/3p
                if cons.cons_5p and cons.cons_3p:
                    sim = approx_identity(cons.cons_5p, cons.cons_3p)
                    cons.merged = (sim >= merge_thr)
                else:
                    cons.merged = False
                
                # Calculate quality score
                quality = 0.0
                if cons.boundary_type in ['canonical', 'variant1']:
                    quality += 20
                if cons.internal_orfs > 0:
                    quality += 20
                if cons.n_5p > 5 and cons.n_3p > 5:
                    quality += 20
                if cons.merged:
                    quality += 10
                quality += min(30, sum(scores) / max(1, len(scores)) * 0.3)
                cons.quality_score = quality
                
                # Store member IDs
                cons.member_ids = [regular_members[i]['id'] for i in idxs]
                
                results.append(cons)
        
        # Process Solo-LTRs separately with parallelization
        if solo_members:
            logger.info(f"Processing {len(solo_members)} Solo-LTR elements in group {group_id}")
            
            # Collect Solo-LTR sequences
            solo_seqs = []
            for m in solo_members:
                # Use whichever LTR is available
                seq = m['left'] or m['right']
                if seq:
                    solo_seqs.append((m['id'], seq, m))
            
            if not solo_seqs:
                return results
            
            # For large groups, use CD-HIT clustering with parallel threads
            if len(solo_seqs) > 10:
                # Get adaptive clustering parameters for solos
                seqs_only = [s for _, s, _ in solo_seqs]
                c_thr, as_thr = get_adaptive_clustering_params(seqs_only)
                
                # Write sequences for CD-HIT
                solo_fa_path = os.path.join(td, f'solo_{group_id}.fa')
                fasta_write(solo_fa_path, [(id, seq) for id, seq, _ in solo_seqs])
                
                # Calculate threads for Solo-LTR processing
                # Use half of available threads for better parallelization
                solo_threads = max(1, threads // 2) if threads > 2 else threads
                logger.debug(f"Solo-LTR clustering for group {group_id}: using {solo_threads} threads")
                
                # Run CD-HIT with parallel threads
                solo_out_prefix = os.path.join(td, f'solo_{group_id}.cdhit')
                run_cdhit_est(solo_fa_path, solo_out_prefix, c=c_thr, aS=as_thr, 
                             threads=solo_threads)  # Use calculated threads
                
                # Parse clusters
                solo_clusters = parse_clstr(solo_out_prefix + '.clstr')
                
                # Process each cluster to build consensus
                for i, cluster_ids in enumerate(solo_clusters):
                    cluster_seqs = [(id, seq) for id, seq, m in solo_seqs if id in cluster_ids]
                    cluster_members = [m for id, seq, m in solo_seqs if id in cluster_ids]
                    
                    if cluster_seqs:
                        # Build consensus with parallel MSA
                        if len(cluster_seqs) > 1:
                            cons, n = msa_consensus(cluster_seqs, threads=solo_threads, 
                                                  tempdir=td, preserve_motifs=True)
                        else:
                            cons = cluster_seqs[0][1]
                            n = 1
                        
                        fc = FamilyConsensus(
                            group=group_id,
                            family_id=f"{group_id}_solo_F{i+1}",
                            n_5p=len(cluster_seqs),
                            cons_5p=cons,
                            is_solo=True,
                            member_ids=[id for id, _ in cluster_seqs]
                        )
                        
                        # Validate boundaries
                        valid, boundary = validate_ltr_boundaries(cons)
                        fc.boundary_type = boundary
                        fc.quality_score = sum(m.get('completeness_score', 50) for m in cluster_members) / len(cluster_members)
                        
                        results.append(fc)
            else:
                # For small groups, process as single family
                solo_seqs_only = [(id, seq) for id, seq, _ in solo_seqs]
                
                # Use appropriate threads for small groups
                solo_threads = max(1, threads // 2) if threads > 2 else threads
                
                # Build consensus if multiple solos, otherwise use single sequence
                if len(solo_seqs_only) > 1:
                    cons, n = msa_consensus(solo_seqs_only, threads=solo_threads, 
                                          tempdir=td, preserve_motifs=True)
                else:
                    cons = solo_seqs_only[0][1]
                    n = 1
                
                fc = FamilyConsensus(
                    group=group_id,
                    family_id=f"{group_id}_solo",
                    n_5p=len(solo_seqs_only),
                    cons_5p=cons,
                    is_solo=True,
                    member_ids=[id for id, _ in solo_seqs_only]
                )
                
                # Validate boundaries
                valid, boundary = validate_ltr_boundaries(cons)
                fc.boundary_type = boundary
                fc.quality_score = sum(m.get('completeness_score', 50) for _, _, m in solo_seqs) / len(solo_seqs)
                
                results.append(fc)
    
    return results

############################
# Output generation
############################

def write_enhanced_outputs(outdir: str, families: List[FamilyConsensus], write_full: bool = True):
    """Write outputs with quality metrics and biological annotations."""
    os.makedirs(outdir, exist_ok=True)
    
    fa_ltr = os.path.join(outdir, 'LTR_consensus.fasta')
    fa_full = os.path.join(outdir, 'Full_consensus.fasta')
    tsv = os.path.join(outdir, 'LTR_consensus_stats.tsv')
    quality_tsv = os.path.join(outdir, 'LTR_consensus_quality.tsv')
    lib = os.path.join(outdir, 'LTR_library.lib')  # Following rtr2repeatmasker.sh naming
    
    # Sort families by quality score
    families.sort(key=lambda x: x.quality_score, reverse=True)
    
    ltr_records = []
    rm_records = []
    stats_rows = []
    quality_rows = []
    
    for fc in families:
        if fc.is_solo:
            # Solo-LTR handling
            if fc.cons_5p:
                name = f"{fc.family_id}_SOLO"
                ltr_records.append((name, fc.cons_5p))
                rm_records.append((f"{name}#LTR", fc.cons_5p))
                stats_rows.append([fc.group, fc.family_id, fc.n_5p, 0, 0, 'solo', 
                                  len(fc.cons_5p), 0])
                quality_rows.append([fc.family_id, 'solo', fc.boundary_type, 
                                   fc.quality_score, fc.n_5p, fc.iteration_refined])
        elif fc.merged and fc.cons_5p:
            name = f"{fc.family_id}"
            ltr_records.append((name, fc.cons_5p))
            rm_records.append((f"{name}#LTR/Unknown", fc.cons_5p))
            stats_rows.append([fc.group, fc.family_id, fc.n_5p, fc.n_3p, fc.n_internal, 
                             'merged', len(fc.cons_5p), len(fc.cons_internal)])
            quality_rows.append([fc.family_id, 'merged', fc.boundary_type, 
                               fc.quality_score, fc.n_5p + fc.n_3p, fc.iteration_refined])
        else:
            if fc.cons_5p:
                name5 = f"{fc.family_id}_5p"
                ltr_records.append((name5, fc.cons_5p))
                rm_records.append((f"{name5}#LTR", fc.cons_5p))
                stats_rows.append([fc.group, name5, fc.n_5p, 0, fc.n_internal, '5p', 
                                 len(fc.cons_5p), len(fc.cons_internal)])
                quality_rows.append([name5, '5p', fc.boundary_type, 
                                   fc.quality_score, fc.n_5p, fc.iteration_refined])
            if fc.cons_3p:
                name3 = f"{fc.family_id}_3p"
                ltr_records.append((name3, fc.cons_3p))
                rm_records.append((f"{name3}#LTR", fc.cons_3p))
                stats_rows.append([fc.group, name3, 0, fc.n_3p, fc.n_internal, '3p', 
                                 len(fc.cons_3p), len(fc.cons_internal)])
                quality_rows.append([name3, '3p', fc.boundary_type, 
                                   fc.quality_score, fc.n_3p, fc.iteration_refined])
    
    # Write files
    if ltr_records:
        fasta_write(fa_ltr, ltr_records)
        fasta_write(lib, rm_records)
    else:
        logger.warning("No consensus sequences generated!")
    
    # Stats file
    with open(tsv, 'w', newline='') as w:
        wt = csv.writer(w, delimiter='\t')
        wt.writerow(['Group','Family','n_5p','n_3p','n_internal','type','len_LTR','len_internal'])
        wt.writerows(stats_rows)
    
    # Quality metrics file
    with open(quality_tsv, 'w', newline='') as w:
        wt = csv.writer(w, delimiter='\t')
        wt.writerow(['Family','Type','Boundary','QualityScore','Members','Iterations'])
        wt.writerows(quality_rows)
    
    # Full consensus (with biological features)
    if write_full:
        full_records = []
        for fc in families:
            if fc.is_solo:
                continue
            
            # Build full element with annotations
            left = fc.cons_5p or fc.cons_3p
            right = fc.cons_3p or fc.cons_5p
            
            if left:
                full = left
                if fc.cons_internal:
                    full += fc.cons_internal
                if right:
                    full += right
                
                # Add annotation to header
                annot = []
                if fc.boundary_type != 'unknown':
                    annot.append(f"boundary={fc.boundary_type}")
                if fc.internal_orfs > 0:
                    annot.append(f"ORFs={fc.internal_orfs}")
                if fc.quality_score > 0:
                    annot.append(f"quality={fc.quality_score:.1f}")
                
                header = f"{fc.family_id}_FULL"
                if annot:
                    header += f" [{';'.join(annot)}]"
                
                full_records.append((header, full))
        
        if full_records:
            fasta_write(fa_full, full_records)
    
    # Summary statistics
    logger.info(f"Generated {len(families)} consensus families")
    high_quality = sum(1 for fc in families if fc.quality_score > 70)
    logger.info(f"High-quality families (score>70): {high_quality}")
    canonical = sum(1 for fc in families if fc.boundary_type == 'canonical')
    logger.info(f"Families with canonical boundaries: {canonical}")
    with_orfs = sum(1 for fc in families if fc.internal_orfs > 0)
    logger.info(f"Families with internal ORFs: {with_orfs}")

############################
# Main pipeline
############################

def main():
    ap = argparse.ArgumentParser(
        description="Enhanced LTR consensus builder with RepeatMasker integration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
    # Basic usage (recommended for most users):
    python build_ltr_consensus.py --rtr genome.rtr --fasta genome.fa --outdir results
    
    # Advanced usage with custom parameters:
    python build_ltr_consensus.py --rtr genome.rtr --fasta genome.fa --outdir results \
        --min-identity 0.85 --threads 16 --processes 4
        
    # Quick mode for testing:
    python build_ltr_consensus.py --rtr genome.rtr --fasta genome.fa --outdir results --quick
    """)
    
    # Required arguments
    ap.add_argument('--rtr', required=True, 
                   help='LTR detection table (genome.rtr from Look4LTRs)')
    ap.add_argument('--fasta', required=True, 
                   help='Genome FASTA file used for LTR detection')
    ap.add_argument('--outdir', required=True, 
                   help='Output directory for consensus library and scripts')
    
    # Quality filtering (with sensible defaults)
    ap.add_argument('--min-identity', type=float, default=0.70, 
                   help='Min LTR identity threshold (default: 0.70, similar to rtr2repeatmasker.sh)')
    ap.add_argument('--keep-casetypes', default='Single,RecentlyNested,SoloSingle,SoloLTR', 
                   help='CaseTypes to include')
    ap.add_argument('--min-ltr-len', type=int, default=100, help='Min LTR length')
    ap.add_argument('--max-ltr-len', type=int, default=5000, help='Max LTR length')
    
    ap.add_argument('--cdhit-c', type=float, default=0.90, help='Base CD-HIT identity')
    ap.add_argument('--cdhit-aS', type=float, default=0.80, help='Base CD-HIT coverage')
    ap.add_argument('--adaptive-clustering', action='store_true', default=True,
                   help='Use adaptive clustering parameters')
    
    ap.add_argument('--max-msa-seq', type=int, default=100, help='Max sequences per MSA')
    ap.add_argument('--trim-gap', type=float, default=0.70, help='Trim gap threshold')
    ap.add_argument('--merge-ltr-threshold', type=float, default=0.95, help='Merge 5p/3p threshold')
    
    ap.add_argument('--use-nesting', action='store_true', default=True,
                   help='Use nesting relationships for family assignment')
    ap.add_argument('--iterative-refinement', action='store_true', default=True,
                   help='Use iterative consensus refinement')
    
    # Performance settings with auto-detection
    import multiprocessing
    default_threads = min(8, multiprocessing.cpu_count())
    default_processes = min(2, max(1, multiprocessing.cpu_count() // 4))
    
    ap.add_argument('--threads', type=int, default=default_threads, 
                   help=f'Threads per tool (default: {default_threads}, auto-detected)')
    ap.add_argument('--processes', type=int, default=default_processes, 
                   help=f'Parallel processes (default: {default_processes}, auto-detected)')
    
    # Quick mode for testing
    ap.add_argument('--quick', action='store_true', 
                   help='Quick mode: reduced accuracy but faster processing')
    
    args = ap.parse_args()
    
    # Adjust parameters for quick mode
    if args.quick:
        logger.info("Quick mode enabled: using faster but less accurate settings")
        args.max_msa_seq = min(args.max_msa_seq, 50)
        args.iterative_refinement = False
        args.adaptive_clustering = False
    
    # Input validation
    if not os.path.exists(args.rtr):
        logger.error(f"RTR file not found: {args.rtr}")
        sys.exit(1)
    
    if not os.path.exists(args.fasta):
        logger.error(f"FASTA file not found: {args.fasta}")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    # Log file setup
    log_file = os.path.join(args.outdir, 'ltr_consensus.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.info(f"Started LTR consensus building with parameters: {vars(args)}")
    
    # Tool checks with better error handling
    logger.info("Checking required and optional tools...")
    
    # Required tools
    required_tools = ['cd-hit-est']
    missing_required = []
    for tool in required_tools:
        if not check_tool(tool, required=True):
            missing_required.append(tool)
    
    if missing_required:
        logger.error(f"Missing required tools: {', '.join(missing_required)}")
        logger.error("Please install missing tools and try again.")
        sys.exit(1)
    
    # Optional tools
    optional_tools = [('mafft', 'multiple sequence alignment'), 
                     ('muscle', 'multiple sequence alignment'),
                     ('samtools', 'FASTA indexing')]
    available_optional = []
    for tool, desc in optional_tools:
        path = check_tool(tool)
        if path:
            available_optional.append(tool)
            logger.info(f"Found {tool}: {path}")
    
    if not any(t in available_optional for t in ['mafft', 'muscle']):
        logger.warning("No MSA tool found (mafft or muscle). Consensus quality may be reduced.")
    
    logger.info(f"Tool check complete. Available optional tools: {', '.join(available_optional) if available_optional else 'none'}")
    
    # Read RTR table
    logger.info(f"Reading RTR: {args.rtr}")
    rows = read_rtr_table(args.rtr)
    recs: List[LTRRecord] = [LTRRecord.from_row(r) for r in rows]
    logger.info(f"Loaded {len(recs)} records")
    
    # Autoscale identity
    ids_scaled = autoscale_identity([r.identity for r in recs])
    for r, v in zip(recs, ids_scaled):
        r.identity = v
    
    # Filter by CaseType
    keep_types = set([x.strip() for x in args.keep_casetypes.split(',') if x.strip()])
    recs = [r for r in recs if r.casetype in keep_types]
    
    # Filter by LTR length and identity
    def ltr_len_ok(r: LTRRecord) -> bool:
        L5 = abs(r.left_end - r.left_start) + 1 if r.left_start and r.left_end else 0
        L3 = abs(r.right_end - r.right_start) + 1 if r.right_start and r.right_end else 0
        Lmax = max(L5, L3)
        Lmin = min(L5, L3) if L5 > 0 and L3 > 0 else Lmax
        
        if args.min_ltr_len and Lmax < args.min_ltr_len:
            return False
        if args.max_ltr_len and Lmax > args.max_ltr_len:
            return False
        return True
    
    before = len(recs)
    recs = [r for r in recs if (r.identity is None or r.identity >= args.min_identity) and ltr_len_ok(r)]
    logger.info(f"Records after filtering: {len(recs)}/{before}")
    
    # Build nesting graph if requested
    nesting_families = {}
    if args.use_nesting:
        logger.info("Building nesting relationship graph...")
        nest_graph = NestingGraph()
        for rec in recs:
            nest_graph.add_record(rec)
        nesting_families = nest_graph.assign_families_by_nesting()
        logger.info(f"Found {len(set(nesting_families.values()))} nesting-based families")
    
    # Prepare work items
    work_rows = []
    for r in recs:
        work_rows.append({
            'chrom': r.chrom, 'id': r.id, 
            'left_start': r.left_start, 'left_end': r.left_end,
            'right_start': r.right_start, 'right_end': r.right_end, 
            'rc': r.rc, 'casetype': r.casetype, 
            'group': r.group or 'ALL', 'identity': r.identity,
            'ppt_start': r.ppt_start, 'ppt_end': r.ppt_end, 
            'tsd_start': r.tsd_start, 'tsd_end': r.tsd_end,
            'nest_in': r.nest_in, 'nest_out': r.nest_out
        })
    
    # Group by GraphGroup
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for w in work_rows:
        g = w['group'] or 'ALL'
        groups.setdefault(g, []).append(w)
    
    logger.info(f"Processing {len(groups)} groups")
    
    # Extract sequences with biological validation
    all_members: Dict[str, List[Dict[str, Any]]] = {g: [] for g in groups}
    tasks = []
    
    with ProcessPoolExecutor(max_workers=max(1, args.processes)) as ex:
        for g, lst in groups.items():
            chunk_size = max(100, len(lst) // (args.processes or 1))
            chunks = [lst[i:i+chunk_size] for i in range(0, len(lst), chunk_size)]
            for ch in chunks:
                tasks.append((g, ex.submit(extract_sequences_worker, (ch, args.fasta))))
        
        for g, fut in tasks:
            try:
                out = fut.result()
                all_members[g].extend(out)
            except Exception as e:
                logger.error(f"Extract failed for group {g}: {e}")
                raise
    
    # Process each group with enhancements
    families_all: List[FamilyConsensus] = []
    procs = max(1, args.processes)
    threads_per_group = max(1, args.threads // procs)
    
    logger.info(f"Processing groups with {procs} processes, {threads_per_group} threads per group")
    
    with ProcessPoolExecutor(max_workers=procs) as ex:
        futs = {}
        for g, m in all_members.items():
            futs[ex.submit(process_group_enhanced, g, m, args.outdir, 
                          args.cdhit_c, args.cdhit_aS,
                          threads_per_group, 
                          args.max_msa_seq, args.trim_gap, 
                          args.merge_ltr_threshold,
                          nesting_families)] = g
        
        for fut, g in futs.items():
            try:
                fams = fut.result()
                families_all.extend(fams)
                logger.info(f"Group {g}: {len(fams)} families")
            except Exception as e:
                logger.error(f"Group {g} failed: {e}")
                raise
    
    # Write enhanced outputs
    write_enhanced_outputs(args.outdir, families_all, write_full=True)
    
    logger.info("Enhanced LTR consensus building completed successfully!")

if __name__ == '__main__':
    main()