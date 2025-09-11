#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simplified LTR consensus builder focused on RepeatMasker compatibility.

This simplified version removes over-engineering and focuses on core functionality:
- Fast and accurate consensus building
- RepeatMasker-compatible output
- Minimal dependencies
- Clear, maintainable code

USAGE:
    python build_ltr_consensus_simplified.py \
        --rtr genome.rtr \
        --fasta genome.fa \
        --outdir out_ltr_consensus

Dependencies: cd-hit-est, mafft (optional)
"""

import os
import sys
import csv
import re
import math
import shutil
import tempfile
import logging
import argparse
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import subprocess

# Setup logging
logger = logging.getLogger("ltr_consensus_simple")
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

# Simple LTR boundary patterns (canonical only)
LTR_BOUNDARIES = {
    'TG': 'CA',  # Canonical
    'CA': 'CA',  # Variant
    'TG': 'TA',  # Common variant
}

@dataclass
class LTRRecord:
    """Simple LTR record with essential fields only."""
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
                return float(x)
            except (ValueError, TypeError):
                return None
        
        return LTRRecord(
            chrom=s(r.get('chrom') or r.get('Chr') or r.get('seqid')),
            id=str(r.get('ID') or r.get('id') or ''),
            left_start=num(r.get('LeftStart')) or 0,
            left_end=num(r.get('LeftEnd')) or 0,
            right_start=num(r.get('RightStart')) or 0,
            right_end=num(r.get('RightEnd')) or 0,
            rc=num(r.get('RC')) or 0,
            casetype=s(r.get('CaseType') or r.get('type') or 'Single'),
            group=str(r.get('GraphGroup') or r.get('group') or 'ALL'),
            identity=f0_1(r.get('LTRIdentity')),
        )
    
    def normalize(self):
        """Normalize coordinates."""
        if self.left_start > self.left_end:
            self.left_start, self.left_end = self.left_end, self.left_start
        if self.right_start > self.right_end:
            self.right_start, self.right_end = self.right_end, self.right_start

@dataclass
class FamilyConsensus:
    """Simple family consensus."""
    group: str
    family_id: str
    n_5p: int = 0
    n_3p: int = 0
    cons_5p: str = ''
    cons_3p: str = ''
    merged: bool = False
    quality_score: float = 0.0
    member_ids: List[str] = None
    
    def __post_init__(self):
        if self.member_ids is None:
            self.member_ids = []

class FastaAccessor:
    """Simple FASTA accessor with minimal dependencies."""
    
    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.sequences = {}
        self._load_sequences()
    
    def _load_sequences(self):
        """Load sequences into memory for fast access."""
        try:
            # Try pyfaidx first
            import pyfaidx
            self.impl = pyfaidx.Fasta(self.fasta_path, as_raw=True, sequence_always_upper=True)
            self.mode = 'pyfaidx'
            logger.info("Using pyfaidx for FASTA access")
        except ImportError:
            # Fallback to Bio.SeqIO
            try:
                from Bio import SeqIO
                self.sequences = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(self.fasta_path, "fasta")}
                self.mode = 'seqio'
                logger.info("Using Bio.SeqIO for FASTA access")
            except ImportError:
                # Manual parsing as last resort
                self._manual_parse()
                self.mode = 'manual'
                logger.info("Using manual parsing for FASTA access")
    
    def _manual_parse(self):
        """Manual FASTA parsing."""
        with open(self.fasta_path) as f:
            seq_id = None
            seq_parts = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if seq_id:
                        self.sequences[seq_id] = ''.join(seq_parts).upper()
                    seq_id = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if seq_id:
                self.sequences[seq_id] = ''.join(seq_parts).upper()
    
    def fetch(self, chrom: str, start1: int, end1: int) -> str:
        """Fetch sequence using 1-based coordinates."""
        if start1 is None or end1 is None or start1 <= 0 or end1 < start1:
            return ""
        
        if self.mode == 'pyfaidx':
            return str(self.impl[chrom][start1-1:end1])
        else:
            seq = self.sequences.get(chrom, '')
            if seq:
                return seq[start1-1:end1]
            return ""

def check_tool(tool_name: str) -> Optional[str]:
    """Check if tool is available."""
    return shutil.which(tool_name)

def rc_seq(seq: str) -> str:
    """Reverse complement."""
    comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(comp)[::-1]

def validate_ltr_boundaries(seq: str) -> bool:
    """Simple LTR boundary validation."""
    if not seq or len(seq) < 4:
        return False
    
    seq_upper = seq.upper()
    start_di = seq_upper[:2]
    end_di = seq_upper[-2:]
    
    # Check canonical boundaries
    return (start_di == 'TG' and end_di == 'CA') or \
           (start_di == 'CA' and end_di == 'CA') or \
           (start_di == 'TG' and end_di == 'TA')

def simple_consensus(sequences: List[str], min_frac_major: float = 0.70) -> str:
    """Build consensus using majority rule with stricter thresholds."""
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
        
        # Find most frequent base
        most_freq_base = max(bases.items(), key=lambda x: x[1])
        if most_freq_base[1] / total >= min_frac_major:
            consensus.append(most_freq_base[0])
        else:
            consensus.append('N')
    
    return ''.join(consensus).rstrip('N')

def run_mafft(in_fa: str, out_fa: str) -> bool:
    """Run MAFFT alignment."""
    exe = check_tool('mafft')
    if not exe:
        return False
    
    if not os.path.exists(in_fa) or os.path.getsize(in_fa) == 0:
        return False
    
    cmd = [exe, '--auto', '--quiet', in_fa]
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, check=True, timeout=300)
        with open(out_fa, 'w') as f:
            f.write(result.stdout)
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False

def msa_consensus(sequences: List[Tuple[str, str]], tempdir: str) -> Tuple[str, int]:
    """MSA-based consensus with fallback to simple consensus."""
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
    success = run_mafft(in_fa, aln_fa)
    
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
            consensus = simple_consensus(aligned_seqs, min_frac_major=0.70)
            return (consensus, len(sequences))
    
    # Fallback to simple consensus without alignment
    seqs_only = [seq for _, seq in sequences]
    consensus = simple_consensus(seqs_only, min_frac_major=0.60)  # Lower threshold for unaligned
    return (consensus, len(sequences))

def run_cdhit_est(fa_path: str, out_prefix: str, c: float = 0.85, threads: int = 1) -> str:
    """Run CD-HIT-EST clustering."""
    exe = check_tool('cd-hit-est')
    if not exe:
        raise RuntimeError("cd-hit-est not found")
    
    cmd = [exe, '-i', fa_path, '-o', out_prefix, '-c', str(c), '-aS', '0.80',
           '-T', str(max(1, threads)), '-M', '0', '-d', '0']
    
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        logger.error(f"CD-HIT-EST timed out for {fa_path}")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"CD-HIT-EST failed for {fa_path}")
        raise
    
    return out_prefix + '.clstr'

def parse_cdhit_clusters(clstr_path: str) -> List[List[str]]:
    """Parse CD-HIT cluster file."""
    clusters = []
    current_cluster = []
    
    with open(clstr_path) as f:
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

def process_group(group_id: str, members: List[Dict], outdir: str, threads: int = 1) -> List[FamilyConsensus]:
    """Process a group of LTR elements."""
    results = []
    
    if not members:
        return results
    
    # Handle single member
    if len(members) == 1:
        m = members[0]
        fc = FamilyConsensus(
            group=group_id,
            family_id=f"{group_id}_single",
            n_5p=1 if m['left'] else 0,
            n_3p=1 if m['right'] else 0,
            cons_5p=m['left'] if m['left'] else '',
            cons_3p=m['right'] if m['right'] else '',
            member_ids=[m['id']]
        )
        
        # Simple quality score
        fc.quality_score = 50 if fc.cons_5p or fc.cons_3p else 0
        if fc.cons_5p and fc.cons_3p:
            # Check if 5p and 3p should be merged
            sim = sum(1 for a, b in zip(fc.cons_5p, fc.cons_3p) if a == b) / max(len(fc.cons_5p), len(fc.cons_3p))
            fc.merged = (sim >= 0.90)
        
        results.append(fc)
        return results
    
    with tempfile.TemporaryDirectory(prefix=f"grp_{group_id}_", dir=outdir) as td:
        # Separate sequences by type
        five_records = [(f"{group_id}_5p_{i}", m['left']) for i, m in enumerate(members) if m['left']]
        three_records = [(f"{group_id}_3p_{i}", m['right']) for i, m in enumerate(members) if m['right']]
        
        def cluster_and_consensus(records: List[Tuple[str, str]], end_name: str) -> List[Tuple[str, str, List[str]]]:
            if not records:
                return []
            
            if len(records) == 1:
                return [(f"{group_id}_{end_name}_F1", records[0][1], [records[0][0]])]
            
            # Write sequences for clustering
            fa_path = os.path.join(td, f'{end_name}.fa')
            with open(fa_path, 'w') as f:
                for header, seq in records:
                    f.write(f'>{header}\n{seq}\n')
            
            # Cluster with CD-HIT
            out_prefix = os.path.join(td, f'{end_name}.cdhit')
            clstr_path = run_cdhit_est(fa_path, out_prefix, c=0.85, threads=threads)
            clusters = parse_cdhit_clusters(clstr_path)
            
            # Build consensus for each cluster
            consensuses = []
            seq_dict = dict(records)
            
            for k, cluster_ids in enumerate(clusters):
                if len(cluster_ids) < 2:  # Skip single-member clusters
                    continue
                
                cluster_seqs = [(cid, seq_dict[cid]) for cid in cluster_ids if cid in seq_dict]
                
                if cluster_seqs:
                    consensus, n = msa_consensus(cluster_seqs, td)
                    if consensus:
                        family_id = f"{group_id}_{end_name}_F{k+1}"
                        consensuses.append((family_id, consensus, cluster_ids))
            
            return consensuses
        
        # Process 5p and 3p separately
        five_families = cluster_and_consensus(five_records, '5p')
        three_families = cluster_and_consensus(three_records, '3p')
        
        # Create family objects
        for family_id, consensus, member_ids in five_families:
            fc = FamilyConsensus(
                group=group_id,
                family_id=family_id,
                n_5p=len(member_ids),
                cons_5p=consensus,
                member_ids=member_ids
            )
            fc.quality_score = 60 + min(40, len(member_ids) * 2)  # Simple scoring
            results.append(fc)
        
        for family_id, consensus, member_ids in three_families:
            fc = FamilyConsensus(
                group=group_id,
                family_id=family_id,
                n_3p=len(member_ids),
                cons_3p=consensus,
                member_ids=member_ids
            )
            fc.quality_score = 60 + min(40, len(member_ids) * 2)
            results.append(fc)
    
    return results

def extract_sequences(records: List[LTRRecord], fasta_path: str) -> Dict[str, List[Dict]]:
    """Extract sequences for all records."""
    fa = FastaAccessor(fasta_path)
    groups = defaultdict(list)
    
    for rec in records:
        rec.normalize()
        
        left = fa.fetch(rec.chrom, rec.left_start, rec.left_end)
        right = fa.fetch(rec.chrom, rec.right_start, rec.right_end)
        
        if rec.rc == 1:
            left, right = rc_seq(right), rc_seq(left)
        
        groups[rec.group].append({
            'id': rec.id,
            'left': left,
            'right': right,
            'casetype': rec.casetype
        })
    
    return dict(groups)

def write_outputs(outdir: str, families: List[FamilyConsensus]):
    """Write simplified outputs."""
    os.makedirs(outdir, exist_ok=True)
    
    # Sort by quality
    families.sort(key=lambda x: x.quality_score, reverse=True)
    
    # LTR consensus file
    ltr_file = os.path.join(outdir, 'LTR_consensus.fasta')
    lib_file = os.path.join(outdir, 'LTR_library.lib')
    stats_file = os.path.join(outdir, 'LTR_stats.tsv')
    
    ltr_records = []
    lib_records = []
    stats_rows = []
    
    for fc in families:
        if fc.merged and fc.cons_5p:
            # Merged LTR
            name = f"{fc.family_id}"
            ltr_records.append((name, fc.cons_5p))
            lib_records.append((f"{name}#LTR/Unknown", fc.cons_5p))
            stats_rows.append([fc.group, fc.family_id, fc.n_5p + fc.n_3p, 'merged', len(fc.cons_5p)])
        else:
            # Separate 5p and 3p
            if fc.cons_5p:
                name5 = f"{fc.family_id}_5p"
                ltr_records.append((name5, fc.cons_5p))
                lib_records.append((f"{name5}#LTR", fc.cons_5p))
                stats_rows.append([fc.group, name5, fc.n_5p, '5p', len(fc.cons_5p)])
            
            if fc.cons_3p:
                name3 = f"{fc.family_id}_3p"
                ltr_records.append((name3, fc.cons_3p))
                lib_records.append((f"{name3}#LTR", fc.cons_3p))
                stats_rows.append([fc.group, name3, fc.n_3p, '3p', len(fc.cons_3p)])
    
    # Write FASTA files
    def write_fasta(path: str, records: List[Tuple[str, str]]):
        with open(path, 'w') as f:
            for header, seq in records:
                f.write(f'>{header}\n')
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')
    
    if ltr_records:
        write_fasta(ltr_file, ltr_records)
        write_fasta(lib_file, lib_records)
    
    # Write stats
    with open(stats_file, 'w') as f:
        f.write('Group\tFamily\tMembers\tType\tLength\n')
        for row in stats_rows:
            f.write('\t'.join(map(str, row)) + '\n')
    
    logger.info(f"Generated {len(families)} consensus families")
    logger.info(f"Output files: {ltr_file}, {lib_file}, {stats_file}")

def read_rtr_table(path: str) -> List[Dict[str, str]]:
    """Read RTR table."""
    with open(path, 'r', newline='') as f:
        # Auto-detect delimiter
        sample = f.read(1024)
        f.seek(0)
        
        if '\t' in sample:
            delimiter = '\t'
        elif ',' in sample:
            delimiter = ','
        else:
            delimiter = '\t'
        
        reader = csv.DictReader(f, delimiter=delimiter)
        return [dict(row) for row in reader]

def autoscale_identity(id_list: List[Optional[float]]) -> List[Optional[float]]:
    """Auto-scale identity values from 0-100 to 0-1 if needed."""
    vals = [v for v in id_list if v is not None and not math.isnan(v)]
    if not vals:
        return id_list
    
    max_val = max(vals)
    if max_val > 1.5:  # Likely 0-100 scale
        logger.info("Scaling identity values from 0-100 to 0-1")
        return [None if v is None else v/100.0 for v in id_list]
    
    return id_list

def main():
    parser = argparse.ArgumentParser(
        description="Simplified LTR consensus builder",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
    # Basic usage:
    python build_ltr_consensus_simplified.py --rtr genome.rtr --fasta genome.fa --outdir results
    
    # With custom parameters:
    python build_ltr_consensus_simplified.py --rtr genome.rtr --fasta genome.fa --outdir results \\
        --min-identity 0.80 --threads 8
        """)
    
    # Required arguments
    parser.add_argument('--rtr', required=True, help='LTR detection table (genome.rtr)')
    parser.add_argument('--fasta', required=True, help='Genome FASTA file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    
    # Optional parameters
    parser.add_argument('--min-identity', type=float, default=0.70, help='Min LTR identity (default: 0.70)')
    parser.add_argument('--min-ltr-len', type=int, default=100, help='Min LTR length (default: 100)')
    parser.add_argument('--max-ltr-len', type=int, default=5000, help='Max LTR length (default: 5000)')
    parser.add_argument('--keep-casetypes', default='Single,RecentlyNested,SoloSingle,SoloLTR',
                       help='CaseTypes to include')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.rtr):
        logger.error(f"RTR file not found: {args.rtr}")
        sys.exit(1)
    
    if not os.path.exists(args.fasta):
        logger.error(f"FASTA file not found: {args.fasta}")
        sys.exit(1)
    
    # Check required tools
    if not check_tool('cd-hit-est'):
        logger.error("cd-hit-est not found. Please install CD-HIT.")
        sys.exit(1)
    
    if not check_tool('mafft'):
        logger.warning("mafft not found. Will use simple consensus method.")
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    # Setup logging to file
    log_file = os.path.join(args.outdir, 'ltr_consensus_simple.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    logger.info(f"Starting LTR consensus building with simplified approach")
    logger.info(f"Parameters: {vars(args)}")
    
    # Read and filter records
    logger.info(f"Reading RTR table: {args.rtr}")
    rows = read_rtr_table(args.rtr)
    records = [LTRRecord.from_row(row) for row in rows]
    logger.info(f"Loaded {len(records)} records")
    
    # Auto-scale identity
    identities = autoscale_identity([r.identity for r in records])
    for r, identity in zip(records, identities):
        r.identity = identity
    
    # Filter by casetype
    keep_types = set(args.keep_casetypes.split(','))
    records = [r for r in records if r.casetype in keep_types]
    
    # Filter by length and identity
    def is_valid_ltr(r: LTRRecord) -> bool:
        # Check identity
        if r.identity is not None and r.identity < args.min_identity:
            return False
        
        # Check LTR lengths
        left_len = abs(r.left_end - r.left_start) + 1 if r.left_start and r.left_end else 0
        right_len = abs(r.right_end - r.right_start) + 1 if r.right_start and r.right_end else 0
        max_len = max(left_len, right_len)
        
        if max_len < args.min_ltr_len or max_len > args.max_ltr_len:
            return False
        
        return True
    
    before_filter = len(records)
    records = [r for r in records if is_valid_ltr(r)]
    logger.info(f"Records after filtering: {len(records)}/{before_filter}")
    
    if not records:
        logger.error("No records remaining after filtering")
        sys.exit(1)
    
    # Extract sequences
    logger.info("Extracting sequences...")
    groups = extract_sequences(records, args.fasta)
    logger.info(f"Found {len(groups)} groups")
    
    # Process groups
    all_families = []
    with ProcessPoolExecutor(max_workers=min(4, len(groups))) as executor:
        futures = {}
        for group_id, members in groups.items():
            future = executor.submit(process_group, group_id, members, args.outdir, args.threads)
            futures[future] = group_id
        
        for future in futures:
            group_id = futures[future]
            try:
                families = future.result()
                all_families.extend(families)
                logger.info(f"Group {group_id}: {len(families)} families")
            except Exception as e:
                logger.error(f"Group {group_id} failed: {e}")
                raise
    
    # Write outputs
    write_outputs(args.outdir, all_families)
    logger.info("LTR consensus building completed successfully!")

if __name__ == '__main__':
    main()