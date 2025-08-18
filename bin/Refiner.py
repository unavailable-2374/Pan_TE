#!/usr/bin/env python3
"""
Enhanced Refiner.py based on bin/Refiner/phase2 logic
Replaces the original Refiner.py with advanced consensus building capabilities
"""

import os
import sys
import logging
import argparse
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict, Counter
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)

class EnhancedTERefiner:
    """Enhanced TE Consensus Builder based on phase2 logic"""
    
    def __init__(self, min_score=150, gap_init=20, gap_ext=5, threads=4):
        self.min_score = min_score
        self.gap_init = gap_init
        self.gap_ext = gap_ext
        self.threads = threads or 4
        
        # Try to find required tools
        self.mafft_path = self._find_tool('mafft')
        if not self.mafft_path:
            logger.warning("MAFFT not found in PATH, will use basic consensus")
        
        logger.info(f"Enhanced TE Refiner initialized: min_score={min_score}, threads={threads}")
    
    def _find_tool(self, tool_name: str) -> Optional[str]:
        """Find tool in PATH"""
        import shutil
        return shutil.which(tool_name)
    
    def build_consensus(self, input_file: str, output_file: str) -> bool:
        """Main consensus building function"""
        try:
            # Read input sequences
            sequences = list(SeqIO.parse(input_file, "fasta"))
            if not sequences:
                logger.error(f"No sequences found in {input_file}")
                return False
            
            logger.info(f"Processing {len(sequences)} sequences from {input_file}")
            
            # Convert to internal format
            copies = []
            for i, seq in enumerate(sequences):
                copies.append({
                    'id': seq.id,
                    'sequence': str(seq.seq),
                    'length': len(seq.seq),
                    'description': seq.description
                })
            
            # Build consensus using phase2 logic
            consensus_result = self._build_enhanced_consensus(copies)
            
            if not consensus_result:
                logger.error("Failed to build consensus")
                return False
            
            # Write output
            self._write_consensus(consensus_result, output_file)
            logger.info(f"Consensus written to {output_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error building consensus: {e}")
            return False
    
    def _build_enhanced_consensus(self, copies: List[Dict]) -> Optional[Dict]:
        """Build consensus using enhanced logic from phase2"""
        
        if not copies:
            return None
        
        if len(copies) == 1:
            # Single sequence - use as is but try to improve
            return self._create_single_consensus(copies[0])
        
        # Analyze sequence characteristics
        characteristics = self._analyze_sequence_characteristics(copies)
        logger.info(f"Sequence characteristics: {characteristics}")
        
        # Check if we should split into subfamilies
        if self._should_split_subfamilies(characteristics, len(copies)):
            logger.info("Splitting into subfamilies")
            subfamilies = self._identify_subfamilies(copies, characteristics)
            
            # Build consensus for largest subfamily
            if subfamilies:
                largest_subfamily = max(subfamilies, key=len)
                logger.info(f"Using largest subfamily with {len(largest_subfamily)} sequences")
                return self._build_msa_consensus(largest_subfamily, characteristics)
        
        # Build single consensus from all sequences
        return self._build_msa_consensus(copies, characteristics)
    
    def _analyze_sequence_characteristics(self, copies: List[Dict]) -> Dict:
        """Analyze sequence characteristics"""
        lengths = [copy['length'] for copy in copies]
        sequences = [copy['sequence'] for copy in copies]
        
        # Calculate pairwise identities (simplified)
        identities = []
        if len(sequences) > 1:
            # Sample pairs for efficiency
            max_pairs = min(100, len(sequences) * (len(sequences) - 1) // 2)
            pairs_sampled = 0
            
            for i in range(len(sequences)):
                for j in range(i + 1, len(sequences)):
                    if pairs_sampled >= max_pairs:
                        break
                    identity = self._calculate_identity(sequences[i], sequences[j])
                    identities.append(identity)
                    pairs_sampled += 1
                if pairs_sampled >= max_pairs:
                    break
        
        characteristics = {
            'count': len(copies),
            'avg_length': np.mean(lengths),
            'length_std': np.std(lengths),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'avg_identity': np.mean(identities) if identities else 100.0,
            'identity_std': np.std(identities) if identities else 0.0
        }
        
        return characteristics
    
    def _calculate_identity(self, seq1: str, seq2: str) -> float:
        """Calculate simple sequence identity"""
        if len(seq1) == 0 or len(seq2) == 0:
            return 0.0
        
        # Align sequences (simple approach)
        min_len = min(len(seq1), len(seq2))
        matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
        return (matches / min_len) * 100.0
    
    def _should_split_subfamilies(self, characteristics: Dict, num_copies: int) -> bool:
        """Decide if sequences should be split into subfamilies"""
        # Split if:
        # 1. Many sequences with high diversity
        # 2. Large length variation
        if num_copies < 4:
            return False
        
        if (num_copies >= 10 and 
            characteristics['avg_identity'] < 85.0 and 
            characteristics['length_std'] / characteristics['avg_length'] > 0.2):
            return True
        
        return False
    
    def _identify_subfamilies(self, copies: List[Dict], characteristics: Dict) -> List[List[Dict]]:
        """Identify subfamilies using clustering"""
        
        if len(copies) < 4:
            return [copies]
        
        try:
            # Build distance matrix
            n = len(copies)
            distances = np.zeros((n, n))
            
            for i in range(n):
                for j in range(i + 1, n):
                    identity = self._calculate_identity(copies[i]['sequence'], copies[j]['sequence'])
                    distance = 100.0 - identity
                    distances[i, j] = distance
                    distances[j, i] = distance
            
            # Simple clustering - group sequences within distance threshold
            threshold = 15.0  # 85% identity threshold
            visited = [False] * n
            subfamilies = []
            
            for i in range(n):
                if visited[i]:
                    continue
                
                subfamily = [copies[i]]
                visited[i] = True
                
                for j in range(i + 1, n):
                    if not visited[j] and distances[i, j] <= threshold:
                        subfamily.append(copies[j])
                        visited[j] = True
                
                if len(subfamily) >= 2:  # Only keep subfamilies with multiple members
                    subfamilies.append(subfamily)
            
            # If no good subfamilies, return all as one
            if not subfamilies:
                return [copies]
            
            logger.info(f"Identified {len(subfamilies)} subfamilies")
            return subfamilies
            
        except Exception as e:
            logger.warning(f"Subfamily identification failed: {e}, using all sequences")
            return [copies]
    
    def _build_msa_consensus(self, copies: List[Dict], characteristics: Dict) -> Optional[Dict]:
        """Build consensus using multiple sequence alignment"""
        
        if not copies:
            return None
        
        if len(copies) == 1:
            return self._create_single_consensus(copies[0])
        
        # Try MAFFT if available
        if self.mafft_path:
            consensus_seq = self._run_mafft_consensus(copies, characteristics)
            if consensus_seq:
                return self._create_consensus_record(consensus_seq, copies, "mafft")
        
        # Fallback to simple consensus
        logger.info("Using simple consensus method")
        consensus_seq = self._build_simple_consensus(copies)
        if consensus_seq:
            return self._create_consensus_record(consensus_seq, copies, "simple")
        
        # Last resort - use longest sequence
        longest = max(copies, key=lambda x: x['length'])
        return self._create_single_consensus(longest)
    
    def _run_mafft_consensus(self, copies: List[Dict], characteristics: Dict) -> Optional[str]:
        """Run MAFFT and build consensus from MSA"""
        
        try:
            # Create temporary input file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_input:
                for i, copy in enumerate(copies):
                    tmp_input.write(f">seq_{i}\n{copy['sequence']}\n")
                tmp_input_path = tmp_input.name
            
            # Choose MAFFT algorithm based on characteristics
            if len(copies) > 100:
                algorithm = ['--auto']
            elif characteristics.get('avg_identity', 70) > 80:
                algorithm = ['--localpair', '--maxiterate', '1000']
            else:
                algorithm = ['--genafpair', '--maxiterate', '1000']
            
            # Run MAFFT
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_output:
                tmp_output_path = tmp_output.name
            
            cmd = [self.mafft_path] + algorithm + ['--quiet', '--thread', str(self.threads), tmp_input_path]
            
            logger.info(f"Running MAFFT: {' '.join(cmd)}")
            
            with open(tmp_output_path, 'w') as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
            
            if result.returncode != 0:
                logger.warning(f"MAFFT failed: {result.stderr}")
                return None
            
            # Build consensus from MSA
            consensus_seq = self._consensus_from_msa(tmp_output_path)
            
            # Cleanup
            os.unlink(tmp_input_path)
            os.unlink(tmp_output_path)
            
            return consensus_seq
            
        except Exception as e:
            logger.warning(f"MAFFT consensus failed: {e}")
            return None
    
    def _consensus_from_msa(self, msa_file: str) -> Optional[str]:
        """Build consensus from MSA file"""
        
        try:
            aligned_seqs = []
            for record in SeqIO.parse(msa_file, "fasta"):
                aligned_seqs.append(str(record.seq))
            
            if not aligned_seqs:
                return None
            
            # Build consensus position by position
            consensus = []
            msa_length = len(aligned_seqs[0])
            
            for pos in range(msa_length):
                column = [seq[pos] for seq in aligned_seqs if pos < len(seq)]
                column_clean = [c for c in column if c not in '-']
                
                if not column_clean:
                    continue  # Skip gap-only columns
                
                # Get most common nucleotide
                counter = Counter(column_clean)
                most_common = counter.most_common(1)[0]
                
                # Use majority rule with minimum coverage
                if most_common[1] >= len(column_clean) * 0.3:  # At least 30% support
                    consensus.append(most_common[0])
            
            consensus_seq = ''.join(consensus)
            
            # Remove excessive gaps and clean up
            consensus_seq = consensus_seq.replace('-', '')
            
            return consensus_seq if len(consensus_seq) >= 50 else None
            
        except Exception as e:
            logger.warning(f"Error building consensus from MSA: {e}")
            return None
    
    def _build_simple_consensus(self, copies: List[Dict]) -> Optional[str]:
        """Build simple consensus without MSA"""
        
        if not copies:
            return None
        
        # For simple consensus, just use the longest sequence
        # In a more sophisticated implementation, we could do position-wise consensus
        longest = max(copies, key=lambda x: x['length'])
        return longest['sequence']
    
    def _create_single_consensus(self, copy: Dict) -> Dict:
        """Create consensus record from single sequence"""
        return {
            'id': copy['id'] + '_consensus',
            'sequence': copy['sequence'],
            'source_id': copy['id'],
            'copy_number': 1,
            'length': copy['length'],
            'method': 'single'
        }
    
    def _create_consensus_record(self, consensus_seq: str, copies: List[Dict], method: str) -> Dict:
        """Create consensus record"""
        
        # Use first sequence ID as base
        base_id = copies[0]['id']
        if '_' in base_id:
            base_id = base_id.split('_')[0]
        
        return {
            'id': f"{base_id}_consensus",
            'sequence': consensus_seq,
            'source_id': base_id,
            'copy_number': len(copies),
            'length': len(consensus_seq),
            'method': method,
            'original_length': np.mean([c['length'] for c in copies]),
            'improvement_ratio': len(consensus_seq) / np.mean([c['length'] for c in copies])
        }
    
    def _write_consensus(self, consensus_result: Dict, output_file: str):
        """Write consensus to output file"""
        
        record = SeqRecord(
            Seq(consensus_result['sequence']),
            id=consensus_result['id'],
            description=f"consensus from {consensus_result['copy_number']} copies, method={consensus_result['method']}"
        )
        
        with open(output_file, 'w') as handle:
            SeqIO.write([record], handle, "fasta")


def main():
    """Command line interface - maintains compatibility with original Refiner.py"""
    
    parser = argparse.ArgumentParser(
        description='Enhanced TE Consensus Builder based on advanced phase2 logic'
    )
    
    parser.add_argument('input', help='Input FASTA file with family sequences')
    parser.add_argument('output', help='Output FASTA file for consensus')
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='Number of threads (default: 4)')
    parser.add_argument('--min-score', type=float, default=150,
                       help='Minimum alignment score (default: 150)')
    parser.add_argument('--gap-init', type=int, default=20,
                       help='Gap initiation penalty (default: 20)')
    parser.add_argument('--gap-ext', type=int, default=5,
                       help='Gap extension penalty (default: 5)')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate input file
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        # Create refiner and build consensus
        refiner = EnhancedTERefiner(
            min_score=args.min_score,
            gap_init=args.gap_init,
            gap_ext=args.gap_ext,
            threads=args.threads
        )
        
        success = refiner.build_consensus(args.input, args.output)
        
        if success:
            logger.info("Consensus building completed successfully")
            sys.exit(0)
        else:
            logger.error("Consensus building failed")
            sys.exit(1)
            
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()