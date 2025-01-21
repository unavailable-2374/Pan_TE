#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from collections import defaultdict, Counter
import subprocess
import os
import logging
from itertools import combinations

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class RMBlastAlignment:
    def __init__(self, query_id, subject_id, score, query_start, query_end,
                 subject_start, subject_end, alignment, orientation):
        self.query_id = query_id
        self.subject_id = subject_id
        self.score = score
        self.query_start = query_start
        self.query_end = query_end
        self.subject_start = subject_start
        self.subject_end = subject_end
        self.alignment = alignment
        self.orientation = orientation

class TEConsensusBuilder:
    def __init__(self, rmblast_dir, makeblastdb_path, matrix_path,
                 min_score=150, gap_init=-25, gap_ext=-5, 
                 threads=None):
        self.rmblast_path = os.path.join(rmblast_dir, "rmblastn")
        self.makeblastdb_path = makeblastdb_path
        self.matrix_path = matrix_path
        self.min_score = min_score
        self.gap_init = gap_init
        self.gap_ext = gap_ext
        self.threads = threads or 1

    def prepare_blast_db(self, fasta_file):
        cmd = [
            self.makeblastdb_path,
            "-in", fasta_file,
            "-dbtype", "nucl"
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"Created BLAST database for {fasta_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create BLAST database: {e.stderr.decode()}")
            raise

    def run_rmblast(self, query_file, subject_file):

        self.prepare_blast_db(subject_file)

        cmd = [
            self.rmblast_path,
            "-query", query_file,
            "-db", subject_file,
            "-outfmt", "6 qseqid sseqid score qstart qend sstart send qseq sseq sstrand",
            "-matrix", self.matrix_path,
            "-gapopen", str(self.gap_init),
            "-gapextend", str(self.gap_ext),
            "-dust", "no",  
            "-soft_masking", "false",
            "-num_threads", str(self.threads),
            "-complexity_adjust",
            "-evalue", "1e-10",
            "-word_size", "7",
            "-window_size", "40",      
            "-xdrop_gap", "50",        
            "-xdrop_gap_final", "100"  
        ]

        logger.info(f"Running RMBlast command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            alignments = []
            for line in result.stdout.split('\n'):
                if line.strip():
                    fields = line.split('\t')
                    if len(fields) >= 9:
                        alignment = RMBlastAlignment(
                            query_id=fields[0],
                            subject_id=fields[1],
                            score=float(fields[2]),
                            query_start=int(fields[3]),
                            query_end=int(fields[4]),
                            subject_start=int(fields[5]),
                            subject_end=int(fields[6]),
                            alignment=(fields[7], fields[8]),
                            orientation=fields[9] if len(fields) > 9 else 'plus'
                        )
                        if alignment.score >= self.min_score:
                            alignments.append(alignment)
            
            return alignments

        except subprocess.CalledProcessError as e:
            logger.error(f"RMBlast failed with error: {e.stderr}")
            raise
        finally:
            for ext in ['.nin', '.nsq', '.nhr']:
                try:
                    os.remove(subject_file + ext)
                except OSError:
                    pass

    def find_best_reference(self, sequences):

        if not os.path.exists('tmp'):
            os.makedirs('tmp')
        
        temp_name = os.path.join('tmp', f'ref_sequences_{os.getpid()}.fa')
        with open(temp_name, 'w') as temp_file:
            SeqIO.write(sequences, temp_file, "fasta")

        try:
            alignments = self.run_rmblast(temp_name, temp_name)
            
            sequence_scores = defaultdict(float)
            for aln in alignments:
                if aln.query_id != aln.subject_id: 
                    sequence_scores[aln.query_id] += aln.score
            
            if sequence_scores:
                best_seq_id = max(sequence_scores.items(), key=lambda x: x[1])[0]
                return next(seq for seq in sequences if seq.id == best_seq_id)
            else:
                return sequences[0]  

        finally:

            try:
                os.remove(temp_name)

                for ext in ['.nin', '.nsq', '.nhr']:
                    db_file = temp_name + ext
                    if os.path.exists(db_file):
                        os.remove(db_file)
            except OSError as e:
                logger.warning(f"Error cleaning up temporary files: {e}")

    def build_multiple_alignment(self, sequences, reference):

        if not os.path.exists('tmp'):
            os.makedirs('tmp')

        ref_name = os.path.join('tmp', f'reference_{os.getpid()}.fa')
        query_name = os.path.join('tmp', f'queries_{os.getpid()}.fa')
        
        with open(ref_name, 'w') as ref_file, open(query_name, 'w') as query_file:

            SeqIO.write([reference], ref_file, "fasta")

            SeqIO.write(sequences, query_file, "fasta")
        
        try:

            alignments = self.run_rmblast(query_name, ref_name)
            
            aligned_sequences = self.process_alignments(alignments, reference)
            return aligned_sequences
            
        finally:

            for fname in [ref_name, query_name]:
                try:
                    os.remove(fname)
                    for ext in ['.nin', '.nsq', '.nhr']:
                        db_file = fname + ext
                        if os.path.exists(db_file):
                            os.remove(db_file)
                except OSError as e:
                    logger.warning(f"Error cleaning up temporary files: {e}")

    def process_alignments(self, alignments, reference):

        ref_length = len(reference.seq)
        query_aln_dict = {}

        query_score_dict = {}

        all_query_ids = set([aln.query_id for aln in alignments])
        for qid in all_query_ids:
            query_aln_dict[qid] = ['-' for _ in range(ref_length)]
            query_score_dict[qid] = [0.0 for _ in range(ref_length)]

        for aln in alignments:
            qid = aln.query_id
            if aln.orientation == 'plus':
                qseq = aln.alignment[0]
            else:
                qseq = str(Seq(aln.alignment[0]).reverse_complement())

            subj_start = aln.subject_start - 1
            subj_end   = aln.subject_end   - 1

            aligned_len = len(qseq)
            expected_len = (subj_end - subj_start + 1)
            if aligned_len != expected_len:
                continue
            for i in range(aligned_len):
                ref_pos = subj_start + i
                base    = qseq[i]
                if aln.score > query_score_dict[qid][ref_pos]:
                    query_aln_dict[qid][ref_pos] = base
                    query_score_dict[qid][ref_pos] = aln.score

        final_seqs = []
        for qid in all_query_ids:
            merged_seq = ''.join(query_aln_dict[qid])
            final_seqs.append(merged_seq)

        return final_seqs
    def build_consensus(self, aligned_sequences):
        if not aligned_sequences:
            return ""
        
        seq_length = len(aligned_sequences[0])
        consensus = []
        
        for i in range(seq_length):
            bases = [seq[i] for seq in aligned_sequences if i < len(seq)]
            base_counts = Counter(base for base in bases if base != '-')
            
            if base_counts:
                consensus.append(max(base_counts.items(), key=lambda x: x[1])[0])
            else:
                consensus.append('-')
                
        return ''.join(consensus).replace('-', '')

    def build_te_consensus(self, input_file, output_file):
        try:
            logger.info("Reading sequences...")
            sequences = list(SeqIO.parse(input_file, "fasta"))
            if not sequences:
                raise ValueError(f"No sequences found in {input_file}")
            logger.info(f"Read {len(sequences)} sequences")

            logger.info("Finding best reference sequence...")
            reference = self.find_best_reference(sequences)
            logger.info(f"Selected {reference.id} as reference sequence")

            logger.info("Building multiple alignment...")
            aligned_seqs = self.build_multiple_alignment(sequences, reference)
            logger.info("Completed alignment")

            logger.info("Building consensus sequence...")
            consensus = self.build_consensus(aligned_seqs)
            logger.info(f"Generated consensus sequence of length {len(consensus)}")

            consensus_record = SeqRecord(
                Seq(consensus),
                id=os.path.splitext(os.path.basename(input_file))[0],
                description=f"consensus from {len(sequences)} sequences"
            )
            SeqIO.write(consensus_record, output_file, "fasta")
            logger.info(f"Written consensus to {output_file}")

            stats_file = f"{output_file}.stats"
            with open(stats_file, "w") as f:
                f.write(f"Original sequences: {len(sequences)}\n")
                f.write(f"Consensus length: {len(consensus)}\n")
                f.write(f"Reference sequence: {reference.id}\n")
            logger.info(f"Written statistics to {stats_file}")

        except Exception as e:
            logger.error(f"Error in consensus building: {str(e)}")
            raise

class Config:
    def __init__(self):
        import shutil
        
        self.rmblastn = shutil.which('rmblastn')
        self.makeblastdb = shutil.which('makeblastdb')
        
        if not self.rmblastn:
            raise FileNotFoundError("rmblastn not found in PATH")
        if not self.makeblastdb:
            raise FileNotFoundError("makeblastdb not found in PATH")

        conda_env_dir = os.path.dirname(os.path.dirname(self.rmblastn))

        possible_matrix_paths = [
            os.path.join(conda_env_dir, 'share/RepeatModeler/Matrices/ncbi/nt/comparison.matrix'),
            os.path.join(conda_env_dir, 'share/RepeatMasker/Libraries/Dfam.hmm'),
            '/usr/share/RepeatMasker/Matrices/nt',
            '/usr/local/share/RepeatMasker/Matrices/nt',
            os.path.join(conda_env_dir, 'share/RepeatMasker/Matrices/BLOSUM62')
        ]
        
        self.matrix_path = None
        for path in possible_matrix_paths:
            if os.path.exists(path):class Config:
    def __init__(self):
        import shutil
        
        conda_env_path = self.get_current_conda_env()
        if conda_env_path:
            self.rmblastn = os.path.join(conda_env_path, 'bin', 'rmblastn')
            self.makeblastdb = os.path.join(conda_env_path, 'bin', 'makeblastdb')
            logger.info(f"Using Conda environment: {conda_env_path}")
        else:
            self.rmblastn = shutil.which('rmblastn')
            self.makeblastdb = shutil.which('makeblastdb')
            logger.info("No Conda environment detected, using system PATH")
        
        if not os.path.exists(self.rmblastn):
            raise FileNotFoundError(f"rmblastn not found at {self.rmblastn}")
        if not os.path.exists(self.makeblastdb):
            raise FileNotFoundError(f"makeblastdb not found at {self.makeblastdb}")

        conda_env_dir = conda_env_path if conda_env_path else os.path.dirname(os.path.dirname(self.rmblastn))

        possible_matrix_paths = [
            os.path.join(conda_env_dir, 'share/RepeatModeler/Matrices/ncbi/nt/comparison.matrix'),
            os.path.join(conda_env_dir, 'share/RepeatMasker/Libraries/Dfam.hmm'),
            os.path.join(conda_env_dir, 'share/RepeatMasker/Matrices/nt'),
            os.path.join(conda_env_dir, 'share/RepeatMasker/Matrices/BLOSUM62')
        ]
        
        self.matrix_path = None
        for path in possible_matrix_paths:
            if os.path.exists(path):
                self.matrix_path = path
                break
                
        if self.matrix_path is None:
            self.matrix_path = 'BLOSUM62'
            
        logger.info(f"Using rmblastn from: {self.rmblastn}")
        logger.info(f"Using scoring matrix: {self.matrix_path}")

    def get_current_conda_env(self):

        conda_prefix = os.environ.get('CONDA_PREFIX')
        if conda_prefix:
            return conda_prefix

        import sys
        if 'conda' in sys.prefix:
            return sys.prefix

        try:
            conda_path = subprocess.check_output(['which', 'conda']).decode().strip()
            if conda_path:
                result = subprocess.check_output(['conda', 'info', '--json']).decode()
                import json
                conda_info = json.loads(result)
                active_prefix = conda_info.get('active_prefix')
                if active_prefix:
                    return active_prefix
        except (subprocess.CalledProcessError, json.JSONDecodeError):
            pass

        return None

if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description='Build consensus sequence for transposable element family using RMBlast')
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('output', help='Output FASTA file')
    parser.add_argument('-t', '--threads', type=int,
                       help='Number of threads')
    parser.add_argument('--min-score', type=float, default=150,
                       help='Minimum alignment score (default: 150)')
    parser.add_argument('--gap-init', type=int, default=20,
                       help='Gap initiation penalty (default: 20)')
    parser.add_argument('--gap-ext', type=int, default=5,
                       help='Gap extension penalty (default: 5)')
    
    args = parser.parse_args()
    
    try:
        config = Config()
        
        builder = TEConsensusBuilder(
            rmblast_dir=os.path.dirname(config.rmblastn),
            makeblastdb_path=config.makeblastdb,
            matrix_path=config.matrix_path,
            min_score=args.min_score,
            gap_init=args.gap_init,
            gap_ext=args.gap_ext,
            threads=args.threads
        )
        builder.build_te_consensus(args.input, args.output)
        
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)