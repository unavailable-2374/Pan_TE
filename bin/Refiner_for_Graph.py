#!/usr/bin/env python3
import os
import sys
import time
import random
import subprocess
import shutil
import tempfile
import argparse
import json
import logging
from collections import defaultdict, Counter
from io import StringIO

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

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


class SequenceClusterer:
    def __init__(self, te_builder, distance_threshold=0.7):
        self.te_builder = te_builder
        self.distance_threshold = distance_threshold

    def calculate_distance_matrix(self, sequences):
        n_seqs = len(sequences)
        unique_id = f"{int(time.time())}_{random.randint(1000, 9999)}"

        if not hasattr(self.te_builder, 'temp_dir'):
            raise ValueError("TEConsensusBuilder temp_dir not initialized")

        temp_name = os.path.join(self.te_builder.temp_dir, f'ref_sequences_{unique_id}.fa')
        with open(temp_name, 'w') as temp_file:
            SeqIO.write(sequences, temp_file, "fasta")

        try:
            alignments = self.te_builder.run_rmblast(temp_name, temp_name)
            seq_lengths = {seq.id: len(seq.seq) for seq in sequences}
            seq_id_to_idx = {seq.id: idx for idx, seq in enumerate(sequences)}

            similarities = np.zeros((n_seqs, n_seqs))
            for aln in alignments:
                if aln.query_id != aln.subject_id:
                    i = seq_id_to_idx[aln.query_id]
                    j = seq_id_to_idx[aln.subject_id]
                    min_len = min(seq_lengths[aln.query_id], seq_lengths[aln.subject_id])
                    norm_score = aln.score / (min_len * 2)
                    similarities[i, j] = max(similarities[i, j], norm_score)
                    similarities[j, i] = similarities[i, j]

            np.fill_diagonal(similarities, 1.0)
            distances = 1.0 - similarities
            distances = np.maximum(distances, 0.0)

            logger.info(f"Distance matrix shape: {distances.shape}")
            logger.info(f"Distance range: [{distances.min()}, {distances.max()}]")

            return distances

        finally:
            try:
                os.remove(temp_name)
                for ext in ['.nin', '.nsq', '.nhr']:
                    db_file = temp_name + ext
                    if os.path.exists(db_file):
                        os.remove(db_file)
            except OSError as e:
                logger.warning(f"Error cleaning up temporary files: {e}")

    def cluster_sequences(self, sequences):
        if len(sequences) == 1:
            return [[sequences[0]]]

        logger.info("Calculating distance matrix...")
        distances = self.calculate_distance_matrix(sequences)

        logger.info("Performing hierarchical clustering...")
        try:
            distances = np.maximum(distances, distances.T)
            condensed_distances = squareform(distances)
            linkage_matrix = linkage(condensed_distances, method='average')
            clusters = fcluster(linkage_matrix, t=self.distance_threshold, criterion='distance')

            cluster_dict = defaultdict(list)
            for seq, cluster_id in zip(sequences, clusters):
                cluster_dict[cluster_id].append(seq)

            return list(cluster_dict.values())

        except Exception as e:
            logger.error(f"Clustering error: {str(e)}")
            logger.error(f"Distance matrix stats - min: {distances.min()}, max: {distances.max()}, mean: {distances.mean()}")
            raise


class TEConsensusBuilder:
    def __init__(self, rmblast_dir, makeblastdb_path, matrix_path,
                 min_score=150, gap_init=20, gap_ext=5, threads=None):
        self.rmblast_path = os.path.join(rmblast_dir, "rmblastn")
        self.makeblastdb_path = makeblastdb_path
        self.matrix_path = matrix_path
        self.min_score = min_score
        self.gap_init = gap_init
        self.gap_ext = gap_ext
        self.threads = threads or 1
        self.temp_dir = None

        self.mafft_path = shutil.which("mafft")
        if not self.mafft_path:
            raise FileNotFoundError("mafft not found in system PATH")

    def prepare_blast_db(self, fasta_file):
        cmd = [
            self.makeblastdb_path,
            "-in", fasta_file,
            "-dbtype", "nucl",
            "-parse_seqids"
        ]
        try:
            logger.info(f"Creating BLAST database: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info("BLAST database creation successful")
            if result.stderr:
                logger.debug(f"makeblastdb stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create BLAST database: {e.stderr}")
            raise

    def run_rmblast(self, query_file, subject_file):
        current_dir = os.getcwd()
        query_path = os.path.abspath(query_file)
        subject_path = os.path.abspath(subject_file)
        work_dir = os.path.dirname(subject_path)

        os.chdir(work_dir)
        try:
            self.prepare_blast_db(subject_path)

            if not os.path.exists(query_path):
                logger.error(f"Query file not found: {query_path}")
                raise FileNotFoundError(f"Query file not found: {query_path}")
            if not os.path.exists(subject_path):
                logger.error(f"Subject file not found: {subject_path}")
                raise FileNotFoundError(f"Subject file not found: {subject_path}")

            self.prepare_blast_db(os.path.basename(subject_path))

            db_files = [os.path.basename(subject_path) + ext for ext in ['.nhr', '.nin', '.nsq']]
            for db_file in db_files:
                if not os.path.exists(db_file):
                    raise FileNotFoundError(f"BLAST database file not found: {db_file}")

            cmd = [
                self.rmblast_path,
                "-query", query_path,
                "-db", subject_path,
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
            os.chdir(current_dir)

    def build_multiple_alignment(self, sequences, reference=None):
        with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fa", dir=self.temp_dir) as temp_in:
            SeqIO.write(sequences, temp_in, "fasta")
            temp_in_name = temp_in.name

        try:
            cmd = [self.mafft_path, "--auto", "--thread", str(self.threads), temp_in_name]
            logger.info("Running MAFFT for multiple sequence alignment: " + " ".join(cmd))
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            aligned_output = result.stdout

            aligned_seqs = list(SeqIO.parse(StringIO(aligned_output), "fasta"))
            return [str(record.seq) for record in aligned_seqs]

        except subprocess.CalledProcessError as e:
            logger.error(f"MAFFT failed: {e.stderr}")
            raise

        finally:
            os.remove(temp_in_name)

    def build_consensus(self, aligned_sequences):
        if not aligned_sequences:
            return ""

        seq_length = len(aligned_sequences[0])
        consensus = []

        for i in range(seq_length):
            column = [seq[i] for seq in aligned_sequences if i < len(seq)]
            base_counts = Counter(b for b in column if b != '-')
            if base_counts:
                consensus.append(base_counts.most_common(1)[0][0])
            else:
                consensus.append('-')

        return ''.join(consensus).replace('-', '')

    def mash_precluster(self, sequences, mash_threshold=0.1):
        temp_fa = None
        sketch_prefix = None
        try:
            with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fa", dir=self.temp_dir) as temp:
                SeqIO.write(sequences, temp, "fasta")
                temp_fa = temp.name
                logger.info(f"Created temporary FASTA file: {temp_fa}")
                temp.flush()
                os.fsync(temp.fileno())

            seq_id_map = {str(seq.id): seq for seq in sequences}
            
            sketch_prefix = temp_fa + ".mash"
            
            mash_path = shutil.which("mash")
            if not mash_path:
                raise FileNotFoundError("mash command not found in system PATH")
            logger.info(f"Found mash at: {mash_path}")
                
            sketch_cmd = [mash_path, "sketch", "-k", "16", "-m", "2", "-o", sketch_prefix, temp_fa]
            logger.info("Running mash sketch: " + " ".join(sketch_cmd))
            result = subprocess.run(sketch_cmd, 
                                check=True, 
                                capture_output=True, 
                                text=True)
            if result.stderr:
                logger.info(f"Mash sketch stderr: {result.stderr}")

            msh_file = sketch_prefix + ".msh"
            if not os.path.exists(msh_file):
                raise FileNotFoundError(f"Mash sketch output file not found: {msh_file}")
            logger.info(f"Mash sketch file created: {msh_file}")

            dist_cmd = [mash_path, "dist", msh_file, temp_fa]
            logger.info("Running mash dist: " + " ".join(dist_cmd))
            result = subprocess.run(dist_cmd, 
                                check=True, 
                                capture_output=True, 
                                text=True)
            
            if result.stderr:
                logger.info(f"Mash dist stderr: {result.stderr}")
            
            mash_output = result.stdout
            if not mash_output.strip():
                raise RuntimeError("Mash dist produced no output")

            try:
                uf = UnionFind([str(seq.id) for seq in sequences])
                logger.info(f"Initialized UnionFind with {len(sequences)} sequences")
                
                cluster_count = 0
                for line in mash_output.strip().split("\n"):
                    parts = line.split("\t") 
                    if len(parts) < 3:
                        continue
                    
                    seq1_path = parts[0]
                    seq2_path = parts[1]
                    seq1_id = os.path.basename(seq1_path)
                    seq2_id = os.path.basename(seq2_path)
                    
                    try:
                        dist = float(parts[2])
                        if dist < mash_threshold and seq1_id in seq_id_map and seq2_id in seq_id_map:
                            uf.union(seq1_id, seq2_id)
                            cluster_count += 1
                    except (ValueError, KeyError) as e:
                        logger.warning(f"Error processing distance for {seq1_id}-{seq2_id}: {e}")
                        continue
                        
                logger.info(f"Created {cluster_count} initial clusters")

                clusters_dict = defaultdict(list)
                for seq in sequences:
                    seq_id = str(seq.id)
                    root_id = uf.find(seq_id)
                    clusters_dict[root_id].append(seq)
                    
                preclusters = list(clusters_dict.values())
                logger.info(f"Formed {len(preclusters)} preclusters")

                final_clusters = []
                clusterer = SequenceClusterer(self)
                for i, group in enumerate(preclusters):
                    logger.info(f"Processing precluster {i+1}/{len(preclusters)} with {len(group)} sequences")
                    if len(group) == 1:
                        final_clusters.append(group)
                    else:
                        refined = clusterer.cluster_sequences(group)
                        final_clusters.extend(refined)
                        
                logger.info(f"Created {len(final_clusters)} final clusters")
                return final_clusters

            except Exception as e:
                logger.error(f"Error during clustering: {str(e)}")
                raise

        except Exception as e:
            logger.error(f"Error in mash_precluster: {str(e)}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
            raise

        finally:
            try:
                if temp_fa and os.path.exists(temp_fa):
                    os.remove(temp_fa)
                    logger.info(f"Removed temporary file: {temp_fa}")
                if sketch_prefix:
                    msh_file = sketch_prefix + ".msh"
                    if os.path.exists(msh_file):
                        os.remove(msh_file)
                        logger.info(f"Removed mash sketch file: {msh_file}")
            except OSError as e:
                logger.warning(f"Error cleaning up temporary files: {e}")

    def hybrid_cluster(self, sequences):
        with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fa", dir=self.temp_dir) as temp_in:
            SeqIO.write(sequences, temp_in, "fasta")
            temp_in_name = temp_in.name

        temp_out_name = temp_in_name + ".cdhit"
        cmd = [
            "cd-hit-est",
            "-i", temp_in_name,
            "-o", temp_out_name,
            "-aS", "0.8",
            "-c", "0.8",
            "-g", "0",
            "-G", "0",
            "-A", "80",
            "-M", "10000",
            "-t", str(self.threads)
        ]
        logger.info("Running cd-hit-est: " + " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info("cd-hit-est finished successfully.")
        except subprocess.CalledProcessError as e:
            logger.error("cd-hit-est failed: " + e.stderr)
            raise

        preclustered_sequences = list(SeqIO.parse(temp_out_name, "fasta"))

        os.remove(temp_in_name)
        os.remove(temp_out_name)

        clusterer = SequenceClusterer(self)
        refined_clusters = clusterer.cluster_sequences(preclustered_sequences)
        return refined_clusters

    def build_clustered_consensus(self, input_file, output_file):
        original_dir = os.getcwd()
        try:
            input_path = os.path.abspath(input_file)
            output_path = os.path.abspath(output_file)
            output_dir = os.path.dirname(output_path)
            logger.info(f"Using base directory for consensus building: {output_dir}")

            timestamp = int(time.time())
            self.temp_dir = os.path.join(output_dir, f'tmp_{timestamp}')
            os.makedirs(self.temp_dir, exist_ok=True)
            os.chdir(self.temp_dir)
            logger.info(f"Created temporary directory: {self.temp_dir}")

            logger.info(f"Reading sequences from: {input_path}")
            sequences = list(SeqIO.parse(input_path, "fasta"))
            if not sequences:
                raise ValueError(f"No sequences found in {input_path}")
            logger.info(f"Read {len(sequences)} sequences")

            if any(len(seq.seq) > 50000 for seq in sequences):
                logger.info("Long sequences (>50k) detected; using k-mer (Mash) pre-clustering.")
                clusters = self.mash_precluster(sequences)
            elif os.path.getsize(input_path) > 2 * 1024 * 1024:
                logger.info("Large input file detected; using hybrid clustering strategy with cd-hit-est pre-clustering.")
                clusters = self.hybrid_cluster(sequences)
            else:
                logger.info("Using standard RMBlast-based clustering.")
                clusterer = SequenceClusterer(self)
                clusters = clusterer.cluster_sequences(sequences)

            logger.info(f"Found {len(clusters)} clusters")

            consensus_records = []
            for i, cluster in enumerate(clusters, 1):
                logger.info(f"Processing cluster {i} with {len(cluster)} sequences")
                if len(cluster) == 1:
                    consensus_seq = str(cluster[0].seq)
                    consensus_desc = f"single sequence from cluster {i}"
                else:
                    aligned_seqs = self.build_multiple_alignment(cluster)
                    consensus_seq = self.build_consensus(aligned_seqs)
                    consensus_desc = f"consensus from {len(cluster)} sequences in cluster {i}"
                consensus_record = SeqRecord(
                    Seq(consensus_seq),
                    id=f"{os.path.splitext(os.path.basename(input_file))[0]}_cluster_{i}",
                    description=consensus_desc
                )
                consensus_records.append(consensus_record)

            SeqIO.write(consensus_records, output_path, "fasta")
            logger.info(f"Written {len(consensus_records)} consensus sequences to {output_path}")

            stats_file = f"{output_path}.stats"
            with open(stats_file, "w") as f:
                f.write(f"Original sequences: {len(sequences)}\n")
                f.write(f"Number of clusters: {len(clusters)}\n")
                for i, cluster in enumerate(clusters, 1):
                    f.write(f"Cluster {i} size: {len(cluster)}\n")
            logger.info(f"Written statistics to {stats_file}")

        except Exception as e:
            logger.error(f"Error in consensus building: {str(e)}")
            raise
        finally:
            os.chdir(original_dir)
            if self.temp_dir and os.path.exists(self.temp_dir):
                logger.info("Temporary directory will be cleaned by calling script")


class UnionFind:
    def __init__(self, elements):
        self.parent = {e: e for e in elements}

    def find(self, a):
        if self.parent[a] != a:
            self.parent[a] = self.find(self.parent[a])
        return self.parent[a]

    def union(self, a, b):
        rootA = self.find(a)
        rootB = self.find(b)
        if rootA != rootB:
            self.parent[rootB] = rootA


class Config:
    def __init__(self):
        conda_env_path = self.get_current_conda_env()
        if conda_env_path:
            self.rmblastn = os.path.join(conda_env_path, 'bin', 'rmblastn')
            self.makeblastdb = os.path.join(conda_env_path, 'bin', 'makeblastdb')
            logger.info(f"Using Conda environment: {conda_env_path}")
        else:
            self.rmblastn = shutil.which('rmblastn')
            self.makeblastdb = shutil.which('makeblastdb')
            logger.info("No Conda environment detected, using system PATH")

        if not self.rmblastn or not os.path.exists(self.rmblastn):
            raise FileNotFoundError(f"rmblastn not found at {self.rmblastn}")
        if not self.makeblastdb or not os.path.exists(self.makeblastdb):
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
                conda_info = json.loads(result)
                active_prefix = conda_info.get('active_prefix')
                if active_prefix:
                    return active_prefix
        except (subprocess.CalledProcessError, json.JSONDecodeError):
            pass

        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Build consensus sequences for clustered transposable element families using RMBlast for clustering and MAFFT for multiple sequence alignment. '
                    'When sequences longer than 50k exist, a k-mer (Mash) pre-clustering is applied first, followed by RMBlast-based fine clustering; '
                    'otherwise, a hybrid clustering (using cd-hit-est pre-clustering) or standard RMBlast-based clustering is used based on input file size.')
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('output', help='Output FASTA file')
    parser.add_argument('-t', '--threads', type=int,
                        help='Number of threads', default=1)
    parser.add_argument('--min-score', type=float, default=150,
                        help='Minimum alignment score (default: 150)')
    parser.add_argument('--gap-init', type=int, default=20,
                        help='Gap initiation penalty (default: 20)')
    parser.add_argument('--gap-ext', type=int, default=5,
                        help='Gap extension penalty (default: 5)')
    parser.add_argument('--distance-threshold', type=float, default=0.7,
                        help='Distance threshold for clustering (default: 0.7)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger(__name__).setLevel(logging.DEBUG)

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
        builder.build_clustered_consensus(args.input, args.output)

    except Exception as e:
        logger.error(f"Error details: {str(e)}")
        logger.error(f"Full traceback: {traceback.format_exc()}")
        sys.exit(1)
