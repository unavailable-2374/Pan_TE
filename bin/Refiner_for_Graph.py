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
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class SequenceClusterer:
    def __init__(self, te_builder, distance_threshold=0.7):
        self.te_builder = te_builder
        self.distance_threshold = distance_threshold
    
    def calculate_distance_matrix(self, sequences):
        n_seqs = len(sequences)
        distances = np.zeros((n_seqs, n_seqs))
        
        # 创建临时文件
        if not os.path.exists('tmp'):
            os.makedirs('tmp')
        temp_name = os.path.join('tmp', f'all_sequences_{os.getpid()}.fa')
        with open(temp_name, 'w') as temp_file:
            SeqIO.write(sequences, temp_file, "fasta")
            
        try:
            # 使用RMBlast获取所有序列之间的比对
            alignments = self.te_builder.run_rmblast(temp_name, temp_name)
            
            # 根据比对结果计算距离
            seq_lengths = {seq.id: len(seq.seq) for seq in sequences}
            seq_id_to_idx = {seq.id: idx for idx, seq in enumerate(sequences)}
            
            # 初始化距离矩阵为1（最大距离）
            distances.fill(1.0)
            np.fill_diagonal(distances, 0.0)  # 对角线为0
            
            # 根据比对结果更新距离
            for aln in alignments:
                if aln.query_id != aln.subject_id:
                    i = seq_id_to_idx[aln.query_id]
                    j = seq_id_to_idx[aln.subject_id]
                    
                    # 计算相似度得分
                    similarity = aln.score / max(seq_lengths[aln.query_id], 
                                              seq_lengths[aln.subject_id])
                    
                    # 将相似度转换为距离（1 - 相似度）
                    distance = 1.0 - similarity
                    
                    # 取最小距离（最大相似度）
                    distances[i,j] = min(distances[i,j], distance)
                    distances[j,i] = distances[i,j]  # 保持对称性
                    
            return distances
            
        finally:
            # 清理临时文件
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
            
        # 计算距离矩阵
        logger.info("Calculating distance matrix...")
        distances = self.calculate_distance_matrix(sequences)
        
        # 执行层次聚类
        logger.info("Performing hierarchical clustering...")
        linkage_matrix = linkage(squareform(distances), method='average')
        
        # 根据距离阈值切分聚类
        clusters = fcluster(linkage_matrix, t=self.distance_threshold, 
                          criterion='distance')
        
        # 将序列分配到对应的簇
        cluster_dict = defaultdict(list)
        for seq, cluster_id in zip(sequences, clusters):
            cluster_dict[cluster_id].append(seq)
            
        return list(cluster_dict.values())

class TEConsensusBuilder:
    # ... (保持原有的TEConsensusBuilder类代码不变) ...

    def build_clustered_consensus(self, input_file, output_file):
        try:
            logger.info("Reading sequences...")
            sequences = list(SeqIO.parse(input_file, "fasta"))
            if not sequences:
                raise ValueError(f"No sequences found in {input_file}")
            logger.info(f"Read {len(sequences)} sequences")

            # 创建聚类器并执行聚类
            clusterer = SequenceClusterer(self)
            clusters = clusterer.cluster_sequences(sequences)
            logger.info(f"Found {len(clusters)} clusters")

            # 为每个簇构建共识序列
            consensus_records = []
            for i, cluster in enumerate(clusters, 1):
                logger.info(f"Processing cluster {i} with {len(cluster)} sequences")
                
                if len(cluster) == 1:
                    # 如果簇中只有一条序列，直接使用该序列
                    consensus_seq = str(cluster[0].seq)
                    consensus_desc = f"single sequence from cluster {i}"
                else:
                    # 对多序列簇构建共识序列
                    reference = self.find_best_reference(cluster)
                    aligned_seqs = self.build_multiple_alignment(cluster, reference)
                    consensus_seq = self.build_consensus(aligned_seqs)
                    consensus_desc = f"consensus from {len(cluster)} sequences in cluster {i}"

                consensus_record = SeqRecord(
                    Seq(consensus_seq),
                    id=f"{os.path.splitext(os.path.basename(input_file))[0]}_cluster_{i}",
                    description=consensus_desc
                )
                consensus_records.append(consensus_record)

            # 写入所有共识序列
            SeqIO.write(consensus_records, output_file, "fasta")
            logger.info(f"Written {len(consensus_records)} consensus sequences to {output_file}")

            # 写入统计信息
            stats_file = f"{output_file}.stats"
            with open(stats_file, "w") as f:
                f.write(f"Original sequences: {len(sequences)}\n")
                f.write(f"Number of clusters: {len(clusters)}\n")
                for i, cluster in enumerate(clusters, 1):
                    f.write(f"Cluster {i} size: {len(cluster)}\n")
            logger.info(f"Written statistics to {stats_file}")

        except Exception as e:
            logger.error(f"Error in consensus building: {str(e)}")
            raise

if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description='Build consensus sequences for clustered transposable element families using RMBlast')
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
    parser.add_argument('--distance-threshold', type=float, default=0.7,
                       help='Distance threshold for clustering (default: 0.7)')
    
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
        # 使用新的聚类功能构建共识序列
        builder.build_clustered_consensus(args.input, args.output)
        
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)