#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline
import numpy as np
import networkx as nx
from collections import defaultdict
import subprocess
import tempfile
import os
import multiprocessing
from itertools import combinations

class ConservativeMSAConsensusBuilder:
    def __init__(self, 
                 min_similarity=0.6,           # 降低相似度阈值
                 min_coverage=0.3,             # 最小覆盖度要求
                 min_base_freq=0.5,            # 最小碱基频率要求
                 threads=None):
        self.min_similarity = min_similarity
        self.min_coverage = min_coverage
        self.min_base_freq = min_base_freq
        self.threads = threads or multiprocessing.cpu_count()

    def read_sequences(self, fasta_file):
        return list(SeqIO.parse(fasta_file, "fasta"))

    def preprocess_sequences(self, sequences):
        """预处理序列:
        1. 计算长度分布
        2. 过滤掉过短的序列
        3. 识别完整序列
        """
        lengths = [len(seq) for seq in sequences]
        median_len = np.median(lengths)
        std_len = np.std(lengths)
        
        # 保留长度在中位数±2个标准差范围内的序列
        filtered_seqs = [
            seq for seq in sequences
            if median_len - 2*std_len <= len(seq) <= median_len + 2*std_len
        ]
        
        print(f"Filtered {len(sequences) - len(filtered_seqs)} sequences")
        return filtered_seqs

    def run_mafft(self, sequences):
        """运行MAFFT进行比对"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fa') as temp_in:
            SeqIO.write(sequences, temp_in.name, "fasta")
            temp_in_name = temp_in.name

        temp_out = tempfile.NamedTemporaryFile(delete=False, suffix='.aln')
        temp_out.close()

        try:
            # 使用 ginsi 算法提高精确度
            mafft_cline = MafftCommandline(
                input=temp_in_name,
                thread=self.threads,
                maxiterate=1000,  # 增加迭代次数
                globalpair=True,  # 使用全局配对
                adjustdirection=True  # 自动调整序列方向
            )
            stdout, stderr = mafft_cline()
            
            with open(temp_out.name, "w") as handle:
                handle.write(stdout)
            
            alignment = list(SeqIO.parse(temp_out.name, "fasta"))
            
        finally:
            os.unlink(temp_in_name)
            os.unlink(temp_out.name)
            
        return alignment

    def analyze_alignment(self, alignment):
        """分析比对结果,包括:
        1. 位置覆盖度
        2. 碱基频率
        3. 保守性
        """
        msa = MultipleSeqAlignment(alignment)
        seq_count = len(msa)
        alignment_length = msa.get_alignment_length()
        
        coverage = []
        conservation = []
        base_freqs = []
        
        for i in range(alignment_length):
            column = msa[:, i]
            
            # 计算非gap的比例
            non_gaps = sum(1 for base in column if base != '-')
            cov = non_gaps / seq_count
            coverage.append(cov)
            
            # 计算碱基频率
            base_counts = defaultdict(int)
            for base in column:
                if base != '-':
                    base_counts[base] += 1
            
            if non_gaps > 0:
                max_freq = max(base_counts.values()) / non_gaps
                max_base = max(base_counts.items(), key=lambda x: x[1])[0]
            else:
                max_freq = 0
                max_base = '-'
                
            base_freqs.append((max_base, max_freq))
            
            # 计算保守性(香农熵)
            freqs = [count/non_gaps for count in base_counts.values()]
            entropy = -sum(f * np.log2(f) for f in freqs) if freqs else 0
            conservation.append(entropy)
            
        return coverage, base_freqs, conservation

    def generate_conservative_consensus(self, alignment):
        """生成保守的共识序列:
        1. 只在覆盖度高的位置生成共识
        2. 要求较高的碱基一致性
        3. 保留高度保守的区域
        """
        msa = MultipleSeqAlignment(alignment)
        coverage, base_freqs, conservation = self.analyze_alignment(alignment)
        
        consensus = []
        for i in range(msa.get_alignment_length()):
            if (coverage[i] >= self.min_coverage and  # 覆盖度要求
                base_freqs[i][1] >= self.min_base_freq):  # 碱基频率要求
                consensus.append(base_freqs[i][0])
            else:
                consensus.append('N')  # 不确定位点用N标记
        
        # 移除连续的N
        clean_consensus = []
        n_count = 0
        for base in consensus:
            if base != 'N':
                if n_count > 0 and n_count < 10:  # 短的N区域保留
                    clean_consensus.extend(['N'] * n_count)
                n_count = 0
                clean_consensus.append(base)
            else:
                n_count += 1
                
        return ''.join(clean_consensus)

    def build_consensus(self, input_file, output_file):
        """主流程"""
        try:
            print("Reading sequences...")
            sequences = self.read_sequences(input_file)
            if not sequences:
                raise ValueError(f"No sequences found in {input_file}")
            print(f"Read {len(sequences)} sequences")

            print("Preprocessing sequences...")
            filtered_seqs = self.preprocess_sequences(sequences)
            print(f"Retained {len(filtered_seqs)} sequences")

            print("Running multiple sequence alignment...")
            aligned_seqs = self.run_mafft(filtered_seqs)
            print("Completed alignment")

            print("Generating consensus sequence...")
            consensus = self.generate_conservative_consensus(aligned_seqs)
            print(f"Generated consensus sequence of length {len(consensus)}")

            # 统计结果
            n_count = consensus.count('N')
            n_percent = (n_count / len(consensus)) * 100
            print(f"Consensus contains {n_count} N's ({n_percent:.2f}%)")

            # 写入结果
            consensus_record = SeqRecord(
                Seq(consensus),
                id="consensus",
                description=f"length={len(consensus)} from_{len(filtered_seqs)}_sequences"
            )
            SeqIO.write(consensus_record, output_file, "fasta")
            print(f"Written consensus to {output_file}")

            # 生成统计报告
            with open(output_file + ".stats", "w") as f:
                f.write(f"Original sequences: {len(sequences)}\n")
                f.write(f"Filtered sequences: {len(filtered_seqs)}\n")
                f.write(f"Consensus length: {len(consensus)}\n")
                f.write(f"N positions: {n_count} ({n_percent:.2f}%)\n")

        except Exception as e:
            print(f"Error: {str(e)}")
            raise

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Build conservative consensus sequence')
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('output', help='Output FASTA file')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads')
    parser.add_argument('-s', '--similarity', type=float, default=0.6, 
                        help='Minimum similarity threshold (default: 0.6)')
    parser.add_argument('-c', '--coverage', type=float, default=0.3,
                        help='Minimum coverage threshold (default: 0.3)')
    parser.add_argument('-f', '--frequency', type=float, default=0.5,
                        help='Minimum base frequency threshold (default: 0.5)')
    
    args = parser.parse_args()
    
    builder = ConservativeMSAConsensusBuilder(
        min_similarity=args.similarity,
        min_coverage=args.coverage,
        min_base_freq=args.frequency,
        threads=args.threads
    )
    builder.build_consensus(args.input, args.output)
