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
        """创建BLAST数据库"""
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
        """运行RMBlast进行序列比对"""
        # 准备BLAST数据库
        self.prepare_blast_db(subject_file)

        # 设置RMBlast命令
        cmd = [
            self.rmblast_path,
            "-query", query_file,
            "-db", subject_file,
            "-outfmt", "6 qseqid sseqid score qstart qend sstart send qseq sseq sstrand",
            "-matrix", self.matrix_path,
            "-gapopen", str(self.gap_init),
            "-gapextend", str(self.gap_ext),
            "-dust", "no",  # 关闭 DUST 过滤
            "-soft_masking", "false",
            "-num_threads", str(self.threads),
            "-complexity_adjust",
            "-evalue", "1e-10",
            "-word_size", "7",
            "-window_size", "40",      # 设置窗口大小
            "-xdrop_gap", "50",        # X-dropoff for preliminary gapped extensions
            "-xdrop_gap_final", "100"  # X-dropoff for final gapped extensions
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
            # 清理临时文件
            for ext in ['.nin', '.nsq', '.nhr']:
                try:
                    os.remove(subject_file + ext)
                except OSError:
                    pass

    def find_best_reference(self, sequences):
        """找到最佳的参考序列"""
        # 确保tmp目录存在
        if not os.path.exists('tmp'):
            os.makedirs('tmp')
        
        # 在当前目录的tmp下创建临时文件
        temp_name = os.path.join('tmp', f'ref_sequences_{os.getpid()}.fa')
        with open(temp_name, 'w') as temp_file:
            SeqIO.write(sequences, temp_file, "fasta")

        try:
            # 运行所有序列对的比对
            alignments = self.run_rmblast(temp_name, temp_name)
            
            # 计算每个序列的总得分
            sequence_scores = defaultdict(float)
            for aln in alignments:
                if aln.query_id != aln.subject_id:  # 排除自比对
                    sequence_scores[aln.query_id] += aln.score
            
            # 选择得分最高的序列
            if sequence_scores:
                best_seq_id = max(sequence_scores.items(), key=lambda x: x[1])[0]
                return next(seq for seq in sequences if seq.id == best_seq_id)
            else:
                return sequences[0]  # 如果没有有效比对，返回第一个序列

        finally:
            # 清理临时文件
            try:
                os.remove(temp_name)
                # 清理BLAST数据库文件
                for ext in ['.nin', '.nsq', '.nhr']:
                    db_file = temp_name + ext
                    if os.path.exists(db_file):
                        os.remove(db_file)
            except OSError as e:
                logger.warning(f"Error cleaning up temporary files: {e}")

    def build_multiple_alignment(self, sequences, reference):
        """基于参考序列构建多序列比对"""
        # 确保tmp目录存在
        if not os.path.exists('tmp'):
            os.makedirs('tmp')
        
        # 在当前目录的tmp下创建临时文件
        ref_name = os.path.join('tmp', f'reference_{os.getpid()}.fa')
        query_name = os.path.join('tmp', f'queries_{os.getpid()}.fa')
        
        with open(ref_name, 'w') as ref_file, open(query_name, 'w') as query_file:
            # 写入参考序列
            SeqIO.write([reference], ref_file, "fasta")
            # 写入查询序列
            SeqIO.write(sequences, query_file, "fasta")
        
        try:
            # 运行比对
            alignments = self.run_rmblast(query_name, ref_name)
            
            # 处理比对结果
            aligned_sequences = self.process_alignments(alignments, reference)
            return aligned_sequences
            
        finally:
            # 清理临时文件和BLAST数据库文件
            for fname in [ref_name, query_name]:
                try:
                    os.remove(fname)
                    # 清理BLAST数据库文件
                    for ext in ['.nin', '.nsq', '.nhr']:
                        db_file = fname + ext
                        if os.path.exists(db_file):
                            os.remove(db_file)
                except OSError as e:
                    logger.warning(f"Error cleaning up temporary files: {e}")

    def process_alignments(self, alignments, reference):

        ref_length = len(reference.seq)
        # 记录：query_id -> [一个字符数组（长 = ref_length），初始全为 '-']
        query_aln_dict = {}

        # 还要记录每个 query、每个 ref-position 上的分数 (便于重叠时优先保留)
        query_score_dict = {}

        # 初始化
        all_query_ids = set([aln.query_id for aln in alignments])
        for qid in all_query_ids:
            query_aln_dict[qid] = ['-' for _ in range(ref_length)]
            query_score_dict[qid] = [0.0 for _ in range(ref_length)]

        # 依次处理每条比对
        for aln in alignments:
            qid = aln.query_id
            # 如果方向是 'minus'，要先反向互补
            if aln.orientation == 'plus':
                qseq = aln.alignment[0]
            else:
                qseq = str(Seq(aln.alignment[0]).reverse_complement())

            # subject_start, subject_end 基于 1-based 坐标
            # Python list 是 0-based，所以要注意减1
            subj_start = aln.subject_start - 1
            subj_end   = aln.subject_end   - 1

            # 注意：如果 qseq 长度和 (subj_end-subj_start+1) 不一致，要做一些保护
            # 这里假设它们相等，即无内部 gap
            aligned_len = len(qseq)
            expected_len = (subj_end - subj_start + 1)
            if aligned_len != expected_len:
                # 如果不匹配，你需要做额外处理或跳过这个对齐
                continue

            # 现在把 qseq 的每个碱基放到 query_aln_dict[qid] 的正确位置上
            # 并比较 alignment score，若当前位置已有碱基且新对齐分数更高，就替换
            for i in range(aligned_len):
                ref_pos = subj_start + i
                base    = qseq[i]
                # 如果当前对齐分数更高，就替换
                if aln.score > query_score_dict[qid][ref_pos]:
                    query_aln_dict[qid][ref_pos] = base
                    query_score_dict[qid][ref_pos] = aln.score

        # 最后，把每个 query_id 的字符数组变成字符串
        final_seqs = []
        for qid in all_query_ids:
            merged_seq = ''.join(query_aln_dict[qid])
            final_seqs.append(merged_seq)

        return final_seqs

    def build_consensus(self, aligned_sequences):
        """从对齐序列构建共识序列"""
        if not aligned_sequences:
            return ""
        
        # 确保所有序列长度相同
        seq_length = len(aligned_sequences[0])
        consensus = []
        
        # 对每个位置计算最常见的碱基
        for i in range(seq_length):
            bases = [seq[i] for seq in aligned_sequences if i < len(seq)]
            base_counts = Counter(base for base in bases if base != '-')
            
            if base_counts:
                # 选择最频繁的碱基
                consensus.append(max(base_counts.items(), key=lambda x: x[1])[0])
            else:
                consensus.append('-')
                
        return ''.join(consensus).replace('-', '')

    def build_te_consensus(self, input_file, output_file):
        """主要的共识序列构建流程"""
        try:
            # 1. 读取序列
            logger.info("Reading sequences...")
            sequences = list(SeqIO.parse(input_file, "fasta"))
            if not sequences:
                raise ValueError(f"No sequences found in {input_file}")
            logger.info(f"Read {len(sequences)} sequences")

            # 2. 找到最佳参考序列
            logger.info("Finding best reference sequence...")
            reference = self.find_best_reference(sequences)
            logger.info(f"Selected {reference.id} as reference sequence")

            # 3. 基于参考序列构建多序列比对
            logger.info("Building multiple alignment...")
            aligned_seqs = self.build_multiple_alignment(sequences, reference)
            logger.info("Completed alignment")

            # 4. 构建共识序列
            logger.info("Building consensus sequence...")
            consensus = self.build_consensus(aligned_seqs)
            logger.info(f"Generated consensus sequence of length {len(consensus)}")

            # 5. 保存结果
            consensus_record = SeqRecord(
                Seq(consensus),
                id=os.path.splitext(os.path.basename(input_file))[0],
                description=f"consensus from {len(sequences)} sequences"
            )
            SeqIO.write(consensus_record, output_file, "fasta")
            logger.info(f"Written consensus to {output_file}")

            # 6. 生成报告
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
        
        # 查找可执行文件
        self.rmblastn = shutil.which('rmblastn')
        self.makeblastdb = shutil.which('makeblastdb')
        
        if not self.rmblastn:
            raise FileNotFoundError("rmblastn not found in PATH")
        if not self.makeblastdb:
            raise FileNotFoundError("makeblastdb not found in PATH")
            
        # 从 rmblastn 路径获取基础目录
        conda_env_dir = os.path.dirname(os.path.dirname(self.rmblastn))
        
        # 设置常见的矩阵文件位置
        possible_matrix_paths = [
            # Conda 环境中的标准位置
            os.path.join(conda_env_dir, 'share/RepeatModeler/Matrices/ncbi/nt/comparison.matrix'),
            os.path.join(conda_env_dir, 'share/RepeatMasker/Libraries/Dfam.hmm'),
            # 系统安装的标准位置
            '/usr/share/RepeatMasker/Matrices/nt',
            '/usr/local/share/RepeatMasker/Matrices/nt',
            # 如果都找不到，使用默认的 BLOSUM62
            os.path.join(conda_env_dir, 'share/RepeatMasker/Matrices/BLOSUM62')
        ]
        
        # 查找矩阵文件
        self.matrix_path = None
        for path in possible_matrix_paths:
            if os.path.exists(path):
                self.matrix_path = path
                break
                
        if self.matrix_path is None:
            # 如果找不到任何矩阵文件，使用默认的 BLOSUM62
            self.matrix_path = 'BLOSUM62'
            
        logger.info(f"Using scoring matrix: {self.matrix_path}")

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
        # 获取配置
        config = Config()
        
        # 注意这里的参数名已修改为正确的形式
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
