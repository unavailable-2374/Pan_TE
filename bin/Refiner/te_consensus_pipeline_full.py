#!/usr/bin/env python3
"""
TE Consensus Pipeline - 完整集成版本
包含所有功能的单脚本版本，便于集成到其他工作流

基于优化的copy-based clustering算法，包含完整的错误处理和恢复机制
"""

import os
import sys
import logging
import argparse
import json
import pickle
import tempfile
import shutil
import subprocess
import gc
import hashlib
import re
import threading
import queue
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional, Set
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import time
from datetime import datetime
from functools import wraps, partial

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# ==================== 配置和常量 ====================

# IUPAC混合碱基编码
IUPAC_CODES = {
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
    'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
    'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

class PipelineConfig:
    """完整的流程配置类"""
    def __init__(self, repeatscout_file: str, genome_file: str, output_dir: str, **kwargs):
        # 基本路径
        self.repeatscout_file = repeatscout_file
        self.genome_file = genome_file
        self.output_dir = output_dir
        
        # Phase 1 筛选参数
        self.min_length = kwargs.get('min_length', 80)
        self.max_length = kwargs.get('max_length', 20000)
        self.max_n_percent = kwargs.get('max_n_percent', 0.2)
        self.gc_min = kwargs.get('gc_min', 0.2)
        self.gc_max = kwargs.get('gc_max', 0.8)
        self.dust_threshold = kwargs.get('dust_threshold', 7)
        self.entropy_threshold = kwargs.get('entropy_threshold', 1.0)
        
        # Phase 2 聚类参数
        self.min_copy_number = kwargs.get('min_copy_number', 5)
        self.identity_threshold = kwargs.get('identity_threshold', 0.85)
        self.coverage_threshold = kwargs.get('coverage_threshold', 0.60)
        self.max_recruits_per_family = kwargs.get('max_recruits_per_family', 30)
        self.column_coverage_threshold = kwargs.get('column_coverage_threshold', 0.3)
        
        # Phase 3 过滤参数
        self.redundancy_threshold_masking = kwargs.get('redundancy_threshold_masking', 0.95)
        self.redundancy_threshold_analysis = kwargs.get('redundancy_threshold_analysis', 0.90)
        self.chimera_threshold = kwargs.get('chimera_threshold', 0.4)
        self.boundary_quality_threshold = kwargs.get('boundary_quality_threshold', 0.6)
        
        # 性能参数
        self.threads = kwargs.get('threads', 8)
        self.batch_size = kwargs.get('batch_size', 100)
        self.max_retries = kwargs.get('max_retries', 3)
        self.retry_delay = kwargs.get('retry_delay', 5)
        self.use_parallel = kwargs.get('use_parallel', True)
        self.memory_limit = kwargs.get('memory_limit', None)
        
        # 外部工具路径
        self.repeatmasker_path = kwargs.get('repeatmasker_path', 'RepeatMasker')
        self.mafft_path = kwargs.get('mafft_path', 'mafft')
        self.blast_path = kwargs.get('blast_path', 'blastn')
        self.cdhit_path = kwargs.get('cdhit_path', 'cd-hit-est')
        
        # 目录设置
        self.temp_dir = os.path.join(output_dir, 'temp_work')
        self.cache_dir = os.path.join(output_dir, 'cache')
        self.checkpoint_dir = os.path.join(output_dir, 'checkpoints')
        
        # 其他选项
        self.keep_temp = kwargs.get('keep_temp', False)
        self.keep_checkpoints = kwargs.get('keep_checkpoints', False)
        self.resume = kwargs.get('resume', False)
        self.clear_cache = kwargs.get('clear_cache', False)
        self.verbose = kwargs.get('verbose', False)
        
        # 创建必要目录
        self._create_directories()
        
        # 清理选项
        if self.clear_cache and os.path.exists(self.cache_dir):
            shutil.rmtree(self.cache_dir)
            os.makedirs(self.cache_dir)
    
    def _create_directories(self):
        """创建必要的目录"""
        for dir_path in [self.output_dir, self.temp_dir, self.cache_dir, self.checkpoint_dir]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    def save(self, filepath: str):
        """保存配置到JSON文件"""
        config_dict = {k: v for k, v in self.__dict__.items() if not k.startswith('_')}
        with open(filepath, 'w') as f:
            json.dump(config_dict, f, indent=2)
    
    @classmethod
    def load(cls, filepath: str):
        """从JSON文件加载配置"""
        with open(filepath, 'r') as f:
            config_dict = json.load(f)
        return cls(**config_dict)

# ==================== 日志和错误处理 ====================

def setup_logging(output_dir: str, verbose: bool = False):
    """设置日志系统"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    # 配置根日志记录器
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(os.path.join(output_dir, 'pipeline.log')),
            logging.StreamHandler()
        ]
    )
    
    # 减少第三方库的日志输出
    logging.getLogger('Bio').setLevel(logging.WARNING)
    
    return logging.getLogger(__name__)

class RobustRunner:
    """错误处理和重试机制"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def run_with_retry(self, func, *args, **kwargs):
        """带重试的函数执行"""
        max_retries = self.config.max_retries
        retry_delay = self.config.retry_delay
        
        for attempt in range(max_retries):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if attempt < max_retries - 1:
                    self.logger.warning(f"Attempt {attempt + 1} failed: {e}. Retrying in {retry_delay} seconds...")
                    time.sleep(retry_delay * (2 ** attempt))  # 指数退避
                else:
                    self.logger.error(f"All {max_retries} attempts failed.")
                    raise
    
    def run_with_checkpoint(self, func, checkpoint_name: str, *args, **kwargs):
        """带检查点的函数执行"""
        checkpoint_file = os.path.join(self.config.checkpoint_dir, f"{checkpoint_name}.pkl")
        
        # 尝试从检查点恢复
        if self.config.resume and os.path.exists(checkpoint_file):
            self.logger.info(f"Resuming from checkpoint: {checkpoint_name}")
            with open(checkpoint_file, 'rb') as f:
                return pickle.load(f)
        
        # 执行函数
        result = self.run_with_retry(func, *args, **kwargs)
        
        # 保存检查点
        if self.config.keep_checkpoints:
            with open(checkpoint_file, 'wb') as f:
                pickle.dump(result, f)
            self.logger.debug(f"Checkpoint saved: {checkpoint_name}")
        
        return result

# ==================== 缓存系统 ====================

def cache_result(cache_dir: str):
    """结果缓存装饰器"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # 生成缓存键
            cache_key = hashlib.md5(
                f"{func.__name__}_{str(args)}_{str(kwargs)}".encode()
            ).hexdigest()
            cache_file = os.path.join(cache_dir, f"{func.__name__}_{cache_key}.pkl")
            
            # 尝试从缓存加载
            if os.path.exists(cache_file):
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)
            
            # 执行函数并缓存结果
            result = func(*args, **kwargs)
            os.makedirs(cache_dir, exist_ok=True)
            with open(cache_file, 'wb') as f:
                pickle.dump(result, f)
            
            return result
        return wrapper
    return decorator

# ==================== 序列处理工具 ====================

def load_sequences(fasta_file: str) -> List[Dict]:
    """加载FASTA序列"""
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append({
            'id': str(record.id),
            'sequence': str(record.seq).upper(),
            'description': str(record.description),
            'length': len(record.seq)
        })
    return sequences

def save_sequences(sequences: List[Dict], output_file: str):
    """保存序列到FASTA文件"""
    records = []
    for seq in sequences:
        # 构建描述信息
        desc_parts = [seq.get('id', 'unknown')]
        if 'copy_number' in seq:
            desc_parts.append(f"copies={seq['copy_number']}")
        if 'avg_identity' in seq:
            desc_parts.append(f"identity={seq['avg_identity']:.1f}")
        if 'quality_score' in seq:
            desc_parts.append(f"quality={seq['quality_score']:.2f}")
        
        record = SeqRecord(
            Seq(seq['sequence']),
            id=seq['id'],
            description=' '.join(desc_parts)
        )
        records.append(record)
    
    SeqIO.write(records, output_file, "fasta")

def extract_sequence_from_genome(genome_file: str, chrom: str, start: int, end: int, strand: str = '+') -> str:
    """从基因组提取序列"""
    try:
        # 使用samtools或其他工具提取序列
        # 这里使用简化版本，实际应该使用索引访问
        genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
        
        if chrom in genome_dict:
            seq = str(genome_dict[chrom].seq[start-1:end]).upper()
            if strand == '-':
                seq = str(Seq(seq).reverse_complement())
            return seq
    except Exception as e:
        logging.getLogger(__name__).error(f"Error extracting sequence: {e}")
    
    return ""

def calculate_identity(seq1: str, seq2: str) -> float:
    """计算两个序列的相似度"""
    if not seq1 or not seq2:
        return 0.0
    
    # 简单的全局比对相似度
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    
    if min_len == 0:
        return 0.0
    
    # 滑动窗口找最佳匹配
    best_identity = 0
    for offset in range(-min_len//2, min_len//2):
        matches = 0
        comparisons = 0
        
        for i in range(max_len):
            pos1 = i
            pos2 = i + offset
            
            if 0 <= pos1 < len(seq1) and 0 <= pos2 < len(seq2):
                comparisons += 1
                if seq1[pos1] == seq2[pos2]:
                    matches += 1
        
        if comparisons > 0:
            identity = matches / comparisons
            best_identity = max(best_identity, identity)
    
    return best_identity

def detect_tsd(genome_file: str, chrom: str, start: int, end: int, window: int = 20) -> Optional[Dict]:
    """检测Target Site Duplication"""
    try:
        # 提取侧翼序列
        upstream = extract_sequence_from_genome(genome_file, chrom, start - window, start)
        downstream = extract_sequence_from_genome(genome_file, chrom, end, end + window)
        
        if not upstream or not downstream:
            return None
        
        # 查找TSD
        for tsd_len in range(2, min(15, window)):
            upstream_tsd = upstream[-tsd_len:]
            downstream_tsd = downstream[:tsd_len]
            
            if upstream_tsd == downstream_tsd:
                return {
                    'sequence': upstream_tsd,
                    'length': tsd_len,
                    'type': 'perfect'
                }
        
        # 查找近似TSD
        for tsd_len in range(4, min(10, window)):
            upstream_tsd = upstream[-tsd_len:]
            downstream_tsd = downstream[:tsd_len]
            
            identity = calculate_identity(upstream_tsd, downstream_tsd)
            if identity >= 0.8:
                return {
                    'sequence': upstream_tsd,
                    'length': tsd_len,
                    'type': 'imperfect',
                    'identity': identity
                }
        
    except Exception as e:
        logging.getLogger(__name__).debug(f"TSD detection failed: {e}")
    
    return None

# ==================== 复杂度计算 ====================

def calculate_dust_score(sequence: str, window_size: int = 64) -> float:
    """计算DUST复杂度分数"""
    if len(sequence) < 3:
        return 0
    
    sequence = sequence.upper()
    total_score = 0
    num_windows = 0
    
    for i in range(0, len(sequence) - window_size + 1, window_size // 2):
        window = sequence[i:i + window_size]
        
        # 计算三联体频率
        triplet_counts = defaultdict(int)
        for j in range(len(window) - 2):
            triplet = window[j:j+3]
            if 'N' not in triplet:
                triplet_counts[triplet] += 1
        
        # 计算窗口分数
        window_score = 0
        total_triplets = len(window) - 2
        
        for count in triplet_counts.values():
            if count > 1:
                window_score += (count * (count - 1)) / 2
        
        if total_triplets > 0:
            window_score = window_score / (total_triplets / 2)
            total_score += min(window_score, 100)
            num_windows += 1
    
    if num_windows > 0:
        return total_score / num_windows
    
    return 0

def calculate_shannon_entropy(sequence: str) -> float:
    """计算Shannon熵"""
    if not sequence:
        return 0
    
    sequence = sequence.upper()
    base_counts = Counter(c for c in sequence if c in 'ACGT')
    total = sum(base_counts.values())
    
    if total == 0:
        return 0
    
    entropy = 0
    for count in base_counts.values():
        if count > 0:
            prob = count / total
            entropy -= prob * np.log2(prob)
    
    # 归一化到[0, 1]
    return entropy / 2.0

def calculate_complexity_scores(sequence: str) -> Dict[str, float]:
    """计算多维复杂度分数"""
    scores = {
        'dust': calculate_dust_score(sequence),
        'entropy': calculate_shannon_entropy(sequence),
        'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence) if sequence else 0,
        'n_percent': sequence.count('N') / len(sequence) if sequence else 1
    }
    
    # 计算综合分数
    scores['combined'] = (
        (100 - scores['dust']) / 100 * 0.4 +
        scores['entropy'] * 0.3 +
        (1 - abs(scores['gc_content'] - 0.5) * 2) * 0.2 +
        (1 - scores['n_percent']) * 0.1
    )
    
    return scores

# ==================== RepeatMasker处理 ====================

def run_repeatmasker_single(seed_seq: Dict, genome_file: str, config: PipelineConfig, params: Dict = None) -> List[Dict]:
    """对单个序列运行RepeatMasker"""
    logger = logging.getLogger(__name__)
    
    # 创建临时文件
    temp_lib = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False, dir=config.temp_dir)
    temp_lib.write(f">{seed_seq['id']}\n{seed_seq['sequence']}\n")
    temp_lib.close()
    
    temp_genome = tempfile.NamedTemporaryFile(suffix='.fa', delete=False, dir=config.temp_dir)
    temp_genome.close()
    shutil.copy(genome_file, temp_genome.name)
    
    try:
        # 构建命令
        cmd = [
            config.repeatmasker_path,
            '-lib', temp_lib.name,
            '-pa', str(max(1, config.threads // 4))
        ]
        
        # 添加参数
        if params:
            if params.get('s'):
                cmd.append('-s')
            if params.get('q'):
                cmd.append('-q')
            if params.get('no_is'):
                cmd.append('-no_is')
            if params.get('nolow'):
                cmd.append('-nolow')
            if params.get('div'):
                cmd.extend(['-div', str(params['div'])])
            if params.get('cutoff'):
                cmd.extend(['-cutoff', str(params['cutoff'])])
        else:
            cmd.extend(['-nolow', '-no_is', '-norna', '-div', '25'])
        
        cmd.append(temp_genome.name)
        
        # 运行RepeatMasker
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=config.temp_dir)
        
        # 解析输出
        hits = []
        out_file = f"{temp_genome.name}.out"
        
        if os.path.exists(out_file):
            with open(out_file, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith(' '):
                        parts = line.strip().split()
                        if len(parts) >= 14:
                            try:
                                hit = {
                                    'score': int(parts[0]) if parts[0].isdigit() else 0,
                                    'divergence': float(parts[1]) if parts[1] != '*' else 0,
                                    'identity': 100 - float(parts[1]) if parts[1] != '*' else 0,
                                    'chrom': parts[4],
                                    'start': int(parts[5]),
                                    'end': int(parts[6]),
                                    'strand': '+' if parts[8] == '+' else '-',
                                    'repeat_class': parts[10] if len(parts) > 10 else 'Unknown'
                                }
                                hits.append(hit)
                            except (ValueError, IndexError):
                                continue
        
        return hits
        
    except Exception as e:
        logger.error(f"RepeatMasker failed: {e}")
        return []
    
    finally:
        # 清理临时文件
        for f in [temp_lib.name, temp_genome.name, f"{temp_genome.name}.out", 
                  f"{temp_genome.name}.cat", f"{temp_genome.name}.masked"]:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except:
                    pass

def run_repeatmasker_batch(sequences: List[Dict], genome_file: str, config: PipelineConfig) -> Dict:
    """批量运行RepeatMasker"""
    logger = logging.getLogger(__name__)
    results = {}
    
    # 创建临时库文件
    temp_lib = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False, dir=config.temp_dir)
    for seq in sequences:
        temp_lib.write(f">{seq['id']}\n{seq['sequence']}\n")
    temp_lib.close()
    
    try:
        # 运行RepeatMasker
        cmd = [
            config.repeatmasker_path,
            '-lib', temp_lib.name,
            '-pa', str(config.threads),
            '-nolow', '-no_is', '-norna',
            '-div', '25',
            genome_file
        ]
        
        logger.info(f"Running RepeatMasker with {len(sequences)} sequences...")
        subprocess.run(cmd, capture_output=True, text=True, cwd=config.temp_dir)
        
        # 解析结果
        out_file = f"{genome_file}.out"
        if os.path.exists(out_file):
            with open(out_file, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith(' '):
                        parts = line.strip().split()
                        if len(parts) >= 14:
                            try:
                                seq_id = parts[9]
                                if seq_id not in results:
                                    results[seq_id] = {'hits': [], 'total_bp': 0}
                                
                                hit = {
                                    'score': int(parts[0]) if parts[0].isdigit() else 0,
                                    'divergence': float(parts[1]) if parts[1] != '*' else 0,
                                    'identity': 100 - float(parts[1]) if parts[1] != '*' else 0,
                                    'chrom': parts[4],
                                    'start': int(parts[5]),
                                    'end': int(parts[6]),
                                    'strand': '+' if parts[8] == '+' else '-'
                                }
                                results[seq_id]['hits'].append(hit)
                                results[seq_id]['total_bp'] += hit['end'] - hit['start']
                            except (ValueError, IndexError):
                                continue
    
    finally:
        if os.path.exists(temp_lib.name):
            os.remove(temp_lib.name)
    
    return results

# ==================== MAFFT处理 ====================

def run_mafft(sequences: List[Dict], config: PipelineConfig, algorithm: str = 'auto') -> Optional[Dict]:
    """运行MAFFT多序列比对"""
    logger = logging.getLogger(__name__)
    
    if len(sequences) < 2:
        return None
    
    # 创建临时文件
    temp_input = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False, dir=config.temp_dir)
    for seq in sequences:
        temp_input.write(f">{seq.get('id', 'seq')}\n{seq.get('sequence', '')}\n")
    temp_input.close()
    
    temp_output = tempfile.NamedTemporaryFile(suffix='.aln', delete=False, dir=config.temp_dir)
    temp_output.close()
    
    try:
        # 选择算法
        if algorithm == 'auto':
            if len(sequences) > 200:
                algo_params = '--retree 1'
            elif len(sequences) > 100:
                algo_params = '--retree 2'
            else:
                algo_params = '--localpair --maxiterate 1000'
        elif algorithm == 'localpair':
            algo_params = '--localpair --maxiterate 1000'
        elif algorithm == 'genafpair':
            algo_params = '--genafpair --maxiterate 1000'
        else:
            algo_params = '--auto'
        
        # 运行MAFFT
        cmd = f"{config.mafft_path} {algo_params} --thread {config.threads} --quiet {temp_input.name} > {temp_output.name}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"MAFFT failed: {result.stderr}")
            return None
        
        # 解析结果
        alignment = {}
        for record in SeqIO.parse(temp_output.name, "fasta"):
            alignment[str(record.id)] = str(record.seq)
        
        return alignment
        
    except Exception as e:
        logger.error(f"MAFFT error: {e}")
        return None
    
    finally:
        for f in [temp_input.name, temp_output.name]:
            if os.path.exists(f):
                os.remove(f)

def build_consensus_from_alignment(alignment: Dict, min_coverage: float = 0.3, use_iupac: bool = True) -> str:
    """从多序列比对构建共识序列"""
    if not alignment:
        return ""
    
    # 获取对齐长度
    aln_length = len(next(iter(alignment.values())))
    consensus = []
    
    for pos in range(aln_length):
        # 统计每个位置的碱基
        bases = []
        for seq in alignment.values():
            if pos < len(seq) and seq[pos] not in '-N':
                bases.append(seq[pos].upper())
        
        # 需要足够的覆盖度
        if len(bases) >= len(alignment) * min_coverage:
            base_counts = Counter(bases)
            total = sum(base_counts.values())
            
            # 选择共识碱基
            if use_iupac and len(base_counts) > 1:
                # 使用IUPAC混合碱基
                consensus_base = get_iupac_consensus(base_counts, total)
            else:
                # 使用最常见的碱基
                consensus_base = base_counts.most_common(1)[0][0]
            
            consensus.append(consensus_base)
    
    return ''.join(consensus)

def get_iupac_consensus(base_counts: Counter, total: int, threshold: float = 0.75) -> str:
    """获取IUPAC混合碱基"""
    # 如果一个碱基占主导地位
    for base, count in base_counts.items():
        if count / total >= threshold:
            return base
    
    # 检查混合碱基
    bases_present = set(base_counts.keys())
    
    for iupac, bases in IUPAC_CODES.items():
        if bases_present == set(bases):
            return iupac
    
    # 默认使用N
    return 'N'

# ==================== Phase 1: 完整的筛选和评分 ====================

class SequenceScreenerFull:
    """完整的Phase 1实现"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.runner = RobustRunner(config)
        self.logger = logging.getLogger(__name__)
    
    def run(self) -> Dict:
        """执行Phase 1"""
        self.logger.info("="*60)
        self.logger.info("Phase 1: Sequence Screening and Scoring")
        self.logger.info("="*60)
        
        # 加载序列
        sequences = load_sequences(self.config.repeatscout_file)
        self.logger.info(f"Loaded {len(sequences)} sequences from RepeatScout")
        
        # 基本过滤
        filtered_sequences = self._basic_filtering(sequences)
        
        # 计算复杂度分数
        self._calculate_complexity_scores(filtered_sequences)
        
        # RepeatMasker评估
        rm_results = self._repeatmasker_evaluation(filtered_sequences)
        
        # 综合评分和分类
        classified = self._classify_sequences(filtered_sequences, rm_results)
        
        # 生成摘要
        summary = (f"Total: {len(sequences)}, "
                  f"Filtered: {len(filtered_sequences)}, "
                  f"A: {len(classified['a_sequences'])}, "
                  f"B: {len(classified['b_sequences'])}, "
                  f"C: {len(classified['c_sequences'])}")
        
        classified['summary'] = summary
        self.logger.info(f"Phase 1 complete: {summary}")
        
        return classified
    
    def _basic_filtering(self, sequences: List[Dict]) -> List[Dict]:
        """基本长度和质量过滤"""
        filtered = []
        
        for seq in sequences:
            seq_len = len(seq['sequence'])
            n_percent = seq['sequence'].count('N') / seq_len if seq_len > 0 else 1
            gc_content = (seq['sequence'].count('G') + seq['sequence'].count('C')) / seq_len if seq_len > 0 else 0
            
            if (self.config.min_length <= seq_len <= self.config.max_length and
                n_percent <= self.config.max_n_percent and
                self.config.gc_min <= gc_content <= self.config.gc_max):
                filtered.append(seq)
        
        self.logger.info(f"Basic filtering: {len(filtered)}/{len(sequences)} sequences passed")
        return filtered
    
    def _calculate_complexity_scores(self, sequences: List[Dict]):
        """计算复杂度分数"""
        self.logger.info("Calculating complexity scores...")
        
        for seq in sequences:
            scores = calculate_complexity_scores(seq['sequence'])
            seq.update(scores)
    
    def _repeatmasker_evaluation(self, sequences: List[Dict]) -> Dict:
        """RepeatMasker评估"""
        self.logger.info("Running RepeatMasker evaluation...")
        
        # 分批处理以避免内存问题
        batch_size = min(100, len(sequences))
        all_results = {}
        
        for i in range(0, len(sequences), batch_size):
            batch = sequences[i:i+batch_size]
            batch_results = run_repeatmasker_batch(batch, self.config.genome_file, self.config)
            all_results.update(batch_results)
        
        return all_results
    
    def _classify_sequences(self, sequences: List[Dict], rm_results: Dict) -> Dict:
        """序列分类"""
        a_sequences = []
        b_sequences = []
        c_sequences = []
        scores = {}
        
        for seq in sequences:
            # 获取各项分数
            complexity_score = seq.get('combined', 0)
            
            # RepeatMasker覆盖度分数
            rm_data = rm_results.get(seq['id'], {})
            total_bp = rm_data.get('total_bp', 0)
            num_hits = len(rm_data.get('hits', []))
            
            # 归一化覆盖度
            coverage_score = min(total_bp / 50000, 1.0)  # 50kb作为参考
            hit_score = min(num_hits / 100, 1.0)  # 100个hits作为参考
            
            # 综合评分
            final_score = (
                complexity_score * 0.4 +
                coverage_score * 0.3 +
                hit_score * 0.3
            )
            
            # 保存分数
            scores[seq['id']] = {
                'complexity': complexity_score,
                'coverage': coverage_score,
                'hits': hit_score,
                'final': final_score
            }
            
            seq['final_score'] = final_score
            seq['rm_hits'] = rm_data.get('hits', [])
            
            # 分类
            if final_score >= 0.65:
                seq['class'] = 'A'
                a_sequences.append(seq)
            elif final_score >= 0.45:
                seq['class'] = 'B'
                b_sequences.append(seq)
            else:
                seq['class'] = 'C'
                c_sequences.append(seq)
        
        self.logger.info(f"Classification complete: A={len(a_sequences)}, B={len(b_sequences)}, C={len(c_sequences)}")
        
        return {
            'a_sequences': a_sequences,
            'b_sequences': b_sequences,
            'c_sequences': c_sequences,
            'rm_detailed_results': rm_results,
            'scores': scores
        }

# ==================== Phase 2: 完整的共识构建 ====================

class ConsensusBuilderFull:
    """完整的Phase 2实现 - 基于拷贝的聚类"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.runner = RobustRunner(config)
        self.logger = logging.getLogger(__name__)
        
        # 设置并行worker数
        cpu_count = mp.cpu_count()
        if config.threads >= 32:
            self.max_workers = min(config.threads, cpu_count * 2)
        elif config.threads >= 16:
            self.max_workers = min(config.threads, max(cpu_count, 16))
        else:
            self.max_workers = max(4, min(config.threads, cpu_count))
        
        self.max_workers = max(4, min(self.max_workers, 128))
        self.logger.info(f"Phase 2 initialized with {self.max_workers} workers")
    
    def run(self, phase1_output: Dict) -> List[Dict]:
        """执行Phase 2"""
        self.logger.info("="*60)
        self.logger.info("Phase 2: Consensus Building with Copy-based Clustering")
        self.logger.info("="*60)
        
        # 选择候选序列
        candidates = self._select_candidates(phase1_output)
        
        # 按复杂度排序（优化并行处理）
        sorted_candidates = self._sort_by_complexity(candidates)
        
        # 并行处理序列
        all_consensus = self._process_sequences_parallel(sorted_candidates)
        
        self.logger.info(f"Phase 2 complete: Generated {len(all_consensus)} consensus sequences")
        return all_consensus
    
    def _select_candidates(self, phase1_output: Dict) -> List[Dict]:
        """选择候选序列"""
        candidates = []
        rm_results = phase1_output.get('rm_detailed_results', {})
        scores = phase1_output.get('scores', {})
        
        # 所有A类序列
        for seq in phase1_output['a_sequences']:
            seq['rm_hits'] = rm_results.get(seq['id'], {}).get('hits', [])
            seq['quality_class'] = 'A'
            candidates.append(seq)
        
        # 高分B类序列
        for seq in phase1_output['b_sequences']:
            if scores.get(seq['id'], {}).get('final', 0) >= 0.50:
                seq['rm_hits'] = rm_results.get(seq['id'], {}).get('hits', [])
                seq['quality_class'] = 'B'
                candidates.append(seq)
        
        self.logger.info(f"Selected {len(candidates)} candidate sequences for consensus building")
        return candidates
    
    def _sort_by_complexity(self, sequences: List[Dict]) -> List[Dict]:
        """按处理复杂度排序"""
        for seq in sequences:
            seq_length = len(seq.get('sequence', ''))
            estimated_copies = len(seq.get('rm_hits', []))
            
            # 复杂度分数
            complexity_score = np.log1p(seq_length) * np.log1p(estimated_copies)
            if seq.get('quality_class') == 'A':
                complexity_score *= 1.2
            
            seq['_complexity_score'] = complexity_score
        
        sorted_seqs = sorted(sequences, key=lambda x: x.get('_complexity_score', 0), reverse=True)
        
        if sorted_seqs:
            self.logger.info(f"Complexity sorting: highest={sorted_seqs[0]['_complexity_score']:.2f}, "
                           f"lowest={sorted_seqs[-1]['_complexity_score']:.2f}")
        
        return sorted_seqs
    
    def _process_sequences_parallel(self, sequences: List[Dict]) -> List[Dict]:
        """并行处理序列"""
        all_consensus = []
        total = len(sequences)
        completed = 0
        
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            # 提交所有任务
            future_to_seq = {}
            for seq in sequences:
                future = executor.submit(
                    process_sequence_full,
                    seq,
                    self.config.genome_file,
                    self.config
                )
                future_to_seq[future] = seq
            
            # 收集结果
            for future in as_completed(future_to_seq):
                completed += 1
                seq = future_to_seq[future]
                
                try:
                    consensus_list = future.result(timeout=300)
                    if consensus_list:
                        all_consensus.extend(consensus_list)
                        self.logger.debug(f"[{completed}/{total}] {seq['id']}: "
                                        f"generated {len(consensus_list)} consensus")
                except Exception as e:
                    self.logger.error(f"[{completed}/{total}] Failed {seq['id']}: {e}")
                
                # 进度报告
                if completed % 10 == 0 or completed == total:
                    self.logger.info(f"Progress: {completed}/{total} sequences processed, "
                                   f"{len(all_consensus)} consensus generated")
        
        return all_consensus

def process_sequence_full(seed_seq: Dict, genome_file: str, config: PipelineConfig) -> List[Dict]:
    """处理单个序列 - 完整版本"""
    logger = logging.getLogger(__name__)
    
    try:
        # 1. 获取基因组拷贝
        if 'rm_hits' in seed_seq and seed_seq['rm_hits']:
            # 使用Phase 1的结果
            genome_copies = extract_genome_copies_from_hits(seed_seq['rm_hits'], genome_file, config)
        else:
            # 重新运行RepeatMasker
            rm_params = {
                's': True,
                'no_is': True,
                'nolow': True,
                'div': 25,
                'cutoff': 200
            }
            rm_hits = run_repeatmasker_single(seed_seq, genome_file, config, params=rm_params)
            genome_copies = extract_genome_copies_from_hits(rm_hits, genome_file, config)
        
        # 检查拷贝数
        is_a_class = seed_seq.get('quality_class') == 'A'
        min_copies = max(2, config.min_copy_number // 2) if is_a_class else config.min_copy_number
        
        if not genome_copies or len(genome_copies) < min_copies:
            logger.debug(f"Insufficient copies for {seed_seq['id']}: {len(genome_copies)}")
            return []
        
        # 2. 聚类拷贝
        subfamilies = cluster_copies_full(genome_copies, config)
        
        if not subfamilies:
            return []
        
        # 3. 构建共识
        consensus_list = []
        for subfamily_id, subfamily_copies in enumerate(subfamilies):
            min_subfamily_size = 1 if is_a_class else max(2, config.min_copy_number // 2)
            
            if len(subfamily_copies) < min_subfamily_size:
                continue
            
            consensus = build_subfamily_consensus_full(
                seed_seq,
                subfamily_copies,
                subfamily_id,
                genome_file,
                config
            )
            
            if consensus:
                consensus_list.append(consensus)
        
        return consensus_list
        
    except Exception as e:
        logger.error(f"Error processing {seed_seq['id']}: {e}")
        return []

def extract_genome_copies_from_hits(rm_hits: List[Dict], genome_file: str, config: PipelineConfig) -> List[Dict]:
    """从RepeatMasker hits提取基因组拷贝"""
    copies = []
    max_hits = min(200, len(rm_hits))
    
    for hit in rm_hits[:max_hits]:
        try:
            # 提取序列
            copy_seq = extract_sequence_from_genome(
                genome_file,
                hit['chrom'],
                hit['start'],
                hit['end'],
                hit.get('strand', '+')
            )
            
            if copy_seq and 30 <= len(copy_seq) <= 50000:
                # 检测TSD
                tsd = detect_tsd(genome_file, hit['chrom'], hit['start'], hit['end'])
                
                copies.append({
                    'sequence': copy_seq,
                    'chrom': hit['chrom'],
                    'start': hit['start'],
                    'end': hit['end'],
                    'strand': hit.get('strand', '+'),
                    'identity': hit.get('identity', 0),
                    'score': hit.get('score', 0),
                    'length': len(copy_seq),
                    'has_tsd': tsd is not None,
                    'tsd': tsd
                })
        except Exception:
            continue
    
    return copies

def cluster_copies_full(copies: List[Dict], config: PipelineConfig) -> List[List[Dict]]:
    """完整的拷贝聚类"""
    if len(copies) < 2:
        return [copies]
    
    try:
        # 小数据集用完整聚类
        if len(copies) <= 20:
            return hierarchical_clustering(copies, config)
        
        # 大数据集用混合策略
        return hybrid_clustering(copies, config)
        
    except Exception as e:
        logging.getLogger(__name__).error(f"Clustering failed: {e}")
        return [copies]

def hierarchical_clustering(copies: List[Dict], config: PipelineConfig) -> List[List[Dict]]:
    """层次聚类"""
    n = len(copies)
    distance_matrix = np.zeros((n, n))
    
    # 计算距离矩阵
    for i in range(n):
        for j in range(i+1, n):
            identity = calculate_identity(copies[i]['sequence'], copies[j]['sequence'])
            distance = 1.0 - identity
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance
    
    # 聚类
    condensed = squareform(distance_matrix)
    linkage_matrix = linkage(condensed, method='average')
    
    threshold = 1.0 - config.identity_threshold
    clusters = fcluster(linkage_matrix, threshold, criterion='distance')
    
    # 组织结果
    subfamilies = defaultdict(list)
    for i, cluster_id in enumerate(clusters):
        subfamilies[cluster_id].append(copies[i])
    
    return list(subfamilies.values())

def hybrid_clustering(copies: List[Dict], config: PipelineConfig) -> List[List[Dict]]:
    """混合聚类策略"""
    # 按长度预分组
    length_groups = defaultdict(list)
    for copy in copies:
        length_key = len(copy['sequence']) // 200 * 200
        length_groups[length_key].append(copy)
    
    all_subfamilies = []
    for group in length_groups.values():
        if len(group) <= 2:
            all_subfamilies.append(group)
        elif len(group) <= 15:
            subfams = hierarchical_clustering(group, config)
            all_subfamilies.extend(subfams)
        else:
            subfams = sampled_clustering(group, config)
            all_subfamilies.extend(subfams)
    
    return [sf for sf in all_subfamilies if len(sf) >= 2]

def sampled_clustering(copies: List[Dict], config: PipelineConfig) -> List[List[Dict]]:
    """采样聚类"""
    # 选择代表序列
    representatives = copies[::5]
    if len(representatives) < 3:
        representatives = copies[:3]
    
    # 聚类代表序列
    rep_subfamilies = hierarchical_clustering(representatives, config)
    
    # 分配其他序列
    final_subfamilies = [[] for _ in rep_subfamilies]
    
    for copy in copies:
        best_family = 0
        best_identity = 0
        
        for fam_id, rep_family in enumerate(rep_subfamilies):
            if rep_family:
                identity = calculate_identity(copy['sequence'], rep_family[0]['sequence'])
                if identity > best_identity:
                    best_identity = identity
                    best_family = fam_id
        
        if best_identity >= config.identity_threshold:
            final_subfamilies[best_family].append(copy)
        else:
            final_subfamilies.append([copy])
    
    return [sf for sf in final_subfamilies if len(sf) >= 2]

def build_subfamily_consensus_full(seed_seq: Dict, subfamily_copies: List[Dict], 
                                  subfamily_id: int, genome_file: str,
                                  config: PipelineConfig) -> Optional[Dict]:
    """构建亚家族共识序列"""
    logger = logging.getLogger(__name__)
    
    try:
        # 限制序列数量
        if len(subfamily_copies) > 30:
            sorted_copies = sorted(subfamily_copies, key=lambda x: x.get('identity', 0), reverse=True)
            subfamily_copies = sorted_copies[:30]
        
        # 准备序列
        sequences = []
        for i, copy in enumerate(subfamily_copies):
            sequences.append({
                'id': f"{seed_seq['id']}_copy_{i}",
                'sequence': copy['sequence']
            })
        
        # 多序列比对
        if len(sequences) > 1:
            alignment = run_mafft(sequences, config)
            if alignment:
                consensus_seq = build_consensus_from_alignment(
                    alignment,
                    min_coverage=config.column_coverage_threshold
                )
            else:
                # MAFFT失败，使用最长序列
                longest = max(subfamily_copies, key=lambda x: len(x['sequence']))
                consensus_seq = longest['sequence']
        else:
            consensus_seq = subfamily_copies[0]['sequence']
        
        if not consensus_seq or len(consensus_seq) < config.min_length:
            return None
        
        # 计算统计信息
        avg_identity = np.mean([c.get('identity', 0) for c in subfamily_copies])
        
        # TSD信息
        tsd_info = {'has_tsd': False}
        tsd_sequences = [c['tsd']['sequence'] for c in subfamily_copies 
                        if c.get('has_tsd') and c.get('tsd')]
        if tsd_sequences:
            tsd_counter = Counter(tsd_sequences)
            if tsd_counter:
                most_common_tsd = tsd_counter.most_common(1)[0][0]
                tsd_info = {
                    'has_tsd': True,
                    'tsd_sequence': most_common_tsd,
                    'tsd_support': tsd_counter[most_common_tsd] / len(subfamily_copies)
                }
        
        return {
            'id': f"{seed_seq['id']}_subfamily_{subfamily_id}",
            'sequence': consensus_seq,
            'seed_id': seed_seq['id'],
            'subfamily_id': subfamily_id,
            'num_copies': len(subfamily_copies),
            'avg_identity': avg_identity,
            'quality_class': seed_seq.get('quality_class', 'A'),
            'has_tsd': tsd_info['has_tsd'],
            'tsd_sequence': tsd_info.get('tsd_sequence', ''),
            'tsd_support': tsd_info.get('tsd_support', 0)
        }
        
    except Exception as e:
        logger.error(f"Failed to build consensus: {e}")
        return None

# ==================== Phase 3: 质量控制和去冗余 ====================

class QualityControllerFull:
    """完整的Phase 3实现"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def run(self, consensus_list: List[Dict]) -> Dict:
        """执行Phase 3"""
        self.logger.info("="*60)
        self.logger.info("Phase 3: Quality Control and Redundancy Removal")
        self.logger.info("="*60)
        
        # 质量过滤
        filtered = self._quality_filtering(consensus_list)
        
        # 嵌合体检测
        non_chimeric = self._chimera_detection(filtered)
        
        # 边界质量评估
        boundary_filtered = self._boundary_quality_filter(non_chimeric)
        
        # 去冗余
        self._remove_redundancy(boundary_filtered)
        
        # 生成统计
        stats = self._generate_statistics(consensus_list, boundary_filtered)
        
        return stats
    
    def _quality_filtering(self, consensus_list: List[Dict]) -> List[Dict]:
        """质量过滤"""
        filtered = []
        
        for consensus in consensus_list:
            # 长度过滤
            if len(consensus['sequence']) < self.config.min_length:
                continue
            
            # 拷贝数过滤
            if consensus.get('num_copies', 0) < 2:
                continue
            
            # 相似度过滤
            if consensus.get('avg_identity', 0) < 70:
                continue
            
            filtered.append(consensus)
        
        self.logger.info(f"Quality filtering: {len(filtered)}/{len(consensus_list)} passed")
        return filtered
    
    def _chimera_detection(self, consensus_list: List[Dict]) -> List[Dict]:
        """嵌合体检测（简化版）"""
        # 这里应该实现完整的嵌合体检测
        # 暂时使用简单的长度异常检测
        non_chimeric = []
        
        for consensus in consensus_list:
            seq_len = len(consensus['sequence'])
            
            # 检查异常长度
            if seq_len > 15000:
                # 可能是嵌合体
                self.logger.debug(f"Potential chimera: {consensus['id']} ({seq_len}bp)")
                continue
            
            non_chimeric.append(consensus)
        
        self.logger.info(f"Chimera detection: {len(non_chimeric)}/{len(consensus_list)} passed")
        return non_chimeric
    
    def _boundary_quality_filter(self, consensus_list: List[Dict]) -> List[Dict]:
        """边界质量过滤"""
        filtered = []
        
        for consensus in consensus_list:
            # TSD支持度
            tsd_support = consensus.get('tsd_support', 0)
            
            # 计算边界质量分数
            boundary_score = 0.5  # 基础分数
            if consensus.get('has_tsd'):
                boundary_score += 0.3 * tsd_support
            
            # 拷贝数支持
            if consensus.get('num_copies', 0) > 10:
                boundary_score += 0.2
            
            consensus['boundary_score'] = boundary_score
            
            if boundary_score >= self.config.boundary_quality_threshold:
                filtered.append(consensus)
        
        self.logger.info(f"Boundary quality filter: {len(filtered)}/{len(consensus_list)} passed")
        return filtered
    
    def _remove_redundancy(self, consensus_list: List[Dict]):
        """去冗余"""
        if not consensus_list:
            return
        
        # 创建临时文件
        temp_input = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False, 
                                                 dir=self.config.temp_dir)
        for consensus in consensus_list:
            temp_input.write(f">{consensus['id']}\n{consensus['sequence']}\n")
        temp_input.close()
        
        try:
            # 95%去冗余（masking库）
            masking_output = os.path.join(self.config.output_dir, 'consensus_masking.fa')
            self._run_cdhit(temp_input.name, masking_output, 
                          self.config.redundancy_threshold_masking)
            
            # 90%去冗余（分析库）
            analysis_output = os.path.join(self.config.output_dir, 'consensus_analysis.fa')
            self._run_cdhit(temp_input.name, analysis_output,
                          self.config.redundancy_threshold_analysis)
            
            # 统计
            masking_count = sum(1 for _ in SeqIO.parse(masking_output, "fasta"))
            analysis_count = sum(1 for _ in SeqIO.parse(analysis_output, "fasta"))
            
            self.logger.info(f"Redundancy removal: {masking_count} (95%), {analysis_count} (90%)")
            
        finally:
            if os.path.exists(temp_input.name):
                os.remove(temp_input.name)
    
    def _run_cdhit(self, input_file: str, output_file: str, threshold: float):
        """运行CD-HIT"""
        try:
            cmd = [
                self.config.cdhit_path,
                '-i', input_file,
                '-o', output_file,
                '-c', str(threshold),
                '-n', '11',
                '-T', str(self.config.threads),
                '-M', '0',
                '-d', '0'
            ]
            
            subprocess.run(cmd, check=True, capture_output=True)
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"CD-HIT failed: {e}")
            # 如果失败，直接复制文件
            shutil.copy(input_file, output_file)
    
    def _generate_statistics(self, input_list: List[Dict], output_list: List[Dict]) -> Dict:
        """生成统计信息"""
        stats = {
            'input_consensus': len(input_list),
            'filtered_consensus': len(output_list),
            'avg_copy_number': np.mean([c.get('num_copies', 0) for c in output_list]) if output_list else 0,
            'avg_identity': np.mean([c.get('avg_identity', 0) for c in output_list]) if output_list else 0,
            'with_tsd': sum(1 for c in output_list if c.get('has_tsd')),
            'summary': f"Input:{len(input_list)}, Output:{len(output_list)}"
        }
        
        # 保存统计
        stats_file = os.path.join(self.config.output_dir, 'statistics.txt')
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        self.logger.info(f"Statistics saved to {stats_file}")
        return stats

# ==================== 主流程控制 ====================

class TEConsensusPipeline:
    """完整的TE共识序列构建流程"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.runner = RobustRunner(config)
        self.logger = setup_logging(config.output_dir, config.verbose)
        
        # 初始化各阶段
        self.phase1 = SequenceScreenerFull(config)
        self.phase2 = ConsensusBuilderFull(config)
        self.phase3 = QualityControllerFull(config)
        
        self.logger.info("TE Consensus Pipeline initialized")
        self.logger.info(f"Configuration: {vars(config)}")
    
    def run(self) -> Dict:
        """执行完整流程"""
        start_time = time.time()
        
        self.logger.info("="*80)
        self.logger.info("TE CONSENSUS PIPELINE - FULL VERSION")
        self.logger.info(f"Input: {self.config.repeatscout_file}")
        self.logger.info(f"Genome: {self.config.genome_file}")
        self.logger.info(f"Output: {self.config.output_dir}")
        self.logger.info(f"Threads: {self.config.threads}")
        self.logger.info("="*80)
        
        try:
            # Phase 1: 筛选和评分
            phase1_output = self.runner.run_with_checkpoint(
                self.phase1.run,
                checkpoint_name="phase1_complete"
            )
            
            # Phase 2: 共识构建
            consensus_list = self.runner.run_with_checkpoint(
                lambda: self.phase2.run(phase1_output),
                checkpoint_name="phase2_complete"
            )
            
            # Phase 3: 质量控制
            final_stats = self.runner.run_with_checkpoint(
                lambda: self.phase3.run(consensus_list),
                checkpoint_name="phase3_complete"
            )
            
            # 保存配置
            config_file = os.path.join(self.config.output_dir, 'pipeline_config.json')
            self.config.save(config_file)
            
            # 完成
            elapsed_time = time.time() - start_time
            self.logger.info("="*80)
            self.logger.info("PIPELINE COMPLETED SUCCESSFULLY")
            self.logger.info(f"Total time: {elapsed_time/60:.1f} minutes")
            self.logger.info(f"Results: {final_stats['summary']}")
            self.logger.info("Output files:")
            self.logger.info(f"  - {os.path.join(self.config.output_dir, 'consensus_masking.fa')}")
            self.logger.info(f"  - {os.path.join(self.config.output_dir, 'consensus_analysis.fa')}")
            self.logger.info(f"  - {os.path.join(self.config.output_dir, 'statistics.txt')}")
            self.logger.info(f"  - {os.path.join(self.config.output_dir, 'pipeline.log')}")
            self.logger.info("="*80)
            
            return final_stats
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}")
            raise
        
        finally:
            # 清理
            self._cleanup()
    
    def _cleanup(self):
        """清理临时文件"""
        if not self.config.keep_temp and os.path.exists(self.config.temp_dir):
            try:
                shutil.rmtree(self.config.temp_dir)
                self.logger.info("Temporary files cleaned up")
            except Exception as e:
                self.logger.warning(f"Failed to clean temp files: {e}")
        
        if not self.config.keep_checkpoints and os.path.exists(self.config.checkpoint_dir):
            try:
                # 只删除本次运行的检查点
                for checkpoint in ['phase1_complete.pkl', 'phase2_complete.pkl', 'phase3_complete.pkl']:
                    checkpoint_file = os.path.join(self.config.checkpoint_dir, checkpoint)
                    if os.path.exists(checkpoint_file):
                        os.remove(checkpoint_file)
                self.logger.info("Checkpoints cleaned up")
            except Exception as e:
                self.logger.warning(f"Failed to clean checkpoints: {e}")

# ==================== API函数 ====================

def run_pipeline(repeatscout_file: str, genome_file: str, output_dir: str, **kwargs) -> Dict:
    """
    运行完整的TE共识序列构建流程
    
    参数:
        repeatscout_file: RepeatScout输出文件路径
        genome_file: 参考基因组文件路径
        output_dir: 输出目录
        **kwargs: 其他配置参数
    
    返回:
        包含统计信息的字典
    """
    config = PipelineConfig(repeatscout_file, genome_file, output_dir, **kwargs)
    pipeline = TEConsensusPipeline(config)
    return pipeline.run()

def run_phase1_only(repeatscout_file: str, genome_file: str, output_dir: str, **kwargs) -> Dict:
    """只运行Phase 1筛选"""
    config = PipelineConfig(repeatscout_file, genome_file, output_dir, **kwargs)
    logger = setup_logging(output_dir, kwargs.get('verbose', False))
    phase1 = SequenceScreenerFull(config)
    return phase1.run()

def run_phase2_only(phase1_output: Dict, genome_file: str, output_dir: str, **kwargs) -> List[Dict]:
    """只运行Phase 2共识构建"""
    config = PipelineConfig("dummy", genome_file, output_dir, **kwargs)
    logger = setup_logging(output_dir, kwargs.get('verbose', False))
    phase2 = ConsensusBuilderFull(config)
    return phase2.run(phase1_output)

def run_phase3_only(consensus_list: List[Dict], output_dir: str, **kwargs) -> Dict:
    """只运行Phase 3质量控制"""
    config = PipelineConfig("dummy", "dummy", output_dir, **kwargs)
    logger = setup_logging(output_dir, kwargs.get('verbose', False))
    phase3 = QualityControllerFull(config)
    return phase3.run(consensus_list)

# ==================== 命令行接口 ====================

def main():
    """命令行接口"""
    parser = argparse.ArgumentParser(
        description='TE Consensus Pipeline - Full Version',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  %(prog)s -r repeatscout.fa -g genome.fa -o output_dir
  
  # With custom parameters
  %(prog)s -r repeatscout.fa -g genome.fa -o output_dir -t 32 --min-copy-number 10
  
  # Resume from checkpoint
  %(prog)s -r repeatscout.fa -g genome.fa -o output_dir --resume
  
  # Keep all intermediate files
  %(prog)s -r repeatscout.fa -g genome.fa -o output_dir --keep-temp --keep-checkpoints
        """
    )
    
    # 必需参数
    parser.add_argument('-r', '--repeatscout', required=True, 
                       help='RepeatScout output file (filtered with filter-stage-1.prl)')
    parser.add_argument('-g', '--genome', required=True,
                       help='Reference genome file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory')
    
    # 性能参数
    perf_group = parser.add_argument_group('Performance options')
    perf_group.add_argument('-t', '--threads', type=int, default=8,
                           help='Number of threads (default: 8)')
    perf_group.add_argument('--batch-size', type=int, default=100,
                           help='Batch size for processing (default: 100)')
    
    # 过滤参数
    filter_group = parser.add_argument_group('Filtering options')
    filter_group.add_argument('--min-length', type=int, default=80,
                             help='Minimum sequence length (default: 80)')
    filter_group.add_argument('--max-length', type=int, default=20000,
                             help='Maximum sequence length (default: 20000)')
    filter_group.add_argument('--max-n-percent', type=float, default=0.2,
                             help='Maximum N content (default: 0.2)')
    filter_group.add_argument('--dust-threshold', type=float, default=7,
                             help='DUST complexity threshold (default: 7)')
    
    # 聚类参数
    cluster_group = parser.add_argument_group('Clustering options')
    cluster_group.add_argument('--min-copy-number', type=int, default=5,
                              help='Minimum copy number (default: 5)')
    cluster_group.add_argument('--identity-threshold', type=float, default=0.85,
                              help='Identity threshold for clustering (default: 0.85)')
    cluster_group.add_argument('--coverage-threshold', type=float, default=0.60,
                              help='Coverage threshold (default: 0.60)')
    
    # 输出参数
    output_group = parser.add_argument_group('Output options')
    output_group.add_argument('--redundancy-masking', type=float, default=0.95,
                             help='Redundancy threshold for masking library (default: 0.95)')
    output_group.add_argument('--redundancy-analysis', type=float, default=0.90,
                             help='Redundancy threshold for analysis library (default: 0.90)')
    
    # 工具路径
    tool_group = parser.add_argument_group('External tool paths')
    tool_group.add_argument('--repeatmasker-path', default='RepeatMasker',
                           help='Path to RepeatMasker (default: RepeatMasker)')
    tool_group.add_argument('--mafft-path', default='mafft',
                           help='Path to MAFFT (default: mafft)')
    tool_group.add_argument('--cdhit-path', default='cd-hit-est',
                           help='Path to CD-HIT-EST (default: cd-hit-est)')
    
    # 其他选项
    other_group = parser.add_argument_group('Other options')
    other_group.add_argument('--resume', action='store_true',
                            help='Resume from checkpoint if available')
    other_group.add_argument('--clear-cache', action='store_true',
                            help='Clear cache before running')
    other_group.add_argument('--keep-temp', action='store_true',
                            help='Keep temporary files')
    other_group.add_argument('--keep-checkpoints', action='store_true',
                            help='Keep checkpoint files')
    other_group.add_argument('-v', '--verbose', action='store_true',
                            help='Verbose output')
    other_group.add_argument('-c', '--config', 
                            help='Load configuration from JSON file')
    
    args = parser.parse_args()
    
    # 准备配置
    if args.config:
        # 从配置文件加载
        config = PipelineConfig.load(args.config)
        # 覆盖命令行参数
        config.repeatscout_file = args.repeatscout
        config.genome_file = args.genome
        config.output_dir = args.output
    else:
        # 从命令行参数创建配置
        kwargs = {
            'threads': args.threads,
            'batch_size': args.batch_size,
            'min_length': args.min_length,
            'max_length': args.max_length,
            'max_n_percent': args.max_n_percent,
            'dust_threshold': args.dust_threshold,
            'min_copy_number': args.min_copy_number,
            'identity_threshold': args.identity_threshold,
            'coverage_threshold': args.coverage_threshold,
            'redundancy_threshold_masking': args.redundancy_masking,
            'redundancy_threshold_analysis': args.redundancy_analysis,
            'repeatmasker_path': args.repeatmasker_path,
            'mafft_path': args.mafft_path,
            'cdhit_path': args.cdhit_path,
            'resume': args.resume,
            'clear_cache': args.clear_cache,
            'keep_temp': args.keep_temp,
            'keep_checkpoints': args.keep_checkpoints,
            'verbose': args.verbose
        }
        config = PipelineConfig(args.repeatscout, args.genome, args.output, **kwargs)
    
    # 运行流程
    try:
        pipeline = TEConsensusPipeline(config)
        stats = pipeline.run()
        print(f"\nPipeline completed successfully!")
        print(f"Results: {stats['summary']}")
        sys.exit(0)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()