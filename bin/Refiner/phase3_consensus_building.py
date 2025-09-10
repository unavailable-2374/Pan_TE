"""
Phase 3: 共识序列构建和优化（重构版本）
专注功能：基于Phase 2处理后的序列构建高质量共识序列
输入：Phase 2输出的拆分序列、非嵌合体序列和未拆分序列
输出：高质量共识序列库，准备用于基因组屏蔽
"""

import logging
import gc
import numpy as np
from typing import Dict, List, Tuple, Any, Optional
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing as mp
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import queue
import threading
import time
import psutil
import os

from config import PipelineConfig
from utils.robust_runner import RobustRunner

logger = logging.getLogger(__name__)


class ConsensusBuilder:
    """Phase 3: 全面的共识序列构建器（重构版本）
    
    主要功能：
    1. 接收Phase 2处理后的所有序列（拆分片段、非嵌合体、未拆分序列）
    2. 对每个序列进行共识构建和优化
    3. 生成高质量的共识序列库
    4. 为Phase 4基因组屏蔽做准备
    """
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.runner = RobustRunner(config)
        
        # 获取系统内存信息
        try:
            mem = psutil.virtual_memory()
            available_memory_gb = mem.available / (1024**3)
            total_memory_gb = mem.total / (1024**3)
            logger.info(f"System memory: {available_memory_gb:.1f}GB available / {total_memory_gb:.1f}GB total")
        except:
            available_memory_gb = 16  # 默认假设16GB可用
            total_memory_gb = 32
            logger.warning("Could not detect system memory, assuming 16GB available")
        
        # 设置并行处理workers - 基于内存和CPU数量
        cpu_count = mp.cpu_count()
        
        # 估算每个worker需要的内存（降低估算值以允许更多并行）
        # 实际上大部分序列处理不需要太多内存
        memory_per_worker = 0.5  # GB (降低到0.5GB per worker)
        max_workers_by_memory = max(1, int(available_memory_gb / memory_per_worker))
        
        # 优化：直接使用config.threads作为并行度
        # 让系统充分利用所有指定的线程数
        max_workers_by_cpu = config.threads
        
        # 取内存和CPU限制的较小值
        self.max_workers = min(max_workers_by_memory, max_workers_by_cpu)
        self.max_workers = max(1, min(self.max_workers, 128))  # 至少1个worker，最多128个
        
        # 如果内存限制了workers数量，记录警告
        if max_workers_by_memory < max_workers_by_cpu:
            logger.warning(f"Memory constraints: reducing workers from {max_workers_by_cpu} to {self.max_workers}")
        
        # 计算每个子进程应该使用的线程数
        # 如果使用多进程，每个进程只使用1个线程（因为进程本身就是并行单位）
        if self.max_workers > 1:
            self.threads_per_worker = 1
            self.mafft_threads = 1  # 多进程模式下，每个MAFFT使用1个线程
        else:
            # 单进程模式，可以使用多个线程
            self.threads_per_worker = config.threads
            self.mafft_threads = min(config.threads, 4)  # 限制MAFFT最多使用4个线程
        
        logger.info(f"Phase 3 Consensus Builder initialized:")
        logger.info(f"  - Workers: {self.max_workers} (memory-aware)")
        logger.info(f"  - Threads per worker: {self.threads_per_worker}")
        logger.info(f"  - MAFFT threads: {self.mafft_threads}")
        logger.info(f"  - Available memory: {available_memory_gb:.1f}GB")
    
    def _sort_sequences_by_complexity(self, sequences: List[Dict]) -> List[Dict]:
        """
        按照预估处理复杂度排序序列
        复杂度 = 序列长度 × 拷贝数
        优先处理高复杂度序列以避免长尾效应
        """
        # 计算每个序列的复杂度分数
        for seq in sequences:
            # 序列长度
            seq_length = len(seq.get('sequence', ''))
            
            # 估算拷贝数（从Phase 1的RepeatMasker结果）
            rm_hits = seq.get('rm_hits', [])
            estimated_copies = len(rm_hits) if rm_hits else 1
            
            # 复杂度分数：长度和拷贝数的乘积
            # 使用对数缩放避免极端值主导
            complexity_score = np.log1p(seq_length) * np.log1p(estimated_copies)
            
            # 添加质量权重：A类序列可能更重要
            if seq.get('quality_class') == 'A':
                complexity_score *= 1.2
            
            seq['_complexity_score'] = complexity_score
        
        # 按复杂度降序排序（最复杂的先处理）
        sorted_sequences = sorted(sequences, key=lambda x: x.get('_complexity_score', 0), reverse=True)
        
        # 记录排序信息
        if sequences:
            top_5 = sorted_sequences[:5]
            bottom_5 = sorted_sequences[-5:]
            logger.info(f"Sequence complexity sorting:")
            
            # 格式化top 5信息
            top_5_info = []
            for s in top_5:
                seq_len = len(s.get('sequence', ''))
                num_copies = len(s.get('rm_hits', []))
                top_5_info.append((s['id'], f"{seq_len}bp", f"{num_copies} copies"))
            logger.info(f"  Top 5 complex: {top_5_info}")
            
            # 格式化bottom 5信息
            bottom_5_info = []
            for s in bottom_5:
                seq_len = len(s.get('sequence', ''))
                num_copies = len(s.get('rm_hits', []))
                bottom_5_info.append((s['id'], f"{seq_len}bp", f"{num_copies} copies"))
            logger.info(f"  Bottom 5 simple: {bottom_5_info}")
        
        return sorted_sequences
    
    def run(self, phase2_output: Dict) -> Dict[str, Any]:
        """执行Phase 3主流程 - 专注于所有序列的共识构建"""
        
        # 获取Phase 2处理后的所有序列
        input_sequences = self._extract_sequences_from_phase2(phase2_output)
        
        logger.info(f"Phase 3 input: {len(input_sequences)} sequences for consensus building")
        
        # 分类输入序列
        sequence_categories = self._categorize_input_sequences(input_sequences)
        logger.info(f"Input sequence categories:")
        logger.info(f"  - Non-chimeric: {len(sequence_categories['non_chimeric'])}")
        logger.info(f"  - Chimera segments: {len(sequence_categories['chimera_segments'])}")
        logger.info(f"  - Unsplit chimeras: {len(sequence_categories['unsplit_chimeras'])}")
        logger.info(f"  - Processing failures: {len(sequence_categories['failed_sequences'])}")
        
        # 统计信息
        stats = {
            'input_sequences': len(input_sequences),
            'processed_sequences': 0,
            'processed': 0,
            'single_consensus': 0,
            'multiple_consensus': 0,
            'failed': 0,
            'failed_sequences': 0,
            'total_consensus': 0,
            'high_quality_consensus': 0,
            'medium_quality_consensus': 0,
            'low_quality_consensus': 0
        }
        
        # 所有输入序列都是共识构建的候选（Phase 1已经筛选过）
        candidate_sequences = input_sequences
        logger.info(f"Processing {len(candidate_sequences)} sequences for consensus building")
        
        # 按复杂度排序优化处理顺序
        sorted_candidates = self._sort_sequences_by_complexity(candidate_sequences)
        
        # 并行构建共识序列
        all_consensus = []
        
        if sorted_candidates:
            all_consensus = self._build_consensus_parallel(sorted_candidates, stats)
        
        # 质量过滤、分级和去重
        filtered_consensus = self._filter_and_grade_consensus(all_consensus, stats)
        deduplicated_consensus = self._deduplicate_consensus_sequences(filtered_consensus, stats)
        
        # 生成Phase 4兼容的输出结果
        result = self._generate_phase4_compatible_output(deduplicated_consensus, stats, sequence_categories)
        
        logger.info(f"Phase 3 complete:")
        logger.info(f"  - Input: {stats['input_sequences']} sequences")
        logger.info(f"  - Processed: {stats['processed_sequences']} candidates")
        logger.info(f"  - Generated: {stats['total_consensus']} consensus sequences")
        logger.info(f"  - Quality distribution: High={stats['high_quality_consensus']}, "
                   f"Medium={stats['medium_quality_consensus']}, Low={stats['low_quality_consensus']}")
        
        return result
    
    
    def _process_sequences_dynamic(self, sequences: List[Dict], stats: Dict) -> List[Dict]:
        """动态并行处理序列 - 内存优化版本"""
        all_consensus = []
        total_sequences = len(sequences)
        
        # 优化：按预估处理时间排序（长度×拷贝数），优先处理耗时任务
        sorted_sequences = self._sort_sequences_by_complexity(sequences)
        
        # 进度追踪
        progress_lock = threading.Lock()
        last_report_time = time.time()
        report_interval = 10  # 每10秒报告一次进度
        
        def update_progress(result: Dict, seq_id: str):
            """更新进度和统计"""
            nonlocal last_report_time
            
            with progress_lock:
                stats['processed'] += 1
                
                if result['consensus_count'] == 0:
                    stats['failed'] += 1
                elif result['consensus_count'] == 1:
                    stats['single_consensus'] += 1
                else:
                    stats['multiple_consensus'] += 1
                stats['total_consensus'] += result['consensus_count']
                
                all_consensus.extend(result['consensus_list'])
                
                # 定期报告进度
                current_time = time.time()
                if current_time - last_report_time > report_interval or stats['processed'] == total_sequences:
                    percent = (stats['processed'] / total_sequences) * 100
                    logger.info(f"Progress: {stats['processed']}/{total_sequences} ({percent:.1f}%) sequences processed")
                    logger.info(f"  Generated {stats['total_consensus']} consensus sequences")
                    logger.info(f"  Single: {stats['single_consensus']}, Multiple: {stats['multiple_consensus']}, Failed: {stats['failed']}")
                    
                    # 估算剩余时间
                    if stats['processed'] > 0 and stats['processed'] < total_sequences:
                        elapsed = current_time - self.start_time
                        rate = stats['processed'] / elapsed
                        remaining = (total_sequences - stats['processed']) / rate
                        logger.info(f"  Estimated time remaining: {remaining/60:.1f} minutes")
                    
                    # 显示内存使用情况
                    try:
                        mem = psutil.virtual_memory()
                        logger.info(f"  Memory usage: {mem.percent:.1f}% ({mem.used/(1024**3):.1f}GB / {mem.total/(1024**3):.1f}GB)")
                    except:
                        pass
                    
                    last_report_time = current_time
        
        # 记录开始时间
        self.start_time = time.time()
        
        # 检测是否在集群环境（通过环境变量判断）
        is_cluster = os.environ.get('SLURM_JOB_ID') or os.environ.get('PBS_JOBID') or os.environ.get('LSB_JOBID')
        
        # 决定是否使用多进程
        use_multiprocessing = self.max_workers > 1
        
        # 在内存受限的环境中，考虑使用串行处理
        try:
            mem = psutil.virtual_memory()
            if mem.available < 4 * (1024**3):  # 少于4GB可用内存
                logger.warning(f"Low memory detected ({mem.available/(1024**3):.1f}GB), using sequential processing")
                use_multiprocessing = False
        except:
            pass
        
        if not use_multiprocessing:
            # 串行处理 - 避免fork()内存问题
            logger.info(f"Using sequential processing for {len(sorted_sequences)} sequences")
            
            for seq in sorted_sequences:
                try:
                    expansion_result = expand_and_improve_sequence(
                        seq,
                        self.config,
                        self.mafft_threads
                    )
                    update_progress(expansion_result, seq['id'])
                    
                    if expansion_result['consensus_count'] > 1:
                        logger.debug(f"{seq['id']}: generated {expansion_result['consensus_count']} consensus sequences")
                    
                except Exception as e:
                    logger.error(f"Failed processing {seq['id']}: {e}")
                    failed_result = {
                        'source_id': seq['id'],
                        'consensus_count': 0,
                        'consensus_list': []
                    }
                    update_progress(failed_result, seq['id'])
                
                # 主动垃圾回收
                if stats['processed'] % 10 == 0:
                    gc.collect()
        else:
            # 并行处理 - 根据任务特性选择合适的执行器
            # 如果主要是I/O密集型（如调用外部工具），使用ThreadPoolExecutor
            # 如果是CPU密集型（如序列处理），使用ProcessPoolExecutor
            
            # 判断是否主要是I/O密集型（基于经验判断）
            cpu_count = mp.cpu_count()
            use_thread_pool = self.max_workers > cpu_count * 2  # 如果workers远超CPU数，说明期望更多并发
            
            if use_thread_pool:
                logger.info(f"Using ThreadPoolExecutor with {self.max_workers} workers (I/O-bound optimization)")
                executor_class = ThreadPoolExecutor
            else:
                logger.info(f"Using ProcessPoolExecutor with {self.max_workers} workers (CPU-bound optimization)")
                executor_class = ProcessPoolExecutor
                
                # 对于ProcessPoolExecutor，尝试使用spawn方法
                try:
                    original_method = mp.get_start_method()
                    if original_method != 'spawn':
                        logger.info(f"Switching multiprocessing from '{original_method}' to 'spawn' method")
                        mp.set_start_method('spawn', force=True)
                except (RuntimeError, ValueError) as e:
                    logger.debug(f"Could not set spawn method: {e}")
            
            # 使用选定的执行器进行并行处理
            with executor_class(max_workers=self.max_workers) as executor:
                # 优化：提交所有任务到线程池，充分利用所有workers
                futures_to_seq = {}
                
                logger.info(f"Submitting {len(sorted_sequences)} sequences to {self.max_workers} workers for parallel processing")
                
                # 一次性提交所有任务，让executor的内部队列管理任务分配
                # 这样可以确保所有workers始终有任务执行，最大化并行效率
                for seq in sorted_sequences:
                    future = executor.submit(
                        expand_and_improve_sequence,
                        seq,
                        self.config,
                        self.mafft_threads
                    )
                    futures_to_seq[future] = seq
                
                logger.info(f"All {len(sorted_sequences)} tasks submitted successfully")
                
                # 使用as_completed处理完成的任务，这样可以立即处理完成的结果
                # 而不需要等待整个批次完成
                completed_count = 0
                memory_check_interval = max(10, len(sorted_sequences) // 20)  # 每5%检查一次内存
                
                for future in as_completed(futures_to_seq):
                    seq = futures_to_seq[future]
                    completed_count += 1
                    
                    try:
                        expansion_result = future.result(timeout=300)
                        update_progress(expansion_result, seq['id'])
                        
                        if expansion_result['consensus_count'] > 1:
                            logger.debug(f"{seq['id']}: generated {expansion_result['consensus_count']} consensus sequences")
                    
                    except Exception as e:
                        logger.error(f"Failed processing {seq['id']}: {e}")
                        failed_result = {
                            'source_id': seq['id'],
                            'consensus_count': 0,
                            'consensus_list': []
                        }
                        update_progress(failed_result, seq['id'])
                    
                    # 定期检查内存和进行垃圾回收
                    if completed_count % memory_check_interval == 0:
                        gc.collect()
                        try:
                            mem = psutil.virtual_memory()
                            if mem.percent > 85:
                                logger.warning(f"High memory usage ({mem.percent:.1f}%), triggering aggressive GC")
                                gc.collect(2)  # Full collection
                                time.sleep(1)  # Brief pause to allow memory reclamation
                        except:
                            pass
        
        # 最终统计
        total_time = time.time() - self.start_time
        logger.info(f"Phase 2 processing completed in {total_time/60:.1f} minutes")
        logger.info(f"Average processing rate: {total_sequences/(total_time+0.001):.1f} sequences/second")
        
        return all_consensus

    def _extract_sequences_from_phase2(self, phase2_output: Dict) -> List[Dict]:
        """从Phase 2输出中提取所有序列"""
        sequences = []
        
        # 获取所有类型的序列
        if isinstance(phase2_output, dict):
            # 首先检查新的统一格式 (Phase 2重构版本)
            if 'sequences' in phase2_output:
                # 新格式：所有序列都在'sequences'键下
                all_sequences = phase2_output.get('sequences', [])
                for seq in all_sequences:
                    if 'source' not in seq:
                        seq['source'] = 'unknown'
                    sequences.append(seq)
            else:
                # 处理旧格式的嵌合体拆分结果
                split_sequences = phase2_output.get('split_sequences', [])
                non_chimeric = phase2_output.get('non_chimeric_sequences', [])
                unsplit_chimeras = phase2_output.get('unsplit_chimeras', [])
                failed_sequences = phase2_output.get('failed_sequences', [])
                
                # 为每个序列标记来源
                for seq in split_sequences:
                    seq['source'] = 'chimera_split'
                    sequences.append(seq)
                
                for seq in non_chimeric:
                    seq['source'] = 'non_chimeric'
                    sequences.append(seq)
                
                for seq in unsplit_chimeras:
                    seq['source'] = 'unsplit_chimera'
                    sequences.append(seq)
                
                for seq in failed_sequences:
                    seq['source'] = 'processing_failed'
                    sequences.append(seq)
                    
        elif isinstance(phase2_output, list):
            # 兼容列表格式的输入
            for seq in phase2_output:
                if 'source' not in seq:
                    seq['source'] = 'unknown'
                sequences.append(seq)
        
        logger.info(f"Extracted {len(sequences)} sequences from Phase 2 output")
        
        return sequences


def expand_and_improve_sequence(seed_seq: Dict, config: PipelineConfig, mafft_threads: int = 1) -> Dict:
    """
    扩展和改进单个候选序列
    目标：提高序列质量和代表性，不是减少数量
    
    Args:
        seed_seq: 种子序列
        config: 配置对象
        mafft_threads: MAFFT使用的线程数（避免线程爆炸）
    """
    
    result = {
        'source_id': seed_seq['id'],
        'consensus_count': 0,
        'consensus_list': []
    }
    
    try:
        # 1. 获取基因组拷贝（优先使用Phase 1的结果）
        copies = get_genome_copies(seed_seq, config)
        
        if not copies or len(copies) < 2:
            # 拷贝太少，直接使用原始序列
            if len(copies) == 1:
                # 只有一个拷贝，使用该拷贝改进原始序列
                improved = improve_with_single_copy(seed_seq, copies[0])
            else:
                # 没有拷贝，使用原始序列
                improved = create_default_consensus(seed_seq)
            
            if improved:
                result['consensus_list'].append(improved)
                result['consensus_count'] = 1
            return result
        
        logger.debug(f"{seed_seq['id']}: Found {len(copies)} copies")
        
        # 2. 分析拷贝特征
        characteristics = analyze_copy_characteristics(copies)
        
        # 3. 根据特征决定是否需要分亚家族
        if should_split_subfamilies(characteristics, len(copies)):
            # 识别亚家族并分别构建共识
            subfamilies = identify_subfamilies_for_expansion(copies, characteristics)
            logger.debug(f"{seed_seq['id']}: Split into {len(subfamilies)} subfamilies")
            
            for sf_idx, subfamily in enumerate(subfamilies):
                if len(subfamily) >= 2:  # 至少2个拷贝才构建共识
                    consensus = build_expanded_consensus(
                        subfamily,
                        seed_seq,
                        subfamily_id=sf_idx,
                        characteristics=characteristics,
                        config=config,
                        mafft_threads=mafft_threads
                    )
                    if consensus:
                        result['consensus_list'].append(consensus)
        else:
            # 不分亚家族，直接构建单一共识
            consensus = build_expanded_consensus(
                copies,
                seed_seq,
                subfamily_id=0,
                characteristics=characteristics,
                config=config,
                mafft_threads=mafft_threads
            )
            if consensus:
                result['consensus_list'].append(consensus)
        
        result['consensus_count'] = len(result['consensus_list'])
        
        # 如果没有生成任何共识，至少保留改进的原始序列
        if result['consensus_count'] == 0:
            improved = improve_original_sequence(seed_seq, copies)
            if improved:
                result['consensus_list'].append(improved)
                result['consensus_count'] = 1
        
    except Exception as e:
        logger.error(f"Error processing {seed_seq['id']}: {e}")
        # 出错时至少保留原始序列
        fallback = create_default_consensus(seed_seq)
        if fallback:
            result['consensus_list'].append(fallback)
            result['consensus_count'] = 1
    
    return result


def get_genome_copies(seed_seq: Dict, config: PipelineConfig) -> List[Dict]:
    """
    获取基因组拷贝，优先使用Phase 1的结果
    """
    from utils.sequence_utils import extract_sequence_from_genome, detect_tsd
    
    # 优先使用Phase 1提供的RepeatMasker结果
    if 'rm_hits' in seed_seq and seed_seq['rm_hits']:
        rm_hits = seed_seq['rm_hits']
        logger.debug(f"{seed_seq['id']}: Using {len(rm_hits)} hits from Phase 1")
    else:
        # Phase 1没有结果，需要重新运行RepeatMasker
        logger.debug(f"{seed_seq['id']}: Running RepeatMasker (no Phase 1 results)")
        from utils.alignment_utils import run_repeatmasker_single
        
        rm_params = {
            's': True,  # 敏感模式
            'no_is': True,
            'nolow': True,
            'div': 25,  # 允许25%的分化
            'cutoff': 200,
            'pa': 1
        }
        rm_hits = run_repeatmasker_single(seed_seq, config.genome_file, params=rm_params, config=config)
    
    if not rm_hits:
        return []
    
    # 应用优化的二轮采样策略：使用外部工具批量提取
    # 这样避免了提取所有45,081个序列的巨大开销，并且大幅提升提取速度
    try:
        from two_stage_sampling_optimized import apply_optimized_two_stage_sampling
        logger.info(f"Using optimized batch extraction for {len(rm_hits)} RM hits of {seed_seq['id']}")
        copies = apply_optimized_two_stage_sampling(rm_hits, seed_seq, config)
    except ImportError:
        # 降级到原始版本
        from two_stage_sampling import apply_two_stage_sampling
        logger.info(f"Falling back to standard extraction for {len(rm_hits)} RM hits of {seed_seq['id']}")
        copies = apply_two_stage_sampling(rm_hits, seed_seq, config)
    
    logger.info(f"Selected {len(copies)} representative copies after two-stage sampling")
    
    return copies


def analyze_copy_characteristics(copies: List[Dict]) -> Dict:
    """分析拷贝特征"""
    if not copies:
        return {}
    
    identities = [c.get('identity', 0) for c in copies]
    lengths = [c['length'] for c in copies]
    
    characteristics = {
        'copy_number': len(copies),
        'avg_identity': np.mean(identities),
        'std_identity': np.std(identities),
        'min_identity': min(identities),
        'max_identity': max(identities),
        'identity_range': max(identities) - min(identities),
        'avg_length': np.mean(lengths),
        'std_length': np.std(lengths),
        'length_cv': np.std(lengths) / np.mean(lengths) if np.mean(lengths) > 0 else 0,
        'has_tsd_ratio': sum(1 for c in copies if c.get('has_tsd', False)) / len(copies)
    }
    
    return characteristics


def should_split_subfamilies(characteristics: Dict, copy_number: int) -> bool:
    """
    决定是否需要分割成亚家族
    只在有明显差异时才分割，目的是提高注释的全面性
    """
    
    # 拷贝数太少，不分割
    if copy_number < 5:
        return False
    
    # identity差异很大，说明可能有不同的亚家族
    if characteristics['identity_range'] > 20:  # 20%的差异
        return True
    
    # 长度变异很大，可能有不同的变体（完整vs截断）
    if characteristics['length_cv'] > 0.3:  # 变异系数>30%
        return True
    
    # 大家族更可能有亚家族结构
    if copy_number > 20 and characteristics['std_identity'] > 5:
        return True
    
    return False


def identify_subfamilies_for_expansion(copies: List[Dict], characteristics: Dict) -> List[List[Dict]]:
    """
    识别亚家族，目的是发现不同的TE变体
    """
    from utils.sequence_utils import calculate_sequence_identity
    
    if len(copies) < 5:
        return [copies]
    
    # 构建简化的距离矩阵（主要基于identity和长度）
    n = len(copies)
    distance_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            # 序列相似度距离
            identity = copies[i].get('identity', 100) / 100.0
            identity_j = copies[j].get('identity', 100) / 100.0
            avg_identity = (identity + identity_j) / 2
            seq_distance = 1 - avg_identity
            
            # 长度差异距离
            len_i, len_j = copies[i]['length'], copies[j]['length']
            if max(len_i, len_j) > 0:
                length_ratio = min(len_i, len_j) / max(len_i, len_j)
                length_distance = 1 - length_ratio
            else:
                length_distance = 1.0
            
            # 综合距离（简单加权）
            total_distance = 0.7 * seq_distance + 0.3 * length_distance
            distance_matrix[i, j] = total_distance
            distance_matrix[j, i] = total_distance
    
    # 层次聚类
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    
    condensed_distances = squareform(distance_matrix)
    linkage_matrix = linkage(condensed_distances, method='average')
    
    # 动态阈值：基于距离分布
    distances = distance_matrix[np.triu_indices_from(distance_matrix, k=1)]
    if len(distances) > 0:
        # 使用距离的第75百分位作为阈值
        threshold = np.percentile(distances, 75)
        threshold = max(0.15, min(0.4, threshold))  # 限制在合理范围
    else:
        threshold = 0.25
    
    # 聚类
    clusters = fcluster(linkage_matrix, threshold, criterion='distance')
    
    # 组织成亚家族
    subfamilies = defaultdict(list)
    for i, cluster_id in enumerate(clusters):
        subfamilies[cluster_id].append(copies[i])
    
    # 过滤太小的亚家族（合并到最近的）
    final_subfamilies = []
    for subfamily in subfamilies.values():
        if len(subfamily) >= 2:
            final_subfamilies.append(subfamily)
        else:
            # 单个拷贝，尝试合并到最近的亚家族
            if final_subfamilies:
                final_subfamilies[0].extend(subfamily)
            else:
                final_subfamilies.append(subfamily)
    
    return final_subfamilies


def build_expanded_consensus(copies: List[Dict],
                            seed_seq: Dict,
                            subfamily_id: int,
                            characteristics: Dict,
                            config: PipelineConfig,
                            mafft_threads: int = 1) -> Optional[Dict]:
    """
    构建扩展的共识序列
    目标：生成更完整、更准确的共识
    """
    from utils.alignment_utils import run_mafft
    from utils.sequence_utils import build_consensus_from_msa, extend_consensus_boundaries
    
    try:
        if len(copies) == 1:
            # 单个拷贝，直接使用但尝试改进
            consensus_seq = copies[0]['sequence']
            # 尝试扩展边界
            if copies[0].get('extended_sequence'):
                consensus_seq = extend_consensus_boundaries(
                    consensus_seq,
                    copies[0]['extended_sequence']
                )
        else:
            # 多序列比对 - 始终使用MAFFT获得高质量共识
            sequences = [{'id': f"copy_{i}", 'sequence': copy['sequence']} for i, copy in enumerate(copies)]
            
            # 基于TE生物学特性选择MAFFT算法 - 优化保留长度
            avg_identity = characteristics.get('avg_identity', 70)
            length_cv = characteristics.get('length_cv', 0)  # 长度变异系数
            
            if len(sequences) > 100:
                # 大数据集使用最快算法
                algorithm = 'auto'
                logger.debug(f"Large dataset: using auto algorithm")
            elif avg_identity > 85 and length_cv < 0.3:
                # 高相似度且长度一致：年轻TE家族或同亚家族
                algorithm = 'localpair'  # L-INS-i: 适合高度相似序列
                logger.debug(f"Young TE family (id={avg_identity:.1f}%, cv={length_cv:.2f}): using localpair")
            elif length_cv > 0.5:
                # 长度变异大：使用E-INS-i以更好处理长度差异
                algorithm = 'genafpair'  # E-INS-i: 处理大插入/缺失，保留长序列
                logger.debug(f"TE with length variation (id={avg_identity:.1f}%, cv={length_cv:.2f}): using genafpair")
            elif avg_identity > 65:
                # 中等相似度，长度相对一致：保守TE类型(如SINE)
                algorithm = 'localpair'  # L-INS-i
                logger.debug(f"Conservative TE (id={avg_identity:.1f}%, cv={length_cv:.2f}): using localpair")
            else:
                # 低相似度或高度分化：古老TE家族
                algorithm = 'genafpair'  # E-INS-i: 处理高分化序列
                logger.debug(f"Ancient TE (id={avg_identity:.1f}%, cv={length_cv:.2f}): using genafpair")
            
            # TE特异的预处理：处理极端情况
            if avg_identity < 50:
                logger.warning(f"Very low identity ({avg_identity:.1f}%) - may not be same TE family")
                # 对于极低相似度，使用最宽松的参数
                algorithm = 'genafpair'
                
            # 检查长度变异是否过大 - 优化：保留长序列
            if length_cv > 1.0:
                logger.debug(f"Extreme length variation (cv={length_cv:.2f}) - intelligent filtering")
                # 智能过滤：保留最长的序列和有代表性的序列
                lengths = [len(seq['sequence']) for seq in sequences]
                max_len = max(lengths)
                median_len = np.median(lengths)
                
                # 策略1：总是保留最长的序列（可能是完整的TE）
                longest_seqs = [seq for seq in sequences if len(seq['sequence']) >= max_len * 0.95]
                
                # 策略2：保留中位数附近的代表性序列
                median_seqs = [seq for seq in sequences 
                              if median_len * 0.8 <= len(seq['sequence']) <= median_len * 1.5]
                
                # 合并两组序列
                filtered_sequences = list({seq['id']: seq for seq in longest_seqs + median_seqs}.values())
                
                if len(filtered_sequences) >= 2:
                    sequences = filtered_sequences
                    logger.debug(f"Retained {len(sequences)} sequences: {len(longest_seqs)} long + {len(median_seqs)} median")
                else:
                    logger.warning("Too few sequences remain after filtering - using all sequences")
            
            logger.debug(f"Running MAFFT with {algorithm} algorithm on {len(sequences)} sequences using {mafft_threads} threads")
            msa_result = run_mafft(sequences, algorithm=algorithm, thread=mafft_threads, config=config)
            
            if not msa_result:
                # MSA失败，使用TE特异的后备策略 - 优先保留长度
                logger.warning(f"MAFFT failed for {seed_seq['id']}, using length-preserving fallback strategy")
                
                # 策略：优先选择最长的高质量序列
                # 计算综合分数：长度权重70%，质量权重30%
                scored_copies = []
                for copy in copies:
                    length_score = copy['length'] / 1000  # 每1kb加1分
                    quality_score = copy.get('quality_score', copy.get('identity', 0) / 100)
                    combined_score = length_score * 0.7 + quality_score * 100 * 0.3
                    scored_copies.append((combined_score, copy))
                
                scored_copies.sort(key=lambda x: x[0], reverse=True)
                best_copy = scored_copies[0][1]
                consensus_seq = best_copy['sequence']
                
                logger.debug(f"Using best length-quality balanced sequence: "
                           f"{best_copy['length']}bp, identity={best_copy.get('identity', 0):.1f}%")
                
                # 如果最长序列比平均长度长很多，记录警告
                avg_length = np.mean([c['length'] for c in copies])
                if best_copy['length'] > avg_length * 1.5:
                    logger.info(f"Selected sequence is significantly longer than average "
                              f"({best_copy['length']}bp vs {avg_length:.0f}bp average) - likely full-length TE")
            else:
                # 构建高质量共识序列，根据TE特性调整参数 - 优化保留长度
                # 降低min_coverage以保留更多边界序列
                min_coverage = 0.2 if avg_identity > 75 else 0.15  # 降低覆盖度要求，保留边界
                quality_threshold = 0.5 if avg_identity > 80 else 0.4  # 适度降低质量阈值
                
                # 如果有长序列，进一步降低覆盖度要求
                seq_lengths = [len(seq['sequence']) for seq in sequences]
                if max(seq_lengths) > np.median(seq_lengths) * 1.5:
                    min_coverage = max(0.1, min_coverage - 0.1)  # 为长序列降低覆盖度要求
                    logger.debug(f"Long sequences detected, reducing min_coverage to {min_coverage}")
                
                consensus_seq = build_consensus_from_msa(
                    msa_result,
                    min_coverage=min_coverage,
                    use_iupac=True,   # 使用IUPAC模糊碱基表示变异位点
                    quality_threshold=quality_threshold
                )
                
                # TE特异的质量检查
                if consensus_seq and len(consensus_seq) > 0:
                    # 检查N含量
                    n_content = consensus_seq.count('N') / len(consensus_seq)
                    
                    # 检查模糊碱基含量（IUPAC codes: R,Y,S,W,K,M,B,D,H,V,N）
                    ambiguous_bases = set('RYSWKMBDHVN')
                    ambiguity_count = sum(1 for base in consensus_seq.upper() if base in ambiguous_bases)
                    ambiguity_ratio = ambiguity_count / len(consensus_seq)
                    
                    # 根据TE年龄/相似度动态设置模糊碱基阈值
                    if avg_identity > 85:
                        # 年轻TE：高相似度，应该有清晰的共识
                        max_ambiguity_ratio = 0.05  # 最多5%模糊碱基
                        logger.debug(f"Young TE family (identity={avg_identity:.1f}%): max ambiguity ratio = 5%")
                    elif avg_identity > 70:
                        # 中等年龄TE：适度变异
                        max_ambiguity_ratio = 0.10  # 最多10%模糊碱基
                        logger.debug(f"Middle-aged TE family (identity={avg_identity:.1f}%): max ambiguity ratio = 10%")
                    else:
                        # 古老TE：高度分化，但仍需保持特异性
                        max_ambiguity_ratio = 0.15  # 最多15%模糊碱基
                        logger.debug(f"Ancient TE family (identity={avg_identity:.1f}%): max ambiguity ratio = 15%")
                    
                    # 如果N含量过高或模糊碱基过多，需要处理
                    if n_content > 0.5 or ambiguity_ratio > max_ambiguity_ratio:
                        logger.warning(f"Poor consensus quality for {seed_seq['id']}_sf{subfamily_id}: "
                                     f"N-content={n_content:.1%}, ambiguity={ambiguity_ratio:.1%} "
                                     f"(threshold={max_ambiguity_ratio:.1%})")
                        
                        # 策略1：如果模糊碱基过多，尝试重建共识（不使用IUPAC）
                        if ambiguity_ratio > max_ambiguity_ratio and len(sequences) > 1:
                            logger.debug(f"Rebuilding consensus without IUPAC codes due to high ambiguity")
                            # 重新运行MAFFT但不使用模糊碱基
                            consensus_seq_strict = build_consensus_from_msa(
                                msa_result,
                                min_coverage=min_coverage,
                                use_iupac=False,  # 不使用模糊碱基
                                quality_threshold=quality_threshold
                            )
                            
                            # 检查新共识是否更好
                            if consensus_seq_strict and len(consensus_seq_strict) >= len(consensus_seq) * 0.9:
                                consensus_seq = consensus_seq_strict
                                logger.debug(f"Using strict consensus without ambiguous bases")
                            else:
                                # 策略2：优先使用最长的高质量序列（保留完整性）
                                # 按长度和质量综合排序
                                sorted_copies = sorted(copies, 
                                                     key=lambda x: (x['length'] * 0.6 + 
                                                                   x.get('quality_score', x.get('identity', 0)) * x['length'] * 0.4),
                                                     reverse=True)
                                longest_good_copy = sorted_copies[0]
                                consensus_seq = longest_good_copy['sequence']
                                logger.debug(f"Using longest high-quality sequence ({longest_good_copy['length']}bp) due to poor MSA quality")
                        else:
                            # N含量过高，使用最长的高质量序列
                            # 优先考虑长度，其次考虑质量
                            sorted_copies = sorted(copies,
                                                 key=lambda x: (x['length'] * 0.7 + 
                                                               x.get('quality_score', x.get('identity', 0)) * 1000 * 0.3),
                                                 reverse=True)
                            longest_copy = sorted_copies[0]
                            consensus_seq = longest_copy['sequence']
                            logger.debug(f"Using longest sequence ({longest_copy['length']}bp) due to high N-content")
        
        # 尝试扩展边界（如果可能）
        if consensus_seq and len(copies) > 0:
            # 找到最长的拷贝作为扩展参考
            longest_copy = max(copies, key=lambda x: x['length'])
            if longest_copy.get('extended_sequence') or len(longest_copy['sequence']) > len(consensus_seq):
                # 使用最长拷贝尝试扩展共识序列
                from utils.sequence_utils import extend_consensus_boundaries
                reference_seq = longest_copy.get('extended_sequence', longest_copy['sequence'])
                extended_consensus = extend_consensus_boundaries(consensus_seq, reference_seq)
                if len(extended_consensus) > len(consensus_seq):
                    logger.debug(f"Extended consensus from {len(consensus_seq)}bp to {len(extended_consensus)}bp")
                    consensus_seq = extended_consensus
        
        if not consensus_seq or len(consensus_seq) < 50:
            logger.debug(f"Filtered out consensus for {seed_seq['id']}_sf{subfamily_id}: sequence too short ({len(consensus_seq) if consensus_seq else 0}bp < 50bp)")
            return None
        
        # 计算改进程度
        original_length = len(seed_seq['sequence'])
        improved_length = len(consensus_seq)
        improvement_ratio = improved_length / original_length if original_length > 0 else 1.0
        
        # 创建共识记录
        consensus_record = {
            'id': f"{seed_seq['id']}_sf{subfamily_id}",
            'sequence': consensus_seq,
            'source_id': seed_seq['id'],
            'quality_class': seed_seq.get('quality_class', 'A'),
            'subfamily_id': subfamily_id,
            'copy_number': len(copies),
            'avg_identity': np.mean([c.get('identity', 0) for c in copies]),
            'length': len(consensus_seq),
            'original_length': original_length,
            'improvement_ratio': improvement_ratio,
            'characteristics': characteristics
        }
        
        # 添加TSD信息
        tsd_sequences = []
        for c in copies:
            if c.get('has_tsd', False) and c.get('tsd'):
                tsd_info = c.get('tsd')
                if isinstance(tsd_info, dict) and 'sequence' in tsd_info:
                    tsd_sequences.append(tsd_info['sequence'])
                elif isinstance(tsd_info, str):
                    tsd_sequences.append(tsd_info)
        
        if tsd_sequences:
            from collections import Counter
            tsd_counter = Counter(tsd_sequences)
            if tsd_counter:
                most_common_tsd = tsd_counter.most_common(1)[0][0]
                consensus_record['tsd'] = most_common_tsd
                consensus_record['tsd_support'] = tsd_counter[most_common_tsd] / len(copies)
        
        # 添加质量评分
        quality_score = calculate_consensus_quality(consensus_record, copies)
        consensus_record['quality_score'] = quality_score
        
        return consensus_record
        
    except Exception as e:
        import traceback
        logger.error(f"Error building consensus: {e}")
        logger.error(f"Full traceback: {traceback.format_exc()}")
        return None


def calculate_consensus_quality(consensus: Dict, copies: List[Dict]) -> float:
    """计算共识序列质量"""
    
    quality_components = []
    
    # 1. 拷贝数支持度
    copy_support = min(consensus['copy_number'] / 10, 1.0)
    quality_components.append(copy_support * 0.3)
    
    # 2. 序列一致性
    identity_score = consensus['avg_identity'] / 100.0
    quality_components.append(identity_score * 0.3)
    
    # 3. 长度改进
    if consensus['improvement_ratio'] > 1.0:
        improvement_score = min((consensus['improvement_ratio'] - 1.0) / 0.5, 1.0)
    else:
        improvement_score = consensus['improvement_ratio']
    quality_components.append(improvement_score * 0.2)
    
    # 4. TSD支持
    tsd_score = consensus.get('tsd_support', 0)
    quality_components.append(tsd_score * 0.1)
    
    # 5. 长度合理性
    length = consensus['length']
    if 100 <= length <= 5000:
        length_score = 1.0
    elif length < 100:
        length_score = length / 100
    else:
        length_score = max(0.5, 1.0 - (length - 5000) / 10000)
    quality_components.append(length_score * 0.1)
    
    return sum(quality_components)


def improve_with_single_copy(seed_seq: Dict, copy: Dict) -> Dict:
    """使用单个拷贝改进原始序列"""
    # 优化：优先选择更长的序列
    if copy['length'] > len(seed_seq['sequence']):
        # 拷贝更长，使用拷贝
        sequence = copy['sequence']
        logger.debug(f"Using longer copy: {copy['length']}bp vs original {len(seed_seq['sequence'])}bp")
    elif copy.get('identity', 0) > 80:
        # 拷贝质量更高
        sequence = copy['sequence']
    else:
        # 保留原始序列
        sequence = seed_seq['sequence']
    
    # 检查序列长度阈值
    if len(sequence) < 50:
        logger.debug(f"Filtered out single-copy improved sequence for {seed_seq['id']}: too short ({len(sequence)}bp < 50bp)")
        return None
    
    return {
        'id': f"{seed_seq['id']}_improved",
        'sequence': sequence,
        'source_id': seed_seq['id'],
        'quality_class': seed_seq.get('quality_class', 'A'),
        'copy_number': 1,
        'avg_identity': copy.get('identity', 0),
        'quality_score': 0.3,  # 低质量分数因为只有一个拷贝
        'note': 'single_copy_improvement'
    }


def improve_original_sequence(seed_seq: Dict, copies: List[Dict]) -> Dict:
    """使用拷贝信息改进原始序列"""
    # 优化：选择最长的高质量拷贝作为参考
    if copies:
        # 综合考虑长度和质量，长度权重更高
        best_copy = max(copies, key=lambda x: x['length'] * 0.7 + x.get('identity', 0) * x['length'] * 0.003)
        sequence = best_copy['sequence']
        avg_identity = np.mean([c.get('identity', 0) for c in copies])
        logger.debug(f"Selected best copy for improvement: {best_copy['length']}bp, {best_copy.get('identity', 0):.1f}% identity")
    else:
        sequence = seed_seq['sequence']
        avg_identity = 0
    
    # 检查序列长度阈值
    if len(sequence) < 50:
        logger.debug(f"Filtered out original improved sequence for {seed_seq['id']}: too short ({len(sequence)}bp < 50bp)")
        return None
    
    return {
        'id': f"{seed_seq['id']}_improved",
        'sequence': sequence,
        'source_id': seed_seq['id'],
        'quality_class': seed_seq.get('quality_class', 'A'),
        'copy_number': len(copies),
        'avg_identity': avg_identity,
        'quality_score': 0.4,
        'note': 'original_improved'
    }


def create_default_consensus(seed_seq: Dict) -> Dict:
    """创建默认共识（原始序列）"""
    # 检查序列长度阈值
    if len(seed_seq['sequence']) < 50:
        logger.debug(f"Filtered out default consensus for {seed_seq['id']}: too short ({len(seed_seq['sequence'])}bp < 50bp)")
        return None
    
    return {
        'id': f"{seed_seq['id']}_default",
        'sequence': seed_seq['sequence'],
        'source_id': seed_seq['id'],
        'quality_class': seed_seq.get('quality_class', 'A'),
        'copy_number': 0,
        'avg_identity': 0,
        'quality_score': 0.2,  # 最低质量
        'note': 'no_copies_found'
    }


# Phase 3 特有方法 - 添加在文件末尾

def _categorize_input_sequences(self, sequences: List[Dict]) -> Dict[str, List[Dict]]:
    """对输入序列进行分类"""
    categories = {
        'non_chimeric': [],
        'chimera_segments': [],
        'unsplit_chimeras': [],
        'failed_sequences': []
    }
    
    for seq in sequences:
        source = seq.get('source', 'unknown')
        
        if source == 'non_chimeric':
            categories['non_chimeric'].append(seq)
        elif source == 'chimera_split':
            if seq.get('split_from'):
                # 这是一个拆分片段
                categories['chimera_segments'].append(seq)
            else:
                # 这是一个未拆分的嵌合体
                categories['unsplit_chimeras'].append(seq)
        elif source == 'unsplit_chimera':
            categories['unsplit_chimeras'].append(seq)
        elif source == 'processing_failed':
            categories['failed_sequences'].append(seq)
        else:
            # 其他来源，默认作为普通序列处理
            categories['non_chimeric'].append(seq)
    
    return categories


def _build_consensus_parallel(self, sequences: List[Dict], stats: Dict) -> List[Dict]:
    """并行构建共识序列"""
    all_consensus = []
    
    if self.max_workers <= 1 or len(sequences) <= 1:
        logger.info("Using sequential consensus building")
        all_consensus = self._build_consensus_sequential(sequences, stats)
    else:
        logger.info(f"Using parallel consensus building with {self.max_workers} workers")
        # 复用原有的并行处理逻辑
        all_consensus = self._process_sequences_dynamic(sequences, stats)
    
    return all_consensus


def _build_consensus_sequential(self, sequences: List[Dict], stats: Dict) -> List[Dict]:
    """串行构建共识序列"""
    all_consensus = []
    
    for i, seq in enumerate(sequences, 1):
        logger.debug(f"Processing sequence {i}/{len(sequences)}: {seq['id']}")
        
        consensus_result = self._process_single_sequence_for_consensus(seq)
        
        if consensus_result:
            all_consensus.extend(consensus_result['consensus_sequences'])
            stats['processed_sequences'] += 1
            stats['single_consensus'] += len([c for c in consensus_result['consensus_sequences'] 
                                            if c.get('subfamily_count', 1) == 1])
            stats['multiple_consensus'] += len([c for c in consensus_result['consensus_sequences'] 
                                              if c.get('subfamily_count', 1) > 1])
        else:
            stats['failed_sequences'] += 1
        
        # 定期报告进度
        if i % 50 == 0 or i == len(sequences):
            logger.info(f"Progress: {i}/{len(sequences)} sequences processed, "
                       f"{len(all_consensus)} consensus generated")
    
    stats['total_consensus'] = len(all_consensus)
    return all_consensus


def _process_single_sequence_for_consensus(self, sequence: Dict) -> Optional[Dict]:
    """处理单个序列进行共识构建"""
    try:
        seq_id = sequence.get('id', 'unknown')
        seq_data = sequence.get('sequence', '')
        
        if not seq_data:
            logger.warning(f"Empty sequence for {seq_id}")
            return None
        
        # 使用扩展和改进方法（复用现有实现）
        result = expand_and_improve_sequence(sequence, self.config, self.mafft_threads)
        
        if result and result['consensus_count'] > 0:
            return {
                'consensus_sequences': result['consensus_list'],
                'method': f"expanded_{result['consensus_count']}_consensus"
            }
        else:
            # 回退到原序列
            logger.debug(f"Expansion failed for {seq_id}, using original sequence")
            consensus = self._create_single_consensus(sequence)
            return {
                'consensus_sequences': [consensus],
                'method': 'fallback_single'
            }
        
    except Exception as e:
        logger.error(f"Error processing sequence {sequence.get('id', 'unknown')} for consensus: {e}")
        return None


def _create_single_consensus(self, sequence: Dict) -> Dict:
    """创建单一共识序列"""
    consensus = {
        'id': f"{sequence['id']}_consensus",
        'sequence': sequence['sequence'],
        'subfamily_id': 1,
        'subfamily_count': 1,
        'avg_identity': 100.0,
        'avg_length': len(sequence['sequence']),
        'source_sequence': sequence['id'],
        'consensus_method': 'single_sequence',
        'quality_score': self._calculate_consensus_quality_score(100.0, 1)
    }
    
    return consensus


def _calculate_consensus_quality_score(self, avg_identity: float, copy_count: int) -> float:
    """计算共识序列质量分数"""
    # 基于平均相似性和拷贝数的综合评分
    identity_score = avg_identity / 100.0
    copy_score = min(1.0, np.log1p(copy_count) / np.log1p(10))  # 10个拷贝为满分
    
    # 加权平均
    quality_score = identity_score * 0.7 + copy_score * 0.3
    
    return quality_score


def _filter_and_grade_consensus(self, consensus_list: List[Dict], stats: Dict) -> List[Dict]:
    """过滤和分级共识序列"""
    filtered_consensus = []
    
    for consensus in consensus_list:
        quality_score = consensus.get('quality_score', 0)
        seq_length = len(consensus.get('sequence', ''))
        
        # 基本过滤
        if seq_length < self.config.min_length:
            continue
        
        # 质量分级
        if quality_score >= 0.8:
            consensus['quality_grade'] = 'high'
            stats['high_quality_consensus'] += 1
        elif quality_score >= 0.6:
            consensus['quality_grade'] = 'medium'
            stats['medium_quality_consensus'] += 1
        else:
            consensus['quality_grade'] = 'low'
            stats['low_quality_consensus'] += 1
        
        filtered_consensus.append(consensus)
    
    # 按质量分数排序
    filtered_consensus.sort(key=lambda x: x.get('quality_score', 0), reverse=True)
    
    logger.info(f"Filtered and graded {len(filtered_consensus)} consensus sequences")
    return filtered_consensus


def _deduplicate_consensus_sequences(self, consensus_sequences: List[Dict], stats: Dict) -> List[Dict]:
    """去除重复的共识序列"""
    if len(consensus_sequences) <= 1:
        return consensus_sequences
    
    logger.info(f"Deduplicating {len(consensus_sequences)} consensus sequences")
    
    # 按质量分数排序（高到低）
    sorted_consensus = sorted(consensus_sequences, 
                            key=lambda x: x.get('quality_score', 0), 
                            reverse=True)
    
    deduplicated = []
    processed_sequences = set()
    
    for consensus in sorted_consensus:
        seq_id = consensus['id']
        sequence = consensus['sequence']
        
        # 检查是否与已有序列重复
        is_duplicate = False
        
        for existing in deduplicated:
            similarity = self._calculate_sequence_similarity(
                sequence, existing['sequence']
            )
            
            # 如果相似度过高，认为是重复
            if similarity > 0.95:  # 95%相似度阈值
                is_duplicate = True
                logger.debug(f"Removing duplicate: {seq_id} (similarity {similarity:.3f} with {existing['id']})")
                break
        
        if not is_duplicate:
            deduplicated.append(consensus)
            processed_sequences.add(seq_id)
    
    removed_count = len(consensus_sequences) - len(deduplicated)
    stats['deduplicated_sequences'] = removed_count
    
    logger.info(f"Deduplication complete: removed {removed_count} duplicate sequences, "
               f"kept {len(deduplicated)} unique sequences")
    
    return deduplicated


def _calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
    """计算序列相似度"""
    if not seq1 or not seq2:
        return 0.0
    
    # 使用较短序列的长度
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    
    # 长度差异太大，不认为是重复
    if max_len / min_len > 1.2:
        return 0.0
    
    # 计算对齐的相似度
    matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
    
    # 考虑长度差异的惩罚
    length_penalty = min_len / max_len
    
    return (matches / min_len) * length_penalty


def _generate_phase4_compatible_output(self, consensus_sequences: List[Dict], 
                                      stats: Dict, sequence_categories: Dict) -> Dict[str, Any]:
    """生成Phase 4兼容的输出格式"""
    
    # 按质量分级分离序列
    high_quality = [seq for seq in consensus_sequences if seq.get('quality_grade') == 'high']
    medium_quality = [seq for seq in consensus_sequences if seq.get('quality_grade') == 'medium']
    low_quality = [seq for seq in consensus_sequences if seq.get('quality_grade') == 'low']
    
    # 生成masking library（95%去重，用于基因组屏蔽）
    masking_library = high_quality + medium_quality[:len(medium_quality)//2]  # 取一半medium质量
    
    # 生成analysis library（90%去重，包含更多序列用于分析）
    analysis_library = consensus_sequences
    
    result = {
        'masking_library': masking_library,
        'analysis_library': analysis_library,
        'statistics': stats,
        'phase3_complete': True,
        'sequence_categories': sequence_categories,  # 传递给Phase 4作参考
        'quality_distribution': {
            'high': len(high_quality),
            'medium': len(medium_quality),
            'low': len(low_quality)
        }
    }
    
    logger.info(f"Generated Phase 4 compatible output:")
    logger.info(f"  - Masking library: {len(masking_library)} sequences")
    logger.info(f"  - Analysis library: {len(analysis_library)} sequences")
    
    return result


# 为ConsensusBuilder类添加这些方法 - 在类定义后绑定
# _extract_sequences_from_phase2 already defined inside the class
ConsensusBuilder._categorize_input_sequences = _categorize_input_sequences
ConsensusBuilder._build_consensus_parallel = _build_consensus_parallel
ConsensusBuilder._build_consensus_sequential = _build_consensus_sequential
ConsensusBuilder._process_single_sequence_for_consensus = _process_single_sequence_for_consensus
ConsensusBuilder._create_single_consensus = _create_single_consensus
ConsensusBuilder._calculate_consensus_quality_score = _calculate_consensus_quality_score
ConsensusBuilder._filter_and_grade_consensus = _filter_and_grade_consensus
ConsensusBuilder._deduplicate_consensus_sequences = _deduplicate_consensus_sequences
ConsensusBuilder._calculate_sequence_similarity = _calculate_sequence_similarity
ConsensusBuilder._generate_phase4_compatible_output = _generate_phase4_compatible_output
