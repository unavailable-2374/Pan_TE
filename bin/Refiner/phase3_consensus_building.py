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
import re

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
        except Exception as e:
            available_memory_gb = 16  # 默认假设16GB可用
            total_memory_gb = 32
            logger.warning(f"Could not detect system memory ({e}), assuming 16GB available")
        
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
        
        # 统计信息 - 增强版本
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
            'low_quality_consensus': 0,
            'rejected_consensus': 0  # 新增：被拒绝的序列数
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
        
        # Log stratified processing path distribution
        if all_consensus:
            path_counts = {'fast_path_boundary_only': 0, 'minimal_path': 0, 'standard': 0, 'other': 0}
            for c in all_consensus:
                method = c.get('consensus_method', c.get('note', ''))
                if 'fast_path' in str(method):
                    path_counts['fast_path_boundary_only'] += 1
                elif 'minimal_path' in str(method):
                    path_counts['minimal_path'] += 1
                elif method in ('', 'no_copies_found', 'original_improved'):
                    path_counts['other'] += 1
                else:
                    path_counts['standard'] += 1
            total_paths = sum(path_counts.values())
            logger.info(f"Stratified processing path distribution ({total_paths} total):")
            for path_name, count in path_counts.items():
                pct = (count / total_paths * 100) if total_paths > 0 else 0
                logger.info(f"  - {path_name}: {count} ({pct:.1f}%)")

        # 质量过滤、分级和去重
        filtered_consensus = self._filter_and_grade_consensus(all_consensus, stats)
        deduplicated_consensus = self._deduplicate_consensus_sequences(filtered_consensus, stats)
        
        # 生成Phase 4兼容的输出结果
        result = self._generate_phase4_compatible_output(deduplicated_consensus, stats, sequence_categories)
        
        logger.info(f"Enhanced Phase 3 complete:")
        logger.info(f"  - Input: {stats['input_sequences']} sequences")
        logger.info(f"  - Processed: {stats['processed_sequences']} candidates")
        logger.info(f"  - Generated: {stats['total_consensus']} consensus sequences")
        logger.info(f"  - Quality distribution: High={stats['high_quality_consensus']}, "
                   f"Medium={stats['medium_quality_consensus']}, Low={stats['low_quality_consensus']}, "
                   f"Rejected={stats['rejected_consensus']}")
        logger.info(f"  - Enhancement features: Adaptive thresholds, Biological validation, Smart selection")
        
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
                    except Exception:
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
        except Exception:
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
                
                # mp.set_start_method('spawn') is called once at pipeline startup in main.py
            
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
                        except Exception:
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

    def _categorize_input_sequences(self, sequences: List[Dict]) -> Dict[str, List[Dict]]:
        """Categorize input sequences by source type."""
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
                    categories['chimera_segments'].append(seq)
                else:
                    categories['unsplit_chimeras'].append(seq)
            elif source == 'unsplit_chimera':
                categories['unsplit_chimeras'].append(seq)
            elif source == 'processing_failed':
                categories['failed_sequences'].append(seq)
            else:
                categories['non_chimeric'].append(seq)

        return categories

    def _build_consensus_parallel(self, sequences: List[Dict], stats: Dict) -> List[Dict]:
        """Build consensus sequences in parallel or sequentially."""
        if self.max_workers <= 1 or len(sequences) <= 1:
            logger.info("Using sequential consensus building")
            return self._build_consensus_sequential(sequences, stats)
        else:
            logger.info(f"Using parallel consensus building with {self.max_workers} workers")
            return self._process_sequences_dynamic(sequences, stats)

    def _build_consensus_sequential(self, sequences: List[Dict], stats: Dict) -> List[Dict]:
        """Build consensus sequences sequentially."""
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
            if i % 50 == 0 or i == len(sequences):
                logger.info(f"Progress: {i}/{len(sequences)} sequences processed, "
                           f"{len(all_consensus)} consensus generated")
        stats['total_consensus'] = len(all_consensus)
        return all_consensus

    def _process_single_sequence_for_consensus(self, sequence: Dict) -> Optional[Dict]:
        """Process single sequence for consensus building."""
        try:
            seq_id = sequence.get('id', 'unknown')
            seq_data = sequence.get('sequence', '')
            if not seq_data:
                logger.warning(f"Empty sequence for {seq_id}")
                return None
            result = expand_and_improve_sequence(sequence, self.config, self.mafft_threads)
            if result and result['consensus_count'] > 0:
                return {
                    'consensus_sequences': result['consensus_list'],
                    'method': f"expanded_{result['consensus_count']}_consensus"
                }
            else:
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
        """Create single consensus sequence."""
        return {
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

    def _calculate_consensus_quality_score(self, avg_identity: float, copy_count: int) -> float:
        """Calculate consensus quality score."""
        identity_score = avg_identity / 100.0
        copy_score = min(1.0, np.log1p(copy_count) / np.log1p(10))
        return identity_score * 0.7 + copy_score * 0.3

    def _filter_and_grade_consensus(self, consensus_list: List[Dict], stats: Dict) -> List[Dict]:
        """Filter and grade consensus sequences - distribution-aware version.

        Uses percentile-based dynamic trust levels instead of fixed thresholds.
        """
        if not consensus_list:
            return []

        copy_numbers = [c.get('copy_number', 0) for c in consensus_list]
        if copy_numbers:
            p25 = np.percentile(copy_numbers, 25)
            p50 = np.percentile(copy_numbers, 50)
            p75 = np.percentile(copy_numbers, 75)
            copy_mean = np.mean(copy_numbers)
            copy_max = max(copy_numbers)
        else:
            p25, p50, p75 = 2, 5, 15
            copy_mean, copy_max = 5, 20

        logger.info(f"Copy number distribution: P25={p25:.1f}, P50={p50:.1f}, P75={p75:.1f}, "
                   f"mean={copy_mean:.1f}, max={copy_max}")

        high_trust_threshold = max(p75, 3)
        medium_high_threshold = max(p50, 2)
        medium_threshold = max(p25, 1)

        logger.info(f"Dynamic trust thresholds: High>={high_trust_threshold:.1f}, "
                   f"MedHigh>={medium_high_threshold:.1f}, Med>={medium_threshold:.1f}")

        validated_consensus = []
        trust_stats = {'high': 0, 'medium_high': 0, 'medium': 0, 'low': 0}

        for consensus in consensus_list:
            seq_length = len(consensus.get('sequence', ''))
            if seq_length < self.config.min_length:
                continue
            copy_number = consensus.get('copy_number', 0)
            base_quality = consensus.get('quality_score', 0)
            quality_class = consensus.get('quality_class', 'B')
            sequence = consensus.get('sequence', '')

            if copy_number >= high_trust_threshold:
                trust_level = 'high'
            elif copy_number >= medium_high_threshold:
                trust_level = 'medium_high'
            elif copy_number >= medium_threshold:
                trust_level = 'medium'
            else:
                trust_level = 'low'
            trust_stats[trust_level] += 1
            consensus['trust_level'] = trust_level

            skip_validation = False
            validation_strictness = 'standard'

            if trust_level in ['high', 'medium_high']:
                skip_validation = True
            elif trust_level == 'medium':
                if quality_class in ['A', 'C_rescued'] or _has_te_structural_features_standalone(sequence):
                    skip_validation = True
                else:
                    validation_strictness = 'light'
            else:
                if quality_class == 'A' or _has_te_structural_features_standalone(sequence):
                    skip_validation = True
                elif seq_length >= 1500:
                    validation_strictness = 'light'

            if skip_validation:
                consensus['validations'] = {}
                consensus['validation_score'] = 1.0
                final_score = base_quality
            else:
                validations = _validate_consensus_quality(consensus)
                consensus['validations'] = validations
                raw_validation_score = sum(validations.values()) / len(validations) if validations else 0
                if validation_strictness == 'light':
                    final_score = base_quality * 0.8 + raw_validation_score * 0.2
                else:
                    final_score = base_quality * 0.7 + raw_validation_score * 0.3
                consensus['validation_score'] = raw_validation_score

            consensus['final_quality_score'] = final_score
            validated_consensus.append(consensus)

        logger.info(f"Trust level distribution: {trust_stats}")

        if not validated_consensus:
            return []

        quality_scores = [c.get('final_quality_score', 0) for c in validated_consensus]
        thresholds = _calculate_adaptive_thresholds(quality_scores)
        logger.info(f"Quality thresholds: High={thresholds['high']:.3f}, "
                   f"Medium={thresholds['medium']:.3f}, Low={thresholds['low']:.3f}")

        filtered_consensus = []
        rejected_count = 0

        for consensus in validated_consensus:
            final_score = consensus.get('final_quality_score', 0)
            validation_score = consensus.get('validation_score', 0)
            quality_class = consensus.get('quality_class', 'B')
            trust_level = consensus.get('trust_level', 'low')
            seq_length = len(consensus.get('sequence', ''))
            should_reject = False

            if trust_level in ['high', 'medium_high']:
                if final_score < 0.1:
                    should_reject = True
            elif trust_level == 'medium':
                if validation_score < 0.2 and final_score < thresholds['low'] * 0.7:
                    if seq_length < 1000 and quality_class not in ['A', 'C_rescued']:
                        should_reject = True
            else:
                if validation_score < 0.15 and final_score < thresholds['low'] * 0.8:
                    if seq_length < 1000 and quality_class not in ['A', 'C_rescued']:
                        should_reject = True

            if should_reject:
                rejected_count += 1
                continue

            if final_score >= thresholds['high']:
                consensus['quality_grade'] = 'high'
                stats['high_quality_consensus'] += 1
            elif final_score >= thresholds['medium']:
                consensus['quality_grade'] = 'medium'
                stats['medium_quality_consensus'] += 1
            else:
                consensus['quality_grade'] = 'low'
                stats['low_quality_consensus'] += 1
            filtered_consensus.append(consensus)

        stats['rejected_consensus'] = rejected_count
        filtered_consensus.sort(key=lambda x: x.get('final_quality_score', 0), reverse=True)
        logger.info(f"Filtering complete: {len(filtered_consensus)} passed, {rejected_count} rejected")
        return filtered_consensus

    def _deduplicate_consensus_sequences(self, consensus_sequences: List[Dict], stats: Dict) -> List[Dict]:
        """Deduplicate consensus sequences."""
        if len(consensus_sequences) <= 1:
            return consensus_sequences

        logger.info(f"Deduplicating {len(consensus_sequences)} consensus sequences")
        sorted_consensus = sorted(consensus_sequences,
                                  key=lambda x: x.get('quality_score', 0),
                                  reverse=True)
        deduplicated = []
        for consensus in sorted_consensus:
            is_duplicate = False
            for existing in deduplicated:
                similarity = self._calculate_sequence_similarity(
                    consensus['sequence'], existing['sequence'])
                if similarity > 0.95:
                    is_duplicate = True
                    break
            if not is_duplicate:
                deduplicated.append(consensus)

        removed_count = len(consensus_sequences) - len(deduplicated)
        stats['deduplicated_sequences'] = removed_count
        logger.info(f"Deduplication: removed {removed_count}, kept {len(deduplicated)}")
        return deduplicated

    def _calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity."""
        if not seq1 or not seq2:
            return 0.0
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))
        if max_len / min_len > 1.2:
            return 0.0
        matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
        length_penalty = min_len / max_len
        return (matches / min_len) * length_penalty

    def _generate_phase4_compatible_output(self, consensus_sequences: List[Dict],
                                           stats: Dict, sequence_categories: Dict) -> Dict[str, Any]:
        """Generate Phase 4 compatible output."""
        high_quality = [seq for seq in consensus_sequences if seq.get('quality_grade') == 'high']
        medium_quality = [seq for seq in consensus_sequences if seq.get('quality_grade') == 'medium']
        low_quality = [seq for seq in consensus_sequences if seq.get('quality_grade') == 'low']

        logger.info(f"Quality distribution: High={len(high_quality)}, "
                   f"Medium={len(medium_quality)}, Low={len(low_quality)}")

        masking_library = high_quality.copy()
        selected_medium = medium_quality
        if medium_quality:
            masking_library.extend(selected_medium)

        analysis_library = high_quality + medium_quality
        self._save_enhanced_analysis_library_fasta(analysis_library, stats)

        result = {
            'masking_library': masking_library,
            'analysis_library': analysis_library,
            'statistics': stats,
            'phase3_complete': True,
            'sequence_categories': sequence_categories,
            'quality_distribution': {
                'high': len(high_quality),
                'medium': len(medium_quality),
                'low': len(low_quality),
                'masking_selected': {
                    'high': len(high_quality),
                    'medium': len(selected_medium),
                    'low': 0
                }
            },
            'enhancement_metrics': {
                'adaptive_thresholds_used': True,
                'biological_validation': True,
                'quality_components_tracked': True
            }
        }

        logger.info(f"Phase 4 output: Masking={len(masking_library)}, Analysis={len(analysis_library)}")
        return result

    def _save_enhanced_analysis_library_fasta(self, analysis_library: List[Dict], stats: Dict):
        """Save analysis library as FASTA file."""
        output_dir = Path(self.config.output_dir)
        analysis_file = output_dir / "phase3_analysis_library.fa"
        stats_file = output_dir / "phase3_quality_statistics.txt"

        try:
            with open(analysis_file, 'w') as f:
                for seq in analysis_library:
                    seq_id = seq.get('id', 'unknown')
                    sequence = seq.get('sequence', '')
                    header_parts = [
                        f">{seq_id}",
                        f"grade={seq.get('quality_grade', 'unknown')}",
                        f"final_score={seq.get('final_quality_score', 0):.3f}",
                        f"base_score={seq.get('quality_score', 0):.3f}",
                        f"validation={seq.get('validation_score', 0):.3f}",
                        f"copies={seq.get('copy_number', 0)}",
                        f"length={len(sequence)}"
                    ]
                    if seq.get('tsd'):
                        header_parts.append(f"tsd={seq['tsd']}")
                    f.write(" ".join(header_parts) + "\n" + sequence + "\n")
            logger.info(f"Saved analysis library to: {analysis_file}")

            with open(stats_file, 'w') as f:
                f.write("Phase 3 Quality Statistics\n" + "=" * 50 + "\n\n")
                f.write("Overall Statistics:\n")
                for key, value in stats.items():
                    f.write(f"  {key}: {value}\n")
                f.write("\nQuality Component Analysis:\n")
                if analysis_library:
                    component_stats = defaultdict(list)
                    for seq in analysis_library:
                        for comp_name, comp_value in seq.get('quality_components', {}).items():
                            component_stats[comp_name].append(comp_value)
                    for comp_name, values in component_stats.items():
                        if values:
                            f.write(f"  {comp_name}: mean={np.mean(values):.3f}, "
                                    f"std={np.std(values):.3f}, "
                                    f"min={np.min(values):.3f}, max={np.max(values):.3f}\n")
            logger.info(f"Saved quality statistics to: {stats_file}")
        except Exception as e:
            logger.error(f"Failed to save analysis files: {e}")

    def _save_analysis_library_fasta(self, analysis_library: List[Dict]):
        """Backward-compatible wrapper for _save_enhanced_analysis_library_fasta."""
        self._save_enhanced_analysis_library_fasta(analysis_library, {})


def determine_processing_path(seed_seq: Dict, rm_hits: list, config: PipelineConfig) -> str:
    """
    Determine stratified processing path based on sequence confidence.

    Returns:
        'fast'     - High-copy, high-identity: skip MAFFT, boundary-extend only
        'minimal'  - Very low copy: use original sequence as-is
        'standard' - Everything else: full MAFFT rebuild
    """
    if not getattr(config, 'enable_stratified_processing', True):
        return 'standard'

    copy_number = len(rm_hits) if rm_hits else 0
    quality_class = seed_seq.get('quality_class', 'B')

    # Minimal path: 0-1 copies
    max_copies_minimal = getattr(config, 'minimal_path_max_copies', 1)
    if copy_number <= max_copies_minimal:
        return 'minimal'

    # Fast path: high copy + high identity
    min_copies_fast = getattr(config, 'fast_path_min_copies', 10)
    min_identity_fast = getattr(config, 'fast_path_min_identity', 75.0)

    if copy_number >= min_copies_fast:
        # Calculate average identity from rm_hits
        identities = [h.get('identity', 0) for h in rm_hits if 'identity' in h]
        avg_identity = sum(identities) / len(identities) if identities else 0

        if avg_identity >= min_identity_fast:
            # Additional check: structural features or quality class A
            sequence = seed_seq.get('sequence', '')
            if quality_class == 'A' or _has_te_structural_features_standalone(sequence):
                return 'fast'

    return 'standard'


def fast_path_consensus(seed_seq: Dict, copies: list, config: PipelineConfig) -> Dict:
    """
    Fast path: keep RepeatScout's original consensus, only extend boundaries.
    Skips expensive MAFFT MSA entirely.

    For high-copy, high-identity sequences where RepeatScout's consensus
    is already reliable.
    """
    from utils.sequence_utils import extend_consensus_boundaries

    sequence = seed_seq['sequence']
    original_length = len(sequence)

    # Extend boundaries using the longest genomic copy
    if copies:
        longest_copy = max(copies, key=lambda x: x.get('length', 0))
        reference_seq = longest_copy.get('extended_sequence', longest_copy.get('sequence', ''))
        if reference_seq and len(reference_seq) > len(sequence):
            extended = extend_consensus_boundaries(sequence, reference_seq)
            if len(extended) > len(sequence):
                logger.debug(f"{seed_seq['id']}: Fast path extended {len(sequence)}bp -> {len(extended)}bp")
                sequence = extended

    # Calculate quality score from copy statistics (no MAFFT needed)
    identities = [c.get('identity', 0) for c in copies]
    avg_identity = sum(identities) / len(identities) if identities else 0
    copy_number = len(copies)

    # Quality scoring for fast path
    identity_score = avg_identity / 100.0
    copy_score = min(1.0, np.log1p(copy_number) / np.log1p(50))
    quality_score = identity_score * 0.5 + copy_score * 0.3 + 0.2  # Base bonus for qualifying

    if len(sequence) < 50:
        return None

    return {
        'id': f"{seed_seq['id']}_sf0",
        'sequence': sequence,
        'source_id': seed_seq['id'],
        'quality_class': seed_seq.get('quality_class', 'A'),
        'subfamily_id': 0,
        'copy_number': copy_number,
        'avg_identity': avg_identity,
        'length': len(sequence),
        'original_length': original_length,
        'improvement_ratio': len(sequence) / original_length if original_length > 0 else 1.0,
        'quality_score': quality_score,
        'consensus_method': 'fast_path_boundary_only',
        'characteristics': {
            'copy_number': copy_number,
            'avg_identity': avg_identity
        },
        'note': f'fast_path:copies={copy_number},identity={avg_identity:.1f}'
    }


def _select_rmblast_matrix(config: PipelineConfig, divergence: int = 20, gc_pct: int = 41) -> str:
    """Select the best RepeatMasker scoring matrix for rmblastn.

    Args:
        config: Pipeline config with rmblast_matrix_dir
        divergence: Expected divergence level (14, 18, 20, 25)
        gc_pct: Genome GC percentage (35, 37, 39, 41, 43, 45, 47, 49, 51, 53)

    Returns:
        Matrix name (e.g. '20p41g.matrix') or empty string if not available.
    """
    import os
    matrix_dir = getattr(config, 'rmblast_matrix_dir', '')
    if not matrix_dir:
        return ""

    # Try exact match first
    matrix_name = f"{divergence}p{gc_pct}g.matrix"
    if os.path.isfile(os.path.join(matrix_dir, matrix_name)):
        return matrix_name

    # Try nearest GC content
    for gc in [41, 43, 39, 45, 37, 47, 35, 49, 51, 53]:
        matrix_name = f"{divergence}p{gc}g.matrix"
        if os.path.isfile(os.path.join(matrix_dir, matrix_name)):
            return matrix_name

    return ""


def get_genome_copies_blastn(seed_seq: Dict, config: PipelineConfig) -> list:
    """
    RMBlast-based copy recruitment with RepeatMasker-level sensitivity.

    Uses rmblastn (if available) with RepeatMasker's scoring matrices and
    gap parameters, matching the 'general_search_parameters' recipe:
      - matrix: 20p##g.matrix (divergence-optimized)
      - word_size: 9 (sensitive seed matching)
      - gapopen: 24, gapextend: 6 (conservative for TE detection)
      - complexity_adjust: on (Shannon entropy weighting)
      - min_raw_gapped_score: 225 (RepeatMasker minimum)

    Falls back to standard blastn with best-effort parameters if rmblastn
    or matrices are not available.

    Returns list of hit dicts compatible with rm_hits format.
    """
    import subprocess
    import tempfile
    import os

    db_path = getattr(config, 'genome_blast_db', '')
    if not db_path:
        return []

    # Determine whether to use rmblastn or standard blastn
    rmblastn_exe = getattr(config, 'rmblastn_exe', 'rmblastn')
    matrix_dir = getattr(config, 'rmblast_matrix_dir', '')
    use_rmblast = bool(matrix_dir) and os.path.isfile(rmblastn_exe if os.path.sep in rmblastn_exe
                                                       else (subprocess.run(['which', rmblastn_exe],
                                                             capture_output=True, text=True).stdout.strip() or '/nonexistent'))

    query_path = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as qf:
            qf.write(f">{seed_seq['id']}\n{seed_seq['sequence']}\n")
            query_path = qf.name

        if use_rmblast:
            # --- RMBlast path: RepeatMasker-level sensitivity ---
            matrix_name = _select_rmblast_matrix(config)
            if not matrix_name:
                matrix_name = "20p41g.matrix"  # default

            # RepeatMasker general_search recipe translation:
            # gap_initValue=-30, ins_gap_extValue=-6
            # rmblastn gapopen = abs(gap_init - ins_gap_ext) = abs(-30-(-6)) = 24
            # rmblastn gapextend = abs(ins_gap_ext) = 6
            # minscore=225, bandwidth=14
            # xdrop_ungap = minscore*2 = 450
            # xdrop_gap_final = minscore = 225
            # xdrop_gap = minscore/2 = 112
            cmd = [
                rmblastn_exe,
                "-query", query_path,
                "-db", db_path,
                "-outfmt", "6 sseqid sstart send pident length score qstart qend",
                "-matrix", matrix_name,
                "-complexity_adjust",
                "-gapopen", "24",
                "-gapextend", "6",
                "-word_size", "9",
                "-min_raw_gapped_score", "225",
                "-xdrop_ungap", "450",
                "-xdrop_gap_final", "225",
                "-xdrop_gap", "112",
                "-dust", "no",
                "-mask_level", "90",
                "-num_threads", "1",
                "-num_alignments", "9999999",
            ]

            # Set BLASTMAT environment so rmblastn can find matrix files
            env = os.environ.copy()
            env['BLASTMAT'] = matrix_dir

            logger.debug(f"{seed_seq['id']}: Using rmblastn with {matrix_name}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120, env=env)
        else:
            # --- Standard blastn fallback: best-effort parameters ---
            # Cannot match RM sensitivity without scoring matrices, but use
            # smaller word size and lower evalue for better sensitivity
            cmd = [
                config.blastn_exe,
                "-query", query_path,
                "-db", db_path,
                "-outfmt", "6 sseqid sstart send pident length score qstart qend",
                "-evalue", "1e-5",
                "-word_size", "9",
                "-gapopen", "2",
                "-gapextend", "1",
                "-reward", "1",
                "-penalty", "-1",
                "-perc_identity", "55",
                "-max_target_seqs", "500",
                "-num_threads", "1",
                "-dust", "yes"
            ]
            logger.debug(f"{seed_seq['id']}: Using standard blastn (lower sensitivity)")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

        if result.returncode != 0:
            stderr_snippet = result.stderr[:300] if result.stderr else ""
            logger.debug(f"Search failed for {seed_seq['id']} (rc={result.returncode}): {stderr_snippet}")
            return []

        hits = []
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) >= 8:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                identity = float(fields[3])
                length = int(fields[4])

                # Ensure start < end
                if start > end:
                    start, end = end, start

                hits.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'identity': identity,
                    'length': length,
                    'score': int(float(fields[5])),
                    'source': 'rmblastn' if use_rmblast else 'blastn'
                })

        if hits:
            logger.debug(f"{seed_seq['id']}: Found {len(hits)} copies via {'rmblastn' if use_rmblast else 'blastn'}")

        return hits

    except subprocess.TimeoutExpired:
        logger.debug(f"Search timeout for {seed_seq['id']}")
        return []
    except Exception as e:
        logger.debug(f"Search error for {seed_seq['id']}: {e}")
        return []
    finally:
        if query_path and os.path.exists(query_path):
            os.unlink(query_path)


def expand_and_improve_sequence(seed_seq: Dict, config: PipelineConfig, mafft_threads: int = 1) -> Dict:
    """
    扩展和改进单个候选序列
    目标：提高序列质量和代表性，不是减少数量

    Enhanced with stratified processing paths:
    - 'fast': High-copy sequences skip MAFFT, only extend boundaries
    - 'minimal': Very low-copy sequences use original as default consensus
    - 'standard': Full MAFFT rebuild (existing behavior)

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

        # Get rm_hits for stratification decision
        rm_hits = seed_seq.get('rm_hits', [])

        # 2. Determine processing path (stratified optimization)
        processing_path = determine_processing_path(seed_seq, rm_hits, config)

        if processing_path == 'minimal':
            # Minimal path: very low copy, use default consensus
            improved = create_default_consensus(seed_seq)
            if improved:
                improved['consensus_method'] = 'minimal_path'
                improved['note'] = 'minimal_path:low_copy'
                result['consensus_list'].append(improved)
                result['consensus_count'] = 1
            return result

        if not copies or len(copies) < 2:
            # 拷贝太少，直接使用原始序列
            if copies and len(copies) == 1:
                improved = improve_with_single_copy(seed_seq, copies[0])
            else:
                improved = create_default_consensus(seed_seq)

            if improved:
                result['consensus_list'].append(improved)
                result['consensus_count'] = 1
            return result

        logger.debug(f"{seed_seq['id']}: Found {len(copies)} copies, path={processing_path}")

        if processing_path == 'fast':
            # Fast path: skip MAFFT, boundary extension only
            consensus = fast_path_consensus(seed_seq, copies, config)
            if consensus:
                result['consensus_list'].append(consensus)
                result['consensus_count'] = 1
            else:
                # Fallback to standard if fast path fails
                improved = improve_original_sequence(seed_seq, copies)
                if improved:
                    result['consensus_list'].append(improved)
                    result['consensus_count'] = 1
            return result

        # Standard path: full MAFFT rebuild (existing behavior)
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

    Priority: Phase 1 rm_hits > BLASTN > RepeatMasker fallback
    """
    from utils.sequence_utils import extract_sequence_from_genome, detect_tsd

    # 优先使用Phase 1提供的RepeatMasker结果
    if 'rm_hits' in seed_seq and seed_seq['rm_hits']:
        rm_hits = seed_seq['rm_hits']
        logger.debug(f"{seed_seq['id']}: Using {len(rm_hits)} hits from Phase 1")
    else:
        # Try BLASTN first (faster than RepeatMasker)
        rm_hits = []
        blast_db = getattr(config, 'genome_blast_db', '')
        if blast_db:
            logger.debug(f"{seed_seq['id']}: Trying BLASTN copy recruitment")
            rm_hits = get_genome_copies_blastn(seed_seq, config)

        if not rm_hits:
            # Fall back to RepeatMasker
            logger.debug(f"{seed_seq['id']}: Running RepeatMasker (no Phase 1 or BLASTN results)")
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

    # Two-stage sampling: RM-based coarse sampling + sequence-based fine sampling
    # extract_sequence_from_genome now uses indexed access (pysam/samtools) for O(1) per call
    from two_stage_sampling import apply_two_stage_sampling
    logger.info(f"Applying two-stage sampling for {len(rm_hits)} RM hits of {seed_seq['id']}")
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


def _compute_pairwise_identity(seq_a: str, seq_b: str) -> float:
    """Compute pairwise sequence identity between two sequences.

    Uses a fast k-mer Jaccard estimate for speed.  For short sequences
    (<200bp) falls back to direct character comparison after simple
    alignment anchoring.

    Returns identity as fraction 0-1.
    """
    if not seq_a or not seq_b:
        return 0.0

    a = seq_a.upper()
    b = seq_b.upper()

    # For very short sequences, use direct comparison
    min_len = min(len(a), len(b))
    max_len = max(len(a), len(b))

    if min_len < 200:
        # Align shorter to longer using sliding window, pick best match
        if len(a) > len(b):
            a, b = b, a  # a is shorter
        best_matches = 0
        scan_range = len(b) - len(a) + 1
        step = max(1, scan_range // 100)
        for offset in range(0, scan_range, step):
            matches = sum(1 for i in range(len(a)) if a[i] == b[offset + i])
            if matches > best_matches:
                best_matches = matches
        return best_matches / len(a) if len(a) > 0 else 0.0

    # K-mer Jaccard for longer sequences (k=8)
    k = 8
    kmers_a = set()
    kmers_b = set()
    for i in range(len(a) - k + 1):
        kmer = a[i:i + k]
        if 'N' not in kmer:
            kmers_a.add(kmer)
    for i in range(len(b) - k + 1):
        kmer = b[i:i + k]
        if 'N' not in kmer:
            kmers_b.add(kmer)

    if not kmers_a or not kmers_b:
        return 0.0

    intersection = len(kmers_a & kmers_b)
    union = len(kmers_a | kmers_b)
    jaccard = intersection / union if union > 0 else 0.0

    # Convert Jaccard similarity to approximate sequence identity
    # Using Mash-like formula: identity ≈ 1 + (1/k) * ln(2*J / (1+J))
    import math
    if jaccard > 0:
        identity = 1.0 + (1.0 / k) * math.log(2.0 * jaccard / (1.0 + jaccard))
        identity = max(0.0, min(1.0, identity))
    else:
        identity = 0.0

    return identity


def identify_subfamilies_for_expansion(copies: List[Dict], characteristics: Dict) -> List[List[Dict]]:
    """
    Identify subfamilies using actual pairwise sequence identity.

    Previous version used RM-identity-to-seed as a proxy, which
    conflates distinct subfamilies with similar divergence from the
    seed. Now computes real pairwise distances between copies.
    """
    if len(copies) < 5:
        return [copies]

    n = len(copies)

    # Build pairwise distance matrix using actual sequence identity
    distance_matrix = np.zeros((n, n))

    for i in range(n):
        seq_i = copies[i].get('sequence', '')
        for j in range(i + 1, n):
            seq_j = copies[j].get('sequence', '')

            # Compute actual pairwise identity
            pairwise_id = _compute_pairwise_identity(seq_i, seq_j)
            seq_distance = 1.0 - pairwise_id

            # Length difference as secondary signal (weight 20%)
            len_i, len_j = copies[i]['length'], copies[j]['length']
            if max(len_i, len_j) > 0:
                length_distance = 1.0 - min(len_i, len_j) / max(len_i, len_j)
            else:
                length_distance = 1.0

            total_distance = 0.8 * seq_distance + 0.2 * length_distance
            distance_matrix[i, j] = total_distance
            distance_matrix[j, i] = total_distance

    logger.debug(f"Built pairwise distance matrix for {n} copies "
                 f"(mean dist={np.mean(distance_matrix[np.triu_indices(n, k=1)]):.3f})")

    # Hierarchical clustering
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    condensed_distances = squareform(distance_matrix)
    linkage_matrix = linkage(condensed_distances, method='average')

    # Dynamic threshold based on distance distribution
    distances = distance_matrix[np.triu_indices_from(distance_matrix, k=1)]
    if len(distances) > 0:
        threshold = np.percentile(distances, 75)
        threshold = max(0.15, min(0.4, threshold))
    else:
        threshold = 0.25

    clusters = fcluster(linkage_matrix, threshold, criterion='distance')

    # Organize into subfamilies
    subfamilies = defaultdict(list)
    for i, cluster_id in enumerate(clusters):
        subfamilies[cluster_id].append(copies[i])

    # Merge small subfamilies (<2 copies) into nearest cluster
    final_subfamilies = []
    for subfamily in subfamilies.values():
        if len(subfamily) >= 2:
            final_subfamilies.append(subfamily)
        else:
            if final_subfamilies:
                final_subfamilies[0].extend(subfamily)
            else:
                final_subfamilies.append(subfamily)

    logger.debug(f"Identified {len(final_subfamilies)} subfamilies "
                 f"(sizes: {[len(sf) for sf in final_subfamilies]})")

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
    """计算共识序列质量 - 增强版本
    
    改进：
    1. 重新平衡权重，提升TSD和生物学特征的重要性
    2. 降低拷贝数权重，使用对数缩放
    3. 增加序列复杂度和边界质量评估
    4. 考虑序列一致性分布
    """
    
    # 1. TSD支持度 (25% - 提升权重)
    tsd_score = consensus.get('tsd_support', 0)
    if tsd_score > 0.8:  # 强TSD支持给予额外奖励
        tsd_score = min(1.0, tsd_score * 1.1)
    
    # 2. 序列一致性与分布 (20% - 考虑分布)
    avg_identity = consensus.get('avg_identity', 0) / 100.0
    if copies and len(copies) > 1:
        identity_values = [c.get('identity', consensus.get('avg_identity', 0)) for c in copies]
        identity_std = np.std(identity_values) / 100.0
        # 惩罚高方差（表示家族内差异大）
        identity_score = avg_identity * (1 - min(identity_std * 0.3, 0.5))
    else:
        identity_score = avg_identity
    
    # 3. 拷贝数支持 (20% - 降低权重，使用对数缩放)
    copy_number = consensus.get('copy_number', 1)
    copy_score = min(1.0, np.log1p(copy_number) / np.log1p(50))  # 50个拷贝达满分
    
    # 4. 序列复杂度 (15% - 新增)
    sequence = consensus.get('sequence', '')
    complexity_score = _calculate_sequence_complexity_score(sequence)
    
    # 5. 长度完整性 (10%)
    length_score = _evaluate_length_completeness(consensus, copies)
    
    # 6. 边界质量 (10% - 新增)
    boundary_score = _evaluate_boundary_quality(sequence)
    
    # 计算加权总分
    final_score = (
        tsd_score * 0.25 +
        identity_score * 0.20 +
        copy_score * 0.20 +
        complexity_score * 0.15 +
        length_score * 0.10 +
        boundary_score * 0.10
    )
    
    # 保存质量组件用于调试和分析
    consensus['quality_components'] = {
        'tsd': tsd_score,
        'identity': identity_score,
        'copy_support': copy_score,
        'complexity': complexity_score,
        'length': length_score,
        'boundary': boundary_score
    }
    
    return final_score


def _has_te_structural_features_standalone(sequence: str) -> bool:
    """
    Standalone function to check for TE structural features.
    Used in Phase 3 for quality scoring of single/low-copy sequences.

    Features checked:
    1. Terminal Inverted Repeats (TIR)
    2. Long Terminal Repeats (LTR)
    3. Poly-A tail
    """
    if not sequence or len(sequence) < 50:
        return False

    sequence = sequence.upper()
    seq_len = len(sequence)

    # 1. Check for TIR (10-40bp, allowing 2 mismatches)
    for tir_len in range(10, min(41, seq_len // 4)):
        left_tir = sequence[:tir_len]
        right_tir = sequence[-tir_len:]
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        right_tir_rc = ''.join(complement.get(b, 'N') for b in reversed(right_tir))
        mismatches = sum(1 for a, b in zip(left_tir, right_tir_rc) if a != b)
        if mismatches <= 2:
            return True

    # 2. Check for LTR (80-500bp, >75% identity)
    if seq_len >= 260:  # min 80*2 + 100
        for ltr_len in range(80, min(501, seq_len // 3)):
            left_ltr = sequence[:ltr_len]
            right_ltr = sequence[-ltr_len:]
            matches = sum(1 for a, b in zip(left_ltr, right_ltr) if a == b)
            if matches / ltr_len >= 0.75:
                return True

    # 3. Check for poly-A tail
    tail_region = sequence[-50:] if seq_len >= 50 else sequence
    if 'AAAAAA' in tail_region:
        return True

    # 4. Check for TG...CA pattern (LTR signature)
    if sequence[:2] == 'TG' and sequence[-2:] == 'CA':
        return True

    return False


def _calculate_sequence_complexity_simple(sequence: str) -> float:
    """
    Simple sequence complexity calculation using Shannon entropy.
    Returns value between 0 (simple) and 1 (complex).
    """
    if not sequence or len(sequence) < 10:
        return 0.0

    sequence = sequence.upper()

    # Calculate Shannon entropy
    from collections import Counter
    import math

    base_counts = Counter(sequence)
    total = sum(base_counts.values())

    entropy = 0.0
    for count in base_counts.values():
        if count > 0:
            freq = count / total
            entropy -= freq * math.log2(freq)

    # Normalize to 0-1 (max entropy for 4 bases is 2.0)
    normalized_entropy = entropy / 2.0

    # Also check for simple repeats
    simple_repeat_ratio = 0.0
    for pattern in ['AAAA', 'TTTT', 'GGGG', 'CCCC', 'ATAT', 'TATA', 'GCGC', 'CGCG']:
        simple_repeat_ratio += sequence.count(pattern) * len(pattern) / len(sequence)

    # Combine entropy and simple repeat penalty
    complexity = normalized_entropy * (1.0 - min(simple_repeat_ratio, 0.5))

    return min(1.0, max(0.0, complexity))


def _calculate_sequence_complexity_score(sequence: str) -> float:
    """计算序列复杂度分数"""
    if not sequence or len(sequence) < 10:
        return 0.0
    
    sequence = sequence.upper()
    
    # 1. 计算低复杂度比例
    low_complexity_ratio = _calculate_low_complexity_ratio(sequence)
    
    # 2. 计算串联重复比例
    tandem_ratio = _calculate_tandem_repeat_ratio(sequence)
    
    # 3. 计算熵值
    entropy = _calculate_sequence_entropy(sequence)
    
    # 综合评分（高复杂度得高分）
    complexity_score = (
        (1.0 - low_complexity_ratio) * 0.4 +
        (1.0 - tandem_ratio) * 0.3 +
        (entropy / 2.0) * 0.3  # 熵值最大约为2.0
    )
    
    return min(1.0, max(0.0, complexity_score))


def _calculate_low_complexity_ratio(sequence: str) -> float:
    """计算低复杂度区域比例"""
    if not sequence:
        return 1.0
    
    # 简单的低复杂度检测：单核苷酸和二核苷酸重复
    patterns = [
        r'A{5,}', r'T{5,}', r'G{5,}', r'C{5,}',  # 单核苷酸
        r'(AT){4,}', r'(TA){4,}', r'(GC){4,}', r'(CG){4,}',  # 二核苷酸
        r'(AG){4,}', r'(GA){4,}', r'(TC){4,}', r'(CT){4,}',
        r'(AC){4,}', r'(CA){4,}', r'(TG){4,}', r'(GT){4,}'
    ]
    
    low_complexity_length = 0
    covered_positions = set()
    
    for pattern in patterns:
        for match in re.finditer(pattern, sequence):
            start, end = match.span()
            for pos in range(start, end):
                if pos not in covered_positions:
                    covered_positions.add(pos)
                    low_complexity_length += 1
    
    return low_complexity_length / len(sequence)


def _calculate_tandem_repeat_ratio(sequence: str) -> float:
    """计算串联重复比例"""
    if not sequence or len(sequence) < 10:
        return 0.0
    
    # 检测3-6bp的串联重复
    tandem_length = 0
    covered = set()
    
    for repeat_len in range(3, 7):
        for i in range(len(sequence) - repeat_len * 2):
            if i in covered:
                continue
            
            unit = sequence[i:i+repeat_len]
            j = i + repeat_len
            repeat_count = 1
            
            while j + repeat_len <= len(sequence) and sequence[j:j+repeat_len] == unit:
                repeat_count += 1
                j += repeat_len
            
            if repeat_count >= 3:  # 至少重复3次
                for pos in range(i, i + repeat_len * repeat_count):
                    covered.add(pos)
                tandem_length += repeat_len * repeat_count
    
    return len(covered) / len(sequence)


def _calculate_sequence_entropy(sequence: str) -> float:
    """计算序列的信息熵"""
    if not sequence:
        return 0.0
    
    # 计算碱基频率
    base_counts = Counter(sequence)
    seq_len = len(sequence)
    
    entropy = 0.0
    for count in base_counts.values():
        if count > 0:
            freq = count / seq_len
            entropy -= freq * np.log2(freq)
    
    return entropy


def _evaluate_length_completeness(consensus: Dict, copies: List[Dict]) -> float:
    """评估长度完整性"""
    cons_length = len(consensus.get('sequence', ''))
    
    if cons_length < 50:
        return 0.0
    elif cons_length > 20000:
        return 0.5  # 过长可能有问题
    
    # 理想长度范围
    if 100 <= cons_length <= 5000:
        base_score = 1.0
    elif cons_length < 100:
        base_score = cons_length / 100
    else:
        base_score = max(0.5, 1.0 - (cons_length - 5000) / 15000)
    
    # 如果有改进比例，额外加分
    improvement_ratio = consensus.get('improvement_ratio', 1.0)
    if improvement_ratio > 1.0:
        improvement_bonus = min(0.2, (improvement_ratio - 1.0) * 0.4)
        base_score = min(1.0, base_score + improvement_bonus)
    
    return base_score


def _evaluate_boundary_quality(sequence: str) -> float:
    """评估序列边界质量"""
    if not sequence or len(sequence) < 20:
        return 0.0
    
    # 检查5'和3'端的质量
    boundary_len = min(20, len(sequence) // 10)
    
    start_seq = sequence[:boundary_len].upper()
    end_seq = sequence[-boundary_len:].upper()
    
    # 检查N含量
    start_n_ratio = start_seq.count('N') / len(start_seq)
    end_n_ratio = end_seq.count('N') / len(end_seq)
    
    # 检查低复杂度
    start_low_complex = _calculate_low_complexity_ratio(start_seq)
    end_low_complex = _calculate_low_complexity_ratio(end_seq)
    
    # 综合评分
    boundary_score = (
        (1.0 - start_n_ratio) * 0.25 +
        (1.0 - end_n_ratio) * 0.25 +
        (1.0 - start_low_complex) * 0.25 +
        (1.0 - end_low_complex) * 0.25
    )
    
    return boundary_score


def _validate_consensus_quality(consensus: Dict) -> Dict[str, bool]:
    """Multi-level quality validation - relaxed for better sensitivity

    Key changes:
    1. Relaxed entropy threshold (1.0 -> 0.8)
    2. Relaxed GC content range (0.15-0.85 -> 0.10-0.90)
    3. Relaxed N content threshold (0.1 -> 0.15)
    4. Relaxed ambiguity threshold (0.15 -> 0.20)
    """
    validations = {}
    sequence = consensus.get('sequence', '')

    if not sequence:
        return {'sequence_exists': False}

    seq_len = len(sequence)

    # 1. Length validation (unchanged)
    validations['length_valid'] = 50 <= seq_len <= 20000

    # 2. Complexity validation - relaxed tandem ratio threshold
    tandem_ratio = _calculate_tandem_repeat_ratio(sequence)
    validations['complexity_valid'] = tandem_ratio < 0.6  # was 0.5

    # 3. GC content validation - relaxed range
    gc_content = _calculate_gc_content(sequence)
    validations['gc_valid'] = 0.10 < gc_content < 0.90  # was 0.15-0.85

    # 4. N content validation - relaxed threshold
    n_ratio = sequence.upper().count('N') / len(sequence)
    validations['n_content_valid'] = n_ratio < 0.15  # was 0.1

    # 5. TSD validation (if claimed) - unchanged
    if consensus.get('tsd'):
        validations['tsd_valid'] = _verify_tsd_pattern(consensus['tsd'])

    # 6. Ambiguous base validation - relaxed threshold
    ambiguous_ratio = _calculate_ambiguous_base_ratio(sequence)
    validations['ambiguity_valid'] = ambiguous_ratio < 0.20  # was 0.15

    # 7. Minimum entropy validation - relaxed threshold
    # Note: DNA entropy max is ~2.0, threshold 0.8 allows more diversity
    entropy = _calculate_sequence_entropy(sequence.upper())
    validations['entropy_valid'] = entropy > 0.8  # was 1.0

    return validations


def _calculate_gc_content(sequence: str) -> float:
    """计算GC含量"""
    if not sequence:
        return 0.0
    
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    at_count = sequence.count('A') + sequence.count('T')
    total = gc_count + at_count
    
    return gc_count / total if total > 0 else 0.0


def _verify_tsd_pattern(tsd: str) -> bool:
    """验证TSD模式"""
    if not tsd or not isinstance(tsd, str):
        return False
    
    # TSD应该是2-20bp的序列
    if not (2 <= len(tsd) <= 20):
        return False
    
    # 不应该全是N或低复杂度
    if tsd.upper().count('N') / len(tsd) > 0.5:
        return False
    
    # 不应该是单一碱基重复
    if len(set(tsd.upper())) == 1:
        return False
    
    return True


def _calculate_ambiguous_base_ratio(sequence: str) -> float:
    """计算模糊碱基比例"""
    if not sequence:
        return 1.0
    
    sequence = sequence.upper()
    ambiguous_bases = set('RYSWKMBDHVN')
    ambiguous_count = sum(1 for base in sequence if base in ambiguous_bases)
    
    return ambiguous_count / len(sequence)


def _calculate_adaptive_thresholds(quality_scores: List[float]) -> Dict[str, float]:
    """基于数据分布的自适应阈值"""
    
    if not quality_scores:
        return {'high': 0.8, 'medium': 0.6, 'low': 0.3}
    
    # 移除异常值
    scores_array = np.array(quality_scores)
    q1 = np.percentile(scores_array, 25)
    q3 = np.percentile(scores_array, 75)
    iqr = q3 - q1
    
    # 过滤异常值
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    filtered_scores = scores_array[(scores_array >= lower_bound) & (scores_array <= upper_bound)]
    
    if len(filtered_scores) < 3:
        filtered_scores = scores_array
    
    # 使用分位数方法
    if len(filtered_scores) >= 10:
        # 足够的数据点，使用三分位
        high_threshold = np.percentile(filtered_scores, 70)
        medium_threshold = np.percentile(filtered_scores, 40)
        low_threshold = np.percentile(filtered_scores, 10)
    else:
        # 数据点较少，使用更保守的方法
        mean_score = np.mean(filtered_scores)
        std_score = np.std(filtered_scores)
        
        high_threshold = mean_score + 0.5 * std_score
        medium_threshold = mean_score
        low_threshold = mean_score - 0.5 * std_score
    
    # 应用合理的边界限制
    high_threshold = max(0.55, min(0.85, high_threshold))
    medium_threshold = max(0.40, min(0.65, medium_threshold))
    low_threshold = max(0.25, min(0.45, low_threshold))
    
    # 确保阈值的顺序正确
    if medium_threshold >= high_threshold:
        medium_threshold = high_threshold - 0.1
    if low_threshold >= medium_threshold:
        low_threshold = medium_threshold - 0.1
    
    return {
        'high': high_threshold,
        'medium': medium_threshold,
        'low': low_threshold
    }


def _select_representative_medium_sequences(medium_sequences: List[Dict]) -> List[Dict]:
    """基于多样性和互补性智能选择medium序列"""
    
    if len(medium_sequences) <= 5:
        return medium_sequences
    
    logger.info(f"Intelligently selecting from {len(medium_sequences)} medium quality sequences")
    
    selected = []
    remaining = medium_sequences.copy()
    
    # 1. 按亚家族分组，每组选择最佳代表
    subfamily_groups = defaultdict(list)
    for seq in remaining:
        subfamily = seq.get('subfamily_id', 'default')
        subfamily_groups[subfamily].append(seq)
    
    for subfamily, seqs in subfamily_groups.items():
        if seqs:
            # 选择最终质量分数×长度最高的
            best = max(seqs, 
                      key=lambda x: x.get('final_quality_score', 0) * len(x.get('sequence', '')))
            selected.append(best)
            remaining.remove(best)
    
    # 2. 基于长度分布选择互补序列
    length_ranges = [
        (0, 200, 'very_short'),
        (200, 500, 'short'),
        (500, 1500, 'medium_length'),
        (1500, 3000, 'long'),
        (3000, 10000, 'very_long'),
        (10000, float('inf'), 'extra_long')
    ]
    
    for min_len, max_len, category in length_ranges:
        range_seqs = [s for s in remaining 
                     if min_len <= len(s.get('sequence', '')) < max_len]
        if range_seqs:
            # 每个长度范围选择质量最好的
            best_in_range = max(range_seqs, 
                               key=lambda x: x.get('final_quality_score', 0))
            selected.append(best_in_range)
            remaining.remove(best_in_range)
            logger.debug(f"Selected {best_in_range['id']} for {category} range")
    
    # 3. 基于序列特征的多样性选择
    target_count = min(len(medium_sequences) // 2, len(medium_sequences))
    
    while len(selected) < target_count and remaining:
        # 计算每个剩余序列与已选序列的最小相似度
        diversity_scores = []
        
        for seq in remaining:
            if not selected:
                # 如果还没有选择任何序列，使用质量分数
                diversity_score = seq.get('final_quality_score', 0)
            else:
                # 计算与已选序列的最小相似度
                min_similarity = min(
                    _calculate_sequence_feature_similarity(seq, sel)
                    for sel in selected
                )
                # 多样性分数 = (1 - 相似度) × 质量分数
                diversity_score = (1 - min_similarity) * seq.get('final_quality_score', 0)
            
            diversity_scores.append((diversity_score, seq))
        
        # 选择多样性分数最高的
        if diversity_scores:
            diversity_scores.sort(key=lambda x: x[0], reverse=True)
            best_diverse = diversity_scores[0][1]
            selected.append(best_diverse)
            remaining.remove(best_diverse)
    
    logger.info(f"Selected {len(selected)} representative medium sequences")
    return selected


def _calculate_sequence_feature_similarity(seq1: Dict, seq2: Dict) -> float:
    """计算两个序列的特征相似度"""
    
    # 长度相似度
    len1 = len(seq1.get('sequence', ''))
    len2 = len(seq2.get('sequence', ''))
    length_sim = min(len1, len2) / max(len1, len2) if max(len1, len2) > 0 else 0
    
    # GC含量相似度
    gc1 = _calculate_gc_content(seq1.get('sequence', ''))
    gc2 = _calculate_gc_content(seq2.get('sequence', ''))
    gc_diff = abs(gc1 - gc2)
    gc_sim = 1 - gc_diff
    
    # 质量分数相似度
    score1 = seq1.get('final_quality_score', 0)
    score2 = seq2.get('final_quality_score', 0)
    score_sim = 1 - abs(score1 - score2)
    
    # 亚家族相似度
    subfamily_sim = 1.0 if seq1.get('subfamily_id') == seq2.get('subfamily_id') else 0.0
    
    # 加权平均
    feature_similarity = (
        length_sim * 0.3 +
        gc_sim * 0.2 +
        score_sim * 0.2 +
        subfamily_sim * 0.3
    )
    
    return feature_similarity


def _select_valuable_low_sequences(low_sequences: List[Dict]) -> List[Dict]:
    """选择有特殊价值的低质量序列"""
    
    if not low_sequences:
        return []
    
    valuable = []
    
    # 选择特别长的序列（可能是完整但质量较低的TE）
    long_sequences = [s for s in low_sequences 
                     if len(s.get('sequence', '')) > 3000]
    
    if long_sequences:
        # 选择最长的前3个
        long_sequences.sort(key=lambda x: len(x.get('sequence', '')), reverse=True)
        valuable.extend(long_sequences[:3])
    
    # 选择有强TSD支持的序列
    tsd_sequences = [s for s in low_sequences 
                    if s.get('tsd_support', 0) > 0.7 and s not in valuable]
    
    if tsd_sequences:
        # 选择TSD支持最强的前2个
        tsd_sequences.sort(key=lambda x: x.get('tsd_support', 0), reverse=True)
        valuable.extend(tsd_sequences[:2])
    
    if valuable:
        logger.info(f"Selected {len(valuable)} valuable low-quality sequences for inclusion")
    
    return valuable


def improve_with_single_copy(seed_seq: Dict, copy: Dict) -> Dict:
    """
    Improve original sequence using single copy.

    Key change: Single-copy sequences are no longer penalized with fixed low score.
    Rationale: RepeatScout already filtered for copy number (≥10 k-mer occurrences).
    A single RepeatMasker hit may indicate:
    - Ancient TE with divergent copies not detected by RM
    - Young TE just starting to amplify
    - TE with structural features that validate it
    """
    # Select best sequence based on length and quality
    if copy['length'] > len(seed_seq['sequence']):
        sequence = copy['sequence']
        logger.debug(f"Using longer copy: {copy['length']}bp vs original {len(seed_seq['sequence'])}bp")
    elif copy.get('identity', 0) > 80:
        sequence = copy['sequence']
    else:
        sequence = seed_seq['sequence']

    # Check minimum length threshold
    if len(sequence) < 50:
        logger.debug(f"Filtered out single-copy improved sequence for {seed_seq['id']}: too short ({len(sequence)}bp < 50bp)")
        return None

    # Dynamic quality scoring for single-copy sequences
    # Base score starts at 0.35 (not 0.3 death sentence)
    base_score = 0.35
    score_reasons = []

    seq_length = len(sequence)

    # Bonus for length - longer sequences more likely to be real TEs
    if seq_length >= 2000:
        base_score += 0.2
        score_reasons.append("very_long")
    elif seq_length >= 1000:
        base_score += 0.15
        score_reasons.append("long")
    elif seq_length >= 500:
        base_score += 0.1
        score_reasons.append("medium_length")

    # Bonus for structural features
    if _has_te_structural_features_standalone(sequence):
        base_score += 0.15
        score_reasons.append("structural_features")

    # Bonus for sequence complexity (non-simple repeat)
    complexity = _calculate_sequence_complexity_simple(sequence)
    if complexity >= 0.7:
        base_score += 0.1
        score_reasons.append("high_complexity")
    elif complexity >= 0.5:
        base_score += 0.05
        score_reasons.append("medium_complexity")

    # Bonus if quality class is good
    quality_class = seed_seq.get('quality_class', 'B')
    if quality_class == 'A':
        base_score += 0.1
        score_reasons.append("A_class")
    elif quality_class == 'C_rescued':
        base_score += 0.05
        score_reasons.append("rescued")

    final_score = min(0.75, base_score)  # Cap at 0.75 for single copy

    return {
        'id': f"{seed_seq['id']}_improved",
        'sequence': sequence,
        'source_id': seed_seq['id'],
        'quality_class': quality_class,
        'copy_number': 1,
        'avg_identity': copy.get('identity', 0),
        'quality_score': final_score,
        'note': f'single_copy_improvement:{",".join(score_reasons)}'
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


# STRATIFIED BIOLOGICAL FILTERING INTEGRATION
# Override the default _filter_and_grade_consensus and _save_enhanced_analysis_library_fasta
# with stratified versions if available.  The methods are now defined inside the class body
# (no more monkey-patching via setattr).
try:
    from .phase3_stratified_integration import (
        _filter_and_grade_consensus_stratified,
        _save_enhanced_analysis_library_fasta_stratified
    )
except ImportError:
    try:
        from phase3_stratified_integration import (
            _filter_and_grade_consensus_stratified,
            _save_enhanced_analysis_library_fasta_stratified
        )
    except ImportError:
        _filter_and_grade_consensus_stratified = None
        _save_enhanced_analysis_library_fasta_stratified = None

if _filter_and_grade_consensus_stratified is not None:
    ConsensusBuilder._filter_and_grade_consensus = _filter_and_grade_consensus_stratified
    ConsensusBuilder._save_enhanced_analysis_library_fasta = _save_enhanced_analysis_library_fasta_stratified
    logger.info("Stratified biological filtering enabled")
