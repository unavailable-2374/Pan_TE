"""
Phase 2: 嵌合体序列拆分（重构版本）
专注功能：识别和拆分嵌合体序列，为Phase 3的共识构建做准备
输入：Phase 1筛选后的序列（包含嵌合体标记）
输出：拆分后的序列片段，准备进行共识构建
"""

import logging
import gc
import numpy as np
from typing import Dict, List, Tuple, Any, Optional
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp
import time
import psutil
import os
import tempfile

from config import PipelineConfig
from utils.robust_runner import RobustRunner

logger = logging.getLogger(__name__)


class ChimeraSplitter:
    """
    Phase 2: 嵌合体检测和拆分处理
    
    主要功能：
    1. 接收Phase 1的所有候选序列
    2. 逐一检测序列是否为嵌合体
    3. 对嵌合体进行保守的拆分
    4. 将处理后的所有序列传递给Phase 3
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
            available_memory_gb = 16
            total_memory_gb = 32
            logger.warning("Could not detect system memory, assuming 16GB available")
        
        # 设置并行处理workers - 嵌合体拆分主要是I/O密集型（RepeatMasker调用）
        cpu_count = mp.cpu_count()
        memory_per_thread = 0.3  # GB per thread for I/O intensive tasks
        max_workers_by_memory = max(1, int(available_memory_gb / memory_per_thread))
        
        # I/O密集型任务可以使用更多线程，直接使用配置的线程数
        # 但需要考虑内存限制
        self.max_workers = min(config.threads, max_workers_by_memory)
        
        # 如果用户配置的线程数较大，提示实际使用的线程数
        if config.threads > self.max_workers:
            logger.warning(f"Requested {config.threads} threads, but limited to {self.max_workers} due to memory constraints")
        
        # 每个线程使用单线程工具避免嵌套并行
        self.threads_per_worker = 1
        
        logger.info(f"Phase 2 Chimera Splitter initialized:")
        logger.info(f"  - Requested threads: {config.threads}")
        logger.info(f"  - Actual workers: {self.max_workers} threads")
        logger.info(f"  - Available memory: {available_memory_gb:.1f}GB")
        logger.info(f"  - Memory constraint: {max_workers_by_memory} threads")
        logger.info(f"  - CPU count: {cpu_count} cores")
    
    def run(self, phase1_output: Dict) -> Dict[str, Any]:
        """执行Phase 2主流程 - 处理所有候选序列的嵌合体检测和拆分"""
        
        # 获取Phase 1识别的候选序列
        candidate_sequences = phase1_output.get('consensus_candidates', [])
        
        # 保存Phase 1的RepeatMasker结果供重用
        self.phase1_rm_results = phase1_output.get('rm_detailed_results', {})
        
        logger.info(f"Phase 2 input: {len(candidate_sequences)} consensus candidates")
        
        # 检查Phase 1 RepeatMasker结果的可用性
        available_rm_results = len(self.phase1_rm_results)
        reusable_sequences = sum(1 for seq in candidate_sequences if seq['id'] in self.phase1_rm_results)
        logger.info(f"Phase 1 RepeatMasker results: {available_rm_results} available, "
                   f"{reusable_sequences} applicable to candidate sequences")
        
        # 统计信息
        stats = {
            'input_candidates': len(candidate_sequences),
            'chimeric_sequences': 0,       # 检测到的嵌合体
            'non_chimeric_sequences': 0,   # 非嵌合体序列
            'split_sequences': 0,          # 成功拆分的序列
            'unsplit_sequences': 0,        # 未拆分的序列
            'output_segments': 0,          # 拆分后的片段数
            'failed_sequences': 0,         # 处理失败的序列
            'phase1_reused': 0,            # 重用Phase 1结果的数量
            'new_analysis': 0,             # 新进行RepeatMasker分析的数量
            'reuse_success_rate': 0.0      # Phase 1结果重用成功率
        }
        
        # 处理所有候选序列，检测嵌合体并拆分
        processed_sequences = []
        if candidate_sequences:
            logger.info(f"Processing {len(candidate_sequences)} candidate sequences for chimera detection")
            processed_sequences = self._process_candidate_sequences(candidate_sequences)
            
            # 统计结果
            for result in processed_sequences:
                if result.get('source') == 'chimera_split':
                    stats['output_segments'] += 1
                    if result.get('split_from'):
                        # 这是一个拆分片段
                        pass
                    else:
                        # 这是一个未拆分的嵌合体
                        stats['unsplit_sequences'] += 1
                else:
                    # 非嵌合体序列
                    stats['non_chimeric_sequences'] += 1
                
                # 统计分析来源
                analysis_source = result.get('_analysis_source', 'unknown')
                if 'phase1_reused' in analysis_source:
                    stats['phase1_reused'] += 1
                elif 'new_analysis' in analysis_source:
                    stats['new_analysis'] += 1
            
            # 计算统计数据
            unique_split_sources = set()
            for result in processed_sequences:
                if result.get('split_from'):
                    unique_split_sources.add(result['split_from'])
            
            stats['chimeric_sequences'] = len([seq for seq in candidate_sequences 
                                             if any(res.get('split_from') == seq['id'] or 
                                                   (res.get('note', '').startswith('unsplit_chimera') and res['id'] == seq['id'])
                                                   for res in processed_sequences)])
            stats['split_sequences'] = len(unique_split_sources)
            
            total_processed = stats['phase1_reused'] + stats['new_analysis']
            if total_processed > 0:
                stats['reuse_success_rate'] = stats['phase1_reused'] / total_processed
            
            logger.info(f"Chimera detection and splitting results:")
            logger.info(f"  - Input candidates: {stats['input_candidates']}")
            logger.info(f"  - Detected chimeric sequences: {stats['chimeric_sequences']}")
            logger.info(f"  - Successfully split sequences: {stats['split_sequences']}")
            logger.info(f"  - Non-chimeric sequences: {stats['non_chimeric_sequences']}")
            logger.info(f"  - Total output sequences: {len(processed_sequences)}")
            logger.info(f"RepeatMasker analysis efficiency:")
            logger.info(f"  - Phase 1 results reused: {stats['phase1_reused']}")
            logger.info(f"  - New analyses performed: {stats['new_analysis']}")
            logger.info(f"  - Reuse rate: {stats['reuse_success_rate']:.1%}")
        
        # 准备输出：所有处理后的序列
        output_sequences = processed_sequences
        stats['output_total'] = len(output_sequences)
        
        # 生成输出结果
        result = {
            'sequences': output_sequences,
            'statistics': stats,
            'phase2_complete': True
        }
        
        logger.info(f"Phase 2 complete:")
        logger.info(f"  - Input: {stats['input_candidates']} candidate sequences")
        logger.info(f"  - Output: {stats['output_total']} sequences for Phase 3")
        logger.info(f"  - Chimera processing: {stats['chimeric_sequences']} detected, {stats['split_sequences']} split")
        logger.info(f"  - Processing rate: {stats['output_total']/max(1,stats['input_candidates']):.2f}x expansion")
        
        return result
    
    def _process_candidate_sequences(self, candidate_sequences: List[Dict]) -> List[Dict]:
        """处理所有候选序列，检测嵌合体并拆分"""
        if not candidate_sequences:
            return []
        
        # 根据线程数和序列数决定是否使用并行处理
        use_parallel = len(candidate_sequences) > 1 and self.max_workers > 1
        
        logger.info(f"Candidate sequence processing strategy:")
        logger.info(f"  - Sequences: {len(candidate_sequences)}")
        logger.info(f"  - Max workers: {self.max_workers}")
        logger.info(f"  - Use parallel: {use_parallel}")
        
        if not use_parallel:
            logger.info("Using sequential processing")
            return self._process_candidates_sequential(candidate_sequences)
        else:
            logger.info(f"Using parallel processing with {self.max_workers} threads")
            return self._process_candidates_parallel(candidate_sequences)
    
    def _process_candidates_sequential(self, candidate_sequences: List[Dict]) -> List[Dict]:
        """串行处理候选序列"""
        processed_sequences = []
        
        for i, candidate in enumerate(candidate_sequences, 1):
            logger.debug(f"Processing candidate {i}/{len(candidate_sequences)}: "
                        f"{candidate['id']} ({len(candidate['sequence'])}bp)")
            
            result = self._process_single_candidate(candidate)
            processed_sequences.extend(result)
            
            # 定期报告进度
            if i % 10 == 0 or i == len(candidate_sequences):
                logger.info(f"Progress: {i}/{len(candidate_sequences)} candidates processed")
        
        return processed_sequences
    
    def _process_candidates_parallel(self, candidate_sequences: List[Dict]) -> List[Dict]:
        """并行处理候选序列 - 使用ThreadPoolExecutor"""
        processed_sequences = []
        processed_count = 0
        total_count = len(candidate_sequences)
        start_time = time.time()
        
        # 进度报告设置
        last_report_time = time.time()
        report_interval = 30  # 每30秒报告一次
        
        def report_progress():
            nonlocal last_report_time
            current_time = time.time()
            if current_time - last_report_time > report_interval or processed_count == total_count:
                percent = (processed_count / total_count) * 100
                total_outputs = len(processed_sequences)
                elapsed = current_time - start_time
                rate = processed_count / elapsed if elapsed > 0 else 0
                logger.info(f"Progress: {processed_count}/{total_count} ({percent:.1f}%) - "
                           f"{total_outputs} outputs - {rate:.1f} seq/s")
                last_report_time = current_time
        
        # 动态调整批大小以更好地利用线程
        if len(candidate_sequences) <= self.max_workers:
            batch_size = len(candidate_sequences)
        else:
            batch_size = max(self.max_workers * 2, min(200, len(candidate_sequences)))
        
        total_batches = (len(candidate_sequences) + batch_size - 1) // batch_size
        
        logger.info(f"Processing {len(candidate_sequences)} candidates in {total_batches} batches")
        logger.info(f"  - Batch size: {batch_size}")
        logger.info(f"  - Using {self.max_workers} parallel threads")
        
        # 分批处理序列
        for batch_idx in range(total_batches):
            start_idx = batch_idx * batch_size
            end_idx = min(start_idx + batch_size, len(candidate_sequences))
            current_batch = candidate_sequences[start_idx:end_idx]
            
            logger.info(f"Processing batch {batch_idx + 1}/{total_batches} ({len(current_batch)} candidates)")
            
            actual_workers = min(self.max_workers, len(current_batch))
            with ThreadPoolExecutor(max_workers=actual_workers) as executor:
                # 提交任务
                future_to_candidate = {
                    executor.submit(self._process_single_candidate, candidate): candidate
                    for candidate in current_batch
                }
                
                # 处理完成的任务
                batch_processed = 0
                for future in as_completed(future_to_candidate):
                    candidate = future_to_candidate[future]
                    try:
                        result = future.result(timeout=300)  # 5分钟超时
                        processed_sequences.extend(result)
                        processed_count += 1
                        batch_processed += 1
                        
                        # 定期报告进度
                        report_progress()
                        
                    except TimeoutError:
                        logger.warning(f"Timeout processing candidate {candidate.get('id', 'unknown')}")
                        # 超时的序列保持原样
                        fallback_result = self._handle_processing_failure(candidate, "timeout")
                        processed_sequences.extend(fallback_result)
                        processed_count += 1
                    
                    except Exception as e:
                        logger.error(f"Error processing candidate {candidate.get('id', 'unknown')}: {e}")
                        # 失败的序列保持原样
                        fallback_result = self._handle_processing_failure(candidate, "error")
                        processed_sequences.extend(fallback_result)
                        processed_count += 1
            
            # 批次完成后清理内存
            logger.debug(f"Completed batch {batch_idx + 1}/{total_batches}, "
                        f"processed {batch_processed}/{len(current_batch)} candidates")
            gc.collect()
        
        # 最终统计
        total_time = time.time() - start_time
        outputs_generated = len(processed_sequences)
        
        logger.info(f"Parallel candidate processing completed in {total_time:.1f}s:")
        logger.info(f"  - Processed {total_count} candidates ({total_count/total_time:.1f} seq/s)")
        logger.info(f"  - Generated {outputs_generated} output sequences ({outputs_generated/total_time:.1f} outputs/s)")
        
        return processed_sequences
    
    def _process_single_candidate(self, candidate: Dict) -> List[Dict]:
        """处理单个候选序列 - 检测嵌合体并决定是否拆分"""
        try:
            candidate_id = candidate.get('id', 'unknown')
            sequence = candidate.get('sequence', '')
            
            if not sequence:
                logger.warning(f"Empty sequence for candidate {candidate_id}")
                return self._handle_processing_failure(candidate, "empty_sequence")
            
            # 检测是否为嵌合体
            chimera_analysis = self._analyze_for_chimera(candidate)
            
            if chimera_analysis['is_chimeric']:
                logger.debug(f"Detected chimeric sequence: {candidate_id}")
                
                if chimera_analysis.get('is_splittable', False):
                    # 执行拆分
                    segments = self._split_sequence_by_breakpoints(candidate, chimera_analysis['breakpoints'])
                    if segments and len(segments) > 1:
                        logger.debug(f"Split {candidate_id} into {len(segments)} segments")
                        # 标记为拆分结果
                        for segment in segments:
                            segment['source'] = 'chimera_split'
                            segment['_analysis_source'] = chimera_analysis.get('analysis_source', 'unknown')
                        return segments
                    else:
                        # 拆分失败，保持原序列
                        return self._handle_unsplit_chimera(candidate, "split_failed")
                else:
                    # 是嵌合体但不可拆分，保持原序列
                    return self._handle_unsplit_chimera(candidate, chimera_analysis.get('reason', 'not_splittable'))
            else:
                # 非嵌合体，直接传递给Phase 3
                logger.debug(f"Non-chimeric sequence: {candidate_id}")
                result = dict(candidate)  # 复制以避免修改原数据
                result['source'] = 'non_chimeric'
                result['split_from'] = None
                result['_analysis_source'] = chimera_analysis.get('analysis_source', 'unknown')
                return [result]
                
        except Exception as e:
            logger.error(f"Error processing candidate {candidate.get('id', 'unknown')}: {e}")
            return self._handle_processing_failure(candidate, f"processing_error: {str(e)}")
    
    def _analyze_for_chimera(self, candidate: Dict) -> Dict[str, Any]:
        """分析序列是否为嵌合体"""
        try:
            candidate_id = candidate.get('id', 'unknown')
            
            # 优先检查是否可以今Phase 1的RepeatMasker结果分析嵌合体
            phase1_result = self._try_analyze_with_phase1_results(candidate)
            if phase1_result is not None:
                logger.debug(f"Using Phase 1 RepeatMasker results for {candidate_id}")
                phase1_result['analysis_source'] = 'phase1_reused'
                return phase1_result
            else:
                # 如果Phase 1结果不可用，进行新的RepeatMasker自比对分析
                logger.debug(f"Running new RepeatMasker analysis for {candidate_id}")
                chimera_result = self._analyze_chimera_with_repeatmasker(candidate)
                chimera_result['analysis_source'] = 'new_analysis'
                return chimera_result
                
        except Exception as e:
            logger.error(f"Error analyzing chimera for {candidate.get('id', 'unknown')}: {e}")
            return {
                'is_chimeric': False,
                'is_splittable': False,
                'breakpoints': [],
                'reason': f'analysis_failed: {str(e)}',
                'analysis_source': 'error'
            }
    
    def _handle_processing_failure(self, candidate: Dict, reason: str) -> List[Dict]:
        """处理失败的候选序列 - 保持原样传递给Phase 3"""
        result = dict(candidate)  # 复制以避免修改原数据
        result['source'] = 'processing_failed'
        result['split_from'] = None
        result['processing_failure_reason'] = reason
        result['note'] = f'processing_failed: {reason}'
        return [result]
    
    def _split_chimeric_sequences(self, chimeric_sequences: List[Dict]) -> List[Dict]:
        """分离嵌合体序列 - 使用线程池并行处理"""
        if not chimeric_sequences:
            return []
        
        # 根据线程数和序列数决定是否使用并行处理
        use_parallel = len(chimeric_sequences) > 1 and self.max_workers > 1
        
        logger.info(f"Chimeric splitting strategy:")
        logger.info(f"  - Sequences: {len(chimeric_sequences)}")
        logger.info(f"  - Max workers: {self.max_workers}")
        logger.info(f"  - Use parallel: {use_parallel}")
        
        if not use_parallel:
            logger.info("Using sequential processing")
            return self._split_chimeric_sequences_sequential(chimeric_sequences)
        else:
            logger.info(f"Using parallel processing with {self.max_workers} threads")
            return self._split_chimeric_sequences_parallel(chimeric_sequences)
    
    def _split_chimeric_sequences_sequential(self, chimeric_sequences: List[Dict]) -> List[Dict]:
        """串行分离嵌合体序列"""
        separated_sequences = []
        
        for i, chimera in enumerate(chimeric_sequences, 1):
            logger.debug(f"Processing chimeric sequence {i}/{len(chimeric_sequences)}: "
                        f"{chimera['id']} ({len(chimera['sequence'])}bp)")
            
            result = self._process_single_chimera(chimera)
            separated_sequences.extend(result)
            
            # 定期报告进度
            if i % 10 == 0 or i == len(chimeric_sequences):
                logger.info(f"Progress: {i}/{len(chimeric_sequences)} sequences processed")
        
        return separated_sequences
    
    def _split_chimeric_sequences_parallel(self, chimeric_sequences: List[Dict]) -> List[Dict]:
        """并行分离嵌合体序列 - 使用ThreadPoolExecutor"""
        separated_sequences = []
        processed_count = 0
        total_count = len(chimeric_sequences)
        start_time = time.time()
        
        # 进度报告设置
        last_report_time = time.time()
        report_interval = 30  # 每30秒报告一次
        
        def report_progress():
            nonlocal last_report_time
            current_time = time.time()
            if current_time - last_report_time > report_interval or processed_count == total_count:
                percent = (processed_count / total_count) * 100
                total_segments = len(separated_sequences)
                avg_segments = total_segments / processed_count if processed_count > 0 else 0
                elapsed = current_time - start_time
                rate = processed_count / elapsed if elapsed > 0 else 0
                logger.info(f"Progress: {processed_count}/{total_count} ({percent:.1f}%) - "
                           f"{total_segments} segments - {rate:.1f} seq/s")
                last_report_time = current_time
        
        # 动态调整批大小以更好地利用线程
        # 如果序列数量少于线程数，不需要分批
        if len(chimeric_sequences) <= self.max_workers:
            batch_size = len(chimeric_sequences)
        else:
            # 确保批大小至少是线程数的2倍，以充分利用并行
            batch_size = max(self.max_workers * 2, min(200, len(chimeric_sequences)))
        
        total_batches = (len(chimeric_sequences) + batch_size - 1) // batch_size
        
        logger.info(f"Processing {len(chimeric_sequences)} sequences in {total_batches} batches")
        logger.info(f"  - Batch size: {batch_size}")
        logger.info(f"  - Using {self.max_workers} parallel threads")
        
        # 分批处理序列
        for batch_idx in range(total_batches):
            start_idx = batch_idx * batch_size
            end_idx = min(start_idx + batch_size, len(chimeric_sequences))
            current_batch = chimeric_sequences[start_idx:end_idx]
            
            logger.info(f"Processing batch {batch_idx + 1}/{total_batches} ({len(current_batch)} sequences)")
            
            # 始终使用配置的最大线程数，除非批次很小
            actual_workers = min(self.max_workers, len(current_batch))
            with ThreadPoolExecutor(max_workers=actual_workers) as executor:
                # 提交任务
                future_to_chimera = {
                    executor.submit(self._process_single_chimera, chimera): chimera
                    for chimera in current_batch
                }
                
                # 处理完成的任务
                batch_processed = 0
                for future in as_completed(future_to_chimera):
                    chimera = future_to_chimera[future]
                    try:
                        result = future.result(timeout=300)  # 5分钟超时
                        separated_sequences.extend(result)
                        processed_count += 1
                        batch_processed += 1
                        
                        # 定期报告进度
                        report_progress()
                        
                    except TimeoutError:
                        logger.warning(f"Timeout processing chimeric sequence {chimera.get('id', 'unknown')}")
                        # 超时的序列标记为未分割
                        classified = self._handle_unsplit_chimera(chimera, "timeout")
                        separated_sequences.extend(classified)
                        processed_count += 1
                    
                    except Exception as e:
                        logger.error(f"Error processing chimeric sequence {chimera.get('id', 'unknown')}: {e}")
                        # 失败的序列标记为未分割
                        classified = self._handle_unsplit_chimera(chimera, "error")
                        separated_sequences.extend(classified)
                        processed_count += 1
                        
                        # 调试信息
                        logger.debug(f"Detailed error for {chimera.get('id', 'unknown')}: {str(e)}")
            
            # 批次完成后清理内存
            logger.debug(f"Completed batch {batch_idx + 1}/{total_batches}, "
                        f"processed {batch_processed}/{len(current_batch)} sequences")
            gc.collect()
        
        # 最终统计
        total_time = time.time() - start_time
        segments_generated = len(separated_sequences)
        
        logger.info(f"Parallel chimeric splitting completed in {total_time:.1f}s:")
        logger.info(f"  - Processed {total_count} chimeric sequences ({total_count/total_time:.1f} seq/s)")
        logger.info(f"  - Generated {segments_generated} segments ({segments_generated/total_time:.1f} segments/s)")
        logger.info(f"  - Average segments per sequence: {segments_generated/total_count:.2f}")
        
        return separated_sequences
    
    def _process_single_chimera(self, chimera: Dict) -> List[Dict]:
        """处理单个嵌合体序列 - 保守策略：只有高度确定时才拆分"""
        try:
            chimera_id = chimera.get('id', 'unknown')
            sequence = chimera.get('sequence', '')
            
            if not sequence:
                logger.warning(f"Empty sequence for chimera {chimera_id}")
                return self._handle_unsplit_chimera(chimera, "empty_sequence")
            
            # 优先检查是否可以从Phase 1的RepeatMasker结果分析嵌合体
            phase1_result = self._try_analyze_with_phase1_results(chimera)
            if phase1_result is not None:
                logger.debug(f"Using Phase 1 RepeatMasker results for {chimera_id}")
                split_result = phase1_result
                # 标记为重用Phase 1结果
                chimera['_analysis_source'] = 'phase1_reused'
            else:
                # 如果Phase 1结果不可用，进行新的RepeatMasker自比对分析
                logger.debug(f"Running new RepeatMasker analysis for {chimera_id}")
                split_result = self._analyze_chimera_with_repeatmasker(chimera)
                # 标记为新分析
                chimera['_analysis_source'] = 'new_analysis'
            
            if split_result['is_splittable']:
                # 执行拆分
                segments = self._split_sequence_by_breakpoints(chimera, split_result['breakpoints'])
                logger.debug(f"Split {chimera_id} into {len(segments)} segments")
                return segments
            else:
                # 不拆分，标记为完整序列
                logger.debug(f"Sequence {chimera_id} not split: {split_result['reason']}")
                return self._handle_unsplit_chimera(chimera, split_result['reason'])
                
        except Exception as e:
            logger.error(f"Error processing chimera {chimera.get('id', 'unknown')}: {e}")
            return self._handle_unsplit_chimera(chimera, f"processing_error: {str(e)}")
    
    def _try_analyze_with_phase1_results(self, chimera: Dict) -> Optional[Dict[str, Any]]:
        """尝试使用Phase 1的RepeatMasker结果分析嵌合体"""
        chimera_id = chimera.get('id', 'unknown')
        
        # 检查是否有Phase 1结果
        if not hasattr(self, 'phase1_rm_results') or chimera_id not in self.phase1_rm_results:
            return None
        
        rm_data = self.phase1_rm_results[chimera_id]
        
        # 检查是否有足够的hit信息用于断点分析
        hits = rm_data.get('hits', [])
        if len(hits) < 2:
            logger.debug(f"Insufficient hits ({len(hits)}) in Phase 1 results for {chimera_id}")
            return None
        
        try:
            # 从Phase 1的RepeatMasker hits中提取断点
            breakpoints = self._extract_breakpoints_from_hits(hits, len(chimera['sequence']))
            
            if len(breakpoints) >= 2:
                return {
                    'is_chimeric': True,
                    'is_splittable': True,
                    'breakpoints': breakpoints,
                    'reason': f'phase1_analysis_{len(breakpoints)}_breakpoints'
                }
            else:
                return {
                    'is_chimeric': False,
                    'is_splittable': False,
                    'breakpoints': [],
                    'reason': 'phase1_insufficient_breakpoints'
                }
                
        except Exception as e:
            logger.debug(f"Error analyzing Phase 1 results for {chimera_id}: {e}")
            return None
    
    def _extract_breakpoints_from_hits(self, hits: List[Dict], sequence_length: int) -> List[int]:
        """从RepeatMasker hits中提取潜在断点"""
        breakpoints = []
        
        try:
            # 按位置排序hits
            sorted_hits = sorted(hits, key=lambda x: x.get('query_start', 0))
            
            # 寻找hits之间的间隙作为断点
            for i in range(len(sorted_hits) - 1):
                current_end = sorted_hits[i].get('query_end', 0)
                next_start = sorted_hits[i + 1].get('query_start', 0)
                
                # 如果两个hit之间有足够的间隙（>50bp），在间隙中点设置断点
                if next_start > current_end + 50:
                    breakpoint = current_end + (next_start - current_end) // 2
                    # 确保断点在有效范围内
                    if 50 <= breakpoint <= sequence_length - 50:
                        breakpoints.append(breakpoint)
            
            return sorted(breakpoints)
            
        except Exception as e:
            logger.debug(f"Error extracting breakpoints from hits: {e}")
            return []
    
    def _analyze_chimera_with_repeatmasker(self, chimera: Dict) -> Dict[str, Any]:
        """使用RepeatMasker分析嵌合体序列的可拆分性"""
        try:
            # 创建临时文件
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as temp_fa:
                temp_fa.write(f">{chimera['id']}\n{chimera['sequence']}\n")
                temp_fa_path = temp_fa.name
            
            # 运行RepeatMasker进行自比对分析
            temp_dir = Path(temp_fa_path).parent
            rm_cmd = [
                "RepeatMasker",
                "-no_is",  # 不检查简单重复
                "-nolow",  # 不屏蔽低复杂度
                "-s",      # 慢速但准确
                "-pa", str(self.threads_per_worker),  # 指定RepeatMasker使用的线程数
                "-dir", str(temp_dir),
                "-lib", temp_fa_path,  # 使用自身作为库进行比对
                temp_fa_path
            ]
            
            result = self.runner.run_command(
                rm_cmd,
                timeout=120,  # 2分钟超时
                description=f"RepeatMasker analysis for chimera {chimera['id']}"
            )
            
            # 解析RepeatMasker输出寻找拆分点
            rm_out_path = temp_fa_path + ".out"
            breakpoints = []
            
            if os.path.exists(rm_out_path):
                breakpoints = self._parse_repeatmasker_for_breakpoints(rm_out_path)
            
            # 清理临时文件
            for temp_file in [temp_fa_path, rm_out_path, temp_fa_path + ".masked", temp_fa_path + ".cat"]:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
            
            # 决定是否拆分
            if len(breakpoints) >= 2:  # 至少需要2个拆分点才能产生有意义的片段
                return {
                    'is_chimeric': True,
                    'is_splittable': True,
                    'breakpoints': breakpoints,
                    'reason': f'found_{len(breakpoints)}_breakpoints'
                }
            else:
                return {
                    'is_chimeric': False,
                    'is_splittable': False,
                    'breakpoints': [],
                    'reason': 'insufficient_breakpoints'
                }
                
        except Exception as e:
            logger.debug(f"RepeatMasker analysis failed for {chimera['id']}: {e}")
            return {
                'is_chimeric': False,
                'is_splittable': False,
                'breakpoints': [],
                'reason': f'analysis_failed: {str(e)}'
            }
    
    def _parse_repeatmasker_for_breakpoints(self, rm_out_path: str) -> List[int]:
        """解析RepeatMasker输出文件寻找潜在的拆分点"""
        breakpoints = []
        
        try:
            with open(rm_out_path, 'r') as f:
                lines = f.readlines()
            
            # 跳过头部信息
            data_start = 0
            for i, line in enumerate(lines):
                if line.strip().startswith('SW') or line.strip().startswith('score'):
                    data_start = i + 1
                    break
            
            hits = []
            for line in lines[data_start:]:
                if line.strip() and not line.startswith('There were no'):
                    parts = line.split()
                    if len(parts) >= 7:
                        try:
                            start = int(parts[5])
                            end = int(parts[6])
                            hits.append((start, end))
                        except (ValueError, IndexError):
                            continue
            
            # 寻找重复区域之间的间隙作为潜在拆分点
            if len(hits) >= 2:
                hits.sort()  # 按起始位置排序
                
                for i in range(len(hits) - 1):
                    current_end = hits[i][1]
                    next_start = hits[i + 1][0]
                    
                    # 如果两个重复区域之间有足够的间隙（>50bp），则在间隙中点拆分
                    gap_size = next_start - current_end
                    if gap_size > 50:
                        breakpoint = current_end + gap_size // 2
                        breakpoints.append(breakpoint)
        
        except Exception as e:
            logger.debug(f"Error parsing RepeatMasker output {rm_out_path}: {e}")
        
        return sorted(breakpoints)
    
    def _split_sequence_by_breakpoints(self, chimera: Dict, breakpoints: List[int]) -> List[Dict]:
        """根据拆分点将序列分割成片段"""
        sequence = chimera['sequence']
        chimera_id = chimera['id']
        
        # 确保拆分点在合理范围内
        valid_breakpoints = [bp for bp in breakpoints if 50 <= bp <= len(sequence) - 50]
        
        if not valid_breakpoints:
            return self._handle_unsplit_chimera(chimera, "no_valid_breakpoints")
        
        # 添加序列起始和结束位置
        positions = [0] + sorted(valid_breakpoints) + [len(sequence)]
        
        segments = []
        for i in range(len(positions) - 1):
            start_pos = positions[i]
            end_pos = positions[i + 1]
            
            # 确保片段有合理的长度
            if end_pos - start_pos < 50:  # 最小片段长度50bp
                continue
            
            segment_seq = sequence[start_pos:end_pos]
            segment_id = f"{chimera_id}_seg{i+1}"
            
            # 创建片段记录
            segment = {
                'id': segment_id,
                'sequence': segment_seq,
                'split_from': chimera_id,
                'segment_index': i + 1,
                'segment_start': start_pos,
                'segment_end': end_pos,
                'original_length': len(sequence),
                'quality_class': chimera.get('quality_class', 'B'),  # 继承或默认为B类
                'note': f'chimera_segment_from_{chimera_id}',
                'rm_hits': [],  # 重置RepeatMasker结果，需要重新分析
                'phase1_scores': {},  # 重置Phase1分数，将在Phase3重新计算
                '_analysis_source': chimera.get('_analysis_source', 'unknown')  # 继承分析来源
            }
            
            segments.append(segment)
        
        # 确保至少产生了有效的片段
        if not segments:
            return self._handle_unsplit_chimera(chimera, "no_valid_segments_generated")
        
        return segments
    
    def _handle_unsplit_chimera(self, chimera: Dict, reason: str) -> List[Dict]:
        """处理无法拆分的嵌合体序列 - 保留为完整序列"""
        unsplit = dict(chimera)  # 复制原序列
        unsplit['source'] = 'chimera_split'  # 标记为嵌合体处理结果
        unsplit['note'] = f'unsplit_chimera: {reason}'
        unsplit['split_from'] = None
        unsplit['segment_index'] = None
        
        return [unsplit]