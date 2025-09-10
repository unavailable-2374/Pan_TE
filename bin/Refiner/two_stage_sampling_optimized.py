"""
优化版二轮采样策略：使用外部工具批量提取序列
"""

import logging
import numpy as np
import random
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
import math

# 导入原始的二轮采样类
from two_stage_sampling import TwoStageSamplingStrategy

# 导入新的快速提取器
from utils.fast_sequence_extractor import FastSequenceExtractor, extract_sequences_batch

logger = logging.getLogger(__name__)


class OptimizedTwoStageSamplingStrategy(TwoStageSamplingStrategy):
    """
    优化的二轮采样策略，使用外部工具批量提取序列
    """
    
    def __init__(self, config):
        super().__init__(config)
        # 初始化快速提取器
        self.extractor = None
        if hasattr(config, 'genome_file'):
            self.extractor = FastSequenceExtractor(config.genome_file)
            # 注意：基因组文件已规范化为lcl|格式，无需ID映射
    
    def _extract_sequences_for_selected_hits(self, selected_hits: List[Dict], 
                                            seed_seq: Dict, config) -> List[Dict]:
        """
        优化版：使用安全的序列提取器
        同时处理ID映射和坐标不匹配问题
        """
        if not selected_hits:
            logger.warning("No hits provided for extraction")
            return []
        
        if not config or not hasattr(config, 'genome_file'):
            logger.error("Config or genome_file is missing!")
            return []
        
        logger.info(f"Starting safe batch extraction for {len(selected_hits)} selected hits")
        
        # 使用安全序列提取器
        from utils.safe_sequence_extractor import create_safe_extractor
        extractor = create_safe_extractor(config.genome_file)
        
        # 准备批量提取的区域列表
        regions = []
        for i, hit in enumerate(selected_hits):
            # 检查数据完整性
            required_fields = ['chrom', 'start', 'end']
            if all(field in hit and hit[field] is not None for field in required_fields):
                regions.append({
                    'chrom': hit['chrom'],
                    'start': hit['start'], 
                    'end': hit['end'],
                    'strand': hit.get('strand', '+'),
                    'original_hit': hit,
                    'hit_index': i,
                    'identity': hit.get('identity', 100),
                    'divergence': hit.get('divergence', 0)
                })
            else:
                logger.warning(f"Hit {i} missing required fields: {hit}")
        
        if not regions:
            logger.warning("No valid regions to extract")
            return []
        
        # 批量安全提取序列（包含20bp flanking区域）
        logger.info(f"Safe batch extracting {len(regions)} sequences")
        extracted_results = extractor.extract_batch_safe(
            regions, 
            flanking=20, 
            fallback_strategy='truncate'  # 截断超出范围的坐标
        )
        
        logger.info(f"Successfully extracted {len(extracted_results)} sequences")
        
        # 处理提取结果
        copies = []
        for result in extracted_results:
            if not result.get('sequence'):
                continue
            
            hit = result['original_hit']
            sequence = result['sequence']
            
            # 获取实际提取的坐标
            actual_start = result['start']
            actual_end = result['end']
            original_start = result['original_start']
            original_end = result['original_end']
            flanking = result.get('flanking', 20)
            
            # 计算flanking区域
            left_flanking = max(0, original_start - actual_start)
            right_flanking = max(0, actual_end - original_end)
            
            # 提取核心序列
            if left_flanking > 0 or right_flanking > 0:
                core_start = left_flanking
                core_end = len(sequence) - right_flanking if right_flanking > 0 else len(sequence)
                core_sequence = sequence[core_start:core_end] if core_start < core_end else sequence
                extended_sequence = sequence
            else:
                core_sequence = sequence
                extended_sequence = None
            
            # 检测TSD
            tsd_info = None
            if left_flanking >= 5 and right_flanking >= 5:
                tsd_info = self._detect_tsd_from_flanking(
                    sequence[:left_flanking] if left_flanking > 0 else '',
                    sequence[-right_flanking:] if right_flanking > 0 else ''
                )
            
            # 创建拷贝记录
            copy_record = {
                'sequence': core_sequence,
                'extended_sequence': extended_sequence,
                'chrom': result['original_chrom'],  # 保持原始ID用于记录
                'mapped_chrom': result['chrom'],    # 实际映射的ID
                'start': original_start,
                'end': original_end,
                'actual_start': actual_start,
                'actual_end': actual_end,
                'strand': hit.get('strand', '+'),
                'length': len(core_sequence),
                'identity': hit.get('identity', 100),
                'divergence': hit.get('divergence', 0),
                'has_tsd': tsd_info is not None,
                'tsd': tsd_info,
                'quality_score': self._calculate_copy_quality(hit, tsd_info is not None),
                'coordinate_issues': result.get('coordinate_issues', []),
                'extraction_notes': self._format_extraction_notes(result)
            }
            
            copies.append(copy_record)
        
        logger.info(f"Processed {len(copies)} valid copies from {len(extracted_results)} extractions")
        
        return copies
    
    def _format_extraction_notes(self, result: Dict) -> str:
        """格式化序列提取的注释信息"""
        notes = []
        
        if result.get('coordinate_issues'):
            notes.extend(result['coordinate_issues'])
        
        if result.get('fallback_strategy'):
            notes.append(f"Used fallback strategy: {result['fallback_strategy']}")
        
        if result['original_chrom'] != result['chrom']:
            notes.append(f"ID mapped: {result['original_chrom']} -> {result['chrom']}")
        
        return '; '.join(notes) if notes else ''
    
    def _detect_tsd_from_flanking(self, left_flank: str, right_flank: str) -> Optional[Dict]:
        """检测目标位点重复（TSD）"""
        if not left_flank or not right_flank:
            return None
        
        # 检查常见的TSD长度 (4-20bp)
        for tsd_len in range(4, min(21, len(left_flank)+1, len(right_flank)+1)):
            left_tsd = left_flank[-tsd_len:]
            right_tsd = right_flank[:tsd_len]
            
            # 计算相似度
            matches = sum(1 for a, b in zip(left_tsd, right_tsd) if a == b)
            identity = matches / tsd_len
            
            if identity >= 0.8:  # 80%相似度阈值
                return {
                    'sequence': left_tsd,
                    'length': tsd_len,
                    'identity': identity
                }
        
        return None
    
    def _calculate_copy_quality(self, hit: Dict, has_tsd: bool) -> float:
        """计算拷贝质量分数"""
        quality_score = 0.0
        
        # 基于identity
        identity = hit.get('identity', 100)
        quality_score += (identity / 100) * 0.5
        
        # 基于长度覆盖
        coverage = hit.get('coverage', 100)
        quality_score += (coverage / 100) * 0.3
        
        # TSD奖励
        if has_tsd:
            quality_score += 0.2
        
        return min(quality_score, 1.0)


def apply_optimized_two_stage_sampling(rm_hits: List[Dict], seed_seq: Dict, config) -> List[Dict]:
    """
    优化的二轮采样入口函数
    """
    if not rm_hits:
        return []
    
    logger.info(f"Applying optimized two-stage sampling to {len(rm_hits)} RM hits")
    
    # 小规模数据直接批量提取
    if len(rm_hits) <= 200:
        logger.info(f"Small dataset ({len(rm_hits)} hits), using direct batch extraction")
        strategy = OptimizedTwoStageSamplingStrategy(config)
        return strategy._extract_sequences_for_selected_hits(rm_hits, seed_seq, config)
    
    # 大规模数据使用优化的二轮采样
    strategy = OptimizedTwoStageSamplingStrategy(config)
    return strategy.two_stage_sampling(rm_hits, seed_seq, config)