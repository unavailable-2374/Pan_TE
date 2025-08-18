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
    
    def _extract_sequences_for_selected_hits(self, selected_hits: List[Dict], 
                                            seed_seq: Dict, config) -> List[Dict]:
        """
        优化版：使用外部工具批量提取序列
        替代原来的逐个BioPython提取
        """
        if not selected_hits:
            logger.warning("No hits provided for extraction")
            return []
        
        if not config or not hasattr(config, 'genome_file'):
            logger.error("Config or genome_file is missing!")
            return []
        
        logger.info(f"Starting optimized batch extraction for {len(selected_hits)} selected hits")
        
        # 准备批量提取的区域列表
        regions = []
        valid_indices = []
        
        for i, hit in enumerate(selected_hits):
            # 检查数据完整性
            required_fields = ['chrom', 'start', 'end']
            if all(field in hit and hit[field] is not None for field in required_fields):
                # 添加flanking区域
                flanking = 20
                extract_start = max(0, hit['start'] - flanking)
                extract_end = hit['end'] + flanking
                
                regions.append({
                    'chrom': hit['chrom'],
                    'start': extract_start,
                    'end': extract_end,
                    'strand': hit.get('strand', '+'),
                    'original_hit': hit,
                    'hit_index': i
                })
                valid_indices.append(i)
            else:
                logger.warning(f"Hit {i} missing required fields")
        
        if not regions:
            logger.warning("No valid regions to extract")
            return []
        
        # 批量提取序列
        logger.info(f"Batch extracting {len(regions)} sequences")
        
        if self.extractor:
            # 使用优化的批量提取
            extracted = self.extractor.extract_batch(regions)
        else:
            # 降级到传统批量提取
            extracted = extract_sequences_batch(config.genome_file, regions)
        
        logger.info(f"Successfully extracted {len(extracted)} sequences")
        
        # 处理提取结果
        copies = []
        for result in extracted:
            if not result.get('sequence'):
                continue
            
            hit = result['original_hit']
            sequence = result['sequence']
            
            # 计算实际的flanking长度
            actual_start = result['start']
            actual_end = result['end']
            left_flanking = hit['start'] - actual_start
            right_flanking = actual_end - hit['end']
            
            # 如果获取了flanking序列，提取核心序列和扩展序列
            if left_flanking > 0 or right_flanking > 0:
                core_start = left_flanking
                core_end = len(sequence) - right_flanking
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
                'chrom': hit['chrom'],
                'start': hit['start'],
                'end': hit['end'],
                'strand': hit.get('strand', '+'),
                'length': len(core_sequence),
                'identity': hit.get('identity', 100),
                'divergence': hit.get('divergence', 0),
                'has_tsd': tsd_info is not None,
                'tsd': tsd_info,
                'quality_score': self._calculate_copy_quality(hit, tsd_info is not None)
            }
            
            copies.append(copy_record)
        
        logger.info(f"Processed {len(copies)} valid copies from {len(regions)} extractions")
        
        return copies
    
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