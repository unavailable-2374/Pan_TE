"""
二轮采样策略：先基于RM信息粗采样，再基于序列信息精采样
兼顾性能和精度的最优方案
"""

import logging
import numpy as np
import random
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
import math

logger = logging.getLogger(__name__)


class TwoStageSamplingStrategy:
    """
    二轮采样策略实现
    第一轮：基于RepeatMasker hits的快速采样
    第二轮：基于实际序列特征的精确采样
    """
    
    def __init__(self, config):
        self.config = config
        
        # 采样参数配置
        self.stage1_params = {
            'small': {'max_copies': 500, 'reduction_ratio': 0.1},    # <5K拷贝
            'medium': {'max_copies': 800, 'reduction_ratio': 0.08},  # 5K-15K拷贝
            'large': {'max_copies': 1000, 'reduction_ratio': 0.06},  # 15K-30K拷贝
            'huge': {'max_copies': 1200, 'reduction_ratio': 0.04},   # >30K拷贝
        }
        
        self.stage2_params = {
            'target_copies': 200,  # 最终目标拷贝数
            'quality_threshold': 0.6,  # 质量门槛
            'diversity_weight': 0.7,   # 多样性权重
        }
    
    def two_stage_sampling(self, rm_hits: List[Dict], seed_seq: Dict, config) -> List[Dict]:
        """
        二轮采样主函数
        """
        if not rm_hits:
            return []
        
        n_total = len(rm_hits)
        logger.info(f"Starting two-stage sampling for {n_total} RM hits of {seed_seq['id']}")
        
        # 确定第一轮采样策略
        stage1_strategy = self._determine_stage1_strategy(n_total)
        stage1_config = self.stage1_params[stage1_strategy]
        
        # 第一轮：基于RM信息的粗采样
        stage1_selected = self._stage1_rm_based_sampling(
            rm_hits, 
            stage1_config['max_copies'],
            seed_seq
        )
        
        logger.info(f"Stage 1: Selected {len(stage1_selected)} hits from {n_total} (strategy: {stage1_strategy})")
        
        # 检查选中的hits是否为空
        if not stage1_selected:
            logger.warning("Stage 1: No hits were selected!")
            return []
        
        # 提取第一轮选中的序列
        stage1_copies = self._extract_sequences_for_selected_hits(
            stage1_selected, 
            seed_seq, 
            config
        )
        
        logger.info(f"Stage 1: Successfully extracted {len(stage1_copies)} sequences")
        
        # 第二轮：基于序列信息的精采样
        if len(stage1_copies) <= self.stage2_params['target_copies']:
            # 如果第一轮结果已经足够少，直接返回
            logger.info("Stage 1 result already meets target, skipping stage 2")
            return stage1_copies
        
        stage2_selected = self._stage2_sequence_based_sampling(
            stage1_copies,
            self.stage2_params['target_copies'],
            seed_seq
        )
        
        logger.info(f"Stage 2: Final selection of {len(stage2_selected)} sequences")
        
        return stage2_selected
    
    def _determine_stage1_strategy(self, n_copies: int) -> str:
        """确定第一轮采样策略"""
        if n_copies < 5000:
            return 'small'
        elif n_copies < 15000:
            return 'medium'
        elif n_copies < 30000:
            return 'large'
        else:
            return 'huge'
    
    def _stage1_rm_based_sampling(self, rm_hits: List[Dict], target_count: int, seed_seq: Dict) -> List[Dict]:
        """
        第一轮：基于RepeatMasker信息的分层采样
        重点：快速缩减规模，保持多样性覆盖
        """
        if len(rm_hits) <= target_count:
            return rm_hits
        
        # 1. 多维度分层
        stratified_hits = self._stratify_rm_hits(rm_hits)
        
        # 2. 按层分配配额
        layer_quotas = self._calculate_layer_quotas(stratified_hits, target_count)
        
        # 3. 从每层采样
        selected_hits = []
        for layer_key, layer_hits in stratified_hits.items():
            quota = layer_quotas.get(layer_key, 0)
            if quota > 0:
                layer_selected = self._sample_from_layer(layer_hits, quota, method='quality_spread')
                selected_hits.extend(layer_selected)
        
        # 4. 确保目标数量
        if len(selected_hits) < target_count:
            # 补充采样
            all_selected_ids = set(hit.get('hit_id', f"{hit['chrom']}_{hit['start']}_{hit['end']}") for hit in selected_hits)
            remaining_hits = [hit for hit in rm_hits 
                            if hit.get('hit_id', f"{hit['chrom']}_{hit['start']}_{hit['end']}") not in all_selected_ids]
            
            additional_needed = target_count - len(selected_hits)
            if remaining_hits and additional_needed > 0:
                additional = random.sample(remaining_hits, min(additional_needed, len(remaining_hits)))
                selected_hits.extend(additional)
        
        return selected_hits[:target_count]
    
    def _stratify_rm_hits(self, rm_hits: List[Dict]) -> Dict[Tuple, List[Dict]]:
        """对RM hits进行多维度分层"""
        stratified = defaultdict(list)
        
        for hit in rm_hits:
            # 维度1：Identity层级 (4层)
            identity = hit.get('identity', 0)
            if identity >= 90:
                identity_tier = 3  # 极高相似度
            elif identity >= 80:
                identity_tier = 2  # 高相似度
            elif identity >= 70:
                identity_tier = 1  # 中等相似度
            else:
                identity_tier = 0  # 低相似度
            
            # 维度2：Score层级 (3层)
            score = hit.get('score', 0)
            score_percentiles = self._calculate_score_percentiles(rm_hits)
            if score >= score_percentiles[75]:
                score_tier = 2  # 高分
            elif score >= score_percentiles[25]:
                score_tier = 1  # 中分
            else:
                score_tier = 0  # 低分
            
            # 维度3：长度层级 (基于start-end)
            hit_length = hit.get('end', 0) - hit.get('start', 0)
            length_tier = self._classify_hit_length(hit_length)
            
            # 维度4：染色体分布 (简化为数值hash)
            chrom_tier = hash(hit.get('chrom', 'chr1')) % 4  # 4个染色体组
            
            # 组合键
            layer_key = (identity_tier, score_tier, length_tier, chrom_tier)
            stratified[layer_key].append(hit)
        
        return stratified
    
    def _calculate_score_percentiles(self, rm_hits: List[Dict]) -> Dict[int, float]:
        """计算score的百分位数"""
        scores = [hit.get('score', 0) for hit in rm_hits]
        if not scores:
            return {25: 0, 50: 0, 75: 0}
        
        return {
            25: np.percentile(scores, 25),
            50: np.percentile(scores, 50),
            75: np.percentile(scores, 75)
        }
    
    def _classify_hit_length(self, hit_length: int) -> int:
        """基于hit长度分类"""
        if hit_length >= 10000:
            return 3  # 超长
        elif hit_length >= 5000:
            return 2  # 长
        elif hit_length >= 1000:
            return 1  # 中等
        else:
            return 0  # 短
    
    def _calculate_layer_quotas(self, stratified_hits: Dict, total_quota: int) -> Dict:
        """计算每层的采样配额"""
        layer_quotas = {}
        
        # 基础配额：按层大小比例分配
        total_hits = sum(len(hits) for hits in stratified_hits.values())
        if total_hits == 0:
            return layer_quotas
        
        base_quotas = {}
        for layer_key, layer_hits in stratified_hits.items():
            proportion = len(layer_hits) / total_hits
            base_quotas[layer_key] = int(total_quota * proportion)
        
        # 确保重要层有最小配额
        min_quota_per_layer = max(1, total_quota // (len(stratified_hits) * 4))
        
        allocated = 0
        for layer_key, base_quota in base_quotas.items():
            final_quota = max(min_quota_per_layer, base_quota)
            layer_quotas[layer_key] = final_quota
            allocated += final_quota
        
        # 调整超额部分
        if allocated > total_quota:
            # 按比例缩减
            scale_factor = total_quota / allocated
            for layer_key in layer_quotas:
                layer_quotas[layer_key] = max(1, int(layer_quotas[layer_key] * scale_factor))
        
        return layer_quotas
    
    def _sample_from_layer(self, layer_hits: List[Dict], quota: int, method: str = 'quality_spread') -> List[Dict]:
        """从层内采样"""
        if len(layer_hits) <= quota:
            return layer_hits
        
        if method == 'quality_spread':
            # 质量分布采样：确保覆盖质量梯度
            # 50%高质量 + 30%中等质量 + 20%多样性
            sorted_hits = sorted(layer_hits, 
                               key=lambda x: x.get('identity', 0) * 0.7 + x.get('score', 0) / 1000 * 0.3, 
                               reverse=True)
            
            high_quota = max(1, quota // 2)
            medium_quota = max(1, quota * 3 // 10)
            diversity_quota = quota - high_quota - medium_quota
            
            selected = []
            
            # 高质量选择
            selected.extend(sorted_hits[:high_quota])
            
            # 中等质量选择
            mid_start = len(sorted_hits) // 4
            mid_end = len(sorted_hits) * 3 // 4
            if mid_end > mid_start:
                medium_pool = sorted_hits[mid_start:mid_end]
                selected.extend(medium_pool[:medium_quota])
            
            # 多样性选择（等间隔采样）
            if diversity_quota > 0:
                used_indices = set(range(high_quota)) | set(range(mid_start, mid_start + medium_quota))
                available_hits = [hit for i, hit in enumerate(sorted_hits) if i not in used_indices]
                
                if available_hits:
                    if len(available_hits) <= diversity_quota:
                        selected.extend(available_hits)
                    else:
                        step = len(available_hits) // diversity_quota
                        for i in range(0, len(available_hits), step):
                            if len(selected) < quota:
                                selected.append(available_hits[i])
            
            return selected[:quota]
        
        else:
            # 简单随机采样
            return random.sample(layer_hits, quota)
    
    def _extract_sequences_for_selected_hits(self, selected_hits: List[Dict], seed_seq: Dict, config) -> List[Dict]:
        """为选中的hits提取序列"""
        # 早期检查
        if not selected_hits:
            logger.warning("_extract_sequences_for_selected_hits: No hits provided for extraction")
            return []
        
        if not config or not hasattr(config, 'genome_file'):
            logger.error("Config or genome_file is missing!")
            return []
            
        from utils.sequence_utils import extract_sequence_from_genome, detect_tsd_from_sequence
        
        logger.info(f"Starting sequence extraction for {len(selected_hits)} selected hits")
        logger.debug(f"Genome file: {config.genome_file}")
        logger.debug(f"Seed sequence length: {len(seed_seq.get('sequence', ''))}")
        
        copies = []
        successful_extractions = 0
        failed_extractions = 0
        
        for i, hit in enumerate(selected_hits):
            try:
                # 记录hit信息
                logger.debug(f"Processing hit {i+1}/{len(selected_hits)}: {hit.get('chrom', 'unknown')}:{hit.get('start', 'unknown')}-{hit.get('end', 'unknown')}")
                
                # 检查hit数据完整性
                required_fields = ['chrom', 'start', 'end']
                missing_fields = [field for field in required_fields if field not in hit or hit[field] is None]
                if missing_fields:
                    logger.warning(f"Hit {i} missing required fields: {missing_fields}")
                    failed_extractions += 1
                    continue
                
                # 提取序列
                flanking = 20
                extract_start = max(0, hit['start'] - flanking)
                extract_end = hit['end'] + flanking
                
                # 记录原始坐标用于计算实际flanking
                original_start = hit['start']
                original_end = hit['end']
                
                logger.debug(f"Extracting {hit['chrom']}:{extract_start}-{extract_end} (strand: {hit.get('strand', '+')})")
                
                extended_seq = extract_sequence_from_genome(
                    config.genome_file,
                    hit['chrom'],
                    extract_start,
                    extract_end,
                    hit.get('strand', '+')
                )
                
                if not extended_seq:
                    logger.warning(f"Failed to extract sequence for hit {i}: {hit['chrom']}:{extract_start}-{extract_end}")
                    failed_extractions += 1
                    continue
                
                logger.debug(f"Successfully extracted sequence of length {len(extended_seq)}")
                
                # 计算实际获取的flanking长度
                actual_left_flanking = original_start - extract_start
                actual_right_flanking = extract_end - original_end
                
                # 提取核心序列，考虑实际的flanking长度
                if actual_left_flanking > 0 and actual_right_flanking > 0:
                    if len(extended_seq) > actual_left_flanking + actual_right_flanking:
                        core_seq = extended_seq[actual_left_flanking:len(extended_seq)-actual_right_flanking]
                        logger.debug(f"Core sequence length after flanking removal: {len(core_seq)}")
                    else:
                        core_seq = extended_seq
                        logger.debug(f"Extended sequence too short for flanking removal, using full sequence")
                else:
                    core_seq = extended_seq
                    logger.debug(f"No flanking regions to remove (boundary hit)")
                
                # 检查核心序列质量
                if len(core_seq) == 0:
                    logger.warning(f"Empty core sequence for hit {i}")
                    failed_extractions += 1
                    continue
                
                # 检测TSD
                try:
                    tsd = detect_tsd_from_sequence(extended_seq, flanking_size=flanking)
                    logger.debug(f"TSD detection result: {tsd is not None}")
                except Exception as tsd_error:
                    logger.debug(f"TSD detection failed for hit {i}: {tsd_error}")
                    tsd = None
                
                # 创建copy对象
                copy = {
                    'id': f"copy_{i}",
                    'sequence': core_seq,
                    'extended_sequence': extended_seq,
                    'chr': hit['chrom'],
                    'start': hit['start'],
                    'end': hit['end'],
                    'strand': hit.get('strand', '+'),
                    'identity': hit.get('identity', 0),
                    'divergence': 100 - hit.get('identity', 0),
                    'score': hit.get('score', 0),
                    'length': len(core_seq),
                    'tsd': tsd,
                    'has_tsd': tsd is not None,
                    'rm_source': True  # 标记来源于第一轮RM采样
                }
                
                copies.append(copy)
                successful_extractions += 1
                
                # 定期报告进度
                if (i + 1) % 100 == 0:
                    logger.info(f"Processed {i+1}/{len(selected_hits)} hits: {successful_extractions} successful, {failed_extractions} failed")
                
            except Exception as e:
                logger.error(f"Error extracting sequence for hit {i}: {e}")
                failed_extractions += 1
                continue
        
        # 最终统计和报告
        logger.info(f"Sequence extraction completed: {successful_extractions} successful, {failed_extractions} failed")
        logger.info(f"Total extracted sequences: {len(copies)}")
        
        if len(copies) == 0:
            logger.warning("No sequences were successfully extracted! Possible issues:")
            logger.warning("  1. Chromosome name mismatch between RepeatMasker output and genome file")
            logger.warning("  2. Coordinate errors in RepeatMasker hits")
            logger.warning("  3. Genome file format or accessibility issues")
            
            # 提供一些调试建议
            if selected_hits:
                sample_hit = selected_hits[0]
                logger.warning(f"Sample hit data: {sample_hit}")
        
        return copies
    
    def _stage2_sequence_based_sampling(self, copies: List[Dict], target_count: int, seed_seq: Dict) -> List[Dict]:
        """
        第二轮：基于序列特征的精确采样
        重点：序列质量评估，精确多样性控制
        """
        if len(copies) <= target_count:
            return copies
        
        logger.info(f"Stage 2: Refining {len(copies)} sequences to {target_count}")
        
        # 1. 序列质量评估
        enriched_copies = self._evaluate_sequence_quality(copies, seed_seq)
        
        # 2. 高质量序列优先保留 (40%配额)
        quality_quota = max(5, target_count * 2 // 5)
        quality_selected = self._select_by_quality(enriched_copies, quality_quota)
        
        # 3. 多样性保护采样 (40%配额)
        diversity_quota = max(5, target_count * 2 // 5)
        remaining_copies = [c for c in enriched_copies if c not in quality_selected]
        diversity_selected = self._select_by_diversity(remaining_copies, diversity_quota)
        
        # 4. 结构特征保护 (20%配额)
        structure_quota = target_count - len(quality_selected) - len(diversity_selected)
        remaining_copies = [c for c in enriched_copies 
                          if c not in quality_selected and c not in diversity_selected]
        structure_selected = self._select_by_structure(remaining_copies, structure_quota)
        
        # 5. 合并结果
        final_selected = quality_selected + diversity_selected + structure_selected
        
        return final_selected[:target_count]
    
    def _evaluate_sequence_quality(self, copies: List[Dict], seed_seq: Dict) -> List[Dict]:
        """评估序列质量"""
        logger.info(f"Evaluating quality for {len(copies)} sequences")
        
        seed_length = len(seed_seq.get('sequence', ''))
        logger.debug(f"Seed sequence length for reference: {seed_length}")
        
        quality_stats = {
            'high_quality': 0,
            'medium_quality': 0, 
            'low_quality': 0,
            'total_n_content': 0,
            'total_gc_content': 0,
            'has_tsd_count': 0
        }
        
        for i, copy in enumerate(copies):
            sequence = copy.get('sequence', '')
            
            if len(sequence) == 0:
                logger.warning(f"Copy {i} has empty sequence")
                continue
            
            # 1. N含量评估
            n_content = sequence.count('N') / len(sequence)
            copy['n_content'] = n_content
            quality_stats['total_n_content'] += n_content
            
            # 2. GC含量评估
            gc_count = sequence.count('G') + sequence.count('C')
            gc_content = gc_count / len(sequence)
            copy['gc_content'] = gc_content
            quality_stats['total_gc_content'] += gc_content
            
            # 3. 低复杂度区域评估（简化版）
            complexity_score = self._calculate_sequence_complexity(sequence)
            copy['complexity_score'] = complexity_score
            
            # 4. 长度完整性
            length_completeness = len(sequence) / seed_length if seed_length > 0 else 0
            copy['length_completeness'] = min(length_completeness, 1.0)
            
            # 5. TSD统计
            if copy.get('has_tsd', False):
                quality_stats['has_tsd_count'] += 1
            
            # 6. 综合质量分数
            quality_components = [
                (1 - n_content) * 0.3,  # N含量越低越好
                copy.get('identity', 0) / 100 * 0.4,  # 相似度
                copy['length_completeness'] * 0.2,  # 长度完整性
                (1 if copy.get('has_tsd', False) else 0.5) * 0.1  # TSD存在性
            ]
            copy['sequence_quality'] = sum(quality_components)
            
            # 质量分级
            if copy['sequence_quality'] >= 0.8:
                quality_stats['high_quality'] += 1
            elif copy['sequence_quality'] >= 0.5:
                quality_stats['medium_quality'] += 1
            else:
                quality_stats['low_quality'] += 1
            
            # 详细日志每50个序列
            if (i + 1) % 50 == 0 or i < 5:
                logger.debug(f"Copy {i}: len={len(sequence)}, N={n_content:.3f}, GC={gc_content:.3f}, "
                           f"complexity={complexity_score:.3f}, quality={copy['sequence_quality']:.3f}")
        
        # 报告质量统计
        total_copies = len(copies)
        if total_copies > 0:
            avg_n_content = quality_stats['total_n_content'] / total_copies
            avg_gc_content = quality_stats['total_gc_content'] / total_copies
            tsd_percentage = quality_stats['has_tsd_count'] / total_copies * 100
            
            logger.info(f"Quality assessment summary:")
            logger.info(f"  High quality (≥0.8): {quality_stats['high_quality']} ({quality_stats['high_quality']/total_copies*100:.1f}%)")
            logger.info(f"  Medium quality (0.5-0.8): {quality_stats['medium_quality']} ({quality_stats['medium_quality']/total_copies*100:.1f}%)")
            logger.info(f"  Low quality (<0.5): {quality_stats['low_quality']} ({quality_stats['low_quality']/total_copies*100:.1f}%)")
            logger.info(f"  Average N content: {avg_n_content:.3f}")
            logger.info(f"  Average GC content: {avg_gc_content:.3f}")
            logger.info(f"  Sequences with TSD: {quality_stats['has_tsd_count']} ({tsd_percentage:.1f}%)")
        
        return copies
    
    def _calculate_sequence_complexity(self, sequence: str) -> float:
        """计算序列复杂度（简化版）"""
        if len(sequence) == 0:
            return 0
        
        # 计算4-mer频率分布的熵
        from collections import Counter
        kmers = [sequence[i:i+4] for i in range(len(sequence)-3)]
        if not kmers:
            return 0
        
        kmer_counts = Counter(kmers)
        total_kmers = len(kmers)
        
        entropy = 0
        for count in kmer_counts.values():
            freq = count / total_kmers
            if freq > 0:
                entropy -= freq * math.log2(freq)
        
        # 归一化到0-1范围
        max_entropy = math.log2(min(256, len(kmers)))  # 理论最大熵
        return entropy / max_entropy if max_entropy > 0 else 0
    
    def _select_by_quality(self, copies: List[Dict], quota: int) -> List[Dict]:
        """基于质量选择"""
        if len(copies) <= quota:
            return copies
        
        return sorted(copies, key=lambda x: x.get('sequence_quality', 0), reverse=True)[:quota]
    
    def _select_by_diversity(self, copies: List[Dict], quota: int) -> List[Dict]:
        """基于多样性选择（简化版聚类）"""
        if len(copies) <= quota:
            return copies
        
        # 基于identity和length进行简单分组
        groups = defaultdict(list)
        for copy in copies:
            identity_group = int(copy.get('identity', 0) // 10)  # 10%间隔分组
            length_group = int(math.log2(max(copy.get('length', 1), 1)))  # 对数长度分组
            group_key = (identity_group, length_group)
            groups[group_key].append(copy)
        
        # 从每组选择最佳代表
        selected = []
        per_group_quota = max(1, quota // len(groups))
        
        for group_copies in groups.values():
            group_quota = min(per_group_quota, len(group_copies))
            group_best = sorted(group_copies, 
                              key=lambda x: x.get('sequence_quality', 0), 
                              reverse=True)[:group_quota]
            selected.extend(group_best)
        
        return selected[:quota]
    
    def _select_by_structure(self, copies: List[Dict], quota: int) -> List[Dict]:
        """基于结构特征选择"""
        if len(copies) <= quota or quota <= 0:
            return copies[:quota] if quota > 0 else []
        
        # 优先级：TSD > 高GC > 长序列 > 低N含量
        scored_copies = []
        for copy in copies:
            structure_score = 0
            
            # TSD加分
            if copy.get('has_tsd', False):
                structure_score += 100
            
            # GC含量合理性加分
            gc_content = copy.get('gc_content', 0.5)
            if 0.3 <= gc_content <= 0.7:  # 合理GC范围
                structure_score += 50
            
            # 长度完整性加分
            structure_score += copy.get('length_completeness', 0) * 30
            
            # 低N含量加分
            structure_score += (1 - copy.get('n_content', 1)) * 20
            
            scored_copies.append((structure_score, copy))
        
        # 按结构分数排序选择
        scored_copies.sort(key=lambda x: x[0], reverse=True)
        return [copy for _, copy in scored_copies[:quota]]


def apply_two_stage_sampling(rm_hits: List[Dict], seed_seq: Dict, config) -> List[Dict]:
    """
    二轮采样的接口函数
    """
    if not rm_hits:
        return []
    
    # 小规模数据直接提取所有序列
    if len(rm_hits) <= 200:
        logger.info(f"Small dataset ({len(rm_hits)} hits), extracting all sequences directly")
        strategy = TwoStageSamplingStrategy(config)
        return strategy._extract_sequences_for_selected_hits(rm_hits, seed_seq, config)
    
    # 大规模数据使用二轮采样
    strategy = TwoStageSamplingStrategy(config)
    return strategy.two_stage_sampling(rm_hits, seed_seq, config)