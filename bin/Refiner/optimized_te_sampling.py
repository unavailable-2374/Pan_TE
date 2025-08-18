"""
优化的TE感知采样策略
专门针对高拷贝数序列设计，避免昂贵的层次聚类
基于TE生物学特性的高效采样算法
"""

import logging
import numpy as np
import random
from typing import List, Dict, Tuple, Optional
from collections import defaultdict, Counter
import math

logger = logging.getLogger(__name__)


class OptimizedTESamplingStrategy:
    """
    针对高拷贝数TE的优化采样策略
    核心思想：基于TE生物学特性进行智能分层，避免全量聚类
    """
    
    def __init__(self, config):
        self.config = config
        self.base_limit = getattr(config, 'max_recruits_per_family', 30)
        
        # 不同拷贝数量级的处理策略
        self.copy_count_strategies = {
            'small': {'threshold': 100, 'max_sample': 80, 'strategy': 'full_analysis'},
            'medium': {'threshold': 500, 'max_sample': 100, 'strategy': 'light_clustering'}, 
            'large': {'threshold': 2000, 'max_sample': 120, 'strategy': 'grid_sampling'},
            'huge': {'threshold': 10000, 'max_sample': 150, 'strategy': 'stratified_sampling'},
            'ultra': {'threshold': float('inf'), 'max_sample': 200, 'strategy': 'fast_stratified'}
        }
        
        # TE类型特异参数
        self.te_profiles = {
            'LINE': {
                'length_bins': [0, 1000, 3000, 6000, float('inf')], 
                'identity_importance': 0.8,
                'length_importance': 0.9,
                'tsd_importance': 0.3
            },
            'SINE': {
                'length_bins': [0, 100, 300, 500, float('inf')],
                'identity_importance': 0.7, 
                'length_importance': 0.6,
                'tsd_importance': 0.8
            },
            'LTR': {
                'length_bins': [0, 2000, 5000, 10000, float('inf')],
                'identity_importance': 0.9,
                'length_importance': 0.8, 
                'tsd_importance': 0.9
            },
            'DNA': {
                'length_bins': [0, 500, 1500, 3000, float('inf')],
                'identity_importance': 0.8,
                'length_importance': 0.7,
                'tsd_importance': 0.7
            }
        }
    
    def select_representative_copies(self, copies: List[Dict], seed_seq: Dict) -> List[Dict]:
        """
        主采样函数 - 根据拷贝数量选择最优策略
        """
        if not copies:
            return []
        
        n_copies = len(copies)
        logger.info(f"Optimized sampling for {n_copies} copies of {seed_seq['id']}")
        
        # 预处理：增强元数据
        enriched_copies = self._enrich_metadata(copies, seed_seq)
        
        # 确定处理策略
        strategy_type = self._determine_strategy(n_copies)
        strategy_config = self.copy_count_strategies[strategy_type]
        
        logger.info(f"Using {strategy_config['strategy']} strategy for {n_copies} copies")
        
        # 根据策略选择采样方法
        if strategy_config['strategy'] == 'full_analysis':
            selected = self._full_analysis_sampling(enriched_copies, strategy_config['max_sample'])
        elif strategy_config['strategy'] == 'light_clustering':
            selected = self._light_clustering_sampling(enriched_copies, strategy_config['max_sample'])
        elif strategy_config['strategy'] == 'grid_sampling':
            selected = self._grid_based_sampling(enriched_copies, strategy_config['max_sample'], seed_seq)
        elif strategy_config['strategy'] == 'stratified_sampling':
            selected = self._stratified_sampling(enriched_copies, strategy_config['max_sample'], seed_seq)
        else:  # fast_stratified
            selected = self._fast_stratified_sampling(enriched_copies, strategy_config['max_sample'], seed_seq)
        
        logger.info(f"Selected {len(selected)} representative copies using {strategy_config['strategy']}")
        return selected
    
    def _determine_strategy(self, n_copies: int) -> str:
        """根据拷贝数确定处理策略"""
        for strategy_name, config in self.copy_count_strategies.items():
            if n_copies <= config['threshold']:
                return strategy_name
        return 'ultra'
    
    def _enrich_metadata(self, copies: List[Dict], seed_seq: Dict) -> List[Dict]:
        """增强拷贝元数据，但避免昂贵的序列比较"""
        te_type = self._fast_te_type_inference(copies, seed_seq)
        te_profile = self.te_profiles.get(te_type, self.te_profiles['DNA'])
        
        enriched = []
        for copy in copies:
            enriched_copy = copy.copy()
            
            # 快速质量评分（不涉及序列比对）
            enriched_copy['quality_score'] = self._fast_quality_score(copy, seed_seq, te_profile)
            enriched_copy['te_type'] = te_type
            enriched_copy['length_bin'] = self._assign_length_bin(copy, te_profile['length_bins'])
            enriched_copy['identity_tier'] = self._assign_identity_tier(copy)
            
            enriched.append(enriched_copy)
        
        return enriched
    
    def _fast_te_type_inference(self, copies: List[Dict], seed_seq: Dict) -> str:
        """快速TE类型推断，基于长度和拷贝特征"""
        seed_length = len(seed_seq.get('sequence', ''))
        
        # 采样少量拷贝进行分析（避免全量计算）
        sample_size = min(50, len(copies))
        sample_indices = np.random.choice(len(copies), sample_size, replace=False)
        sample_copies = [copies[i] for i in sample_indices]
        
        avg_length = np.mean([c.get('length', 0) for c in sample_copies])
        length_cv = np.std([c.get('length', 0) for c in sample_copies]) / avg_length if avg_length > 0 else 0
        tsd_ratio = sum(1 for c in sample_copies if c.get('has_tsd', False)) / len(sample_copies)
        
        # 基于统计特征分类
        if avg_length < 500:
            return 'SINE'
        elif avg_length > 5000:
            return 'LINE' 
        elif tsd_ratio > 0.5 and 1000 <= avg_length <= 5000:
            return 'LTR'
        elif length_cv < 0.3:  # 长度一致性高
            return 'DNA'
        else:
            return 'DNA'  # 默认
    
    def _fast_quality_score(self, copy: Dict, seed_seq: Dict, te_profile: Dict) -> float:
        """快速质量评分，避免复杂计算"""
        score_components = []
        
        # Identity权重
        identity = copy.get('identity', 0) / 100.0
        score_components.append(identity * te_profile['identity_importance'])
        
        # Length相对评分  
        seed_length = len(seed_seq.get('sequence', ''))
        copy_length = copy.get('length', 0)
        if seed_length > 0:
            length_ratio = copy_length / seed_length
            # 长度在合理范围内得高分
            if 0.8 <= length_ratio <= 1.2:
                length_score = 1.0
            elif 0.5 <= length_ratio <= 1.5:
                length_score = 0.8
            else:
                length_score = 0.5
        else:
            length_score = 0.5
        score_components.append(length_score * te_profile['length_importance'])
        
        # TSD加分
        tsd_score = 1.0 if copy.get('has_tsd', False) else 0.3
        score_components.append(tsd_score * te_profile['tsd_importance'])
        
        # 基本质量门槛
        basic_quality = 1.0 if copy_length >= 50 and identity >= 50 else 0.2
        score_components.append(basic_quality * 0.1)
        
        return sum(score_components) / (te_profile['identity_importance'] + 
                                      te_profile['length_importance'] + 
                                      te_profile['tsd_importance'] + 0.1)
    
    def _assign_length_bin(self, copy: Dict, length_bins: List[float]) -> int:
        """分配长度分箱"""
        length = copy.get('length', 0)
        for i, threshold in enumerate(length_bins[1:]):
            if length <= threshold:
                return i
        return len(length_bins) - 2
    
    def _assign_identity_tier(self, copy: Dict) -> int:
        """分配相似度层级"""
        identity = copy.get('identity', 0)
        if identity >= 95:
            return 4  # 极高相似度
        elif identity >= 85: 
            return 3  # 高相似度
        elif identity >= 70:
            return 2  # 中等相似度
        elif identity >= 50:
            return 1  # 低相似度
        else:
            return 0  # 极低相似度
    
    def _full_analysis_sampling(self, copies: List[Dict], max_sample: int) -> List[Dict]:
        """小规模数据：可以进行相对完整的分析"""
        if len(copies) <= max_sample:
            return copies
        
        # 简单的质量加多样性采样
        # 1. 高质量拷贝 (50%)
        high_quality_quota = max_sample // 2
        high_quality = sorted(copies, key=lambda x: x['quality_score'], reverse=True)[:high_quality_quota]
        
        # 2. 多样性拷贝 (50%)
        remaining = [c for c in copies if c not in high_quality]
        diversity_quota = max_sample - len(high_quality)
        diversity_selected = self._simple_diversity_sampling(remaining, diversity_quota)
        
        return high_quality + diversity_selected
    
    def _light_clustering_sampling(self, copies: List[Dict], max_sample: int) -> List[Dict]:
        """中等规模数据：轻量级聚类"""
        if len(copies) <= max_sample:
            return copies
        
        # 基于identity和length进行简单分组（避免复杂聚类）
        identity_groups = defaultdict(list)
        for copy in copies:
            identity_tier = copy['identity_tier'] 
            identity_groups[identity_tier].append(copy)
        
        # 从每个组按比例选择
        selected = []
        total_groups = len(identity_groups)
        per_group_quota = max_sample // total_groups
        remainder = max_sample % total_groups
        
        for i, (tier, group_copies) in enumerate(identity_groups.items()):
            group_quota = per_group_quota + (1 if i < remainder else 0)
            if len(group_copies) <= group_quota:
                selected.extend(group_copies)
            else:
                # 组内按质量选择
                group_selected = sorted(group_copies, key=lambda x: x['quality_score'], reverse=True)[:group_quota]
                selected.extend(group_selected)
        
        return selected[:max_sample]
    
    def _grid_based_sampling(self, copies: List[Dict], max_sample: int, seed_seq: Dict) -> List[Dict]:
        """大规模数据：网格采样"""
        if len(copies) <= max_sample:
            return copies
        
        # 构建二维网格：identity × length_bin
        grid = defaultdict(list)
        for copy in copies:
            grid_key = (copy['identity_tier'], copy['length_bin'])
            grid[grid_key].append(copy)
        
        # 计算每个网格的配额
        total_cells = len(grid)
        if total_cells == 0:
            return copies[:max_sample]
        
        base_quota = max_sample // total_cells
        remainder = max_sample % total_cells
        
        selected = []
        sorted_cells = sorted(grid.items(), key=lambda x: len(x[1]), reverse=True)  # 优先大细胞
        
        for i, ((identity_tier, length_bin), cell_copies) in enumerate(sorted_cells):
            cell_quota = base_quota + (1 if i < remainder else 0)
            
            if len(cell_copies) <= cell_quota:
                selected.extend(cell_copies)
            else:
                # 细胞内采样：质量优先 + 随机
                quality_quota = max(1, cell_quota // 2)
                random_quota = cell_quota - quality_quota
                
                # 质量最佳
                quality_selected = sorted(cell_copies, key=lambda x: x['quality_score'], reverse=True)[:quality_quota]
                selected.extend(quality_selected)
                
                # 随机选择剩余
                remaining_copies = [c for c in cell_copies if c not in quality_selected]
                if remaining_copies and random_quota > 0:
                    random_selected = random.sample(remaining_copies, min(random_quota, len(remaining_copies)))
                    selected.extend(random_selected)
        
        return selected[:max_sample]
    
    def _stratified_sampling(self, copies: List[Dict], max_sample: int, seed_seq: Dict) -> List[Dict]:
        """超大规模数据：分层抽样"""
        if len(copies) <= max_sample:
            return copies
        
        # 三层分层：quality_tier × identity_tier × length_bin
        strata = defaultdict(list)
        
        # 添加质量分层
        for copy in copies:
            quality_tier = self._get_quality_tier(copy['quality_score'])
            stratum_key = (quality_tier, copy['identity_tier'], copy['length_bin'])
            strata[stratum_key].append(copy)
        
        # 分层采样策略：优先重要层
        selected = []
        remaining_quota = max_sample
        
        # 第一轮：确保每个重要层都有代表
        important_strata = [(k, v) for k, v in strata.items() if k[0] >= 2]  # 高质量层
        if important_strata:
            min_per_stratum = max(1, remaining_quota // (len(important_strata) * 2))
            for stratum_key, stratum_copies in important_strata:
                take_count = min(min_per_stratum, len(stratum_copies), remaining_quota)
                if take_count > 0:
                    # 层内最优选择
                    best_copies = sorted(stratum_copies, key=lambda x: x['quality_score'], reverse=True)[:take_count]
                    selected.extend(best_copies)
                    remaining_quota -= take_count
                    if remaining_quota <= 0:
                        break
        
        # 第二轮：按层大小比例分配剩余配额
        if remaining_quota > 0:
            all_unselected = [c for c in copies if c not in selected]
            if all_unselected:
                # 简单随机采样填充剩余配额
                additional_count = min(remaining_quota, len(all_unselected))
                additional = random.sample(all_unselected, additional_count)
                selected.extend(additional)
        
        return selected[:max_sample]
    
    def _fast_stratified_sampling(self, copies: List[Dict], max_sample: int, seed_seq: Dict) -> List[Dict]:
        """超大规模数据：极速分层抽样"""
        if len(copies) <= max_sample:
            return copies
        
        logger.info(f"Ultra-high copy count: {len(copies)}, using fast stratified sampling")
        
        # 快速三层采样：不进行复杂计算
        # 1. 顶级质量 (20%)
        top_quota = max(5, max_sample // 5)
        top_copies = sorted(copies, key=lambda x: x['quality_score'], reverse=True)[:top_quota]
        
        # 2. 多样性代表 (60%) - 简化多样性采样
        diversity_quota = max_sample * 3 // 5
        top_ids = set(c['id'] for c in top_copies)
        remaining_copies = [c for c in copies if c['id'] not in top_ids]
        diversity_copies = self._ultra_fast_diversity_sampling(remaining_copies, diversity_quota)
        
        # 3. 随机补充 (20%)
        random_quota = max_sample - len(top_copies) - len(diversity_copies)
        if random_quota > 0:
            # 使用ID集合而不是直接比较字典对象
            selected_ids = set(c['id'] for c in top_copies + diversity_copies)
            unselected = [c for c in copies if c['id'] not in selected_ids]
            if unselected:
                random_supplement = random.sample(unselected, min(random_quota, len(unselected)))
            else:
                random_supplement = []
        else:
            random_supplement = []
        
        final_selected = top_copies + diversity_copies + random_supplement
        return final_selected[:max_sample]
    
    def _ultra_fast_diversity_sampling(self, copies: List[Dict], quota: int) -> List[Dict]:
        """超快速多样性采样 - 避免任何聚类计算"""
        if len(copies) <= quota:
            return copies
        
        # 基于预计算的分箱进行采样
        bins = defaultdict(list)
        for copy in copies:
            # 使用identity_tier和length_bin创建复合键
            bin_key = (copy['identity_tier'], copy['length_bin'])
            bins[bin_key].append(copy)
        
        # 从每个bin均匀采样
        selected = []
        bins_list = list(bins.items())
        per_bin_quota = max(1, quota // len(bins_list))
        
        for bin_key, bin_copies in bins_list:
            if len(selected) >= quota:
                break
            
            take_count = min(per_bin_quota, len(bin_copies), quota - len(selected))
            if take_count > 0:
                # 简单的分布式采样：每隔一定间隔取一个
                if len(bin_copies) <= take_count:
                    selected.extend(bin_copies)
                else:
                    step = len(bin_copies) // take_count
                    for i in range(0, len(bin_copies), step):
                        if len(selected) < quota:
                            selected.append(bin_copies[i])
        
        # 如果还有配额，随机补充
        if len(selected) < quota:
            unselected = [c for c in copies if c not in selected]
            if unselected:
                additional_count = min(quota - len(selected), len(unselected))
                additional = random.sample(unselected, additional_count) 
                selected.extend(additional)
        
        return selected[:quota]
    
    def _simple_diversity_sampling(self, copies: List[Dict], quota: int) -> List[Dict]:
        """简化的多样性采样"""
        if len(copies) <= quota:
            return copies
        
        # 基于identity分布采样
        identity_sorted = sorted(copies, key=lambda x: x.get('identity', 0))
        selected = []
        
        # 等间隔采样
        if quota > 0:
            step = len(identity_sorted) / quota
            for i in range(quota):
                index = int(i * step)
                if index < len(identity_sorted):
                    selected.append(identity_sorted[index])
        
        return selected
    
    def _get_quality_tier(self, quality_score: float) -> int:
        """获取质量层级"""
        if quality_score >= 0.8:
            return 3  # 顶级
        elif quality_score >= 0.6:
            return 2  # 高质量
        elif quality_score >= 0.4:
            return 1  # 中等质量
        else:
            return 0  # 低质量


def apply_optimized_te_sampling(copies: List[Dict], seed_seq: Dict, config) -> List[Dict]:
    """
    优化TE采样的接口函数
    自动根据拷贝数量选择最适合的算法
    """
    if not copies:
        return []
    
    # 小规模数据使用原有策略，大规模使用优化策略
    if len(copies) <= 100:
        # 小规模：使用现有的完整分析
        from te_aware_sampling_strategy import apply_te_aware_sampling
        return apply_te_aware_sampling(copies, seed_seq, config)
    else:
        # 大规模：使用优化策略
        strategy = OptimizedTESamplingStrategy(config)
        return strategy.select_representative_copies(copies, seed_seq)