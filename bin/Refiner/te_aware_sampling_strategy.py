"""
TE感知的拷贝采样策略
基于TE生物学特性的智能拷贝选择算法
"""

import logging
import numpy as np
from typing import List, Dict, Tuple, Optional
from collections import defaultdict, Counter
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
import math

logger = logging.getLogger(__name__)

class TEAwareSamplingStrategy:
    """基于TE生物学特性的拷贝采样策略"""
    
    def __init__(self, config):
        self.config = config
        self.base_limit = getattr(config, 'max_recruits_per_family', 30)
        
        # TE类型特征参数（基于充分性分析优化）
        self.te_type_params = {
            'LINE': {'max_copies': 60, 'diversity_weight': 0.7, 'length_weight': 0.8},
            'SINE': {'max_copies': 50, 'diversity_weight': 0.6, 'length_weight': 0.6},
            'LTR': {'max_copies': 55, 'diversity_weight': 0.8, 'length_weight': 0.9},
            'DNA': {'max_copies': 45, 'diversity_weight': 0.9, 'length_weight': 0.7},
            'Unknown': {'max_copies': 50, 'diversity_weight': 0.7, 'length_weight': 0.7}
        }
    
    def select_representative_copies(self, copies: List[Dict], seed_seq: Dict) -> List[Dict]:
        """
        选择代表性拷贝的主函数
        基于TE特性的多层次采样策略
        """
        if not copies:
            return []
        
        if len(copies) <= self.base_limit:
            return copies
        
        logger.info(f"Selecting representative copies from {len(copies)} total copies for {seed_seq['id']}")
        
        # 1. 预处理：质量评分和TE类型识别
        enriched_copies = self._enrich_copy_metadata(copies, seed_seq)
        
        # 2. TE类型识别
        te_type = self._infer_te_type(enriched_copies, seed_seq)
        te_params = self.te_type_params.get(te_type, self.te_type_params['Unknown'])
        
        # 3. 动态调整采样数量
        target_count = self._calculate_dynamic_target(enriched_copies, te_params)
        
        logger.info(f"TE type: {te_type}, target copies: {target_count}")
        
        # 4. 多层次采样策略
        if len(enriched_copies) <= target_count * 1.2:
            # 接近目标数量，简单质量过滤
            selected = self._quality_based_selection(enriched_copies, target_count)
        else:
            # 大量拷贝，使用分层采样
            selected = self._hierarchical_sampling(enriched_copies, target_count, te_params)
        
        logger.info(f"Selected {len(selected)} representative copies")
        return selected
    
    def _enrich_copy_metadata(self, copies: List[Dict], seed_seq: Dict) -> List[Dict]:
        """增强拷贝元数据"""
        enriched = []
        
        for copy in copies:
            enriched_copy = copy.copy()
            
            # 计算综合质量分数
            quality_score = self._calculate_copy_quality(copy, seed_seq)
            enriched_copy['quality_score'] = quality_score
            
            # 计算相对长度（与种子序列比较）
            seed_length = len(seed_seq.get('sequence', ''))
            copy_length = copy.get('length', len(copy.get('sequence', '')))
            enriched_copy['length_ratio'] = copy_length / seed_length if seed_length > 0 else 0
            
            # 计算完整性评分
            enriched_copy['completeness'] = self._assess_completeness(copy, seed_seq)
            
            # 检测结构特征
            enriched_copy['structural_features'] = self._detect_structural_features(copy)
            
            enriched.append(enriched_copy)
        
        return enriched
    
    def _calculate_copy_quality(self, copy: Dict, seed_seq: Dict) -> float:
        """计算拷贝的综合质量分数"""
        components = []
        
        # 1. 序列相似度 (30%)
        identity = copy.get('identity', 0) / 100.0
        components.append(identity * 0.30)
        
        # 2. 比对分数归一化 (25%)
        score = copy.get('score', 0)
        # 简单归一化：假设好的分数在200-1000范围
        normalized_score = min(score / 500.0, 1.0) if score > 0 else 0
        components.append(normalized_score * 0.25)
        
        # 3. 长度完整性 (20%)
        seed_length = len(seed_seq.get('sequence', ''))
        copy_length = copy.get('length', len(copy.get('sequence', '')))
        length_ratio = copy_length / seed_length if seed_length > 0 else 0
        # 长度在80%-120%范围内最优
        length_score = 1.0 - abs(length_ratio - 1.0) if 0.5 <= length_ratio <= 1.5 else 0.3
        components.append(length_score * 0.20)
        
        # 4. TSD存在性 (15%)
        tsd_score = 1.0 if copy.get('has_tsd', False) else 0.3
        components.append(tsd_score * 0.15)
        
        # 5. 基本质量过滤 (10%)
        basic_quality = 1.0 if copy_length >= 50 and identity >= 50 else 0.2
        components.append(basic_quality * 0.10)
        
        return sum(components)
    
    def _infer_te_type(self, copies: List[Dict], seed_seq: Dict) -> str:
        """推断TE类型"""
        # 基于种子序列长度和拷贝特征推断
        seed_length = len(seed_seq.get('sequence', ''))
        avg_length = np.mean([c.get('length', 0) for c in copies])
        
        # 长度特征判断
        if avg_length < 500:
            return 'SINE'  # 短序列，可能是SINE
        elif avg_length > 5000:
            return 'LINE'  # 长序列，可能是LINE
        elif 1000 <= avg_length <= 5000:
            # 检查是否有LTR特征
            has_ltr_features = any(c.get('has_tsd', False) for c in copies[:10])
            if has_ltr_features:
                return 'LTR'
            else:
                return 'DNA'
        else:
            return 'Unknown'
    
    def _calculate_dynamic_target(self, copies: List[Dict], te_params: Dict) -> int:
        """动态计算目标拷贝数"""
        base_target = te_params['max_copies']
        
        # 根据分化程度调整
        identities = [c.get('identity', 0) for c in copies]
        diversity_score = np.std(identities) / 100.0  # 标准化到0-1
        
        # 高分化家族需要更多代表
        diversity_adjustment = 1.0 + (diversity_score * te_params['diversity_weight'])
        
        # 根据总拷贝数调整（对数增长）
        total_copies = len(copies)
        size_adjustment = 1.0 + math.log10(total_copies / 100.0) * 0.3 if total_copies > 100 else 1.0
        
        # 综合调整
        target = int(base_target * diversity_adjustment * size_adjustment)
        
        # 基于分化程度的额外调整
        if diversity_score > 0.15:  # 高度分化家族
            target = int(target * 1.3)  # 增加30%拷贝数
            logger.debug(f"High diversity detected, increased target to {target}")
        elif diversity_score < 0.05:  # 低分化年轻家族
            target = int(target * 0.8)  # 减少20%拷贝数
            logger.debug(f"Low diversity detected, reduced target to {target}")
        
        # 限制在合理范围内（基于计算能力和生物学需求的平衡）
        if total_copies > 20000:
            # 超高拷贝数：为了计算效率，适度限制但保证多样性覆盖
            target = max(30, min(target, 100))
            logger.debug(f"Ultra-high copy count ({total_copies}), target: {target}")
        elif total_copies > 10000:
            # 高拷贝数：保持较高采样数以捕获多样性
            target = max(25, min(target, 120))
            logger.debug(f"High copy count ({total_copies}), target: {target}")
        elif total_copies > 5000:
            # 中高拷贝数：标准采样（与phase2处理能力匹配）
            target = max(20, min(target, 120))
        else:
            # 普通序列：允许更多采样但考虑处理效率
            target = max(15, min(target, 120))
        
        return target
    
    def _hierarchical_sampling(self, copies: List[Dict], target_count: int, te_params: Dict) -> List[Dict]:
        """分层采样策略"""
        
        # 1. 强制选择高质量拷贝 (30%配额)
        high_quality_quota = max(3, target_count // 3)
        high_quality = sorted(copies, key=lambda x: x['quality_score'], reverse=True)[:high_quality_quota]
        
        # 2. 基于相似度分层采样 (50%配额)
        diversity_quota = target_count // 2
        high_quality_ids = {c['id'] for c in high_quality}
        diversity_samples = self._diversity_based_sampling(
            [c for c in copies if c['id'] not in high_quality_ids], 
            diversity_quota
        )
        
        # 3. 结构特征保护采样 (20%配额)
        structure_quota = target_count - len(high_quality) - len(diversity_samples)
        selected_ids = high_quality_ids | {c['id'] for c in diversity_samples}
        structure_samples = self._structure_aware_sampling(
            [c for c in copies if c['id'] not in selected_ids],
            structure_quota
        )
        
        # 合并结果
        selected = high_quality + diversity_samples + structure_samples
        
        # 如果还不够，随机补充
        if len(selected) < target_count:
            all_selected_ids = selected_ids | {c['id'] for c in structure_samples}
            remaining = [c for c in copies if c['id'] not in all_selected_ids]
            if remaining:
                additional_count = min(target_count - len(selected), len(remaining))
                # 使用索引来避免字典比较问题
                indices = np.random.choice(len(remaining), additional_count, replace=False)
                additional = [remaining[i] for i in indices]
                selected.extend(additional)
        
        return selected[:target_count]
    
    def _diversity_based_sampling(self, copies: List[Dict], quota: int) -> List[Dict]:
        """基于多样性的采样"""
        if len(copies) <= quota:
            return copies
        
        # 基于相似度聚类
        try:
            # 使用identity作为距离度量
            identities = np.array([c.get('identity', 0) for c in copies])
            
            # 如果identity差异很小，使用长度作为补充
            if np.std(identities) < 5:  # identity标准差小于5%
                lengths = np.array([c.get('length', 0) for c in copies])
                # 组合identity和length特征
                features = np.column_stack([identities / 100.0, lengths / np.max(lengths)])
            else:
                features = identities.reshape(-1, 1)
            
            # 层次聚类
            distances = pdist(features, metric='euclidean')
            linkage_matrix = linkage(distances, method='ward')
            
            # 动态确定聚类数
            n_clusters = min(quota, max(3, len(copies) // 5))
            cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
            
            # 从每个聚类选择代表
            selected = []
            clusters = defaultdict(list)
            for i, label in enumerate(cluster_labels):
                clusters[label].append(copies[i])
            
            for cluster_copies in clusters.values():
                # 从每个聚类选择最佳代表
                best = max(cluster_copies, key=lambda x: x['quality_score'])
                selected.append(best)
            
            return selected[:quota]
            
        except Exception as e:
            logger.warning(f"Diversity sampling failed: {e}, using quality-based fallback")
            return sorted(copies, key=lambda x: x['quality_score'], reverse=True)[:quota]
    
    def _structure_aware_sampling(self, copies: List[Dict], quota: int) -> List[Dict]:
        """结构感知采样"""
        if not copies or quota <= 0:
            return []
        
        # 优先级：TSD > 完整长度 > 高分数
        priorities = []
        
        for copy in copies:
            priority = 0
            
            # TSD加分
            if copy.get('has_tsd', False):
                priority += 100
            
            # 完整性加分
            completeness = copy.get('completeness', 0)
            priority += completeness * 50
            
            # 质量分数
            priority += copy.get('quality_score', 0) * 20
            
            priorities.append((priority, copy))
        
        # 按优先级排序选择 (稳定排序，避免字典比较)
        priorities.sort(key=lambda x: x[0], reverse=True)
        return [copy for _, copy in priorities[:quota]]
    
    def _quality_based_selection(self, copies: List[Dict], target_count: int) -> List[Dict]:
        """简单的质量优先选择"""
        return sorted(copies, key=lambda x: x['quality_score'], reverse=True)[:target_count]
    
    def _assess_completeness(self, copy: Dict, seed_seq: Dict) -> float:
        """评估拷贝完整性"""
        seed_length = len(seed_seq.get('sequence', ''))
        copy_length = copy.get('length', 0)
        
        if seed_length == 0:
            return 0.5
        
        # 长度完整性
        length_completeness = min(copy_length / seed_length, 1.0)
        
        # 相似度作为质量指标
        identity_factor = copy.get('identity', 0) / 100.0
        
        # 综合评分
        return (length_completeness * 0.7 + identity_factor * 0.3)
    
    def _detect_structural_features(self, copy: Dict) -> Dict:
        """检测结构特征"""
        features = {
            'has_tsd': copy.get('has_tsd', False),
            'tsd_length': len(copy.get('tsd', '')) if copy.get('tsd') else 0,
            'is_full_length': copy.get('length_ratio', 0) > 0.8,
            'high_identity': copy.get('identity', 0) > 80
        }
        return features


def apply_te_aware_sampling(copies: List[Dict], seed_seq: Dict, config) -> List[Dict]:
    """
    应用TE感知采样策略的接口函数
    可以直接在phase2中调用
    """
    if not copies:
        return []
    
    strategy = TEAwareSamplingStrategy(config)
    return strategy.select_representative_copies(copies, seed_seq)