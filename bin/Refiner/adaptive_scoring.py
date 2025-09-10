"""
自适应TE评分系统
考虑序列长度、基因组大小和TE类型特征的动态评分
"""

import logging
import numpy as np
from typing import Dict, Tuple

logger = logging.getLogger(__name__)

class AdaptiveTEScorer:
    """
    自适应TE评分器，根据序列特征和基因组背景动态调整评分标准
    """
    
    def __init__(self, genome_size: int = None):
        """
        Args:
            genome_size: 基因组大小（bp），用于调整拷贝数期望值
        """
        self.genome_size = genome_size or 1e9  # 默认1Gb
        self.genome_size_gb = self.genome_size / 1e9
        
        # 基因组大小对拷贝数期望的影响因子
        self.genome_factor = np.log10(self.genome_size_gb + 1) + 1
        
    def calculate_adaptive_copy_score(self, copy_number: int, seq_length: int) -> float:
        """
        根据序列长度自适应计算拷贝数分数
        
        原理：
        - 短序列（SINE/MITE）天然倾向有更多拷贝
        - 长序列（LINE/LTR）通常拷贝数较少
        - 考虑基因组大小的影响
        - 拷贝数超过15个直接给满分
        """
        # 拷贝数超过15个直接给满分
        if copy_number >= 15:
            return 1.0
        
        # 基于序列长度的期望拷贝数，限制最大期望为15
        # 所有序列类型的最小拷贝数要求都是3（小于3个认为是噪音）
        if seq_length < 300:  # SINE类
            expected_copies = min(15, 50 / self.genome_factor)
            min_copies = 3
        elif seq_length < 800:  # MITE类
            expected_copies = min(15, 30 / self.genome_factor)
            min_copies = 3
        elif seq_length < 3000:  # 中等长度TE
            expected_copies = min(15, 20 / self.genome_factor)
            min_copies = 3
        elif seq_length < 7000:  # LINE/LTR类
            expected_copies = min(15, 10 / self.genome_factor)
            min_copies = 3
        else:  # 超长TE
            expected_copies = min(15, 5 / self.genome_factor)
            min_copies = 3
        
        # 确保expected_copies至少等于min_copies
        expected_copies = max(expected_copies, min_copies)
        
        # 计算相对拷贝数
        relative_copies = copy_number / expected_copies
        
        # 非线性评分函数
        if copy_number < min_copies:
            # 低于最小阈值（3个拷贝），视为噪音
            score = 0.01  # 极低分数，基本排除
        elif copy_number >= expected_copies:
            # 达到或超过期望值，给高分
            if copy_number >= expected_copies * 1.5:
                score = 1.0  # 远超期望值
            else:
                score = 0.8 + 0.2 * (copy_number - expected_copies) / (expected_copies * 0.5)
        elif relative_copies >= 0.5:
            # 接近期望值
            score = 0.6 + 0.2 * (relative_copies - 0.5) / 0.5
        elif relative_copies >= 0.25:
            # 低于期望但可接受
            score = 0.4 + 0.2 * (relative_copies - 0.25) / 0.25
        else:
            # 明显低于期望
            score = 0.4 * relative_copies / 0.25
        
        return min(1.0, max(0.0, score))
    
    def calculate_adaptive_identity_score(self, avg_identity: float, seq_length: int,
                                         copy_number: int) -> float:
        """
        根据序列特征自适应计算一致性分数
        
        原理：
        - 短序列通常进化更快，可以容忍更低的一致性
        - 高拷贝数的家族可能有更多的序列分化
        """
        # 基于序列长度的一致性期望
        if seq_length < 500:
            # 短序列，容忍更多变异
            thresholds = [65, 55, 45, 35]
        elif seq_length < 3000:
            # 中等长度
            thresholds = [70, 60, 50, 40]
        else:
            # 长序列，期望更高的保守性
            thresholds = [75, 65, 55, 45]
        
        # 根据拷贝数调整阈值（高拷贝数容忍更多变异）
        if copy_number > 100:
            thresholds = [t - 5 for t in thresholds]
        elif copy_number > 50:
            thresholds = [t - 3 for t in thresholds]
        
        # 计算分数
        if avg_identity >= thresholds[0]:
            score = 1.0
        elif avg_identity >= thresholds[1]:
            score = 0.8 + 0.2 * (avg_identity - thresholds[1]) / (thresholds[0] - thresholds[1])
        elif avg_identity >= thresholds[2]:
            score = 0.6 + 0.2 * (avg_identity - thresholds[2]) / (thresholds[1] - thresholds[2])
        elif avg_identity >= thresholds[3]:
            score = 0.3 + 0.3 * (avg_identity - thresholds[3]) / (thresholds[2] - thresholds[3])
        else:
            score = 0.3 * avg_identity / thresholds[3]
        
        return min(1.0, max(0.0, score))
    
    def calculate_length_score(self, seq_length: int, copy_number: int) -> float:
        """
        计算长度合理性分数，考虑拷贝数的影响
        
        原理：
        - 极短序列如果有很多拷贝，可能是简单重复而非TE
        - 极长序列如果只有单拷贝，可能是基因组片段而非TE
        """
        if seq_length < 50:
            # 太短，可能是简单重复
            if copy_number > 100:
                score = 0.1  # 高拷贝的短序列很可能是简单重复
            else:
                score = 0.3
        elif seq_length < 80:
            # 边缘长度
            score = 0.5 if copy_number > 10 else 0.4
        elif seq_length < 100:
            # 最短的TE（如MIR）
            score = 0.7
        elif seq_length <= 15000:
            # 正常TE长度范围
            if seq_length <= 500:
                score = 1.0  # SINE的理想长度
            elif seq_length <= 3000:
                score = 1.0  # 大多数TE的理想长度
            elif seq_length <= 7000:
                score = 0.95  # LINE/LTR的正常长度
            else:
                score = 0.85  # 较长但仍合理
        elif seq_length <= 20000:
            # 超长但可能是复合TE
            score = 0.7 if copy_number >= 2 else 0.5
        else:
            # 极长，可能有问题
            score = 0.4 if copy_number >= 3 else 0.2
        
        return score
    
    def calculate_adaptive_final_score(self, seq_data: Dict, scores: Dict) -> Tuple[float, str]:
        """
        计算自适应最终分数和分类
        
        Returns:
            (final_score, classification_reason)
        """
        seq_length = len(seq_data['sequence'])
        copy_number = scores.get('copy_number', 0)
        avg_identity = scores.get('avg_identity', 0)
        complexity_score = scores.get('complexity', 0)
        
        # 使用自适应评分
        copy_score = self.calculate_adaptive_copy_score(copy_number, seq_length)
        identity_score = self.calculate_adaptive_identity_score(avg_identity, seq_length, copy_number)
        length_score = self.calculate_length_score(seq_length, copy_number)
        
        # 边界质量（如果有）
        boundary_quality = scores.get('boundary_quality', 0.5)
        
        # 动态权重分配
        if seq_length < 500:
            # 短序列：拷贝数最重要
            weights = {
                'copy': 0.45,
                'identity': 0.20,
                'complexity': 0.20,
                'length': 0.10,
                'boundary': 0.05
            }
        elif seq_length < 3000:
            # 中等长度：平衡各因素
            weights = {
                'copy': 0.35,
                'identity': 0.25,
                'complexity': 0.20,
                'length': 0.10,
                'boundary': 0.10
            }
        else:
            # 长序列：结构和一致性更重要
            weights = {
                'copy': 0.30,
                'identity': 0.30,
                'complexity': 0.15,
                'length': 0.10,
                'boundary': 0.15
            }
        
        # 计算加权分数
        final_score = (
            weights['copy'] * copy_score +
            weights['identity'] * identity_score +
            weights['complexity'] * complexity_score +
            weights['length'] * length_score +
            weights['boundary'] * boundary_quality
        )
        
        # 生成分类理由
        reasons = []
        if copy_score >= 0.8:
            reasons.append(f"high_copy({copy_number})")
        elif copy_score <= 0.3:
            reasons.append(f"low_copy({copy_number})")
        
        if identity_score >= 0.8:
            reasons.append(f"high_identity({avg_identity:.1f}%)")
        elif identity_score <= 0.3:
            reasons.append(f"low_identity({avg_identity:.1f}%)")
        
        if seq_length < 100:
            reasons.append("very_short")
        elif seq_length > 10000:
            reasons.append("very_long")
        
        classification_reason = "_".join(reasons) if reasons else "standard"
        
        # 存储详细分数
        scores['adaptive_copy_score'] = copy_score
        scores['adaptive_identity_score'] = identity_score
        scores['adaptive_length_score'] = length_score
        scores['adaptive_final'] = final_score
        scores['classification_reason'] = classification_reason
        
        return final_score, classification_reason
    
    def get_adaptive_thresholds(self, sequences: list) -> Tuple[float, float]:
        """
        基于序列分布动态确定ABC分级阈值
        
        Returns:
            (A_threshold, B_threshold)
        """
        if not sequences:
            return 0.65, 0.45
        
        # 计算序列长度分布
        lengths = [len(seq['sequence']) for seq in sequences]
        median_length = np.median(lengths)
        
        # 基于中位长度调整阈值
        if median_length < 500:
            # 短序列为主，降低阈值
            A_threshold = 0.60
            B_threshold = 0.40
        elif median_length < 2000:
            # 标准长度
            A_threshold = 0.65
            B_threshold = 0.45
        else:
            # 长序列为主，可以稍微提高阈值
            A_threshold = 0.70
            B_threshold = 0.50
        
        # 基于基因组大小微调
        if self.genome_size_gb < 0.1:
            # 小基因组，降低拷贝数要求
            A_threshold -= 0.05
            B_threshold -= 0.05
        elif self.genome_size_gb > 3:
            # 大基因组，提高拷贝数要求
            A_threshold += 0.05
            B_threshold += 0.05
        
        return max(0.5, A_threshold), max(0.3, B_threshold)