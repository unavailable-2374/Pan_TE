import logging
import numpy as np
from collections import Counter
from typing import Dict

logger = logging.getLogger(__name__)

def calculate_dust_score(sequence: str, window: int = 64, threshold: int = 7) -> float:
    """
    计算序列的DUST复杂度分数
    返回值范围：0（高复杂度）到1（低复杂度）
    """
    if not sequence or len(sequence) < window:
        return 0.0
    
    sequence = sequence.upper()
    total_score = 0
    num_windows = 0
    
    # 滑动窗口计算
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        
        # 计算三联体频率
        triplets = Counter()
        for j in range(len(window_seq) - 2):
            triplet = window_seq[j:j + 3]
            if 'N' not in triplet:  # 忽略含N的三联体
                triplets[triplet] += 1
        
        # 计算窗口得分
        window_score = 0
        for count in triplets.values():
            if count > 1:
                window_score += count * (count - 1) / 2
        
        # 标准化得分
        max_score = (window - 2) * (window - 3) / 2
        if max_score > 0:
            normalized_score = window_score / max_score
            total_score += normalized_score
            num_windows += 1
    
    # 返回平均DUST分数
    if num_windows > 0:
        avg_score = total_score / num_windows
        return min(avg_score * threshold, 1.0)
    
    return 0.0

def calculate_shannon_entropy(sequence: str, k: int = 1) -> float:
    """
    计算序列的Shannon熵
    k: k-mer大小（1=单核苷酸，2=二核苷酸）
    返回值：熵值（0=低复杂度，高值=高复杂度）
    """
    if not sequence or len(sequence) < k:
        return 0.0
    
    sequence = sequence.upper()
    
    # 计算k-mer频率
    kmers = Counter()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if 'N' not in kmer:  # 忽略含N的k-mer
            kmers[kmer] += 1
    
    # 计算总数
    total = sum(kmers.values())
    if total == 0:
        return 0.0
    
    # 计算Shannon熵
    entropy = 0.0
    for count in kmers.values():
        if count > 0:
            prob = count / total
            entropy -= prob * np.log2(prob)
    
    return entropy

def calculate_sequence_complexity(sequence: str) -> Dict[str, float]:
    """
    计算序列的多维度复杂度指标
    """
    if not sequence:
        return {
            'dust_score': 0.0,
            'shannon_entropy_1': 0.0,
            'shannon_entropy_2': 0.0,
            'gc_content': 0.0,
            'n_content': 0.0,
            'complexity_score': 0.0
        }
    
    sequence = sequence.upper()
    seq_len = len(sequence)
    
    # 计算基本统计
    gc_count = sequence.count('G') + sequence.count('C')
    n_count = sequence.count('N')
    gc_content = gc_count / seq_len if seq_len > 0 else 0
    n_content = n_count / seq_len if seq_len > 0 else 0
    
    # 计算复杂度指标
    dust_score = calculate_dust_score(sequence)
    shannon_1 = calculate_shannon_entropy(sequence, k=1)
    shannon_2 = calculate_shannon_entropy(sequence, k=2)
    
    # 综合复杂度分数
    complexity_score = (
        0.4 * (1 - dust_score) +  # DUST分数（反向）
        0.3 * (shannon_1 / 2) +    # 单核苷酸熵（标准化）
        0.3 * (shannon_2 / 4)      # 二核苷酸熵（标准化）
    )
    
    return {
        'dust_score': dust_score,
        'shannon_entropy_1': shannon_1,
        'shannon_entropy_2': shannon_2,
        'gc_content': gc_content,
        'n_content': n_content,
        'complexity_score': complexity_score
    }