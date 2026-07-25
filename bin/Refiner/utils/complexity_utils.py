import logging
import numpy as np
from collections import Counter
from typing import Dict

logger = logging.getLogger(__name__)

def calculate_lowcomplexity_fraction(sequence: str, window: int = 64,
                                     dust_cut: float = 0.7) -> float:
    """Fraction of non-overlapping windows that are low-complexity (DUST > dust_cut).

    Complements the average DUST score: a consensus that is mostly simple-repeat with
    a complex patch can pass the global DUST gate yet be >50% junk by window. Returns
    0.0 (no low-complexity windows) .. 1.0 (entirely low-complexity).
    """
    s = sequence.upper()
    if len(s) < window:
        return 1.0 if calculate_dust_score(s) > dust_cut else 0.0
    n = lc = 0
    for i in range(0, len(s) - window + 1, window):
        n += 1
        if calculate_dust_score(s[i:i + window]) > dust_cut:
            lc += 1
    return lc / n if n else 0.0

def calculate_dust_score(sequence: str, window: int = 64,
                         threshold: int = 7) -> float:
    """Return the legacy DUST score using an O(sequence length) scan.

    The previous implementation rebuilt all triplet counts for every
    overlapping window.  This rolling form is numerically equivalent: adding a
    triplet with current count ``c`` increases c*(c-1)/2 by ``c``; removing a
    triplet with count ``c`` decreases it by ``c-1``.
    """
    if not sequence or len(sequence) < window or window < 4:
        return 0.0

    sequence = sequence.upper()
    triplets = Counter()
    for start in range(window - 2):
        triplet = sequence[start:start + 3]
        if 'N' not in triplet:
            triplets[triplet] += 1

    window_score = sum(
        count * (count - 1) / 2
        for count in triplets.values() if count > 1)
    total_score = window_score
    num_windows = len(sequence) - window + 1

    for window_start in range(1, num_windows):
        outgoing = sequence[window_start - 1:window_start + 2]
        if 'N' not in outgoing:
            count = triplets[outgoing]
            window_score -= count - 1
            if count == 1:
                del triplets[outgoing]
            else:
                triplets[outgoing] = count - 1

        incoming_start = window_start + window - 3
        incoming = sequence[incoming_start:incoming_start + 3]
        if 'N' not in incoming:
            count = triplets[incoming]
            window_score += count
            triplets[incoming] = count + 1

        total_score += window_score

    max_score = (window - 2) * (window - 3) / 2
    average_score = total_score / num_windows
    return min(average_score / max_score * threshold, 1.0)

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