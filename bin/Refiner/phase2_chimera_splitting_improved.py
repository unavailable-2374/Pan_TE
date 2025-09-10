"""
Phase 2: 嵌合体序列拆分（改进版本）
基于深度评估结果的优化实现
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Any, Optional
from pathlib import Path
from collections import defaultdict
import tempfile

from config import PipelineConfig
from utils.robust_runner import RobustRunner

logger = logging.getLogger(__name__)


class ImprovedChimeraSplitter:
    """改进的Phase 2嵌合体拆分器"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.runner = RobustRunner(config)
        
        # 自适应参数设置
        self.adaptive_gap_threshold = self._calculate_adaptive_threshold()
        self.min_segment_length = 30  # 降低最小片段长度要求
        self.edge_buffer = 30  # 降低边缘缓冲区要求
        
        logger.info(f"Improved Chimera Splitter initialized:")
        logger.info(f"  - Adaptive gap threshold: {self.adaptive_gap_threshold}bp")
        logger.info(f"  - Min segment length: {self.min_segment_length}bp")
        logger.info(f"  - Edge buffer: {self.edge_buffer}bp")
    
    def _calculate_adaptive_threshold(self) -> int:
        """计算自适应间隙阈值"""
        # 基于基因组大小和配置动态调整
        try:
            import os
            genome_size = os.path.getsize(self.config.genome_file)
            
            # 大基因组使用更严格的阈值
            if genome_size > 1e9:  # >1Gb
                base_threshold = 50
            elif genome_size > 100e6:  # >100Mb
                base_threshold = 40
            else:  # 小基因组
                base_threshold = 30
                
        except:
            base_threshold = 40  # 默认中等阈值
        
        # 允许通过配置覆盖
        return getattr(self.config, 'chimera_gap_threshold', base_threshold)
    
    def analyze_chimera_advanced(self, chimera: Dict) -> Dict[str, Any]:
        """高级嵌合体分析，整合多种策略"""
        
        chimera_id = chimera.get('id', 'unknown')
        sequence = chimera.get('sequence', '')
        
        if not sequence:
            return {'is_splittable': False, 'breakpoints': [], 'reason': 'empty_sequence'}
        
        # 策略1：基于RepeatMasker hits的间隙分析
        rm_breakpoints = self._analyze_repeatmasker_gaps(chimera)
        
        # 策略2：基于序列组成的突变点检测
        composition_breakpoints = self._detect_composition_changes(sequence)
        
        # 策略3：基于k-mer频率的边界检测
        kmer_breakpoints = self._detect_kmer_boundaries(sequence)
        
        # 策略4：TSD（目标位点重复）检测
        tsd_breakpoints = self._detect_tsd_boundaries(sequence)
        
        # 整合所有断点候选
        all_breakpoints = self._merge_breakpoint_candidates(
            rm_breakpoints, 
            composition_breakpoints,
            kmer_breakpoints,
            tsd_breakpoints,
            len(sequence)
        )
        
        # 评估拆分可行性
        if len(all_breakpoints) >= 1:  # 降低要求：1个高置信度断点即可
            # 验证是否为串联重复
            if self._is_tandem_repeat(sequence, all_breakpoints):
                return {
                    'is_splittable': False,
                    'breakpoints': [],
                    'reason': 'tandem_repeat_detected'
                }
            
            return {
                'is_splittable': True,
                'breakpoints': all_breakpoints,
                'reason': f'multi_strategy_{len(all_breakpoints)}_breakpoints',
                'confidence': self._calculate_confidence(all_breakpoints, sequence)
            }
        else:
            return {
                'is_splittable': False,
                'breakpoints': [],
                'reason': 'insufficient_evidence'
            }
    
    def _analyze_repeatmasker_gaps(self, chimera: Dict) -> List[int]:
        """改进的RepeatMasker间隙分析"""
        breakpoints = []
        
        # 从Phase 1结果或新分析获取hits
        hits = self._get_repeatmasker_hits(chimera)
        
        if len(hits) < 2:
            return []
        
        # 按位置排序
        sorted_hits = sorted(hits, key=lambda x: x.get('query_start', 0))
        
        # 检查序列开头的潜在断点
        if sorted_hits[0]['query_start'] > self.edge_buffer * 2:
            breakpoints.append(sorted_hits[0]['query_start'] // 2)
        
        # 检查hits之间的间隙
        for i in range(len(sorted_hits) - 1):
            current_end = sorted_hits[i]['query_end']
            next_start = sorted_hits[i + 1]['query_start']
            
            gap_size = next_start - current_end
            
            # 使用自适应阈值
            if gap_size > self.adaptive_gap_threshold:
                # 智能断点定位：考虑序列特征而非简单中点
                breakpoint = self._find_optimal_breakpoint(
                    chimera['sequence'][current_end:next_start],
                    current_end
                )
                breakpoints.append(breakpoint)
        
        # 检查序列末尾的潜在断点
        seq_length = len(chimera['sequence'])
        if sorted_hits[-1]['query_end'] < seq_length - self.edge_buffer * 2:
            breakpoints.append((sorted_hits[-1]['query_end'] + seq_length) // 2)
        
        return breakpoints
    
    def _detect_composition_changes(self, sequence: str) -> List[int]:
        """检测序列组成的突变点（GC含量、复杂度等）"""
        breakpoints = []
        window_size = 100
        step_size = 20
        
        if len(sequence) < window_size * 2:
            return []
        
        # 计算滑动窗口的GC含量
        gc_values = []
        positions = []
        
        for i in range(0, len(sequence) - window_size, step_size):
            window = sequence[i:i + window_size]
            gc_content = (window.count('G') + window.count('C')) / len(window)
            gc_values.append(gc_content)
            positions.append(i + window_size // 2)
        
        # 检测突变点
        if len(gc_values) > 3:
            gc_diffs = np.diff(gc_values)
            threshold = np.std(gc_diffs) * 2  # 2个标准差
            
            for i, diff in enumerate(gc_diffs):
                if abs(diff) > threshold and threshold > 0.05:  # 显著变化
                    potential_breakpoint = positions[i]
                    if self.edge_buffer <= potential_breakpoint <= len(sequence) - self.edge_buffer:
                        breakpoints.append(potential_breakpoint)
        
        return breakpoints
    
    def _detect_kmer_boundaries(self, sequence: str, k: int = 6) -> List[int]:
        """基于k-mer频率变化检测边界"""
        breakpoints = []
        
        if len(sequence) < 200:  # 序列太短
            return []
        
        window_size = 100
        step_size = 20
        
        # 计算每个窗口的k-mer频谱
        kmer_signatures = []
        positions = []
        
        for i in range(0, len(sequence) - window_size, step_size):
            window = sequence[i:i + window_size]
            signature = self._calculate_kmer_signature(window, k)
            kmer_signatures.append(signature)
            positions.append(i + window_size // 2)
        
        # 比较相邻窗口的k-mer频谱差异
        for i in range(len(kmer_signatures) - 1):
            distance = self._kmer_distance(kmer_signatures[i], kmer_signatures[i + 1])
            
            if distance > 0.3:  # 显著差异阈值
                potential_breakpoint = positions[i]
                if self.edge_buffer <= potential_breakpoint <= len(sequence) - self.edge_buffer:
                    breakpoints.append(potential_breakpoint)
        
        return breakpoints
    
    def _detect_tsd_boundaries(self, sequence: str) -> List[int]:
        """检测目标位点重复（TSD）以识别TE插入边界"""
        breakpoints = []
        
        # TSD通常为4-10bp
        for tsd_len in range(4, 11):
            for i in range(len(sequence) - tsd_len * 2 - 50):  # 至少50bp的TE
                potential_tsd1 = sequence[i:i + tsd_len]
                
                # 在合理距离内寻找匹配的TSD
                for j in range(i + 50, min(i + 20000, len(sequence) - tsd_len)):
                    potential_tsd2 = sequence[j:j + tsd_len]
                    
                    if potential_tsd1 == potential_tsd2:
                        # 找到潜在的TSD对
                        # 断点在第一个TSD之后
                        breakpoint = i + tsd_len
                        if self.edge_buffer <= breakpoint <= len(sequence) - self.edge_buffer:
                            breakpoints.append(breakpoint)
                        break
        
        return list(set(breakpoints))  # 去重
    
    def _merge_breakpoint_candidates(self, *breakpoint_lists, seq_length: int) -> List[int]:
        """整合多个策略的断点候选，提高置信度"""
        all_breakpoints = []
        for bp_list in breakpoint_lists:
            all_breakpoints.extend(bp_list)
        
        if not all_breakpoints:
            return []
        
        # 聚类相近的断点（容差20bp）
        clustered = self._cluster_breakpoints(all_breakpoints, tolerance=20)
        
        # 选择每个簇的代表（支持度最高的位置）
        final_breakpoints = []
        for cluster in clustered:
            if len(cluster) >= 2:  # 至少2个策略支持
                # 使用中位数作为代表
                representative = int(np.median(cluster))
                if self.edge_buffer <= representative <= seq_length - self.edge_buffer:
                    final_breakpoints.append(representative)
        
        return sorted(final_breakpoints)
    
    def _is_tandem_repeat(self, sequence: str, breakpoints: List[int]) -> bool:
        """判断是否为串联重复而非真实嵌合体"""
        if len(breakpoints) < 2:
            return False
        
        # 检查片段之间的相似性
        positions = [0] + sorted(breakpoints) + [len(sequence)]
        segments = []
        
        for i in range(len(positions) - 1):
            segment = sequence[positions[i]:positions[i + 1]]
            if len(segment) >= self.min_segment_length:
                segments.append(segment)
        
        if len(segments) < 2:
            return False
        
        # 计算片段间的相似度
        high_similarity_pairs = 0
        total_pairs = 0
        
        for i in range(len(segments)):
            for j in range(i + 1, len(segments)):
                similarity = self._sequence_similarity(segments[i], segments[j])
                total_pairs += 1
                if similarity > 0.8:  # 80%相似度
                    high_similarity_pairs += 1
        
        # 如果大部分片段对都高度相似，可能是串联重复
        if total_pairs > 0 and high_similarity_pairs / total_pairs > 0.6:
            return True
        
        return False
    
    def _find_optimal_breakpoint(self, gap_sequence: str, offset: int) -> int:
        """在间隙中找到最优断点位置"""
        # 简单实现：寻找最低复杂度的位置
        if len(gap_sequence) < 10:
            return offset + len(gap_sequence) // 2
        
        min_complexity = float('inf')
        best_pos = len(gap_sequence) // 2
        
        window = 10
        for i in range(window, len(gap_sequence) - window):
            subseq = gap_sequence[i-window:i+window]
            complexity = len(set(subseq))  # 简单的复杂度度量
            if complexity < min_complexity:
                min_complexity = complexity
                best_pos = i
        
        return offset + best_pos
    
    def _calculate_kmer_signature(self, sequence: str, k: int) -> Dict[str, float]:
        """计算序列的k-mer频率签名"""
        kmer_counts = defaultdict(int)
        total = 0
        
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            kmer_counts[kmer] += 1
            total += 1
        
        # 归一化
        signature = {}
        for kmer, count in kmer_counts.items():
            signature[kmer] = count / total if total > 0 else 0
        
        return signature
    
    def _kmer_distance(self, sig1: Dict, sig2: Dict) -> float:
        """计算两个k-mer签名之间的距离"""
        all_kmers = set(sig1.keys()) | set(sig2.keys())
        
        if not all_kmers:
            return 0
        
        distance = 0
        for kmer in all_kmers:
            freq1 = sig1.get(kmer, 0)
            freq2 = sig2.get(kmer, 0)
            distance += abs(freq1 - freq2)
        
        return distance / 2  # 归一化到[0,1]
    
    def _cluster_breakpoints(self, breakpoints: List[int], tolerance: int) -> List[List[int]]:
        """聚类相近的断点"""
        if not breakpoints:
            return []
        
        sorted_bp = sorted(breakpoints)
        clusters = [[sorted_bp[0]]]
        
        for bp in sorted_bp[1:]:
            if bp - clusters[-1][-1] <= tolerance:
                clusters[-1].append(bp)
            else:
                clusters.append([bp])
        
        return clusters
    
    def _sequence_similarity(self, seq1: str, seq2: str) -> float:
        """计算两个序列的相似度（简化版）"""
        if not seq1 or not seq2:
            return 0
        
        # 使用较短序列的长度
        min_len = min(len(seq1), len(seq2))
        matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
        
        return matches / min_len
    
    def _calculate_confidence(self, breakpoints: List[int], sequence: str) -> float:
        """计算拆分决策的置信度"""
        if not breakpoints:
            return 0
        
        # 基于多个因素计算置信度
        confidence = 0.5  # 基础置信度
        
        # 断点数量
        if len(breakpoints) >= 2:
            confidence += 0.2
        
        # 片段长度合理性
        positions = [0] + sorted(breakpoints) + [len(sequence)]
        segment_lengths = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
        
        if all(length >= 50 for length in segment_lengths):
            confidence += 0.2
        
        # 片段数量合理性（2-4个片段最可信）
        if 2 <= len(segment_lengths) <= 4:
            confidence += 0.1
        
        return min(1.0, confidence)
    
    def _get_repeatmasker_hits(self, chimera: Dict) -> List[Dict]:
        """获取RepeatMasker hits（从Phase 1或新分析）"""
        # 这里简化处理，实际应该整合Phase 1结果
        return chimera.get('rm_hits', [])