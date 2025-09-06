import logging
from typing import Dict, List, Any
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import multiprocessing

from config import PipelineConfig
from utils.cache_utils import cache_result
from utils.robust_runner import RobustRunner

logger = logging.getLogger(__name__)

def _calculate_single_complexity_parallel(seq_data, dust_window, dust_threshold):
    """独立的进程级函数，用于并行计算单个序列的复杂度分数"""
    from utils.complexity_utils import calculate_dust_score, calculate_shannon_entropy
    
    try:
        seq = seq_data['sequence']
        seq_id = seq_data['id']
        
        # DUST分数
        dust_score = calculate_dust_score(
            seq, 
            window=dust_window,
            threshold=dust_threshold
        )
        
        # Shannon熵
        h_mono = calculate_shannon_entropy(seq, k=1)
        h_di = calculate_shannon_entropy(seq, k=2)
        
        # 基于TE特征的复杂度评分
        if dust_score > 0.9:
            complexity_score = 0.3  # 过于简单，可能是tandem repeat
        elif dust_score > 0.6:
            complexity_score = 0.8  # 典型的TE重复模式
        elif dust_score > 0.3:
            complexity_score = 1.0  # 理想的复杂度
        else:
            complexity_score = 0.7  # 过于复杂
        
        # 使用熵值进行微调
        entropy_factor = (h_mono / 2.0 + h_di / 4.0) / 2
        
        if entropy_factor > 0.5:
            complexity_score = min(1.0, complexity_score * 1.1)
        elif entropy_factor < 0.3:
            complexity_score = max(0.2, complexity_score * 0.9)
        
        return {
            'seq_id': seq_id,
            'complexity': complexity_score,
            'dust_score': dust_score
        }
    except Exception as e:
        print(f"Error processing sequence {seq_data.get('id', 'unknown')}: {e}")
        return None

class SequenceScreenerOptimized:
    """Phase 1: 智能筛选与多维度评分"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.sequences = []
        self.scores = {}
        self.rm_detailed_results = {}  # 保存详细的RepeatMasker结果
        self.runner = RobustRunner(config)
        
        # 基于基因组大小自动设置RepeatMasker快速模式
        self._configure_repeatmasker_mode()
        
    def _configure_repeatmasker_mode(self):
        """基于基因组大小自动配置RepeatMasker模式"""
        import os
        
        try:
            # 获取基因组文件大小
            genome_size = os.path.getsize(self.config.genome_file)
            
            # 如果基因组大于1GB，启用快速模式
            if genome_size > self.config.large_genome_threshold:
                self.config.repeatmasker_quick = True
                logger.info(f"Large genome detected ({genome_size / 1024**3:.2f} GB > "
                           f"{self.config.large_genome_threshold / 1024**3:.1f} GB). "
                           f"Enabling RepeatMasker quick mode for faster processing.")
            else:
                logger.info(f"Standard genome size ({genome_size / 1024**3:.2f} GB). "
                           f"Using RepeatMasker sensitive mode for higher accuracy.")
                
        except (OSError, IOError) as e:
            logger.warning(f"Could not determine genome file size: {e}. "
                          f"Using default RepeatMasker settings.")
    
    def _get_repeatmasker_params(self) -> Dict:
        """获取基于基因组大小优化的RepeatMasker参数"""
        if self.config.repeatmasker_quick:
            # 大基因组快速模式：去掉敏感模式，使用更宽松的参数
            params = {
                'no_is': True,
                'nolow': True,
                'cutoff': 250,  # 稍微提高cutoff，减少低质量匹配
                'pa': max(1, self.config.threads // 2)
            }
            logger.info("Using RepeatMasker quick mode parameters (no -s, higher cutoff)")
        else:
            # 小基因组敏感模式：保持原有的高精度参数
            params = {
                's': True,  # 敏感模式，适合Phase 2重用
                'no_is': True,
                'nolow': True,
                'cutoff': 200,  # 降低cutoff以获得更多hits
                'pa': max(1, self.config.threads // 2)
            }
            logger.info("Using RepeatMasker sensitive mode parameters (-s, lower cutoff)")
        
        return params
    
    def run(self) -> Dict[str, Any]:
        """执行Phase 1"""
        # 使用lambda包装整个Phase 1流程，实现完整的checkpoint
        def run_phase1():
            # 步骤1: 加载并基础过滤
            self.sequences = self.runner.run_with_checkpoint(
                self.load_and_filter,
                checkpoint_name="phase1_sequences"
            )
            
            # 步骤2: RepeatMasker前去冗余（使用phase3标准）
            self.sequences = self.runner.run_with_checkpoint(
                self.pre_repeatmasker_deduplication,
                checkpoint_name="phase1_deduplication"
            )
            logger.info(f"Phase 1 sequences after pre-RepeatMasker deduplication: {len(self.sequences)}")
            
            # 步骤3: 计算复杂度分数
            self.calculate_complexity_scores()
            
            # 步骤4: 计算覆盖度分数（RepeatMasker）
            self.calculate_coverage_scores()
            
            # 步骤5: 综合评分
            self.calculate_final_scores()
            
            # 步骤6: 分类为A/B/C类
            return self.categorize_sequences()
        
        # 使用checkpoint运行整个phase1
        return self.runner.run_with_checkpoint(
            run_phase1,
            checkpoint_name="phase1_complete"
        )
    
    def load_and_filter(self):
        """加载序列并进行基础过滤"""
        from Bio import SeqIO
        sequences = []
        
        logger.info(f"Loading sequences from {self.config.repeatscout_file}")
        
        for record in SeqIO.parse(self.config.repeatscout_file, "fasta"):
            seq_str = str(record.seq).upper()
            seq_len = len(seq_str)
            
            # 长度过滤
            if seq_len < self.config.min_length or seq_len > self.config.max_length:
                continue
            
            # N含量过滤
            n_count = seq_str.count('N')
            if n_count / seq_len > self.config.max_n_percent:
                continue
            
            # 极端GC含量过滤
            gc_count = seq_str.count('G') + seq_str.count('C')
            gc_content = gc_count / seq_len
            if gc_content < 0.1 or gc_content > 0.9:
                continue
            
            sequences.append({
                'id': record.id,
                'sequence': seq_str,
                'length': seq_len,
                'gc_content': gc_content
            })
        
        logger.info(f"Loaded {len(sequences)} sequences after filtering")
        return sequences
    
    def calculate_complexity_scores(self):
        """计算复杂度分数 - 基于TE生物学特征（多进程并行版本）"""
        logger.info(f"Calculating complexity scores for {len(self.sequences)} sequences using {self.config.threads} processes...")
        
        # 准备参数
        dust_window = self.config.dust_window
        dust_threshold = self.config.dust_threshold
        
        # 使用ProcessPoolExecutor进行真正的并行处理
        with ProcessPoolExecutor(max_workers=self.config.threads) as executor:
            # 使用partial固定配置参数
            calc_func = partial(_calculate_single_complexity_parallel, 
                               dust_window=dust_window,
                               dust_threshold=dust_threshold)
            
            # 批量提交并获取结果
            results = list(executor.map(calc_func, self.sequences))
        
        # 收集结果到scores字典
        for result in results:
            if result:  # 确保结果有效
                seq_id = result['seq_id']
                if seq_id not in self.scores:
                    self.scores[seq_id] = {}
                self.scores[seq_id]['complexity'] = result['complexity']
                self.scores[seq_id]['dust_score'] = result['dust_score']
        
        # 输出复杂度分数统计
        complexity_scores = [self.scores[seq['id']]['complexity'] for seq in self.sequences if seq['id'] in self.scores]
        if complexity_scores:
            logger.info(f"Complexity scores - Min: {min(complexity_scores):.3f}, "
                       f"Max: {max(complexity_scores):.3f}, "
                       f"Mean: {np.mean(complexity_scores):.3f}")
    
    def calculate_coverage_scores(self):
        """使用RepeatMasker计算覆盖度分数 - 基于TE生物学特征"""
        from utils.alignment_utils import run_repeatmasker_batch_detailed
        
        # 基于基因组大小动态调整RepeatMasker参数
        rm_params = self._get_repeatmasker_params()
        
        # 批量运行RepeatMasker，获取详细结果
        rm_results = self.runner.run_with_retry(
            lambda: run_repeatmasker_batch_detailed(
                self.sequences,
                self.config.genome_file,
                params=rm_params,
                config=self.config
            )
        )
        
        # 保存详细结果供Phase 2使用
        self.rm_detailed_results = rm_results
        
        # 基于TE生物学特征的覆盖度评分
        for seq_id, rm_data in rm_results.items():
            copy_number = rm_data.get('copy_number', 0)
            avg_identity = rm_data.get('avg_identity', 0)
            
            # 拷贝数评分 - TE的核心特征是多拷贝
            if copy_number >= 10:
                copy_score = 1.0      # 明确的TE家族
            elif copy_number >= 5:
                copy_score = 0.8      # 良好的TE家族
            elif copy_number >= 3:
                copy_score = 0.6      # 最低可接受的TE
            elif copy_number >= 2:
                copy_score = 0.3      # 可疑，可能是segmental duplication
            else:
                copy_score = 0.1      # 不太可能是真正的TE
            
            # 序列一致性评分 - 反映TE家族的进化年龄和活性
            if avg_identity >= 75:
                identity_score = 1.0      # 年轻/活跃的TE
            elif avg_identity >= 65:
                identity_score = 0.8      # 中等年龄的TE
            elif avg_identity >= 55:
                identity_score = 0.6      # 古老但仍可识别的TE
            elif avg_identity >= 45:
                identity_score = 0.3      # 高度退化的TE
            else:
                identity_score = 0.1      # 可能不是同一TE家族
            
            if seq_id not in self.scores:
                self.scores[seq_id] = {}
            
            # 分别保存拷贝数和一致性分数，用于最终评分
            self.scores[seq_id]['copy_score'] = copy_score
            self.scores[seq_id]['identity_score'] = identity_score
            self.scores[seq_id]['copy_number'] = copy_number
            self.scores[seq_id]['avg_identity'] = avg_identity
        
        # 输出拷贝数和一致性统计
        copy_numbers = [self.scores[seq_id]['copy_number'] for seq_id in self.scores]
        identities = [self.scores[seq_id]['avg_identity'] for seq_id in self.scores]
        if copy_numbers:
            logger.info(f"Copy numbers - Min: {min(copy_numbers)}, Max: {max(copy_numbers)}, "
                       f"Mean: {np.mean(copy_numbers):.1f}")
            logger.info(f"Identities - Min: {min(identities):.1f}%, Max: {max(identities):.1f}%, "
                       f"Mean: {np.mean(identities):.1f}%")
    
    def calculate_final_scores(self):
        """计算综合分数 - 基于TE生物学特征的固定权重系统"""
        logger.info("Calculating final scores based on TE biological features...")
        
        for seq_data in self.sequences:
            seq_id = seq_data['id']
            sequence_length = seq_data['length']
            
            # 获取各项分数
            copy_score = self.scores[seq_id].get('copy_score', 0)
            identity_score = self.scores[seq_id].get('identity_score', 0)
            complexity_score = self.scores[seq_id].get('complexity', 0)
            
            # 长度合理性评分（不同TE类型有不同长度特征）
            if 100 <= sequence_length <= 500:
                # SINE/MITE的典型长度
                length_score = 1.0
            elif 500 < sequence_length <= 3000:
                # 大多数TE的理想长度范围
                length_score = 1.0
            elif 3000 < sequence_length <= 7000:
                # LINE/LTR的典型长度
                length_score = 0.9
            elif 7000 < sequence_length <= 15000:
                # 较长但仍合理
                length_score = 0.7
            else:
                # 过短（<100）或过长（>15000）
                length_score = 0.5
            
            # 边界质量（简化版）
            boundary_quality = self._calculate_boundary_quality(seq_data)
            
            # 基于TE生物学特征的权重分配
            # 拷贝数是最重要的TE特征（40%）
            # 一致性反映家族保守性（25%）
            # 复杂度排除简单重复（20%）
            # 长度合理性（10%）
            # 边界质量（5%）
            final_score = (
                0.40 * copy_score +
                0.25 * identity_score +
                0.20 * complexity_score +
                0.10 * length_score +
                0.05 * boundary_quality
            )
            
            # 保存所有分数用于调试
            self.scores[seq_id]['final'] = final_score
            self.scores[seq_id]['length_score'] = length_score
            self.scores[seq_id]['boundary_quality'] = boundary_quality
            
        # 输出最终分数统计
        final_scores = [self.scores[seq['id']]['final'] for seq in self.sequences]
        logger.info(f"Final scores - Min: {min(final_scores):.3f}, "
                   f"Max: {max(final_scores):.3f}, "
                   f"Mean: {np.mean(final_scores):.3f}, "
                   f"Median: {np.median(final_scores):.3f}")
    
    def categorize_sequences(self):
        """将序列分类为A/B/C三类 - 基于TE生物学特征的固定阈值"""
        a_sequences = []
        b_sequences = []
        c_sequences = []
        
        # 基于TE生物学特征的固定阈值
        # 这些阈值基于TE序列的质量标准，而不是数据分布
        # 调整阈值以获得更合理的分布
        A_THRESHOLD = 0.65  # 明确的高质量TE（多拷贝、高一致性、合理复杂度）
        B_THRESHOLD = 0.45  # 可能的TE，需要进一步验证
        
        for seq_data in self.sequences:
            seq_id = seq_data['id']
            final_score = self.scores[seq_id]['final']
            
            if final_score >= A_THRESHOLD:
                a_sequences.append(seq_data)
            elif final_score >= B_THRESHOLD:
                b_sequences.append(seq_data)
            else:
                c_sequences.append(seq_data)
        
        # 计算分数分布用于日志
        all_final_scores = [self.scores[seq['id']]['final'] for seq in self.sequences]
        score_p25 = np.percentile(all_final_scores, 25)
        score_p50 = np.percentile(all_final_scores, 50)
        score_p75 = np.percentile(all_final_scores, 75)
        
        # 输出分类统计
        logger.info(f"Score distribution - P25: {score_p25:.3f}, P50: {score_p50:.3f}, P75: {score_p75:.3f}")
        logger.info(f"Classification thresholds (fixed) - A: >={A_THRESHOLD:.2f}, B: >={B_THRESHOLD:.2f}")
        logger.info(f"Score range: {min(all_final_scores):.3f} - {max(all_final_scores):.3f}")
        
        # 输出前几个高分序列的详细信息
        top_sequences = sorted(self.sequences, 
                              key=lambda x: self.scores[x['id']]['final'], 
                              reverse=True)[:5]
        logger.info("Top 5 sequences (by final score):")
        for seq in top_sequences:
            sid = seq['id']
            s = self.scores[sid]
            logger.info(f"  {sid}: final={s['final']:.3f}, "
                       f"copy_score={s.get('copy_score', 0):.2f}(n={s.get('copy_number', 0)}), "
                       f"identity={s.get('identity_score', 0):.2f}({s.get('avg_identity', 0):.1f}%), "
                       f"complexity={s['complexity']:.2f}")
        
        logger.info(f"Categorization: A={len(a_sequences)}, B={len(b_sequences)}, C={len(c_sequences)}")
        
        # 分类质量检查
        total = len(self.sequences)
        a_pct = len(a_sequences) / total * 100 if total > 0 else 0
        b_pct = len(b_sequences) / total * 100 if total > 0 else 0
        c_pct = len(c_sequences) / total * 100 if total > 0 else 0
        
        logger.info(f"Distribution: A={a_pct:.1f}%, B={b_pct:.1f}%, C={c_pct:.1f}%")
        
        # 基于TE生物学的质量评估
        if len(a_sequences) < 10:
            logger.warning(f"Very few A-class sequences ({len(a_sequences)}). "
                          f"This might indicate low-quality input or stringent filtering.")
        elif len(a_sequences) > total * 0.5:
            logger.info(f"High proportion of A-class sequences ({a_pct:.1f}%). "
                       f"Input data appears to be high quality.")
        
        # 输出A类序列的特征统计
        if a_sequences:
            a_copy_numbers = [self.scores[seq['id']]['copy_number'] for seq in a_sequences]
            a_identities = [self.scores[seq['id']]['avg_identity'] for seq in a_sequences]
            logger.info(f"A-class stats: Copy number {np.mean(a_copy_numbers):.1f}±{np.std(a_copy_numbers):.1f}, "
                       f"Identity {np.mean(a_identities):.1f}±{np.std(a_identities):.1f}%")
        
        return {
            'a_sequences': a_sequences,
            'b_sequences': b_sequences,
            'c_sequences': c_sequences,
            'scores': self.scores,
            'rm_detailed_results': self.rm_detailed_results,  # 传递给Phase 2
            'summary': f"A:{len(a_sequences)}, B:{len(b_sequences)}, C:{len(c_sequences)}"
        }
    
    def _calculate_boundary_quality(self, seq_data):
        """计算序列边界质量指标"""
        sequence = seq_data['sequence']
        seq_len = len(sequence)
        
        if seq_len < 100:
            return 0.5  # 短序列默认中等质量
        
        # 计算两端和中间部分的GC含量差异
        edge_len = min(50, seq_len // 10)  # 取序列长度的10%或最多50bp
        
        left_edge = sequence[:edge_len]
        right_edge = sequence[-edge_len:]
        middle = sequence[seq_len//4:3*seq_len//4]
        
        def gc_content(seq):
            return (seq.count('G') + seq.count('C')) / len(seq) if len(seq) > 0 else 0
        
        gc_left = gc_content(left_edge)
        gc_right = gc_content(right_edge)
        gc_middle = gc_content(middle)
        
        # 计算边界与中间的GC差异（边界质量好的序列应该有相对稳定的GC分布）
        gc_variance = abs(gc_left - gc_middle) + abs(gc_right - gc_middle)
        
        # 计算N含量在边界的分布（好的边界应该N比较少）
        n_left = left_edge.count('N') / len(left_edge) if len(left_edge) > 0 else 0
        n_right = right_edge.count('N') / len(right_edge) if len(right_edge) > 0 else 0
        n_middle = middle.count('N') / len(middle) if len(middle) > 0 else 0
        
        # 边界质量分数（低的GC差异和低的N含量表示高质量）
        gc_score = max(0, 1.0 - gc_variance * 5)  # GC差异越小越好
        n_score = max(0, 1.0 - (n_left + n_right - n_middle))  # 边界N比中间少越好
        
        # 组合分数
        boundary_quality = 0.6 * gc_score + 0.4 * n_score
        
        return min(1.0, max(0.0, boundary_quality))
    
    def pre_repeatmasker_deduplication(self) -> List[Dict]:
        """RepeatMasker前去冗余处理，使用phase3标准（95%阈值）"""
        logger.info("Starting pre-RepeatMasker deduplication using phase3 standards")
        
        if not self.sequences:
            return []
        
        # 使用phase3的去冗余函数（95%阈值）
        from phase3_finalization_relaxed import run_cdhit_optimized
        
        # 转换序列格式为phase3期望的格式
        formatted_sequences = []
        for seq in self.sequences:
            formatted_sequences.append({
                'id': seq['id'],
                'sequence': seq['sequence'],
                'num_copies': 1,  # phase1阶段还没有copy信息
                'length': len(seq['sequence'])
            })
        
        # 按长度排序，确保长序列优先保留
        sorted_sequences = sorted(
            formatted_sequences,
            key=lambda x: x['length'],
            reverse=True
        )
        
        logger.info(f"Running CD-HIT with 90% threshold on {len(sorted_sequences)} sequences")
        
        # 运行CD-HIT去冗余（80%阈值）
        deduplicated = run_cdhit_optimized(
            sorted_sequences,
            threshold=0.80,  # 使用phase3的masking阈值
            threads=self.config.threads,
            config=self.config,
            label="phase1_prefilter"
        )
        
        # 转换回phase1的序列格式
        result_sequences = []
        for seq_data in deduplicated:
            # 找到原始序列
            original_seq = None
            for orig in self.sequences:
                if orig['id'] == seq_data['id']:
                    original_seq = orig
                    break
            
            if original_seq:
                result_sequences.append(original_seq)
            else:
                # 如果找不到原始序列，创建新的
                result_sequences.append({
                    'id': seq_data['id'],
                    'sequence': seq_data['sequence'],
                    'length': len(seq_data['sequence'])
                })
        
        logger.info(f"Pre-RepeatMasker deduplication: {len(self.sequences)} -> {len(result_sequences)} "
                   f"(removed {len(self.sequences) - len(result_sequences)} redundant sequences)")
        
        return result_sequences
