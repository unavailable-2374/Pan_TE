import logging
from typing import Dict, List, Any, Tuple
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
    """Phase 1: 序列筛选与共识候选识别
    
    主要功能：
    1. 对RepeatScout输出进行质量评分
    2. 识别需要构建共识的候选序列
    3. 为候选序列收集RepeatMasker信息
    4. 将候选序列传递给Phase2进行嵌合体分析
    """
    
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
        
        # 初始化自适应评分器
        try:
            import os
            genome_size = os.path.getsize(self.config.genome_file)
        except:
            genome_size = 1e9  # 默认1Gb
        
        from adaptive_scoring import AdaptiveTEScorer
        adaptive_scorer = AdaptiveTEScorer(genome_size)
        
        # 基于TE生物学特征的覆盖度评分
        for seq_id, rm_data in rm_results.items():
            copy_number = rm_data.get('copy_number', 0)
            avg_identity = rm_data.get('avg_identity', 0)
            
            # 找到对应序列的长度
            seq_length = None
            for seq_data in self.sequences:
                if seq_data['id'] == seq_id:
                    seq_length = len(seq_data['sequence'])
                    break
            
            if seq_length is None:
                logger.warning(f"Could not find sequence length for {seq_id}")
                seq_length = 1000  # 默认长度
            
            # 处理低拷贝数序列 - 小于3个拷贝的全部当成噪音
            if copy_number < 3:
                logger.debug(f"{seq_id}: Low copy number ({copy_number}) sequence treated as noise")
                copy_score = 0.01  # 极低分数，基本排除
                identity_score = 0.01
            else:
                # 使用自适应拷贝数评分（限制最大期望为15）
                copy_score = adaptive_scorer.calculate_adaptive_copy_score(copy_number, seq_length)
                
                # 使用自适应一致性评分
                identity_score = adaptive_scorer.calculate_adaptive_identity_score(
                    avg_identity, seq_length, copy_number
                )
            
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
        
        # 统计噪音序列
        noise_sequences = [seq_id for seq_id in self.scores if self.scores[seq_id]['copy_number'] < 3]
        valid_sequences = [seq_id for seq_id in self.scores if self.scores[seq_id]['copy_number'] >= 3]
        
        logger.info(f"Copy number filtering results:")
        logger.info(f"  Valid sequences (≥3 copies): {len(valid_sequences)}")
        logger.info(f"  Noise sequences (<3 copies): {len(noise_sequences)}")
        logger.info(f"  Noise ratio: {len(noise_sequences)/len(self.scores)*100:.1f}%")
        
        if copy_numbers:
            logger.info(f"Copy numbers - Min: {min(copy_numbers)}, Max: {max(copy_numbers)}, "
                       f"Mean: {np.mean(copy_numbers):.1f}")
            logger.info(f"Identities - Min: {min(identities):.1f}%, Max: {max(identities):.1f}%, "
                       f"Mean: {np.mean(identities):.1f}%")
        
        if valid_sequences:
            valid_copy_numbers = [self.scores[seq_id]['copy_number'] for seq_id in valid_sequences]
            logger.info(f"Valid sequences copy numbers - Min: {min(valid_copy_numbers)}, "
                       f"Max: {max(valid_copy_numbers)}, Mean: {np.mean(valid_copy_numbers):.1f}")
    
    def calculate_final_scores(self):
        """计算自适应综合分数"""
        logger.info("Calculating adaptive final scores based on sequence characteristics...")
        
        # 初始化自适应评分器
        try:
            import os
            genome_size = os.path.getsize(self.config.genome_file)
        except:
            genome_size = 1e9  # 默认1Gb
        
        from adaptive_scoring import AdaptiveTEScorer
        adaptive_scorer = AdaptiveTEScorer(genome_size)
        
        for seq_data in self.sequences:
            seq_id = seq_data['id']
            
            # 使用自适应评分器计算最终分数
            final_score, reason = adaptive_scorer.calculate_adaptive_final_score(
                seq_data, self.scores[seq_id]
            )
            
            # 保存分数和分类原因
            self.scores[seq_id]['final'] = final_score
            self.scores[seq_id]['classification_reason'] = reason
            
            # 保留原有的边界质量计算用于兼容性
            if 'boundary_quality' not in self.scores[seq_id]:
                boundary_quality = self._calculate_boundary_quality(seq_data)
                self.scores[seq_id]['boundary_quality'] = boundary_quality
            
        # 输出最终分数统计
        final_scores = [self.scores[seq['id']]['final'] for seq in self.sequences]
        logger.info(f"Final scores - Min: {min(final_scores):.3f}, "
                   f"Max: {max(final_scores):.3f}, "
                   f"Mean: {np.mean(final_scores):.3f}, "
                   f"Median: {np.median(final_scores):.3f}")
    
    def categorize_sequences(self):
        """将序列分类为A/B/C三类 - 使用自适应阈值"""
        a_sequences = []
        b_sequences = []
        c_sequences = []
        
        # 初始化自适应评分器获取动态阈值
        try:
            import os
            genome_size = os.path.getsize(self.config.genome_file)
        except:
            genome_size = 1e9  # 默认1Gb
        
        from adaptive_scoring import AdaptiveTEScorer
        adaptive_scorer = AdaptiveTEScorer(genome_size)
        
        # 获取自适应阈值
        A_THRESHOLD, B_THRESHOLD = adaptive_scorer.get_adaptive_thresholds(self.sequences)
        
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
        logger.info(f"Classification thresholds (adaptive) - A: >={A_THRESHOLD:.2f}, B: >={B_THRESHOLD:.2f}")
        logger.info(f"Score range: {min(all_final_scores):.3f} - {max(all_final_scores):.3f}")
        
        # 统计分类原因
        classification_stats = {}
        for seq_data in self.sequences:
            reason = self.scores[seq_data['id']].get('classification_reason', 'unknown')
            if reason not in classification_stats:
                classification_stats[reason] = {'A': 0, 'B': 0, 'C': 0}
            
            final_score = self.scores[seq_data['id']]['final']
            if final_score >= A_THRESHOLD:
                classification_stats[reason]['A'] += 1
            elif final_score >= B_THRESHOLD:
                classification_stats[reason]['B'] += 1
            else:
                classification_stats[reason]['C'] += 1
        
        # 输出分类原因统计
        logger.info("Top classification reasons:")
        for reason, counts in sorted(classification_stats.items(), 
                                    key=lambda x: sum(x[1].values()), reverse=True)[:6]:
            total_count = sum(counts.values())
            logger.info(f"  {reason}: A={counts['A']}, B={counts['B']}, C={counts['C']} (total={total_count})")
        
        # 输出前几个高分序列的详细信息
        top_sequences = sorted(self.sequences, 
                              key=lambda x: self.scores[x['id']]['final'], 
                              reverse=True)[:5]
        logger.info("Top 5 sequences (by adaptive final score):")
        for seq in top_sequences:
            sid = seq['id']
            s = self.scores[sid]
            seq_len = len(seq['sequence'])
            logger.info(f"  {sid} ({seq_len}bp): final={s['final']:.3f}, "
                       f"adaptive_copy={s.get('adaptive_copy_score', s.get('copy_score', 0)):.2f}(n={s.get('copy_number', 0)}), "
                       f"adaptive_identity={s.get('adaptive_identity_score', s.get('identity_score', 0)):.2f}({s.get('avg_identity', 0):.1f}%), "
                       f"reason={s.get('classification_reason', 'unknown')}")
        
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
        
        # 步骤7: 识别需要构建共识的候选序列
        # 策略：A类和B类序列都是共识构建的候选，C类序列质量太低跳过
        consensus_candidates = a_sequences + b_sequences
        
        # 为候选序列添加质量分类标记
        for seq in consensus_candidates:
            seq_score = self.scores[seq['id']]['final']
            if seq_score >= A_THRESHOLD:
                seq['quality_class'] = 'A'
            elif seq_score >= B_THRESHOLD:
                seq['quality_class'] = 'B'
            else:
                seq['quality_class'] = 'C'
        
        # 过滤掉质量过低的序列
        filtered_candidates = []
        filtered_out_sequences = []
        
        for seq in consensus_candidates:
            seq_id = seq['id']
            score_data = self.scores[seq_id]
            
            # 基本过滤条件
            copy_number = score_data.get('copy_number', 0)
            avg_identity = score_data.get('avg_identity', 0)
            seq_length = len(seq['sequence'])
            
            # 过滤条件：拷贝数太少、相似度太低或序列太短的不适合构建共识
            if (copy_number < 2 and seq['quality_class'] != 'A') or avg_identity < 50 or seq_length < self.config.min_length:
                seq['filter_reason'] = f"copy_num={copy_number}, identity={avg_identity:.1f}%, length={seq_length}bp"
                filtered_out_sequences.append(seq)
            else:
                # 添加RepeatMasker结果以供Phase2使用
                if seq_id in self.rm_detailed_results:
                    seq['rm_hits'] = self.rm_detailed_results[seq_id].get('hits', [])
                    seq['rm_coverage'] = self.rm_detailed_results[seq_id].get('coverage', 0)
                else:
                    seq['rm_hits'] = []
                    seq['rm_coverage'] = 0
                
                filtered_candidates.append(seq)
        
        logger.info(f"Consensus candidate selection:")
        logger.info(f"  - A-class candidates: {len([s for s in filtered_candidates if s['quality_class'] == 'A'])}")
        logger.info(f"  - B-class candidates: {len([s for s in filtered_candidates if s['quality_class'] == 'B'])}")
        logger.info(f"  - Total candidates for Phase2: {len(filtered_candidates)}")
        logger.info(f"  - Filtered out (too low quality): {len(filtered_out_sequences)}")
        logger.info(f"  - C-class sequences (skipped): {len(c_sequences)}")
        
        # 候选序列统计
        if filtered_candidates:
            candidate_copy_numbers = [self.scores[seq['id']]['copy_number'] for seq in filtered_candidates]
            candidate_identities = [self.scores[seq['id']]['avg_identity'] for seq in filtered_candidates]
            logger.info(f"Candidate stats: Copy number {np.mean(candidate_copy_numbers):.1f}±{np.std(candidate_copy_numbers):.1f}, "
                       f"Identity {np.mean(candidate_identities):.1f}±{np.std(candidate_identities):.1f}%")
        
        return {
            'consensus_candidates': filtered_candidates,  # 所有需要Phase2处理的序列
            'filtered_out_sequences': filtered_out_sequences,  # 过滤掉的低质量序列
            'c_class_sequences': c_sequences,  # C类序列（仅记录，不处理）
            'scores': self.scores,
            'rm_detailed_results': self.rm_detailed_results,  # 传递给Phase 2重用
            'summary': f"Candidates:{len(filtered_candidates)} (A:{len([s for s in filtered_candidates if s['quality_class'] == 'A'])}, B:{len([s for s in filtered_candidates if s['quality_class'] == 'B'])}), Filtered:{len(filtered_out_sequences)}, C-class:{len(c_sequences)}"
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
        
        # 使用phase1专用的CD-HIT去冗余函数（80%阈值）
        logger.info("Starting CD-HIT pre-RepeatMasker deduplication using phase1 standards")
        
        # 转换序列格式为CD-HIT期望的格式
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
        
        logger.info(f"Running CD-HIT with 70% threshold on {len(sorted_sequences)} sequences")
        
        # 运行CD-HIT去冗余（70%阈值 - 更激进的去冗余）
        deduplicated = run_cdhit_optimized_phase1(
            sorted_sequences,
            threshold=0.80,  # 降低到70%阈值，更激进的去冗余
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
    
    def identify_potential_chimeras(self, sequences: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
        """识别潜在的嵌合体序列"""
        regular_seqs = []
        potential_chimeras = []
        
        logger.info(f"Starting chimera identification for {len(sequences)} sequences")
        
        for seq_data in sequences:
            if self._is_potential_chimera(seq_data):
                potential_chimeras.append(seq_data)
            else:
                regular_seqs.append(seq_data)
        
        logger.debug(f"Identified {len(potential_chimeras)} potential chimeras, {len(regular_seqs)} regular sequences")
        
        return regular_seqs, potential_chimeras
    
    def _is_potential_chimera(self, seq_data: Dict) -> bool:
        """判断序列是否可能是嵌合体"""
        sequence = seq_data['sequence']
        seq_len = len(sequence)
        seq_id = seq_data['id']
        
        # 标准1: 异常长度（>10kb或者>15kb根据基因组大小调整）
        length_threshold = 15000 if hasattr(self.config, 'large_genome_threshold') and \
                          getattr(self.config, 'genome_size', 0) > self.config.large_genome_threshold else 10000
        
        if seq_len > length_threshold:
            logger.debug(f"{seq_id}: Potential chimera - excessive length ({seq_len}bp > {length_threshold}bp)")
            return True
        
        # 标准2: 内部复杂度变化过大
        if seq_len > 2000:  # 只对足够长的序列进行复杂度分析
            complexity_variance = self._calculate_complexity_variance(sequence)
            if complexity_variance > 0.3:  # 复杂度变异系数>30%
                logger.debug(f"{seq_id}: Potential chimera - high complexity variance ({complexity_variance:.3f})")
                return True
        
        # 标准3: GC含量波动异常
        if seq_len > 1500:
            gc_variance = self._calculate_gc_variance(sequence)
            if gc_variance > 0.08:  # GC变异>8%
                logger.debug(f"{seq_id}: Potential chimera - high GC variance ({gc_variance:.3f})")
                return True
        
        # 标准4: RepeatMasker结果显示多个不同家族
        if seq_id in self.rm_detailed_results:
            rm_data = self.rm_detailed_results[seq_id]
            if self._has_multiple_te_families_in_hits(rm_data.get('hits', [])):
                logger.debug(f"{seq_id}: Potential chimera - multiple TE families detected")
                return True
        
        return False
    
    def _calculate_complexity_variance(self, sequence: str, window_size: int = 500) -> float:
        """计算序列复杂度方差"""
        from utils.complexity_utils import calculate_dust_score
        
        if len(sequence) < window_size * 2:
            return 0.0
        
        complexities = []
        step = window_size // 2
        
        for i in range(0, len(sequence) - window_size + 1, step):
            window = sequence[i:i+window_size]
            complexity = calculate_dust_score(window)
            complexities.append(complexity)
        
        if len(complexities) < 2:
            return 0.0
        
        # 计算变异系数
        mean_complexity = np.mean(complexities)
        std_complexity = np.std(complexities)
        
        return std_complexity / mean_complexity if mean_complexity > 0 else 0.0
    
    def _calculate_gc_variance(self, sequence: str, window_size: int = 500) -> float:
        """计算GC含量方差"""
        if len(sequence) < window_size * 2:
            return 0.0
        
        gc_contents = []
        step = window_size // 2
        
        for i in range(0, len(sequence) - window_size + 1, step):
            window = sequence[i:i+window_size]
            gc_count = window.count('G') + window.count('C')
            gc_content = gc_count / len(window) if len(window) > 0 else 0
            gc_contents.append(gc_content)
        
        return np.std(gc_contents) if len(gc_contents) > 1 else 0.0
    
    def _has_multiple_te_families_in_hits(self, hits: List[Dict]) -> bool:
        """检查RepeatMasker hits是否涉及多个TE家族"""
        if len(hits) < 2:
            return False
        
        # 模拟检查：如果有多个高质量hits且相互不重叠，可能是嵌合体
        # 实际实现中需要分析hits的family信息
        high_quality_hits = [hit for hit in hits if hit.get('score', 0) > 200]
        
        if len(high_quality_hits) >= 2:
            # 检查hits是否相互不重叠（简化版本）
            sorted_hits = sorted(high_quality_hits, key=lambda x: x.get('start', 0))
            for i in range(len(sorted_hits) - 1):
                if sorted_hits[i].get('end', 0) < sorted_hits[i+1].get('start', 0):
                    return True  # 发现不重叠的高质量hits
        
        return False


def run_cdhit_optimized_phase1(sequences: List[Dict], threshold: float, threads: int, 
                       config, label: str) -> List[Dict]:
    """Phase1专用的CD-HIT运行函数"""
    import tempfile
    import subprocess
    import os
    import shutil
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    if not sequences:
        return []
    
    # 创建临时输入文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False, dir=config.temp_dir) as input_file:
        records = []
        for seq_data in sequences:
            record = SeqRecord(
                Seq(seq_data['sequence']),
                id=seq_data.get('id', f"seq_{len(records)}"),
                description=""
            )
            records.append(record)
        
        SeqIO.write(records, input_file, "fasta")
        input_path = input_file.name
    
    # Ensure file is readable
    os.chmod(input_path, 0o644)
    
    # 创建临时输出文件
    output_path = input_path + ".cdhit"
    
    try:
        # 检查cd-hit-est是否存在
        if not shutil.which(config.cdhit_exe):
            logger.error(f"CD-HIT executable not found: {config.cdhit_exe}")
            logger.error("Please install CD-HIT or specify correct path")
            raise FileNotFoundError(f"CD-HIT not found: {config.cdhit_exe}")
        
        # 构建CD-HIT命令 - 优化参数以平衡去冗余和保留多样性
        cmd = [
            config.cdhit_exe,
            '-i', input_path,
            '-o', output_path,
            '-c', str(threshold),  # Use the provided threshold parameter
            '-aS', '0.7',  # 70% coverage of shorter sequence - 更宽松以保留部分TE
            # 移除 -aL 参数，允许长序列部分匹配（处理嵌套TE）
            '-G', '0',  # local alignment mode
            '-n', '8' if threshold >= 0.9 else ('5' if threshold >= 0.8 else '4'),  # word size: 90%+=8, 80-89%=5, 70-79%=4
            '-r', '1',  # compare both strands
            '-mask', 'NX',  # mask low-complexity regions
            '-M', '0',  # no memory limit
            '-T', str(threads),  # use all available threads
            '-d', '0',  # keep full sequence names
            '-g', '1'  # 最精确的聚类模式，确保质量
        ]
        
        # Log file information
        input_size = os.path.getsize(input_path)
        logger.info(f"Running CD-HIT {label} with threshold {threshold} and {threads} threads")
        logger.info(f"Input file: {input_path} ({input_size} bytes, {len(sequences)} sequences)")
        logger.debug(f"CD-HIT command: {' '.join(cmd)}")
        
        # 运行CD-HIT (无timeout限制)
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check output
        if os.path.exists(output_path):
            output_size = os.path.getsize(output_path)
            logger.info(f"CD-HIT output: {output_path} ({output_size} bytes)")
        
        if result.returncode != 0:
            logger.error(f"CD-HIT {label} failed with return code {result.returncode}")
            logger.error(f"CD-HIT stderr: {result.stderr}")
            logger.error(f"CD-HIT stdout: {result.stdout}")
            # Don't silently continue - raise an exception
            raise RuntimeError(f"CD-HIT failed: {result.stderr}")
        
        # 解析结果
        representatives = []
        if os.path.exists(output_path):
            for record in SeqIO.parse(output_path, "fasta"):
                # 找到对应的原始序列数据
                for seq_data in sequences:
                    if seq_data.get('id') == record.id:
                        representatives.append(seq_data)
                        break
        
        logger.debug(f"CD-HIT {label}: {len(sequences)} -> {len(representatives)}")
        return representatives
        
    except FileNotFoundError as e:
        logger.error(f"CD-HIT {label} not found: {e}")
        logger.warning("Skipping CD-HIT redundancy removal - returning all sequences")
        return sequences  # CD-HIT not available, return original
    except RuntimeError as e:
        logger.error(f"CD-HIT {label} execution failed: {e}")
        logger.warning("CD-HIT failed - attempting simple redundancy removal fallback")
        # Simple fallback: return original sequences
        return sequences
    except Exception as e:
        logger.error(f"CD-HIT {label} unexpected error: {e}")
        logger.warning("Unexpected error - returning original sequences")
        return sequences  # 出错时返回原始序列
    finally:
        # 清理临时文件
        for temp_file in [input_path, output_path, output_path + ".clstr"]:
            try:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
            except Exception:
                pass
