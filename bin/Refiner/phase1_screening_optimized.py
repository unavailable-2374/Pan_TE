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
        import re
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

            # CRITICAL FIX: Use full unique ID to distinguish sequences from different RepeatScout runs
            # record.id only returns part before first space (e.g., "R=0")
            # record.description contains full header with run suffix (e.g., "R=0 (RR=1...)_l20_t3_set0")
            # Extract unique ID: "R=X_suffix" format
            full_desc = record.description
            suffix_match = re.search(r'_(l\d+_t\d+_set\d+)$', full_desc)
            if suffix_match:
                unique_id = f"{record.id}_{suffix_match.group(1)}"
            else:
                unique_id = record.id  # Fallback if no suffix found

            sequences.append({
                'id': unique_id,
                'original_id': record.id,  # Keep original R= ID for reference
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
        except (OSError, AttributeError) as e:
            logger.debug(f"Could not get genome size: {e}")
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

            # Processing low copy number sequences
            # Key insight: RepeatScout already filtered sequences with <10 k-mer occurrences
            # Low RepeatMasker copy count may indicate:
            # 1. Ancient/divergent TE (RM can't detect all copies)
            # 2. Young TE with few copies
            # 3. RM parameters not sensitive enough
            # Therefore, we should NOT use copy number as a veto mechanism
            if copy_number < 2:
                # Single copy: give base score, let biological features decide
                logger.debug(f"{seq_id}: Single copy - using base score with biological validation")
                copy_score = 0.35  # Base score for single copy (not a death sentence)
                identity_score = adaptive_scorer.calculate_adaptive_identity_score(
                    avg_identity, seq_length, copy_number
                )
            elif copy_number < 3:
                # 2 copies: moderate score, RepeatScout found it so it's likely real
                logger.debug(f"{seq_id}: Low copy number ({copy_number}) - trusting RepeatScout detection")
                copy_score = 0.45  # Reasonable base score
                identity_score = adaptive_scorer.calculate_adaptive_identity_score(
                    avg_identity, seq_length, copy_number
                )
            else:
                # 3+ copies: use adaptive scoring
                copy_score = adaptive_scorer.calculate_adaptive_copy_score(copy_number, seq_length)
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

        # Populate default scores for sequences not in RM results
        # (e.g., RepeatMasker timed out or returned partial results)
        missing_count = 0
        for seq_data in self.sequences:
            seq_id = seq_data['id']
            if seq_id not in self.scores:
                self.scores[seq_id] = {}
            if 'copy_number' not in self.scores[seq_id]:
                missing_count += 1
                # RepeatScout already validated these as repetitive — trust that
                self.scores[seq_id]['copy_score'] = 0.35
                self.scores[seq_id]['identity_score'] = 0.3
                self.scores[seq_id]['copy_number'] = 0
                self.scores[seq_id]['avg_identity'] = 0
        if missing_count > 0:
            logger.warning(f"RepeatMasker returned no results for {missing_count}/{len(self.sequences)} sequences — using default scores")

        # 输出拷贝数和一致性统计
        copy_numbers = [self.scores[seq_id]['copy_number'] for seq_id in self.scores]
        identities = [self.scores[seq_id]['avg_identity'] for seq_id in self.scores]
        
        # Statistics by copy number tiers (no longer treating low-copy as noise)
        single_copy_seqs = [seq_id for seq_id in self.scores if self.scores[seq_id]['copy_number'] < 2]
        low_copy_seqs = [seq_id for seq_id in self.scores if 2 <= self.scores[seq_id]['copy_number'] < 5]
        medium_copy_seqs = [seq_id for seq_id in self.scores if 5 <= self.scores[seq_id]['copy_number'] < 20]
        high_copy_seqs = [seq_id for seq_id in self.scores if self.scores[seq_id]['copy_number'] >= 20]

        logger.info(f"Copy number distribution (RepeatMasker detection):")
        logger.info(f"  Single copy (<2): {len(single_copy_seqs)} - potential young/divergent TEs")
        logger.info(f"  Low copy (2-4): {len(low_copy_seqs)} - may need biological validation")
        logger.info(f"  Medium copy (5-19): {len(medium_copy_seqs)} - typical TEs")
        logger.info(f"  High copy (≥20): {len(high_copy_seqs)} - well-supported TEs")
        
        if copy_numbers:
            logger.info(f"Copy numbers - Min: {min(copy_numbers)}, Max: {max(copy_numbers)}, "
                       f"Mean: {np.mean(copy_numbers):.1f}")
            logger.info(f"Identities - Min: {min(identities):.1f}%, Max: {max(identities):.1f}%, "
                       f"Mean: {np.mean(identities):.1f}%")
    
    def calculate_final_scores(self):
        """计算自适应综合分数"""
        logger.info("Calculating adaptive final scores based on sequence characteristics...")
        
        # 初始化自适应评分器
        try:
            import os
            genome_size = os.path.getsize(self.config.genome_file)
        except (OSError, AttributeError) as e:
            logger.debug(f"Could not get genome size: {e}")
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
        except (OSError, AttributeError) as e:
            logger.debug(f"Could not get genome size: {e}")
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
        # Strategy change: Rescue C-class sequences that have biological evidence
        # Rationale: RepeatScout already filtered by copy number (≥10 k-mer occurrences)
        # C-class may contain real TEs that:
        # - Are ancient/divergent (low RM detection)
        # - Have structural features (TIR, LTR, poly-A)
        # - Are long enough to be functional

        rescued_c_sequences = []
        unrescued_c_sequences = []

        for seq in c_sequences:
            seq_id = seq['id']
            score_data = self.scores[seq_id]
            sequence = seq['sequence']
            seq_length = len(sequence)
            copy_number = score_data.get('copy_number', 0)

            # C-class rescue conditions - AGGRESSIVE STRATEGY
            # Rationale: RepeatScout already requires ≥10 k-mer occurrences
            # These sequences ARE repetitive, just not detected well by RepeatMasker
            # Default: RESCUE unless obviously problematic

            should_rescue = True  # Default to rescue
            rescue_reason = ["repeatscout_trusted"]
            should_reject = False
            reject_reason = []

            complexity = score_data.get('complexity', 0)

            # Only reject if BOTH conditions are met:
            # 1. Very short (<80bp) - might be artifact
            # 2. Low complexity (<0.3) - might be simple repeat
            if seq_length < 80 and complexity < 0.3:
                should_reject = True
                reject_reason.append(f"short_low_complexity({seq_length}bp,{complexity:.2f})")

            # Additional rescue evidence (for logging)
            if seq_length >= 200:
                rescue_reason.append(f"moderate_length({seq_length}bp)")
            if complexity >= 0.4:
                rescue_reason.append(f"reasonable_complexity({complexity:.2f})")
            if copy_number >= 1:
                rescue_reason.append(f"has_copies({copy_number})")
            if self._has_te_structural_features(sequence):
                rescue_reason.append("structural_features")
            if seq_length >= 500:
                rescue_reason.append(f"long_seq({seq_length}bp)")

            # Final decision
            should_rescue = should_rescue and not should_reject

            if should_rescue:
                seq['quality_class'] = 'C_rescued'
                seq['rescue_reason'] = ','.join(rescue_reason)
                rescued_c_sequences.append(seq)
                logger.debug(f"Rescued C-class sequence {seq_id}: {seq['rescue_reason']}")
            else:
                unrescued_c_sequences.append(seq)

        if rescued_c_sequences:
            logger.info(f"C-class rescue: {len(rescued_c_sequences)} sequences rescued, "
                       f"{len(unrescued_c_sequences)} remain filtered")

        # Combine A, B, and rescued C sequences
        consensus_candidates = a_sequences + b_sequences + rescued_c_sequences

        # 为候选序列添加质量分类标记
        for seq in consensus_candidates:
            seq_score = self.scores[seq['id']]['final']
            if seq_score >= A_THRESHOLD:
                seq['quality_class'] = 'A'
            elif seq_score >= B_THRESHOLD:
                seq['quality_class'] = 'B'
            else:
                seq['quality_class'] = 'C'
        
        # Filter sequences - with relaxed criteria to improve sensitivity
        filtered_candidates = []
        filtered_out_sequences = []

        # Relaxed identity threshold (from 50 to 40)
        # Rationale: Ancient TEs can have <50% identity but still be real
        MIN_IDENTITY_THRESHOLD = 40.0

        for seq in consensus_candidates:
            seq_id = seq['id']
            score_data = self.scores[seq_id]
            sequence = seq['sequence']

            copy_number = score_data.get('copy_number', 0)
            avg_identity = score_data.get('avg_identity', 0)
            seq_length = len(sequence)
            quality_class = seq.get('quality_class', 'C')

            # Determine if sequence should pass filter
            should_filter = False
            filter_reason = []

            # Check minimum length
            if seq_length < self.config.min_length:
                should_filter = True
                filter_reason.append(f"too_short({seq_length}bp)")

            # Check identity - but with biological evidence override
            # RELAXED: Trust RepeatScout - if it found the sequence with ≥10 k-mers, it's likely real
            if avg_identity < MIN_IDENTITY_THRESHOLD:
                # Low identity is acceptable if:
                # 1. Sequence is long (>500bp) - likely real TE even if divergent
                # 2. Has structural features
                # 3. Is A-class (high score from other factors)
                # 4. Is C_rescued (already passed RepeatScout's filter)
                has_override = (
                    seq_length >= 500 or
                    quality_class in ['A', 'C_rescued'] or
                    self._has_te_structural_features(sequence)
                )
                if not has_override:
                    should_filter = True
                    filter_reason.append(f"low_identity({avg_identity:.1f}%)")

            # Check copy number - VERY lenient
            # Trust RepeatScout: only filter extremely short sequences with no evidence
            # C_rescued are always allowed (RepeatScout already validated them)
            if copy_number < 1 and quality_class not in ['A', 'B', 'C_rescued'] and seq_length < 200:
                should_filter = True
                filter_reason.append(f"no_copies_and_short")

            if should_filter:
                seq['filter_reason'] = ','.join(filter_reason)
                filtered_out_sequences.append(seq)
            else:
                # Add RepeatMasker results for Phase2
                if seq_id in self.rm_detailed_results:
                    seq['rm_hits'] = self.rm_detailed_results[seq_id].get('hits', [])
                    seq['rm_coverage'] = self.rm_detailed_results[seq_id].get('coverage', 0)
                else:
                    seq['rm_hits'] = []
                    seq['rm_coverage'] = 0

                filtered_candidates.append(seq)

        # Updated logging to reflect new categories
        a_count = len([s for s in filtered_candidates if s.get('quality_class') == 'A'])
        b_count = len([s for s in filtered_candidates if s.get('quality_class') == 'B'])
        c_rescued_count = len([s for s in filtered_candidates if s.get('quality_class') == 'C_rescued'])

        logger.info(f"Consensus candidate selection:")
        logger.info(f"  - A-class candidates: {a_count}")
        logger.info(f"  - B-class candidates: {b_count}")
        logger.info(f"  - C-class rescued: {c_rescued_count}")
        logger.info(f"  - Total candidates for Phase2: {len(filtered_candidates)}")
        logger.info(f"  - Filtered out: {len(filtered_out_sequences)}")
        logger.info(f"  - C-class not rescued: {len(unrescued_c_sequences)}")
        
        # 候选序列统计
        if filtered_candidates:
            candidate_copy_numbers = [self.scores[seq['id']]['copy_number'] for seq in filtered_candidates]
            candidate_identities = [self.scores[seq['id']]['avg_identity'] for seq in filtered_candidates]
            logger.info(f"Candidate stats: Copy number {np.mean(candidate_copy_numbers):.1f}±{np.std(candidate_copy_numbers):.1f}, "
                       f"Identity {np.mean(candidate_identities):.1f}±{np.std(candidate_identities):.1f}%")
        
        return {
            'consensus_candidates': filtered_candidates,  # All sequences for Phase2
            'filtered_out_sequences': filtered_out_sequences,  # Filtered low-quality sequences
            'c_class_sequences': unrescued_c_sequences,  # C-class not rescued (for reference)
            'rescued_c_sequences': rescued_c_sequences,  # C-class rescued sequences
            'scores': self.scores,
            'rm_detailed_results': self.rm_detailed_results,  # Pass to Phase 2 for reuse
            'summary': f"Candidates:{len(filtered_candidates)} (A:{a_count}, B:{b_count}, C_rescued:{c_rescued_count}), Filtered:{len(filtered_out_sequences)}, C_unrescued:{len(unrescued_c_sequences)}"
        }

    def _has_te_structural_features(self, sequence: str) -> bool:
        """
        Check if sequence has TE structural features.

        Features checked:
        1. Terminal Inverted Repeats (TIR) - DNA transposon signature
        2. Long Terminal Repeats (LTR) - Retrotransposon signature
        3. Poly-A tail - LINE/SINE signature
        4. Target Site Duplication patterns

        Returns:
            True if any structural feature is detected
        """
        if not sequence or len(sequence) < 50:
            return False

        sequence = sequence.upper()
        seq_len = len(sequence)

        # 1. Check for Terminal Inverted Repeats (TIR)
        # TIR typically 10-40bp at both ends, allowing 2 mismatches
        if self._has_terminal_inverted_repeat(sequence, min_len=10, max_mismatch=2):
            return True

        # 2. Check for Long Terminal Repeats (LTR)
        # LTR typically 100-500bp, >80% identity
        if seq_len >= 500 and self._has_long_terminal_repeat(sequence, min_len=80, min_identity=0.75):
            return True

        # 3. Check for poly-A tail (LINE/SINE signature)
        # Look for poly-A (≥6 A's) at 3' end
        tail_region = sequence[-50:] if seq_len >= 50 else sequence
        if 'AAAAAA' in tail_region:
            return True

        # 4. Check for common TSD patterns at boundaries
        # Many TEs have characteristic boundary sequences
        if self._has_characteristic_boundaries(sequence):
            return True

        return False

    def _has_terminal_inverted_repeat(self, sequence: str, min_len: int = 10, max_mismatch: int = 2) -> bool:
        """Check for Terminal Inverted Repeats (TIR)"""
        seq_len = len(sequence)
        if seq_len < min_len * 2:
            return False

        # Check TIR lengths from 10 to 40bp
        for tir_len in range(min_len, min(41, seq_len // 4)):
            left_tir = sequence[:tir_len]
            right_tir = sequence[-tir_len:]

            # Reverse complement of right TIR
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
            right_tir_rc = ''.join(complement.get(b, 'N') for b in reversed(right_tir))

            # Count mismatches
            mismatches = sum(1 for a, b in zip(left_tir, right_tir_rc) if a != b)

            if mismatches <= max_mismatch:
                return True

        return False

    def _has_long_terminal_repeat(self, sequence: str, min_len: int = 80, min_identity: float = 0.75) -> bool:
        """Check for Long Terminal Repeats (LTR)"""
        seq_len = len(sequence)
        if seq_len < min_len * 2 + 100:  # Need space for internal region
            return False

        # Check LTR lengths from 80 to 500bp
        for ltr_len in range(min_len, min(501, seq_len // 3)):
            left_ltr = sequence[:ltr_len]
            right_ltr = sequence[-ltr_len:]

            # Calculate identity
            matches = sum(1 for a, b in zip(left_ltr, right_ltr) if a == b)
            identity = matches / ltr_len

            if identity >= min_identity:
                return True

        return False

    def _has_characteristic_boundaries(self, sequence: str) -> bool:
        """Check for characteristic TE boundary patterns"""
        if len(sequence) < 20:
            return False

        # Common TE boundary patterns
        # TA dinucleotide (common TSD for many DNA transposons)
        if sequence[:2] == 'TA' and sequence[-2:] == 'TA':
            return True

        # TG...CA pattern (LTR retrotransposons often have TG...CA)
        if sequence[:2] == 'TG' and sequence[-2:] == 'CA':
            return True

        # Check for simple direct repeats at boundaries (potential TSD)
        left_10bp = sequence[:10]
        right_10bp = sequence[-10:]

        # Look for 4-8bp direct repeats
        for repeat_len in range(4, 9):
            left_motif = sequence[:repeat_len]
            # Check if same motif appears near the end
            if left_motif in sequence[-20:]:
                return True

        return False

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
        
        logger.info(f"Running CD-HIT with 80% identity threshold on {len(sorted_sequences)} sequences")

        # 运行CD-HIT去冗余（80%相似度阈值）
        deduplicated = run_cdhit_optimized_phase1(
            sorted_sequences,
            threshold=0.80,  # 80% sequence identity threshold
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
        
        # 运行CD-HIT
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        except subprocess.TimeoutExpired:
            logger.error(f"CD-HIT {label} timed out after 30 minutes")
            raise RuntimeError(f"CD-HIT {label} timed out")
        
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
