"""
Phase 3 热修复版本 - 充分利用多线程
主要改进：
1. 充分利用用户设置的线程数
2. 并行化去冗余处理
3. 优化CD-HIT并发策略
"""

import logging
import gc
from typing import List, Dict, Any, Optional, Tuple
import numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from functools import partial
import multiprocessing as mp

from config import PipelineConfig
from utils.robust_runner import RobustRunner

logger = logging.getLogger(__name__)

class ConsensusFilterRelaxed:
    """Phase 3: 优化的质量控制与输出生成"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.runner = RobustRunner(config)
        
        # 充分利用线程数，不受CPU核心数限制
        cpu_count = mp.cpu_count()
        
        # 为Phase 3分配合理的线程数
        if config.threads >= 32:
            self.n_threads = min(config.threads, max(cpu_count, 32))
        elif config.threads >= 16:
            self.n_threads = min(config.threads, max(cpu_count, 16))
        else:
            self.n_threads = max(4, config.threads)
        
        # CD-HIT并发策略
        self.cdhit_workers = min(8, max(2, self.n_threads // 4))  # 2-8个CD-HIT进程
        self.threads_per_cdhit = max(1, self.n_threads // self.cdhit_workers)
        
        # 动态调整嵌合检测参数 - 修改为50bp最小长度
        self.chimera_distance_threshold = getattr(config, 'chimera_segment_distance', 20000)
        self.chimera_min_length = 50  # 统一使用50bp最小长度
        self.chimera_strict_mode = getattr(config, 'chimera_detection_strict', False)
        
        logger.info(f"Phase 3 optimized: {self.n_threads} threads total, "
                   f"{self.cdhit_workers} CD-HIT workers with {self.threads_per_cdhit} threads each")
        logger.info(f"Chimera detection: distance={self.chimera_distance_threshold}, "
                   f"min_length={self.chimera_min_length}, strict={self.chimera_strict_mode}")
        
    def run(self, consensus_list: List[Dict]) -> Dict[str, Any]:
        """执行Phase 3完整流程"""
        logger.info(f"Phase 3 starting with {len(consensus_list)} consensus sequences")
        
        # 步骤1: 嵌合检测（宽松模式）
        if self.chimera_strict_mode:
            logger.info("Running strict chimera detection")
            non_chimeric = self.detect_and_remove_chimeras_parallel(consensus_list)
        else:
            logger.info("Skipping chimera detection (relaxed mode)")
            non_chimeric = consensus_list
        
        logger.info(f"After chimera filtering: {len(non_chimeric)} sequences")
        
        # 步骤2: 去冗余 - 并行处理两个阈值
        logger.info("Running parallel redundancy removal")
        masking_lib, analysis_lib = self.remove_redundancy_parallel_optimized(
            non_chimeric,
            self.config.redundancy_threshold_masking,
            self.config.redundancy_threshold_analysis
        )
        
        # 步骤3: 生成输出文件
        self.write_consensus_libraries(masking_lib, analysis_lib)
        
        # 步骤4: 生成统计报告
        stats = self.generate_statistics(consensus_list, masking_lib, analysis_lib)
        
        return {
            'masking_library': masking_lib,
            'analysis_library': analysis_lib,
            'statistics': stats,
            'summary': f"Final output: {len(masking_lib)} masking, {len(analysis_lib)} analysis"
        }
    
    def detect_and_remove_chimeras_parallel(self, consensus_list: List[Dict]) -> List[Dict]:
        """并行嵌合检测 - 简化版本"""
        logger.info(f"Running simple chimera detection with {self.n_threads} threads")
        
        # 简化的嵌合检测：基于长度和复杂度过滤
        filtered = []
        for seq_data in consensus_list:
            sequence = seq_data.get('sequence', '')
            
            # 基本过滤条件 - 修改为50bp最小长度
            if (len(sequence) >= 50 and  # 统一使用50bp最小长度
                len(sequence) <= 100000 and  # 增大最大长度限制到100kb
                sequence.count('N') / len(sequence) < 0.3):  # N含量过滤
                filtered.append(seq_data)
        
        logger.info(f"Chimera detection: {len(consensus_list)} -> {len(filtered)} sequences")
        return filtered
    
    def adaptive_filter_by_copy_number(self, consensus_list: List[Dict]) -> List[Dict]:
        """自适应拷贝数过滤"""
        if not consensus_list:
            return []
        
        copy_numbers = [seq.get('num_copies', 1) for seq in consensus_list]
        
        # 使用分位数动态确定阈值
        p50 = np.percentile(copy_numbers, 50)
        p25 = np.percentile(copy_numbers, 25)
        
        # 固定阈值：最低5个拷贝要求，不再使用动态阈值
        dynamic_threshold = 5  # 固定为5个拷贝
        
        logger.info(f"Copy number statistics: P25={p25:.1f}, P50={p50:.1f}, "
                   f"threshold={dynamic_threshold}")
        
        filtered = []
        with ProcessPoolExecutor(max_workers=self.n_threads) as executor:
            batch_size = max(10, len(consensus_list) // self.n_threads)
            futures = []
            
            for i in range(0, len(consensus_list), batch_size):
                batch = consensus_list[i:i+batch_size]
                future = executor.submit(
                    filter_batch_by_copy_number,
                    batch,
                    dynamic_threshold
                )
                futures.append(future)
            
            for future in as_completed(futures):
                try:
                    batch_result = future.result()
                    filtered.extend(batch_result)
                except Exception as e:
                    logger.error(f"Copy number filtering batch failed: {e}")
        
        logger.info(f"Copy number filtering: {len(consensus_list)} -> {len(filtered)}")
        return filtered
    
    def remove_redundancy_parallel_optimized(self, consensus_list: List[Dict], 
                                           threshold1: float, threshold2: float) -> Tuple[List[Dict], List[Dict]]:
        """优化的并行去冗余处理"""
        logger.info(f"Running redundancy removal with thresholds {threshold1:.2f} and {threshold2:.2f}")
        
        if not consensus_list:
            return [], []
        
        # 按质量/拷贝数排序，确保高质量序列优先保留
        sorted_list = sorted(
            consensus_list,
            key=lambda x: (x.get('num_copies', 1), len(x.get('sequence', ''))),
            reverse=True
        )
        
        # 策略选择：小数据集用串行CD-HIT，大数据集用分批并行
        if len(sorted_list) <= 1000:
            # 小数据集：直接并行运行两个CD-HIT
            logger.info("Small dataset: running 2 parallel CD-HIT processes")
            with ProcessPoolExecutor(max_workers=2) as executor:
                future_masking = executor.submit(
                    run_cdhit_optimized,
                    sorted_list,
                    threshold1,
                    self.threads_per_cdhit,
                    self.config,
                    "masking"
                )
                
                future_analysis = executor.submit(
                    run_cdhit_optimized,
                    sorted_list,
                    threshold2,
                    self.threads_per_cdhit,
                    self.config,
                    "analysis"
                )
                
                masking_lib = future_masking.result()
                analysis_lib = future_analysis.result()
        else:
            # 大数据集：分批并行处理
            logger.info(f"Large dataset: running batch-parallel CD-HIT with {self.cdhit_workers} workers")
            masking_lib = self.run_cdhit_batch_parallel(sorted_list, threshold1, "masking")
            analysis_lib = self.run_cdhit_batch_parallel(sorted_list, threshold2, "analysis")
        
        logger.info(f"Redundancy removal complete: masking={len(masking_lib)}, analysis={len(analysis_lib)}")
        return masking_lib, analysis_lib
    
    def run_cdhit_batch_parallel(self, sequences: List[Dict], threshold: float, 
                                lib_type: str) -> List[Dict]:
        """分批并行CD-HIT处理"""
        if not sequences:
            return []
        
        # 将序列分批
        batch_size = max(100, len(sequences) // self.cdhit_workers)
        batches = []
        for i in range(0, len(sequences), batch_size):
            batches.append(sequences[i:i+batch_size])
        
        logger.info(f"Running {lib_type} CD-HIT: {len(batches)} batches with {self.cdhit_workers} workers")
        
        # 并行处理各批次
        all_representatives = []
        try:
            with ProcessPoolExecutor(max_workers=self.cdhit_workers) as executor:
                future_to_batch = {}
                for i, batch in enumerate(batches):
                    future = executor.submit(
                        run_cdhit_optimized,
                        batch,
                        threshold,
                        self.threads_per_cdhit,
                        self.config,
                        f"{lib_type}_batch_{i}"
                    )
                    future_to_batch[future] = i
                
                # 收集结果
                for future in as_completed(future_to_batch):
                    batch_id = future_to_batch[future]
                    try:
                        batch_result = future.result(timeout=600)
                        all_representatives.extend(batch_result)
                        logger.debug(f"{lib_type} batch {batch_id} completed: {len(batch_result)} representatives")
                    except Exception as e:
                        logger.error(f"{lib_type} batch {batch_id} failed: {e}")
        
        except Exception as e:
            logger.error(f"Batch parallel CD-HIT failed: {e}")
            # 降级为单个CD-HIT
            logger.warning("Falling back to single CD-HIT process")
            all_representatives = run_cdhit_optimized(
                sequences, threshold, self.n_threads, self.config, lib_type
            )
        
        return all_representatives
    
    def write_consensus_libraries(self, masking_lib: List[Dict], analysis_lib: List[Dict]):
        """写入共识库文件"""
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        
        # 输出文件路径
        masking_file = Path(self.config.output_dir) / "consensus_masking.fa"
        analysis_file = Path(self.config.output_dir) / "consensus_analysis.fa"
        
        # 写入masking library
        masking_records = []
        for i, seq_data in enumerate(masking_lib):
            record = SeqRecord(
                Seq(seq_data['sequence']),
                id=seq_data.get('id', f"consensus_masking_{i}"),
                description=f"copies={seq_data.get('num_copies', 'unknown')}"
            )
            masking_records.append(record)
        
        with open(masking_file, 'w') as f:
            SeqIO.write(masking_records, f, "fasta")
        
        # 写入analysis library (带详细信息)
        analysis_records = []
        for i, seq_data in enumerate(analysis_lib):
            metadata = []
            metadata.append(f"copies={seq_data.get('num_copies', 'unknown')}")
            metadata.append(f"length={len(seq_data.get('sequence', ''))}")
            if 'seed_id' in seq_data:
                metadata.append(f"seed={seq_data['seed_id']}")
            if seq_data.get('has_tsd', False):
                metadata.append(f"TSD={seq_data.get('tsd_sequence', 'unknown')}")
            
            record = SeqRecord(
                Seq(seq_data['sequence']),
                id=seq_data.get('id', f"consensus_analysis_{i}"),
                description=" ".join(metadata)
            )
            analysis_records.append(record)
        
        with open(analysis_file, 'w') as f:
            SeqIO.write(analysis_records, f, "fasta")
        
        logger.info(f"Output files written: {masking_file}, {analysis_file}")
    
    def generate_statistics(self, original: List[Dict], masking: List[Dict], 
                          analysis: List[Dict]) -> Dict[str, Any]:
        """生成统计报告"""
        stats = {
            'input_sequences': len(original),
            'masking_library_size': len(masking),
            'analysis_library_size': len(analysis),
            'total_length_masking': sum(len(seq.get('sequence', '')) for seq in masking),
            'total_length_analysis': sum(len(seq.get('sequence', '')) for seq in analysis),
        }
        
        # 长度分布
        if analysis:
            lengths = [len(seq.get('sequence', '')) for seq in analysis]
            stats['length_distribution'] = {
                'min': min(lengths),
                'max': max(lengths),
                'mean': np.mean(lengths),
                'median': np.median(lengths)
            }
        
        # 拷贝数分布
        if analysis:
            copy_numbers = [seq.get('num_copies', 1) for seq in analysis if seq.get('num_copies')]
            if copy_numbers:
                stats['copy_distribution'] = {
                    'min': min(copy_numbers),
                    'max': max(copy_numbers),
                    'mean': np.mean(copy_numbers),
                    'median': np.median(copy_numbers)
                }
        
        # 写入统计文件
        stats_file = Path(self.config.output_dir) / "statistics.txt"
        with open(stats_file, 'w') as f:
            f.write("TE Consensus Pipeline Statistics\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Input sequences: {stats['input_sequences']}\n")
            f.write(f"Masking library: {stats['masking_library_size']} sequences\n")
            f.write(f"Analysis library: {stats['analysis_library_size']} sequences\n")
            f.write(f"Total length (masking): {stats['total_length_masking']:,} bp\n")
            f.write(f"Total length (analysis): {stats['total_length_analysis']:,} bp\n")
            
            if 'length_distribution' in stats:
                f.write(f"\nLength distribution:\n")
                f.write(f"  Min: {stats['length_distribution']['min']} bp\n")
                f.write(f"  Max: {stats['length_distribution']['max']} bp\n")
                f.write(f"  Mean: {stats['length_distribution']['mean']:.1f} bp\n")
                f.write(f"  Median: {stats['length_distribution']['median']:.1f} bp\n")
            
            if 'copy_distribution' in stats:
                f.write(f"\nCopy number distribution:\n")
                f.write(f"  Min: {stats['copy_distribution']['min']}\n")
                f.write(f"  Max: {stats['copy_distribution']['max']}\n")
                f.write(f"  Mean: {stats['copy_distribution']['mean']:.1f}\n")
                f.write(f"  Median: {stats['copy_distribution']['median']:.1f}\n")
        
        logger.info(f"Statistics written to {stats_file}")
        return stats


def filter_batch_by_copy_number(batch: List[Dict], threshold: int) -> List[Dict]:
    """批处理拷贝数过滤"""
    return [seq for seq in batch if seq.get('num_copies', 1) >= threshold]


def simple_redundancy_removal(sequences: List[Dict], threshold: float) -> List[Dict]:
    """简单的冗余去除备用方案"""
    if not sequences:
        return []
    
    # Sort by length (longest first)
    sorted_seqs = sorted(sequences, key=lambda x: len(x['sequence']), reverse=True)
    
    # Keep representatives
    representatives = []
    for seq in sorted_seqs:
        is_redundant = False
        seq_str = seq['sequence'].upper()
        
        # Check against existing representatives
        for rep in representatives:
            rep_str = rep['sequence'].upper()
            
            # Simple similarity check (containment)
            if len(seq_str) <= len(rep_str):
                if seq_str in rep_str:
                    is_redundant = True
                    break
            else:
                if rep_str in seq_str:
                    # Replace shorter representative with longer sequence
                    representatives.remove(rep)
                    representatives.append(seq)
                    is_redundant = True
                    break
        
        if not is_redundant:
            representatives.append(seq)
    
    logger.info(f"Simple redundancy removal: {len(sequences)} -> {len(representatives)}")
    return representatives


def run_cdhit_optimized(sequences: List[Dict], threshold: float, threads: int, 
                       config, label: str) -> List[Dict]:
    """优化的CD-HIT运行"""
    import tempfile
    import subprocess
    import os
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
        import shutil
        if not shutil.which(config.cdhit_exe):
            logger.error(f"CD-HIT executable not found: {config.cdhit_exe}")
            logger.error("Please install CD-HIT or specify correct path")
            raise FileNotFoundError(f"CD-HIT not found: {config.cdhit_exe}")
        
        # 构建CD-HIT命令
        cmd = [
            config.cdhit_exe,
            '-i', input_path,
            '-o', output_path,
            '-c', str(threshold),  # Use the provided threshold parameter
            '-aS', '0.8',  # 80% coverage of shorter sequence
            '-aL', '0.8',  # 80% coverage of longer sequence
            '-G', '0',  # local alignment mode
            '-n', '8' if threshold >= 0.8 else '5',  # word size based on threshold
            '-r', '1',  # compare both strands
            '-mask', 'NX',  # mask low-complexity regions
            '-M', '0',  # no memory limit
            '-T', str(threads),  # use all available threads
            '-d', '0'  # keep full sequence names
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
        # Simple fallback: sort by length and remove very similar sequences
        return simple_redundancy_removal(sequences, threshold)
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
