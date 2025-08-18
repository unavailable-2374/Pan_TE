#!/usr/bin/env python3
"""
测试Phase 2的线程控制是否正常工作
监控MAFFT进程的线程使用情况
"""

import os
import time
import subprocess
import threading
from pathlib import Path
import psutil
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def monitor_mafft_threads(duration=60, interval=1):
    """
    监控MAFFT进程的线程使用情况
    
    Args:
        duration: 监控持续时间（秒）
        interval: 采样间隔（秒）
    """
    logger.info(f"Starting MAFFT thread monitoring for {duration} seconds...")
    
    mafft_stats = {
        'max_processes': 0,
        'max_threads_per_process': 0,
        'total_max_threads': 0,
        'samples': []
    }
    
    start_time = time.time()
    
    while time.time() - start_time < duration:
        try:
            # 查找所有MAFFT进程
            mafft_processes = []
            for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
                try:
                    if 'mafft' in proc.info['name'].lower():
                        mafft_processes.append(proc)
                    elif proc.info['cmdline'] and any('mafft' in str(arg).lower() for arg in proc.info['cmdline']):
                        mafft_processes.append(proc)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
            
            if mafft_processes:
                sample = {
                    'timestamp': time.time() - start_time,
                    'num_processes': len(mafft_processes),
                    'threads_per_process': [],
                    'total_threads': 0
                }
                
                for proc in mafft_processes:
                    try:
                        num_threads = proc.num_threads()
                        sample['threads_per_process'].append(num_threads)
                        sample['total_threads'] += num_threads
                        
                        # 更新最大值
                        mafft_stats['max_threads_per_process'] = max(
                            mafft_stats['max_threads_per_process'], 
                            num_threads
                        )
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        continue
                
                # 更新统计
                mafft_stats['max_processes'] = max(mafft_stats['max_processes'], sample['num_processes'])
                mafft_stats['total_max_threads'] = max(mafft_stats['total_max_threads'], sample['total_threads'])
                mafft_stats['samples'].append(sample)
                
                # 实时输出
                logger.info(f"[{sample['timestamp']:.1f}s] MAFFT processes: {sample['num_processes']}, "
                          f"Threads: {sample['threads_per_process']}, "
                          f"Total: {sample['total_threads']}")
            
            time.sleep(interval)
            
        except Exception as e:
            logger.error(f"Monitoring error: {e}")
            continue
    
    # 输出统计结果
    logger.info("\n" + "="*60)
    logger.info("MAFFT Thread Usage Summary:")
    logger.info(f"  Max concurrent MAFFT processes: {mafft_stats['max_processes']}")
    logger.info(f"  Max threads per MAFFT process: {mafft_stats['max_threads_per_process']}")
    logger.info(f"  Max total threads at once: {mafft_stats['total_max_threads']}")
    
    if mafft_stats['samples']:
        avg_processes = sum(s['num_processes'] for s in mafft_stats['samples']) / len(mafft_stats['samples'])
        avg_total_threads = sum(s['total_threads'] for s in mafft_stats['samples']) / len(mafft_stats['samples'])
        logger.info(f"  Average MAFFT processes: {avg_processes:.1f}")
        logger.info(f"  Average total threads: {avg_total_threads:.1f}")
    
    return mafft_stats


def test_phase2_threading():
    """
    测试Phase 2的线程控制
    创建一个小型测试案例并监控线程使用
    """
    from config import PipelineConfig
    from phase2_consensus_expansion import ConsensusExpansionBuilder
    
    logger.info("Setting up test environment...")
    
    # 创建测试配置
    config = PipelineConfig(
        repeatscout_file="dummy.fa",
        genome_file="dummy.fa",
        output_dir="test_output",
        threads=32  # 设置较高的线程数来测试控制效果
    )
    
    # 创建Phase 2实例
    phase2 = ConsensusExpansionBuilder(config)
    
    logger.info(f"Phase 2 Configuration:")
    logger.info(f"  Total threads: {config.threads}")
    logger.info(f"  Max workers: {phase2.max_workers}")
    logger.info(f"  Threads per worker: {phase2.threads_per_worker}")
    logger.info(f"  MAFFT threads: {phase2.mafft_threads}")
    
    # 计算预期的最大线程数
    expected_max_threads = phase2.max_workers * phase2.mafft_threads
    logger.info(f"  Expected max total threads: {expected_max_threads}")
    
    # 创建模拟数据
    test_sequences = []
    for i in range(10):  # 创建10个测试序列
        test_sequences.append({
            'id': f'test_seq_{i}',
            'sequence': 'ATGC' * 250,  # 1000bp序列
            'quality_class': 'A',
            'rm_hits': [{'chrom': 'chr1', 'start': j*1000, 'end': j*1000+900} 
                       for j in range(50)]  # 每个序列50个hits
        })
    
    phase1_output = {
        'a_sequences': test_sequences,
        'b_sequences': [],
        'rm_detailed_results': {},
        'scores': {seq['id']: {'final': 0.9} for seq in test_sequences}
    }
    
    # 启动监控线程
    monitor_thread = threading.Thread(
        target=monitor_mafft_threads,
        args=(30,),  # 监控30秒
        daemon=True
    )
    monitor_thread.start()
    
    # 等待监控启动
    time.sleep(2)
    
    logger.info("\nStarting Phase 2 processing...")
    
    try:
        # 运行Phase 2（这会触发MAFFT调用）
        results = phase2.run(phase1_output)
        logger.info(f"Phase 2 completed: {len(results)} consensus sequences generated")
    except Exception as e:
        logger.error(f"Phase 2 failed: {e}")
    
    # 等待监控完成
    monitor_thread.join(timeout=35)
    
    logger.info("\nTest completed!")


if __name__ == "__main__":
    # 检查系统信息
    logger.info("System Information:")
    logger.info(f"  CPU cores: {psutil.cpu_count(logical=False)}")
    logger.info(f"  Logical CPUs: {psutil.cpu_count(logical=True)}")
    logger.info(f"  Memory: {psutil.virtual_memory().total / (1024**3):.1f} GB")
    
    test_phase2_threading()