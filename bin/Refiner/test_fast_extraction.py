#!/usr/bin/env python3
"""
测试快速序列提取性能
比较BioPython vs 外部工具的性能差异
"""

import time
import random
import tempfile
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 导入要测试的模块
from utils.fast_sequence_extractor import FastSequenceExtractor
from utils.sequence_utils import extract_sequence_from_genome


def create_test_genome(size_mb=10, num_chromosomes=5):
    """创建测试基因组"""
    print(f"Creating test genome: {size_mb}MB with {num_chromosomes} chromosomes")
    
    bases = ['A', 'T', 'G', 'C']
    records = []
    
    chr_size = (size_mb * 1000000) // num_chromosomes
    
    for i in range(num_chromosomes):
        # 生成随机序列
        seq = ''.join(random.choices(bases, k=chr_size))
        record = SeqRecord(
            Seq(seq),
            id=f"chr{i+1}",
            description=f"Test chromosome {i+1}"
        )
        records.append(record)
    
    # 写入临时文件
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
    SeqIO.write(records, temp_file.name, 'fasta')
    temp_file.close()
    
    return temp_file.name, chr_size


def generate_random_regions(num_regions, num_chromosomes, chr_size):
    """生成随机区域"""
    regions = []
    
    for _ in range(num_regions):
        chrom = f"chr{random.randint(1, num_chromosomes)}"
        length = random.randint(100, 2000)
        start = random.randint(0, max(0, chr_size - length - 1))
        end = start + length
        strand = random.choice(['+', '-'])
        
        regions.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand
        })
    
    return regions


def test_biopython_extraction(genome_file, regions):
    """测试BioPython提取性能"""
    print(f"\nTesting BioPython extraction for {len(regions)} regions...")
    
    start_time = time.time()
    sequences = []
    
    for region in regions:
        seq = extract_sequence_from_genome(
            genome_file,
            region['chrom'],
            region['start'],
            region['end'],
            region['strand']
        )
        sequences.append(seq)
    
    elapsed = time.time() - start_time
    
    print(f"  BioPython: {elapsed:.2f} seconds")
    print(f"  Rate: {len(regions)/elapsed:.1f} sequences/second")
    
    return sequences, elapsed


def test_fast_extraction(genome_file, regions, tool='auto'):
    """测试快速提取性能"""
    print(f"\nTesting fast extraction ({tool}) for {len(regions)} regions...")
    
    start_time = time.time()
    
    extractor = FastSequenceExtractor(genome_file, tool)
    results = extractor.extract_batch(regions)
    sequences = [r['sequence'] for r in results]
    
    elapsed = time.time() - start_time
    
    print(f"  {extractor.tool}: {elapsed:.2f} seconds")
    print(f"  Rate: {len(regions)/elapsed:.1f} sequences/second")
    
    return sequences, elapsed


def compare_results(seq1_list, seq2_list):
    """比较两种方法的结果"""
    if len(seq1_list) != len(seq2_list):
        print(f"  WARNING: Different number of sequences: {len(seq1_list)} vs {len(seq2_list)}")
        return False
    
    mismatches = 0
    for i, (s1, s2) in enumerate(zip(seq1_list, seq2_list)):
        if s1 != s2:
            mismatches += 1
            if mismatches <= 3:  # 只显示前3个不匹配
                print(f"  Mismatch at position {i}: len={len(s1)} vs {len(s2)}")
    
    if mismatches == 0:
        print(f"  ✓ All {len(seq1_list)} sequences match perfectly")
        return True
    else:
        print(f"  ✗ {mismatches}/{len(seq1_list)} sequences don't match")
        return False


def run_performance_test():
    """运行完整性能测试"""
    print("="*60)
    print("Fast Sequence Extraction Performance Test")
    print("="*60)
    
    # 创建测试基因组
    genome_file, chr_size = create_test_genome(size_mb=50, num_chromosomes=10)
    print(f"Test genome created: {genome_file}")
    
    # 测试不同规模
    test_sizes = [10, 100, 500, 1000, 5000]
    
    results = {}
    
    for num_regions in test_sizes:
        print(f"\n{'='*40}")
        print(f"Testing with {num_regions} regions")
        print('='*40)
        
        # 生成随机区域
        regions = generate_random_regions(num_regions, 10, chr_size)
        
        # 测试BioPython（仅对小数据集）
        if num_regions <= 100:
            bio_seqs, bio_time = test_biopython_extraction(genome_file, regions)
        else:
            print(f"\nSkipping BioPython for {num_regions} regions (too slow)")
            bio_time = None
            bio_seqs = []
        
        # 测试快速提取
        fast_seqs, fast_time = test_fast_extraction(genome_file, regions)
        
        # 比较结果（如果都运行了）
        if bio_seqs:
            print("\nVerifying results:")
            compare_results(bio_seqs, fast_seqs)
        
        # 计算加速比
        if bio_time:
            speedup = bio_time / fast_time
            print(f"\n⚡ Speedup: {speedup:.1f}x faster than BioPython")
        
        results[num_regions] = {
            'biopython': bio_time,
            'fast': fast_time,
            'speedup': bio_time / fast_time if bio_time else None
        }
    
    # 打印总结
    print("\n" + "="*60)
    print("PERFORMANCE SUMMARY")
    print("="*60)
    print(f"{'Regions':<10} {'BioPython':<12} {'Fast Tool':<12} {'Speedup':<10}")
    print("-"*46)
    
    for num_regions, data in results.items():
        bio_str = f"{data['biopython']:.2f}s" if data['biopython'] else "N/A"
        fast_str = f"{data['fast']:.2f}s"
        speedup_str = f"{data['speedup']:.1f}x" if data['speedup'] else "N/A"
        print(f"{num_regions:<10} {bio_str:<12} {fast_str:<12} {speedup_str:<10}")
    
    # 清理
    os.unlink(genome_file)
    if os.path.exists(f"{genome_file}.fai"):
        os.unlink(f"{genome_file}.fai")
    
    print("\nTest completed successfully!")


if __name__ == "__main__":
    run_performance_test()