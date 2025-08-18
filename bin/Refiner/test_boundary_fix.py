#!/usr/bin/env python3
"""
测试序列提取的边界检查修复
"""

import tempfile
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from utils.fast_sequence_extractor_fixed import FastSequenceExtractor


def create_test_genome():
    """创建测试基因组"""
    records = []
    
    # 创建不同长度的染色体
    chr_sizes = {
        'chr1': 10000,
        'chr2': 5000,
        'chr3': 2000,
        'chr4': 1000
    }
    
    for chr_name, size in chr_sizes.items():
        # 生成随机序列
        bases = ['A', 'T', 'G', 'C']
        seq = ''.join(random.choices(bases, k=size))
        record = SeqRecord(
            Seq(seq),
            id=chr_name,
            description=f"Test chromosome {chr_name}"
        )
        records.append(record)
    
    # 写入临时文件
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
    SeqIO.write(records, temp_file.name, 'fasta')
    temp_file.close()
    
    return temp_file.name, chr_sizes


def test_boundary_cases():
    """测试边界情况"""
    print("Testing boundary cases for sequence extraction...")
    
    # 创建测试基因组
    genome_file, chr_sizes = create_test_genome()
    print(f"Created test genome: {genome_file}")
    print(f"Chromosome sizes: {chr_sizes}")
    
    # 初始化提取器
    extractor = FastSequenceExtractor(genome_file, tool_preference='auto')
    print(f"Using extraction tool: {extractor.tool}")
    
    # 测试用例：各种边界情况
    test_cases = [
        # 正常情况
        {'name': 'Normal case', 'chrom': 'chr1', 'start': 100, 'end': 200},
        
        # 起始位置为负
        {'name': 'Negative start', 'chrom': 'chr1', 'start': -50, 'end': 100},
        
        # 结束位置超出染色体长度
        {'name': 'End beyond chr', 'chrom': 'chr2', 'start': 4900, 'end': 5200},
        
        # 起始和结束都超出
        {'name': 'Both beyond', 'chrom': 'chr3', 'start': 1900, 'end': 2500},
        
        # 极端情况：完全超出
        {'name': 'Completely beyond', 'chrom': 'chr4', 'start': 2000, 'end': 3000},
        
        # 反向坐标
        {'name': 'Reversed coords', 'chrom': 'chr1', 'start': 500, 'end': 400},
        
        # 包含flanking的情况
        {'name': 'With flanking', 'chrom': 'chr2', 'start': -20, 'end': 5020},
        
        # 染色体不存在
        {'name': 'Non-existent chr', 'chrom': 'chr99', 'start': 100, 'end': 200},
    ]
    
    print("\n" + "="*60)
    print("Running test cases:")
    print("="*60)
    
    for test_case in test_cases:
        print(f"\nTest: {test_case['name']}")
        print(f"  Input: {test_case['chrom']}:{test_case['start']}-{test_case['end']}")
        
        regions = [{
            'chrom': test_case['chrom'],
            'start': test_case['start'],
            'end': test_case['end'],
            'strand': '+'
        }]
        
        try:
            results = extractor.extract_batch(regions)
            
            if results and results[0]:
                result = results[0]
                seq_len = len(result.get('sequence', ''))
                adjusted_start = result.get('start', test_case['start'])
                adjusted_end = result.get('end', test_case['end'])
                
                print(f"  Adjusted: {result['chrom']}:{adjusted_start}-{adjusted_end}")
                print(f"  Sequence length: {seq_len}")
                
                # 验证长度是否匹配
                expected_len = adjusted_end - adjusted_start
                if seq_len != expected_len:
                    print(f"  WARNING: Length mismatch! Expected {expected_len}, got {seq_len}")
                else:
                    print(f"  ✓ Length matches")
            else:
                print(f"  No result returned")
                
        except Exception as e:
            print(f"  ERROR: {e}")
    
    # 批量测试
    print("\n" + "="*60)
    print("Batch extraction test:")
    print("="*60)
    
    # 生成100个随机区域，包括一些边界情况
    batch_regions = []
    for i in range(100):
        chr_name = random.choice(list(chr_sizes.keys()))
        chr_size = chr_sizes[chr_name]
        
        # 25%的概率生成边界情况
        if random.random() < 0.25:
            # 边界情况
            if random.random() < 0.5:
                # 起始为负
                start = random.randint(-100, 0)
                end = random.randint(100, 500)
            else:
                # 结束超出
                start = chr_size - 100
                end = chr_size + 100
        else:
            # 正常情况
            start = random.randint(0, max(0, chr_size - 200))
            end = start + random.randint(50, 200)
        
        batch_regions.append({
            'chrom': chr_name,
            'start': start,
            'end': end,
            'strand': random.choice(['+', '-'])
        })
    
    print(f"Extracting {len(batch_regions)} regions...")
    
    try:
        results = extractor.extract_batch(batch_regions)
        
        success_count = sum(1 for r in results if r.get('sequence'))
        empty_count = sum(1 for r in results if not r.get('sequence'))
        
        print(f"Results:")
        print(f"  Successful extractions: {success_count}")
        print(f"  Empty sequences: {empty_count}")
        
        # 检查是否有异常长度
        for i, result in enumerate(results):
            if result.get('sequence'):
                seq_len = len(result['sequence'])
                expected_len = result['end'] - result['start']
                if abs(seq_len - expected_len) > 1:  # 允许1bp的误差
                    print(f"  Region {i}: Length mismatch - expected {expected_len}, got {seq_len}")
        
        print(f"\n✓ Batch extraction completed successfully!")
        
    except Exception as e:
        print(f"ERROR in batch extraction: {e}")
        import traceback
        print(traceback.format_exc())
    
    # 清理
    import os
    os.unlink(genome_file)
    if os.path.exists(f"{genome_file}.fai"):
        os.unlink(f"{genome_file}.fai")
    
    print("\nTest completed!")


if __name__ == "__main__":
    test_boundary_cases()