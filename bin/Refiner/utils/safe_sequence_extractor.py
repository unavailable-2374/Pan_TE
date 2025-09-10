#!/usr/bin/env python3
"""
安全的序列提取器
同时处理ID映射和坐标调整问题
"""

import logging
from typing import Optional, Dict, List, Tuple
from Bio import SeqIO
from .id_mapper import ChromosomeIDMapper

logger = logging.getLogger(__name__)

class SafeSequenceExtractor:
    """安全的序列提取器，处理ID映射和坐标不匹配问题"""
    
    def __init__(self, genome_file: str, id_mapper: ChromosomeIDMapper = None):
        self.genome_file = genome_file
        self.id_mapper = id_mapper or ChromosomeIDMapper(genome_file)
        self.sequence_cache = {}  # 缓存序列信息以提升性能
        
    def _load_sequence_info(self, chrom_id: str) -> Optional[Dict]:
        """加载序列信息（长度等）"""
        if chrom_id in self.sequence_cache:
            return self.sequence_cache[chrom_id]
        
        try:
            for record in SeqIO.parse(self.genome_file, "fasta"):
                if record.id == chrom_id:
                    seq_info = {
                        'id': record.id,
                        'length': len(record.seq),
                        'sequence': record.seq
                    }
                    self.sequence_cache[chrom_id] = seq_info
                    return seq_info
        except Exception as e:
            logger.error(f"Failed to load sequence {chrom_id}: {e}")
        
        return None
    
    def safe_extract_sequence(self, chrom_id: str, start: int, end: int, 
                            flanking: int = 0, fallback_strategy: str = 'truncate') -> Optional[Dict]:
        """
        安全的序列提取，处理各种异常情况
        
        Args:
            chrom_id: 染色体ID（可能是lcl|格式）
            start: 起始位置（1-based）
            end: 结束位置（1-based）
            flanking: 额外的flanking区域
            fallback_strategy: 坐标超出范围时的处理策略
                - 'truncate': 截断到有效范围
                - 'skip': 跳过无效的提取
                - 'whole_sequence': 使用整个序列
        
        Returns:
            提取结果字典或None
        """
        # 1. 映射染色体ID
        real_chrom = self.id_mapper.map_chromosome_id(chrom_id)
        
        if not real_chrom:
            logger.warning(f"Cannot map chromosome ID: {chrom_id}")
            return None
        
        # 2. 加载序列信息
        seq_info = self._load_sequence_info(real_chrom)
        if not seq_info:
            logger.error(f"Sequence {real_chrom} not found in genome file")
            return None
        
        seq_length = seq_info['length']
        sequence = seq_info['sequence']
        
        # 3. 坐标有效性检查和调整
        original_start, original_end = start, end
        
        # 添加flanking区域
        extract_start = max(1, start - flanking)
        extract_end = min(seq_length, end + flanking)
        
        # 处理坐标超出范围的情况
        coordinate_issues = []
        
        if start > seq_length:
            coordinate_issues.append(f"start {start} > sequence length {seq_length}")
            if fallback_strategy == 'skip':
                logger.warning(f"Skipping extraction due to invalid start coordinate: {chrom_id}:{start}-{end}")
                return None
            elif fallback_strategy == 'whole_sequence':
                extract_start, extract_end = 1, seq_length
            # truncate策略会在下面处理
        
        if end > seq_length:
            coordinate_issues.append(f"end {end} > sequence length {seq_length}")
            if fallback_strategy == 'truncate':
                extract_end = seq_length
                logger.info(f"Truncating end coordinate from {end} to {seq_length}")
            elif fallback_strategy == 'whole_sequence':
                extract_start, extract_end = 1, seq_length
        
        if start > end:
            coordinate_issues.append(f"start {start} > end {end}")
            logger.error(f"Invalid coordinate range: {chrom_id}:{start}-{end}")
            return None
        
        # 4. 提取序列
        try:
            # 转换为0-based坐标
            extracted_seq = sequence[extract_start-1:extract_end]
            
            result = {
                'sequence': str(extracted_seq),
                'chrom': real_chrom,
                'original_chrom': chrom_id,
                'start': extract_start,
                'end': extract_end,
                'original_start': original_start,
                'original_end': original_end,
                'length': len(extracted_seq),
                'sequence_length': seq_length,
                'flanking': flanking,
                'coordinate_issues': coordinate_issues,
                'fallback_strategy': fallback_strategy if coordinate_issues else None
            }
            
            if coordinate_issues:
                logger.warning(f"Coordinate issues for {chrom_id}:{original_start}-{original_end}: {'; '.join(coordinate_issues)}")
            else:
                logger.debug(f"Successfully extracted {len(extracted_seq)} bp from {real_chrom}:{extract_start}-{extract_end}")
            
            return result
            
        except Exception as e:
            logger.error(f"Failed to extract sequence from {real_chrom}:{extract_start}-{extract_end}: {e}")
            return None
    
    def extract_batch_safe(self, regions: List[Dict], flanking: int = 0, 
                          fallback_strategy: str = 'truncate') -> List[Dict]:
        """
        批量安全序列提取
        
        Args:
            regions: 区域列表，每个元素包含 {'chrom', 'start', 'end'}
            flanking: flanking区域大小
            fallback_strategy: 处理策略
            
        Returns:
            提取结果列表
        """
        results = []
        stats = {
            'total': len(regions),
            'successful': 0,
            'id_mapping_failed': 0,
            'coordinate_issues': 0,
            'extraction_failed': 0
        }
        
        for i, region in enumerate(regions):
            chrom_id = region.get('chrom')
            start = region.get('start')
            end = region.get('end')
            
            if not all([chrom_id, start, end]):
                logger.warning(f"Region {i} missing required fields: {region}")
                stats['extraction_failed'] += 1
                continue
            
            result = self.safe_extract_sequence(
                chrom_id, start, end, flanking, fallback_strategy
            )
            
            if result:
                # 保留原始region的其他信息
                result.update({k: v for k, v in region.items() 
                             if k not in ['chrom', 'start', 'end']})
                results.append(result)
                stats['successful'] += 1
                
                if result.get('coordinate_issues'):
                    stats['coordinate_issues'] += 1
            else:
                if not self.id_mapper.map_chromosome_id(chrom_id):
                    stats['id_mapping_failed'] += 1
                else:
                    stats['extraction_failed'] += 1
        
        logger.info(f"Batch extraction stats: {stats}")
        return results
    
    def diagnose_region(self, chrom_id: str, start: int, end: int) -> Dict:
        """诊断特定区域的问题"""
        diagnosis = {
            'chrom_id': chrom_id,
            'coordinates': f"{start}-{end}",
            'issues': [],
            'suggestions': []
        }
        
        # ID映射诊断
        real_chrom = self.id_mapper.map_chromosome_id(chrom_id)
        if not real_chrom:
            diagnosis['issues'].append(f"Cannot map chromosome ID: {chrom_id}")
            diagnosis['suggestions'].append("Check if RepeatMasker used a different genome file")
            return diagnosis
        
        diagnosis['mapped_chrom'] = real_chrom
        
        # 序列存在性诊断
        seq_info = self._load_sequence_info(real_chrom)
        if not seq_info:
            diagnosis['issues'].append(f"Sequence {real_chrom} not found in genome file")
            diagnosis['suggestions'].append("Verify genome file completeness")
            return diagnosis
        
        seq_length = seq_info['length']
        diagnosis['sequence_length'] = seq_length
        
        # 坐标范围诊断
        if start > seq_length:
            diagnosis['issues'].append(f"Start coordinate {start} exceeds sequence length {seq_length}")
        if end > seq_length:
            diagnosis['issues'].append(f"End coordinate {end} exceeds sequence length {seq_length}")
        if start > end:
            diagnosis['issues'].append(f"Invalid range: start {start} > end {end}")
        
        if diagnosis['issues']:
            diagnosis['suggestions'].extend([
                "Use fallback_strategy='truncate' to handle coordinate overflow",
                "Consider using the whole sequence instead of specific coordinates",
                "Check if coordinates are relative to a different reference"
            ])
        else:
            diagnosis['status'] = 'OK'
        
        return diagnosis

def create_safe_extractor(genome_file: str, repeatmasker_dir: str = None) -> SafeSequenceExtractor:
    """创建安全序列提取器的便捷函数"""
    id_mapper = ChromosomeIDMapper(genome_file, repeatmasker_dir)
    return SafeSequenceExtractor(genome_file, id_mapper)

if __name__ == "__main__":
    # 测试代码
    import sys
    if len(sys.argv) > 1:
        genome_file = sys.argv[1]
        extractor = create_safe_extractor(genome_file)
        
        # 测试问题区域
        test_regions = [
            {'chrom': 'lcl|30_3', 'start': 458290, 'end': 458597},
            {'chrom': 'lcl|17_3', 'start': 301361, 'end': 301512},
            {'chrom': 'chr1', 'start': 1000, 'end': 2000}
        ]
        
        for region in test_regions:
            print(f"\n诊断: {region}")
            diagnosis = extractor.diagnose_region(**region)
            print(f"结果: {diagnosis}")
            
            extraction = extractor.safe_extract_sequence(**region)
            if extraction:
                print(f"提取成功: {extraction['length']} bp")
            else:
                print("提取失败")