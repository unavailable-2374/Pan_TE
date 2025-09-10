#!/usr/bin/env python3
"""
RepeatMasker输出解析器
从RepeatMasker的.out文件中提取真实的序列名映射
"""

import logging
import re
from typing import Dict, List, Optional
from pathlib import Path

logger = logging.getLogger(__name__)

def parse_repeatmasker_out_file(out_file: str) -> Dict[str, str]:
    """
    解析RepeatMasker的.out文件，提取lcl|ID到真实序列名的映射
    
    Args:
        out_file: RepeatMasker输出的.out文件路径
        
    Returns:
        映射字典: {lcl_id: real_sequence_name}
    """
    mappings = {}
    
    if not Path(out_file).exists():
        logger.warning(f"RepeatMasker out file not found: {out_file}")
        return mappings
    
    try:
        with open(out_file, 'r') as f:
            # 跳过前3行的头部信息
            for _ in range(3):
                next(f, None)
            
            for line_num, line in enumerate(f, 4):
                line = line.strip()
                if not line:
                    continue
                
                # RepeatMasker .out格式解析
                # SW  perc perc perc  query      position in query           matching  repeat
                # score div. del. ins.  sequence    begin     end    (left)    repeat   class/family
                parts = line.split()
                if len(parts) >= 5:
                    sequence_name = parts[4]  # 第5列是序列名
                    
                    # 如果序列名包含lcl|格式，这可能是我们需要的映射
                    if 'lcl|' in sequence_name:
                        # 提取lcl|部分
                        lcl_match = re.search(r'lcl\|\d+_\d+', sequence_name)
                        if lcl_match:
                            lcl_id = lcl_match.group()
                            mappings[lcl_id] = sequence_name
                            logger.debug(f"Found mapping: {lcl_id} -> {sequence_name}")
    
    except Exception as e:
        logger.error(f"Failed to parse RepeatMasker out file {out_file}: {e}")
    
    logger.info(f"Extracted {len(mappings)} sequence mappings from {out_file}")
    return mappings

def create_blast_db_mapping(genome_file: str) -> Dict[str, str]:
    """
    通过创建BLAST数据库来获取lcl|ID映射
    这是最可靠的方法，因为我们复现了RepeatMasker的处理过程
    
    Args:
        genome_file: 基因组文件路径
        
    Returns:
        映射字典: {lcl_id: real_sequence_name}
    """
    import subprocess
    import tempfile
    
    mappings = {}
    
    try:
        # 创建临时BLAST数据库
        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = f"{temp_dir}/temp_db"
            
            # 创建BLAST数据库
            cmd = ["makeblastdb", "-in", genome_file, "-dbtype", "nucl", "-out", db_path, "-parse_seqids"]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.error(f"makeblastdb failed: {result.stderr}")
                return mappings
            
            # 使用blastdbcmd获取序列列表和映射
            cmd = ["blastdbcmd", "-db", db_path, "-entry", "all", "-outfmt", "%a %t"]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if line:
                        parts = line.split(' ', 1)
                        if len(parts) >= 2:
                            blast_id = parts[0]
                            original_id = parts[1]
                            
                            # 如果BLAST生成了lcl|格式的ID
                            if blast_id.startswith('lcl|'):
                                mappings[blast_id] = original_id
                                logger.debug(f"BLAST mapping: {blast_id} -> {original_id}")
            else:
                logger.error(f"blastdbcmd failed: {result.stderr}")
    
    except Exception as e:
        logger.error(f"Failed to create BLAST database mapping: {e}")
    
    # 对于已规范化的lcl|格式基因组，通常不需要BLAST映射
    if len(mappings) > 0:
        logger.info(f"Created {len(mappings)} BLAST database mappings")
    else:
        logger.debug(f"Created {len(mappings)} BLAST database mappings")  # 降级为debug，减少噪音
    return mappings

def find_real_sequence_name(lcl_id: str, genome_file: str, repeatmasker_dir: str = None) -> Optional[str]:
    """
    尝试找到lcl|ID对应的真实序列名称
    
    Args:
        lcl_id: lcl|格式的ID，如 "lcl|30_3"
        genome_file: 基因组文件路径
        repeatmasker_dir: RepeatMasker工作目录（如果有）
        
    Returns:
        真实的序列名称，如果找不到则返回None
    """
    # 方法1: 如果有RepeatMasker输出目录，尝试解析.out文件
    if repeatmasker_dir:
        out_files = list(Path(repeatmasker_dir).glob("*.out"))
        for out_file in out_files:
            mappings = parse_repeatmasker_out_file(str(out_file))
            if lcl_id in mappings:
                return mappings[lcl_id]
    
    # 方法2: 创建BLAST数据库映射
    blast_mappings = create_blast_db_mapping(genome_file)
    if lcl_id in blast_mappings:
        return blast_mappings[lcl_id]
    
    # 方法3: 回退到原来的数字映射（作为最后的尝试）
    lcl_match = re.match(r'lcl\|(\d+)_(\d+)', lcl_id)
    if lcl_match:
        index_num = int(lcl_match.group(1))
        logger.warning(f"Using fallback mapping for {lcl_id}, index {index_num}")
        
        # 读取基因组文件，按索引查找
        try:
            from Bio import SeqIO
            for i, record in enumerate(SeqIO.parse(genome_file, "fasta"), 1):
                if i == index_num:
                    logger.info(f"Fallback mapping: {lcl_id} -> {record.id}")
                    return record.id
        except Exception as e:
            logger.error(f"Fallback mapping failed: {e}")
    
    return None

class RepeatMaskerIDResolver:
    """RepeatMasker ID解析器"""
    
    def __init__(self, genome_file: str, repeatmasker_dir: str = None):
        self.genome_file = genome_file
        self.repeatmasker_dir = repeatmasker_dir
        self.mappings = {}
        self._load_mappings()
    
    def _load_mappings(self):
        """加载所有可能的映射"""
        # 从RepeatMasker输出加载
        if self.repeatmasker_dir:
            out_files = list(Path(self.repeatmasker_dir).glob("*.out"))
            for out_file in out_files:
                file_mappings = parse_repeatmasker_out_file(str(out_file))
                self.mappings.update(file_mappings)
        
        # 从BLAST数据库加载
        blast_mappings = create_blast_db_mapping(self.genome_file)
        self.mappings.update(blast_mappings)
        
        # 对于已规范化基因组，通常映射数为0是正常的
        if len(self.mappings) > 0:
            logger.info(f"Loaded {len(self.mappings)} total ID mappings")
        else:
            logger.debug(f"Loaded {len(self.mappings)} total ID mappings")  # 降级为debug
    
    def resolve_id(self, lcl_id: str) -> Optional[str]:
        """解析lcl|ID到真实序列名"""
        # 直接查找
        if lcl_id in self.mappings:
            return self.mappings[lcl_id]
        
        # 动态查找
        return find_real_sequence_name(lcl_id, self.genome_file, self.repeatmasker_dir)
    
    def get_stats(self) -> Dict:
        """获取映射统计信息"""
        return {
            'total_mappings': len(self.mappings),
            'sample_mappings': dict(list(self.mappings.items())[:5])
        }

if __name__ == "__main__":
    # 测试代码
    import sys
    if len(sys.argv) > 1:
        genome_file = sys.argv[1]
        resolver = RepeatMaskerIDResolver(genome_file)
        print(f"Mappings: {resolver.get_stats()}")
        
        # 测试解析
        test_ids = ["lcl|30_3", "lcl|17_3", "lcl|1_1"]
        for test_id in test_ids:
            result = resolver.resolve_id(test_id)
            print(f"{test_id} -> {result}")