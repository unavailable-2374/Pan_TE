"""
高性能序列提取器：使用外部工具替代BioPython的逐个提取
支持samtools、bedtools等工具的批量提取
"""

import os
import tempfile
import subprocess
import logging
from typing import List, Dict, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)


class FastSequenceExtractor:
    """使用外部工具进行高效批量序列提取"""
    
    def __init__(self, genome_file: str, tool_preference: str = 'auto'):
        """
        初始化提取器
        
        Args:
            genome_file: 基因组FASTA文件路径
            tool_preference: 工具偏好 ('samtools', 'bedtools', 'seqtk', 'auto')
        """
        self.genome_file = genome_file
        self.tool = self._select_tool(tool_preference)
        self._ensure_index()
        
    def _select_tool(self, preference: str) -> str:
        """选择可用的提取工具"""
        tools = {
            'samtools': 'samtools',
            'bedtools': 'bedtools', 
            'seqtk': 'seqtk'
        }
        
        if preference != 'auto' and preference in tools:
            if self._check_tool(tools[preference]):
                logger.info(f"Using {preference} for sequence extraction")
                return preference
            else:
                logger.warning(f"{preference} not available, falling back to auto-detection")
        
        # 自动检测可用工具（优先级顺序）
        for tool_name, tool_cmd in [
            ('bedtools', 'bedtools'),  # 最快的批量提取
            ('samtools', 'samtools'),  # 其次快，广泛可用
            ('seqtk', 'seqtk')         # 也很快
        ]:
            if self._check_tool(tool_cmd):
                logger.info(f"Auto-selected {tool_name} for sequence extraction")
                return tool_name
        
        logger.warning("No external extraction tool found, falling back to BioPython")
        return 'biopython'
    
    def _check_tool(self, tool_cmd: str) -> bool:
        """检查工具是否可用"""
        try:
            result = subprocess.run(
                [tool_cmd, '--version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    
    def _ensure_index(self):
        """确保基因组索引存在"""
        if self.tool == 'samtools':
            index_file = f"{self.genome_file}.fai"
            if not os.path.exists(index_file):
                logger.info(f"Creating samtools index for {self.genome_file}")
                subprocess.run(
                    ['samtools', 'faidx', self.genome_file],
                    check=True
                )
        elif self.tool == 'bedtools':
            # bedtools也使用samtools的fai索引
            index_file = f"{self.genome_file}.fai"
            if not os.path.exists(index_file):
                logger.info(f"Creating fai index for bedtools")
                subprocess.run(
                    ['samtools', 'faidx', self.genome_file],
                    check=True
                )
    
    def extract_batch(self, regions: List[Dict]) -> List[Dict]:
        """
        批量提取序列
        
        Args:
            regions: 区域列表，每个包含 'chrom', 'start', 'end', 'strand' 等字段
            
        Returns:
            序列列表，每个包含 'sequence' 字段和原始区域信息
        """
        if not regions:
            return []
        
        if self.tool == 'bedtools':
            return self._extract_with_bedtools(regions)
        elif self.tool == 'samtools':
            return self._extract_with_samtools(regions)
        elif self.tool == 'seqtk':
            return self._extract_with_seqtk(regions)
        else:
            return self._extract_with_biopython(regions)
    
    def _extract_with_bedtools(self, regions: List[Dict]) -> List[Dict]:
        """使用bedtools getfasta批量提取（最快）"""
        try:
            # 创建临时BED文件
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_file:
                bed_path = bed_file.name
                for i, region in enumerate(regions):
                    # BED格式: chrom start end name score strand
                    strand = region.get('strand', '+')
                    bed_file.write(f"{region['chrom']}\t{region['start']}\t{region['end']}\t"
                                 f"region_{i}\t0\t{strand}\n")
            
            # 运行bedtools getfasta
            with tempfile.NamedTemporaryFile(mode='r', suffix='.fa', delete=False) as out_file:
                out_path = out_file.name
            
            cmd = [
                'bedtools', 'getfasta',
                '-fi', self.genome_file,
                '-bed', bed_path,
                '-fo', out_path,
                '-s'  # 考虑链方向
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"bedtools failed: {result.stderr}")
                return self._extract_with_samtools(regions)  # 降级到samtools
            
            # 解析结果
            sequences = []
            with open(out_path, 'r') as f:
                current_seq = []
                current_id = None
                
                for line in f:
                    if line.startswith('>'):
                        if current_id is not None:
                            seq_idx = int(current_id.split('_')[1])
                            sequences.append({
                                **regions[seq_idx],
                                'sequence': ''.join(current_seq)
                            })
                        current_id = line.strip()[1:].split('::')[0]
                        current_seq = []
                    else:
                        current_seq.append(line.strip())
                
                # 处理最后一个序列
                if current_id is not None:
                    seq_idx = int(current_id.split('_')[1])
                    sequences.append({
                        **regions[seq_idx],
                        'sequence': ''.join(current_seq)
                    })
            
            # 清理临时文件
            os.unlink(bed_path)
            os.unlink(out_path)
            
            return sequences
            
        except Exception as e:
            logger.error(f"bedtools extraction failed: {e}")
            return self._extract_with_samtools(regions)
    
    def _extract_with_samtools(self, regions: List[Dict]) -> List[Dict]:
        """使用samtools faidx批量提取"""
        sequences = []
        
        # samtools可以一次提取多个区域
        region_specs = []
        for i, region in enumerate(regions):
            region_spec = f"{region['chrom']}:{region['start']+1}-{region['end']}"
            region_specs.append(region_spec)
        
        try:
            # 批量提取（分批以避免命令行过长）
            batch_size = 100
            for batch_start in range(0, len(region_specs), batch_size):
                batch_specs = region_specs[batch_start:batch_start + batch_size]
                batch_regions = regions[batch_start:batch_start + batch_size]
                
                cmd = ['samtools', 'faidx', self.genome_file] + batch_specs
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    logger.error(f"samtools failed: {result.stderr}")
                    continue
                
                # 解析FASTA输出
                current_seq = []
                current_idx = -1
                
                for line in result.stdout.split('\n'):
                    if line.startswith('>'):
                        if current_idx >= 0:
                            region = batch_regions[current_idx]
                            seq = ''.join(current_seq)
                            
                            # 处理负链
                            if region.get('strand', '+') == '-':
                                seq = self._reverse_complement(seq)
                            
                            sequences.append({
                                **region,
                                'sequence': seq
                            })
                        
                        current_idx += 1
                        current_seq = []
                    elif line:
                        current_seq.append(line.strip())
                
                # 处理最后一个序列
                if current_idx >= 0 and current_idx < len(batch_regions):
                    region = batch_regions[current_idx]
                    seq = ''.join(current_seq)
                    
                    if region.get('strand', '+') == '-':
                        seq = self._reverse_complement(seq)
                    
                    sequences.append({
                        **region,
                        'sequence': seq
                    })
            
            return sequences
            
        except Exception as e:
            logger.error(f"samtools extraction failed: {e}")
            return self._extract_with_biopython(regions)
    
    def _extract_with_seqtk(self, regions: List[Dict]) -> List[Dict]:
        """使用seqtk subseq提取"""
        # seqtk主要用于按名称提取，对于坐标提取不如其他工具方便
        # 这里提供降级到samtools
        return self._extract_with_samtools(regions)
    
    def _extract_with_biopython(self, regions: List[Dict]) -> List[Dict]:
        """使用BioPython提取（最慢，作为后备）"""
        from Bio import SeqIO
        
        sequences = []
        
        # 加载基因组到内存（对大基因组可能有问题）
        logger.warning("Using slow BioPython extraction as fallback")
        genome_dict = {}
        
        with open(self.genome_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                genome_dict[record.id] = str(record.seq)
        
        for region in regions:
            chrom = region['chrom']
            if chrom in genome_dict:
                seq = genome_dict[chrom][region['start']:region['end']]
                
                if region.get('strand', '+') == '-':
                    seq = self._reverse_complement(seq)
                
                sequences.append({
                    **region,
                    'sequence': seq
                })
            else:
                logger.warning(f"Chromosome {chrom} not found in genome")
                sequences.append({
                    **region,
                    'sequence': ''
                })
        
        return sequences
    
    def _reverse_complement(self, seq: str) -> str:
        """获取反向互补序列"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                     'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                     'N': 'N', 'n': 'n'}
        
        return ''.join(complement.get(base, base) for base in reversed(seq))


def extract_sequences_batch(genome_file: str, regions: List[Dict], 
                           tool_preference: str = 'auto') -> List[Dict]:
    """
    便捷函数：批量提取序列
    
    Args:
        genome_file: 基因组文件路径
        regions: 区域列表
        tool_preference: 工具偏好
        
    Returns:
        包含序列的区域列表
    """
    extractor = FastSequenceExtractor(genome_file, tool_preference)
    return extractor.extract_batch(regions)