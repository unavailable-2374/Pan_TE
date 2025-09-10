"""
高性能序列提取器：使用外部工具替代BioPython的逐个提取
支持samtools、bedtools等工具的批量提取
修复版：包含坐标边界检查和更健壮的解析
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
        self.chrom_sizes = self._load_chrom_sizes()
        
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
            ('samtools', 'samtools'),  # 最稳定，优先使用
            ('bedtools', 'bedtools'),  # 快速但解析复杂
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
        index_file = f"{self.genome_file}.fai"
        if not os.path.exists(index_file):
            logger.debug(f"Index not found for {self.genome_file}, attempting to create")
            try:
                subprocess.run(
                    ['samtools', 'faidx', self.genome_file],
                    check=True,
                    capture_output=True
                )
                logger.debug(f"Successfully created index: {index_file}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to create index: {e}")
    
    def _load_chrom_sizes(self) -> Dict[str, int]:
        """从fai索引文件加载染色体大小"""
        chrom_sizes = {}
        index_file = f"{self.genome_file}.fai"
        
        if os.path.exists(index_file):
            try:
                with open(index_file, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            chrom_sizes[parts[0]] = int(parts[1])
                logger.info(f"Loaded {len(chrom_sizes)} chromosome sizes from index")
            except Exception as e:
                logger.warning(f"Failed to load chromosome sizes: {e}")
        
        return chrom_sizes
    
    def _adjust_coordinates(self, chrom: str, start: int, end: int) -> Tuple[int, int]:
        """
        调整坐标以确保不超出染色体边界
        
        Returns:
            调整后的 (start, end)
        """
        # 确保起始位置非负
        start = max(0, start)
        
        # 如果知道染色体大小，确保不超出边界
        if chrom in self.chrom_sizes:
            chrom_size = self.chrom_sizes[chrom]
            end = min(end, chrom_size)
            
            # 如果调整后start >= end，返回有效的小区域
            if start >= end:
                start = max(0, end - 100)
        
        # 最终检查
        if start >= end:
            logger.warning(f"Invalid coordinates after adjustment: {chrom}:{start}-{end}")
            return 0, min(100, end)
        
        return start, end
    
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
        
        # 调整所有区域的坐标
        adjusted_regions = []
        for region in regions:
            adj_start, adj_end = self._adjust_coordinates(
                region['chrom'], 
                region['start'], 
                region['end']
            )
            adjusted_region = {
                **region,
                'start': adj_start,
                'end': adj_end,
                'original_start': region['start'],
                'original_end': region['end']
            }
            adjusted_regions.append(adjusted_region)
        
        if self.tool == 'bedtools':
            return self._extract_with_bedtools(adjusted_regions)
        elif self.tool == 'samtools':
            return self._extract_with_samtools(adjusted_regions)
        elif self.tool == 'seqtk':
            return self._extract_with_seqtk(adjusted_regions)
        else:
            return self._extract_with_biopython(adjusted_regions)
    
    def _extract_with_bedtools(self, regions: List[Dict]) -> List[Dict]:
        """使用bedtools getfasta批量提取"""
        try:
            # 创建临时BED文件
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_file:
                bed_path = bed_file.name
                region_map = {}  # 用于映射region名称到索引
                
                for i, region in enumerate(regions):
                    # 使用更独特的名称避免解析冲突
                    region_name = f"seq{i:06d}"
                    region_map[region_name] = i
                    
                    # BED格式: chrom start end name score strand
                    strand = region.get('strand', '+')
                    bed_file.write(f"{region['chrom']}\t{region['start']}\t{region['end']}\t"
                                 f"{region_name}\t0\t{strand}\n")
            
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
            sequences = {}
            with open(out_path, 'r') as f:
                current_seq = []
                current_name = None
                
                for line in f:
                    if line.startswith('>'):
                        # 保存上一个序列
                        if current_name is not None and current_name in region_map:
                            idx = region_map[current_name]
                            sequences[idx] = ''.join(current_seq)
                        
                        # 解析新序列名称
                        header = line.strip()[1:]
                        # bedtools可能会添加坐标信息，如 seq000001::chr1:100-200(+)
                        # 提取我们的序列名称
                        if '::' in header:
                            current_name = header.split('::')[0]
                        elif ':' in header:
                            current_name = header.split(':')[0]
                        else:
                            current_name = header.split()[0]
                        
                        current_seq = []
                    else:
                        current_seq.append(line.strip())
                
                # 处理最后一个序列
                if current_name is not None and current_name in region_map:
                    idx = region_map[current_name]
                    sequences[idx] = ''.join(current_seq)
            
            # 组装结果
            results = []
            for i, region in enumerate(regions):
                if i in sequences:
                    results.append({
                        **region,
                        'sequence': sequences[i]
                    })
                else:
                    logger.warning(f"Missing sequence for region {i}: {region['chrom']}:{region['start']}-{region['end']}")
                    results.append({
                        **region,
                        'sequence': ''
                    })
            
            # 清理临时文件
            os.unlink(bed_path)
            os.unlink(out_path)
            
            return results
            
        except Exception as e:
            logger.error(f"bedtools extraction failed: {e}")
            import traceback
            logger.debug(f"Traceback: {traceback.format_exc()}")
            return self._extract_with_samtools(regions)
    
    def _extract_with_samtools(self, regions: List[Dict]) -> List[Dict]:
        """使用samtools faidx批量提取（更稳定）"""
        sequences = []
        
        try:
            for region in regions:
                # samtools使用1-based坐标
                region_spec = f"{region['chrom']}:{region['start']+1}-{region['end']}"
                
                cmd = ['samtools', 'faidx', self.genome_file, region_spec]
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    logger.warning(f"samtools failed for {region_spec}: {result.stderr}")
                    sequences.append({
                        **region,
                        'sequence': ''
                    })
                    continue
                
                # 解析FASTA输出
                seq_lines = []
                for line in result.stdout.split('\n'):
                    if not line.startswith('>') and line:
                        seq_lines.append(line.strip())
                
                seq = ''.join(seq_lines)
                
                # 处理负链
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
        
        try:
            with open(self.genome_file, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    genome_dict[record.id] = str(record.seq)
        except Exception as e:
            logger.error(f"Failed to load genome: {e}")
            return [{**r, 'sequence': ''} for r in regions]
        
        for region in regions:
            chrom = region['chrom']
            if chrom in genome_dict:
                # 调整坐标确保不超出序列长度
                chrom_seq = genome_dict[chrom]
                start = max(0, region['start'])
                end = min(len(chrom_seq), region['end'])
                
                if start < end:
                    seq = chrom_seq[start:end]
                else:
                    seq = ''
                
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