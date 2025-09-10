"""
染色体ID映射工具
处理不同格式的染色体名称之间的转换
"""

import logging
import re
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

class ChromosomeIDMapper:
    """染色体ID映射器，处理不同格式间的转换"""
    
    def __init__(self, genome_file: str = None, repeatmasker_dir: str = None):
        """
        Args:
            genome_file: 基因组文件路径，用于自动检测可用的染色体ID
            repeatmasker_dir: RepeatMasker工作目录，用于解析真实的ID映射
        """
        self.genome_file = genome_file
        self.repeatmasker_dir = repeatmasker_dir
        self.available_ids = set()
        self.id_mappings = {}
        self.repeatmasker_resolver = None
        
        if genome_file:
            self._load_available_ids()
            # 注意：对于已规范化的lcl|格式基因组文件，通常不需要RepeatMasker ID解析器
            # 但保留此功能以支持特殊情况
            try:
                from .repeatmasker_parser import RepeatMaskerIDResolver
                self.repeatmasker_resolver = RepeatMaskerIDResolver(genome_file, repeatmasker_dir)
                logger.debug("RepeatMasker ID resolver initialized")  # 降级为debug
            except Exception as e:
                logger.debug(f"RepeatMasker resolver not initialized: {e}")  # 降级为debug
    
    def _load_available_ids(self):
        """从基因组文件加载可用的染色体ID"""
        try:
            from Bio import SeqIO
            
            for record in SeqIO.parse(self.genome_file, "fasta"):
                self.available_ids.add(record.id)
            
            logger.debug(f"Loaded {len(self.available_ids)} chromosome IDs from {self.genome_file}")
            
            # 构建映射关系
            self._build_mappings()
            
        except Exception as e:
            logger.warning(f"Failed to load chromosome IDs from {self.genome_file}: {e}")
    
    def _build_mappings(self):
        """构建不同格式之间的映射关系"""
        for chrom_id in self.available_ids:
            # 为每个ID建立各种可能的映射
            mappings = self._generate_possible_mappings(chrom_id)
            for mapping in mappings:
                self.id_mappings[mapping] = chrom_id
    
    def _generate_possible_mappings(self, chrom_id: str) -> List[str]:
        """为给定的染色体ID生成可能的映射格式"""
        mappings = [chrom_id]  # 原始ID
        
        # 移除常见前缀
        for prefix in ['lcl|', 'ref|', 'gi|']:
            if chrom_id.startswith(prefix):
                base_id = chrom_id[len(prefix):]
                mappings.append(base_id)
                
                # 进一步处理base_id
                mappings.extend(self._process_base_id(base_id))
        
        # 处理其他格式
        mappings.extend(self._process_base_id(chrom_id))
        
        return list(set(mappings))  # 去重
    
    def _process_base_id(self, base_id: str) -> List[str]:
        """处理基础ID，生成可能的变体"""
        mappings = []
        
        # 处理带下划线的ID（如 1_1, 2_3）
        if '_' in base_id:
            # 尝试提取染色体编号
            parts = base_id.split('_')
            if len(parts) >= 2 and parts[0].isdigit():
                # 生成chr格式
                chr_num = parts[0]
                mappings.append(f"chr{chr_num}")
                mappings.append(chr_num)
                
                # 也尝试匹配到带坐标的格式
                for available_id in self.available_ids:
                    if available_id.startswith(f"chr{chr_num}:"):
                        mappings.append(available_id)
        
        # 处理chr格式
        if base_id.startswith('chr'):
            # 移除chr前缀
            num = base_id[3:]
            mappings.append(num)
        elif base_id.isdigit():
            # 数字格式，添加chr前缀
            mappings.append(f"chr{base_id}")
            
            # 也尝试匹配到带坐标的格式
            for available_id in self.available_ids:
                if available_id.startswith(f"chr{base_id}:"):
                    mappings.append(available_id)
        
        # 处理坐标格式（如 chr1:100-200）
        coord_match = re.match(r'(.+):(\d+)-(\d+)', base_id)
        if coord_match:
            chrom_part = coord_match.group(1)
            mappings.append(chrom_part)
            mappings.extend(self._process_base_id(chrom_part))
        
        return mappings
    
    def map_chromosome_id(self, input_id: str) -> Optional[str]:
        """
        将输入的染色体ID映射到基因组文件中实际存在的ID
        
        Args:
            input_id: 需要映射的染色体ID
            
        Returns:
            映射后的染色体ID，如果找不到则返回None
        """
        # 直接匹配
        if input_id in self.available_ids:
            return input_id
        
        # 使用预建的映射
        if input_id in self.id_mappings:
            return self.id_mappings[input_id]
        
        # 对于lcl|格式，优先使用RepeatMasker解析器
        if input_id.startswith('lcl|') and self.repeatmasker_resolver:
            resolved_id = self.repeatmasker_resolver.resolve_id(input_id)
            if resolved_id and resolved_id in self.available_ids:
                self.id_mappings[input_id] = resolved_id
                logger.debug(f"RepeatMasker resolved: {input_id} -> {resolved_id}")
                return resolved_id
        
        # 动态映射尝试
        possible_mappings = self._generate_possible_mappings(input_id)
        for mapping in possible_mappings:
            if mapping in self.available_ids:
                # 缓存这个映射
                self.id_mappings[input_id] = mapping
                return mapping
        
        # 模糊匹配
        fuzzy_match = self._fuzzy_match(input_id)
        if fuzzy_match:
            self.id_mappings[input_id] = fuzzy_match
            return fuzzy_match
        
        logger.warning(f"Cannot map chromosome ID: {input_id}")
        return None
    
    def _fuzzy_match(self, input_id: str) -> Optional[str]:
        """模糊匹配染色体ID"""
        # 提取数字部分进行匹配
        input_nums = re.findall(r'\d+', input_id)
        if not input_nums:
            return None
        
        for available_id in self.available_ids:
            available_nums = re.findall(r'\d+', available_id)
            if input_nums == available_nums:
                logger.debug(f"Fuzzy matched {input_id} -> {available_id}")
                return available_id
        
        return None
    
    def get_mapping_stats(self) -> Dict:
        """获取映射统计信息"""
        return {
            'available_ids': len(self.available_ids),
            'cached_mappings': len(self.id_mappings),
            'sample_available_ids': list(self.available_ids)[:10],
            'sample_mappings': dict(list(self.id_mappings.items())[:10])
        }


def create_id_mapper(genome_file: str) -> ChromosomeIDMapper:
    """创建ID映射器的便捷函数"""
    return ChromosomeIDMapper(genome_file)


def map_repeatmasker_hits(hits: List[Dict], mapper: ChromosomeIDMapper) -> List[Dict]:
    """
    映射RepeatMasker hits中的染色体ID
    
    Args:
        hits: RepeatMasker hits列表
        mapper: ID映射器
        
    Returns:
        映射后的hits列表
    """
    mapped_hits = []
    failed_mappings = 0
    
    for hit in hits:
        original_chrom = hit.get('chrom', '')
        
        # 处理可能包含坐标的染色体名称
        if ':' in original_chrom:
            chrom_part = original_chrom.split(':')[0]
        else:
            chrom_part = original_chrom
        
        # 尝试映射
        mapped_chrom = mapper.map_chromosome_id(chrom_part)
        
        if mapped_chrom:
            # 创建映射后的hit
            mapped_hit = hit.copy()
            mapped_hit['chrom'] = mapped_chrom
            mapped_hit['original_chrom'] = original_chrom
            mapped_hits.append(mapped_hit)
        else:
            failed_mappings += 1
            logger.debug(f"Failed to map chromosome: {original_chrom}")
    
    if failed_mappings > 0:
        logger.warning(f"Failed to map {failed_mappings}/{len(hits)} chromosome IDs")
    
    return mapped_hits