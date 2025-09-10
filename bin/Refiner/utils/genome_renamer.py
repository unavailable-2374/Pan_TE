#!/usr/bin/env python3
"""
基因组序列重命名工具
在RepeatMasker前重命名为简单ID (chr1-chrN)，完成后还原
"""

import logging
import json
import os
from pathlib import Path
from typing import Dict, Tuple, List
from Bio import SeqIO

logger = logging.getLogger(__name__)

class GenomeRenamer:
    """基因组序列重命名器"""
    
    def __init__(self, output_dir: str = None):
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        self.mapping_file = self.output_dir / "genome_id_mapping.json"
        self.original_genome = None
        self.renamed_genome = None
        self.id_mapping = {}  # {new_id: original_id}
        self.reverse_mapping = {}  # {original_id: new_id}
        
    def rename_genome_for_repeatmasker(self, input_genome: str) -> Tuple[str, str]:
        """
        重命名基因组文件为RepeatMasker友好的ID格式
        
        Args:
            input_genome: 原始基因组文件路径
            
        Returns:
            (重命名后的基因组文件路径, 映射文件路径)
        """
        self.original_genome = Path(input_genome)
        self.renamed_genome = self.output_dir / f"renamed_{self.original_genome.name}"
        
        logger.info(f"Renaming genome sequences for RepeatMasker compatibility")
        logger.info(f"Input: {input_genome}")
        logger.info(f"Output: {self.renamed_genome}")
        
        # 读取原始基因组并重命名
        renamed_count = 0
        with open(self.renamed_genome, 'w') as out_file:
            for i, record in enumerate(SeqIO.parse(input_genome, "fasta"), 1):
                original_id = record.id
                new_id = f"chr{i}"
                
                # 记录映射关系
                self.id_mapping[new_id] = original_id
                self.reverse_mapping[original_id] = new_id
                
                # 写入重命名的序列
                out_file.write(f">{new_id}\n")
                
                # 分行写入序列（每行80个字符）
                sequence = str(record.seq)
                for j in range(0, len(sequence), 80):
                    out_file.write(sequence[j:j+80] + "\n")
                
                renamed_count += 1
                
                if renamed_count % 100 == 0:
                    logger.info(f"Renamed {renamed_count} sequences...")
        
        # 保存映射关系
        self._save_mapping()
        
        logger.info(f"Successfully renamed {renamed_count} sequences")
        logger.info(f"Mapping saved to: {self.mapping_file}")
        
        return str(self.renamed_genome), str(self.mapping_file)
    
    def _save_mapping(self):
        """保存ID映射关系到文件"""
        mapping_data = {
            'original_genome': str(self.original_genome),
            'renamed_genome': str(self.renamed_genome),
            'total_sequences': len(self.id_mapping),
            'id_mapping': self.id_mapping,  # {new_id: original_id}
            'reverse_mapping': self.reverse_mapping  # {original_id: new_id}
        }
        
        with open(self.mapping_file, 'w') as f:
            json.dump(mapping_data, f, indent=2)
        
        logger.debug(f"Saved mapping for {len(self.id_mapping)} sequences")
    
    def load_mapping(self, mapping_file: str = None):
        """加载ID映射关系"""
        mapping_path = Path(mapping_file) if mapping_file else self.mapping_file
        
        if not mapping_path.exists():
            raise FileNotFoundError(f"Mapping file not found: {mapping_path}")
        
        with open(mapping_path, 'r') as f:
            mapping_data = json.load(f)
        
        self.original_genome = Path(mapping_data['original_genome'])
        self.renamed_genome = Path(mapping_data['renamed_genome'])
        self.id_mapping = mapping_data['id_mapping']
        self.reverse_mapping = mapping_data['reverse_mapping']
        
        logger.info(f"Loaded mapping for {len(self.id_mapping)} sequences")
    
    def restore_genome_ids(self, masked_genome_file: str, output_file: str = None) -> str:
        """
        将masked基因组的ID还原为原始ID
        
        Args:
            masked_genome_file: RepeatMasker处理后的基因组文件
            output_file: 输出文件路径（可选）
            
        Returns:
            还原后的基因组文件路径
        """
        if not self.id_mapping:
            raise ValueError("No ID mapping loaded. Call load_mapping() first.")
        
        masked_path = Path(masked_genome_file)
        if not output_file:
            output_file = masked_path.parent / f"restored_{masked_path.name}"
        
        output_path = Path(output_file)
        
        logger.info(f"Restoring genome IDs from RepeatMasker output")
        logger.info(f"Input: {masked_genome_file}")
        logger.info(f"Output: {output_path}")
        
        restored_count = 0
        with open(output_path, 'w') as out_file:
            for record in SeqIO.parse(masked_genome_file, "fasta"):
                current_id = record.id
                
                # 查找原始ID
                original_id = self.id_mapping.get(current_id)
                if not original_id:
                    logger.warning(f"No mapping found for ID: {current_id}, keeping as is")
                    original_id = current_id
                
                # 写入还原的序列
                out_file.write(f">{original_id}\n")
                
                # 分行写入序列
                sequence = str(record.seq)
                for i in range(0, len(sequence), 80):
                    out_file.write(sequence[i:i+80] + "\n")
                
                restored_count += 1
                
                if restored_count % 100 == 0:
                    logger.info(f"Restored {restored_count} sequences...")
        
        logger.info(f"Successfully restored {restored_count} sequences to original IDs")
        return str(output_path)
    
    def get_mapping_stats(self) -> Dict:
        """获取映射统计信息"""
        return {
            'total_sequences': len(self.id_mapping),
            'original_genome': str(self.original_genome) if self.original_genome else None,
            'renamed_genome': str(self.renamed_genome) if self.renamed_genome else None,
            'mapping_file': str(self.mapping_file),
            'sample_mappings': dict(list(self.id_mapping.items())[:5])
        }

class MaskedGenomeRestorer:
    """专门处理masked基因组的ID还原"""
    
    def __init__(self, mapping_file: str):
        self.renamer = GenomeRenamer()
        self.renamer.load_mapping(mapping_file)
    
    def restore_masked_files(self, masked_dir: str, output_dir: str = None) -> Dict[str, str]:
        """
        还原目录中所有masked文件的ID
        
        Args:
            masked_dir: 包含masked文件的目录
            output_dir: 输出目录
            
        Returns:
            {original_file: restored_file} 映射
        """
        masked_path = Path(masked_dir)
        if not output_dir:
            output_dir = masked_path.parent / "restored_genomes"
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        restored_files = {}
        
        # 查找所有可能的masked文件
        masked_patterns = [
            "*.fa.masked", "*.fasta.masked", "*_masked.fa", "*_masked.fasta",
            "*.fa.out", "genome_final_masked.fa", "consensus_masking.fa"
        ]
        
        all_files = []
        for pattern in masked_patterns:
            all_files.extend(masked_path.glob(pattern))
        
        logger.info(f"Found {len(all_files)} masked files to restore")
        
        for masked_file in all_files:
            try:
                output_file = output_path / f"restored_{masked_file.name}"
                restored_file = self.renamer.restore_genome_ids(
                    str(masked_file), str(output_file)
                )
                restored_files[str(masked_file)] = restored_file
                logger.info(f"Restored: {masked_file.name} → {Path(restored_file).name}")
            except Exception as e:
                logger.error(f"Failed to restore {masked_file}: {e}")
        
        return restored_files

def integrate_with_repeatmasker_workflow(genome_file: str, work_dir: str) -> Tuple[str, str]:
    """
    集成到RepeatMasker工作流的便捷函数
    
    Args:
        genome_file: 原始基因组文件
        work_dir: 工作目录
        
    Returns:
        (重命名后的基因组文件, 映射文件路径)
    """
    renamer = GenomeRenamer(work_dir)
    return renamer.rename_genome_for_repeatmasker(genome_file)

def restore_after_phase4(mapping_file: str, phase4_output_dir: str, 
                        final_output_dir: str = None) -> Dict[str, str]:
    """
    Phase4完成后还原所有输出文件的ID
    
    Args:
        mapping_file: ID映射文件路径
        phase4_output_dir: Phase4输出目录
        final_output_dir: 最终输出目录
        
    Returns:
        还原后的文件映射
    """
    restorer = MaskedGenomeRestorer(mapping_file)
    return restorer.restore_masked_files(phase4_output_dir, final_output_dir)

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python genome_renamer.py rename <genome_file> [output_dir]")
        print("  python genome_renamer.py restore <masked_file> <mapping_file> [output_file]")
        sys.exit(1)
    
    action = sys.argv[1]
    
    if action == "rename":
        genome_file = sys.argv[2]
        output_dir = sys.argv[3] if len(sys.argv) > 3 else "."
        
        renamer = GenomeRenamer(output_dir)
        renamed_file, mapping_file = renamer.rename_genome_for_repeatmasker(genome_file)
        
        print(f"Renamed genome: {renamed_file}")
        print(f"Mapping file: {mapping_file}")
        print(f"Stats: {renamer.get_mapping_stats()}")
    
    elif action == "restore":
        if len(sys.argv) < 4:
            print("Error: restore requires <masked_file> <mapping_file>")
            sys.exit(1)
        
        masked_file = sys.argv[2]
        mapping_file = sys.argv[3]
        output_file = sys.argv[4] if len(sys.argv) > 4 else None
        
        renamer = GenomeRenamer()
        renamer.load_mapping(mapping_file)
        restored_file = renamer.restore_genome_ids(masked_file, output_file)
        
        print(f"Restored genome: {restored_file}")
    
    else:
        print(f"Unknown action: {action}")
        sys.exit(1)