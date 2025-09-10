#!/usr/bin/env python3
"""
Legacy Data Adapter: 将旧版Phase1结果转换为新架构兼容格式
用于从旧的phase1结果无缝迁移到新的数据流设计
"""

import logging
import pickle
import json
from pathlib import Path
from typing import Dict, List, Any, Optional

logger = logging.getLogger(__name__)


class LegacyPhase1Adapter:
    """将旧版Phase1结果转换为新架构格式的适配器"""
    
    def __init__(self, old_phase1_result_path: str):
        """
        初始化适配器
        
        Args:
            old_phase1_result_path: 旧Phase1结果文件路径（.pkl或.json）
        """
        self.old_result_path = Path(old_phase1_result_path)
        self.old_data = None
        
    def load_old_results(self) -> bool:
        """加载旧的Phase1结果"""
        try:
            if self.old_result_path.suffix == '.pkl':
                with open(self.old_result_path, 'rb') as f:
                    self.old_data = pickle.load(f)
                logger.info(f"Loaded legacy Phase1 results from pickle file: {self.old_result_path}")
            elif self.old_result_path.suffix == '.json':
                with open(self.old_result_path, 'r') as f:
                    self.old_data = json.load(f)
                logger.info(f"Loaded legacy Phase1 results from JSON file: {self.old_result_path}")
            else:
                logger.error(f"Unsupported file format: {self.old_result_path.suffix}")
                return False
                
            return True
            
        except Exception as e:
            logger.error(f"Failed to load legacy Phase1 results: {e}")
            return False
    
    def convert_to_new_format(self) -> Dict[str, Any]:
        """将旧格式转换为新的Phase1输出格式"""
        if not self.old_data:
            raise ValueError("No old data loaded. Call load_old_results() first.")
        
        logger.info("Converting legacy Phase1 results to new format...")
        
        # 提取旧格式的数据
        old_format = self._detect_old_format()
        
        if old_format == "classic":
            return self._convert_classic_format()
        elif old_format == "chimera_split":
            return self._convert_chimera_split_format()
        else:
            raise ValueError(f"Unknown legacy format: {old_format}")
    
    def _detect_old_format(self) -> str:
        """检测旧数据的格式类型"""
        if isinstance(self.old_data, dict):
            # 检查是否有嵌合体分割的旧格式
            if 'chimeric_sequences' in self.old_data:
                return "chimera_split"
            # 检查是否有A/B/C分类的经典格式
            elif 'a_sequences' in self.old_data:
                return "classic"
        
        logger.warning("Unknown legacy format, assuming classic format")
        return "classic"
    
    def _convert_classic_format(self) -> Dict[str, Any]:
        """转换经典的A/B/C分类格式"""
        logger.info("Converting classic A/B/C classification format")
        
        # 提取旧格式数据
        a_sequences = self.old_data.get('a_sequences', [])
        b_sequences = self.old_data.get('b_sequences', [])
        c_sequences = self.old_data.get('c_sequences', [])
        scores = self.old_data.get('scores', {})
        rm_results = self.old_data.get('rm_detailed_results', {})
        
        # 准备候选序列（A类+B类）
        consensus_candidates = []
        filtered_out_sequences = []
        
        # 处理A类序列
        for seq in a_sequences:
            candidate = self._prepare_candidate_sequence(seq, 'A', scores, rm_results)
            if candidate:
                consensus_candidates.append(candidate)
        
        # 处理B类序列
        for seq in b_sequences:
            candidate = self._prepare_candidate_sequence(seq, 'B', scores, rm_results)
            if candidate:
                consensus_candidates.append(candidate)
        
        # C类序列作为过滤掉的低质量序列
        for seq in c_sequences:
            seq_copy = dict(seq)
            seq_copy['quality_class'] = 'C'
            seq_copy['filter_reason'] = 'low_quality_c_class'
            filtered_out_sequences.append(seq_copy)
        
        logger.info(f"Conversion complete:")
        logger.info(f"  - Consensus candidates: {len(consensus_candidates)} (A:{len(a_sequences)}, B:{len(b_sequences)})")
        logger.info(f"  - Filtered out: {len(filtered_out_sequences)} (C-class)")
        
        # 生成新格式的输出
        return {
            'consensus_candidates': consensus_candidates,
            'filtered_out_sequences': filtered_out_sequences,
            'c_class_sequences': c_sequences,
            'scores': scores,
            'rm_detailed_results': rm_results,
            'summary': f"Converted from legacy: Candidates:{len(consensus_candidates)}, Filtered:{len(filtered_out_sequences)}, C-class:{len(c_sequences)}",
            'conversion_source': 'legacy_classic_format'
        }
    
    def _convert_chimera_split_format(self) -> Dict[str, Any]:
        """转换包含嵌合体信息的旧格式"""
        logger.info("Converting chimera-aware legacy format")
        
        # 提取数据
        a_sequences = self.old_data.get('a_sequences', [])
        b_sequences = self.old_data.get('b_sequences', [])
        c_sequences = self.old_data.get('c_sequences', [])
        chimeric_sequences = self.old_data.get('chimeric_sequences', [])
        scores = self.old_data.get('scores', {})
        rm_results = self.old_data.get('rm_detailed_results', {})
        
        # 合并所有需要处理的序列
        all_candidates = []
        
        # A/B类正常序列
        for seq in a_sequences + b_sequences:
            quality_class = 'A' if seq in a_sequences else 'B'
            candidate = self._prepare_candidate_sequence(seq, quality_class, scores, rm_results)
            if candidate:
                all_candidates.append(candidate)
        
        # 嵌合体序列也作为候选序列
        for seq in chimeric_sequences:
            candidate = self._prepare_candidate_sequence(seq, 'B', scores, rm_results)  # 嵌合体默认B类
            if candidate:
                candidate['is_potential_chimera'] = True  # 标记为潜在嵌合体
                all_candidates.append(candidate)
        
        # C类作为过滤序列
        filtered_out_sequences = []
        for seq in c_sequences:
            seq_copy = dict(seq)
            seq_copy['quality_class'] = 'C'
            seq_copy['filter_reason'] = 'low_quality_c_class'
            filtered_out_sequences.append(seq_copy)
        
        logger.info(f"Conversion complete:")
        logger.info(f"  - Total candidates: {len(all_candidates)} (including {len(chimeric_sequences)} potential chimeras)")
        logger.info(f"  - Filtered out: {len(filtered_out_sequences)} (C-class)")
        
        return {
            'consensus_candidates': all_candidates,
            'filtered_out_sequences': filtered_out_sequences,
            'c_class_sequences': c_sequences,
            'scores': scores,
            'rm_detailed_results': rm_results,
            'summary': f"Converted from legacy: Candidates:{len(all_candidates)}, Filtered:{len(filtered_out_sequences)}, C-class:{len(c_sequences)}",
            'conversion_source': 'legacy_chimera_format'
        }
    
    def _prepare_candidate_sequence(self, seq: Dict, quality_class: str, 
                                  scores: Dict, rm_results: Dict) -> Optional[Dict]:
        """准备候选序列，添加必要的字段"""
        try:
            seq_id = seq.get('id', 'unknown')
            sequence = seq.get('sequence', '')
            
            if not sequence:
                logger.warning(f"Skipping sequence with empty content: {seq_id}")
                return None
            
            # 创建候选序列
            candidate = dict(seq)  # 复制原始数据
            
            # 添加质量分类
            candidate['quality_class'] = quality_class
            
            # 添加RepeatMasker结果
            if seq_id in rm_results:
                rm_data = rm_results[seq_id]
                candidate['rm_hits'] = rm_data.get('hits', [])
                candidate['rm_coverage'] = rm_data.get('coverage', 0)
            else:
                candidate['rm_hits'] = []
                candidate['rm_coverage'] = 0
            
            # 添加评分信息
            if seq_id in scores:
                score_data = scores[seq_id]
                candidate['copy_number'] = score_data.get('copy_number', 0)
                candidate['avg_identity'] = score_data.get('avg_identity', 0)
                candidate['final_score'] = score_data.get('final', 0)
            
            return candidate
            
        except Exception as e:
            logger.warning(f"Failed to prepare candidate sequence {seq.get('id', 'unknown')}: {e}")
            return None
    
    def save_converted_results(self, output_path: str, converted_data: Dict[str, Any]):
        """保存转换后的结果"""
        output_file = Path(output_path)
        
        try:
            if output_file.suffix == '.pkl':
                with open(output_file, 'wb') as f:
                    pickle.dump(converted_data, f)
            elif output_file.suffix == '.json':
                with open(output_file, 'w') as f:
                    json.dump(converted_data, f, indent=2, default=str)
            else:
                # 默认保存为pickle
                output_file = output_file.with_suffix('.pkl')
                with open(output_file, 'wb') as f:
                    pickle.dump(converted_data, f)
            
            logger.info(f"Saved converted results to: {output_file}")
            return str(output_file)
            
        except Exception as e:
            logger.error(f"Failed to save converted results: {e}")
            raise


def convert_legacy_phase1(input_path: str, output_path: str = None) -> str:
    """
    便捷函数：转换旧的Phase1结果
    
    Args:
        input_path: 旧Phase1结果文件路径
        output_path: 输出文件路径（可选）
    
    Returns:
        转换后的文件路径
    """
    adapter = LegacyPhase1Adapter(input_path)
    
    if not adapter.load_old_results():
        raise RuntimeError("Failed to load legacy Phase1 results")
    
    converted_data = adapter.convert_to_new_format()
    
    if not output_path:
        input_file = Path(input_path)
        output_path = input_file.parent / f"{input_file.stem}_converted_for_new_pipeline{input_file.suffix}"
    
    return adapter.save_converted_results(output_path, converted_data)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Convert legacy Phase1 results to new pipeline format")
    parser.add_argument("input", help="Legacy Phase1 result file (.pkl or .json)")
    parser.add_argument("-o", "--output", help="Output file path (optional)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    
    try:
        output_file = convert_legacy_phase1(args.input, args.output)
        print(f"✅ Conversion successful!")
        print(f"📁 Converted file: {output_file}")
        print(f"🚀 You can now use this file with the new pipeline")
        
    except Exception as e:
        print(f"❌ Conversion failed: {e}")
        exit(1)