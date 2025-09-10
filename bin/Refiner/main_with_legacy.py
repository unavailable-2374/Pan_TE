#!/usr/bin/env python3
import logging
import argparse
import json
import pickle
from pathlib import Path
from typing import Dict, Any

from config import PipelineConfig
from phase1_screening_optimized import SequenceScreenerOptimized
from phase2_chimera_splitting import ChimeraSplitter
from phase3_consensus_building import ConsensusBuilder
from phase4_genome_masking import GenomeMasker
from utils.robust_runner import RobustRunner
from legacy_adapter import LegacyPhase1Adapter, convert_legacy_phase1

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pipeline.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class TEConsensusPipelineWithLegacy:
    """
    支持加载旧Phase1结果的TE共识序列构建流水线
    
    用法:
    1. 从头开始运行: python main_with_legacy.py -r repeatscout.fa -g genome.fa -o output
    2. 从旧Phase1结果继续: python main_with_legacy.py --load-phase1 old_phase1.pkl -g genome.fa -o output
    """
    
    def __init__(self, config: PipelineConfig, phase1_result_path: str = None):
        self.config = config
        self.runner = RobustRunner(config)
        self.phase1_result_path = phase1_result_path
        
        logger.info("Initializing TE Consensus Pipeline with Legacy Support:")
        logger.info("  Phase 1: Candidate identification (can load legacy results)")
        logger.info("  Phase 2: Chimera detection & splitting")
        logger.info("  Phase 3: Comprehensive consensus building")
        logger.info("  Phase 4: Genome masking preparation")
        
        # 只有在不加载旧结果时才初始化Phase1
        if not phase1_result_path:
            self.phase1 = SequenceScreenerOptimized(config)
        else:
            self.phase1 = None
            logger.info(f"Will load Phase1 results from: {phase1_result_path}")
        
        self.phase2 = ChimeraSplitter(config)
        self.phase3 = ConsensusBuilder(config)
        self.phase4 = GenomeMasker(config)
        
        # 创建必要的目录
        for dir_path in [config.output_dir, config.temp_dir, 
                        config.cache_dir, config.checkpoint_dir]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    def run(self) -> Dict[str, Any]:
        """执行完整流程"""
        logger.info("="*60)
        logger.info("Starting TE Consensus Pipeline with Legacy Support")
        if self.phase1_result_path:
            logger.info(f"Loading Phase1 results from: {self.phase1_result_path}")
        else:
            logger.info(f"Input: {self.config.repeatscout_file}")
        logger.info(f"Genome: {self.config.genome_file}")
        logger.info(f"Output: {self.config.output_dir}")
        logger.info("="*60)
        
        try:
            # Phase 1: 获取候选序列（加载旧结果或重新运行）
            if self.phase1_result_path:
                logger.info("\n>>> Phase 1: Loading Legacy Results")
                phase1_output = self._load_and_convert_phase1_results()
            else:
                logger.info("\n>>> Phase 1: Consensus Candidate Identification")
                phase1_output = self.runner.run_with_checkpoint(
                    self.phase1.run,
                    checkpoint_name="phase1_complete"
                )
            
            logger.info(f"Phase 1 complete: {phase1_output['summary']}")
            
            # Phase 2: 嵌合体检测和拆分
            logger.info("\n>>> Phase 2: Chimera Detection and Splitting")
            phase2_output = self.runner.run_with_checkpoint(
                lambda: self.phase2.run(phase1_output),
                checkpoint_name="phase2_complete"
            )
            # Phase2统计信息
            phase2_stats = phase2_output.get('statistics', {})
            logger.info(f"Phase 2 complete: Processed {phase2_stats.get('input_candidates', 0)} candidate sequences")
            logger.info(f"  - Chimeric sequences detected: {phase2_stats.get('chimeric_sequences', 0)}")
            logger.info(f"  - Successfully split: {phase2_stats.get('split_sequences', 0)}")
            logger.info(f"  - Output sequences for Phase 3: {phase2_stats.get('output_total', 0)}")
            
            # Phase 3: 全面共识构建
            logger.info("\n>>> Phase 3: Comprehensive Consensus Building")
            phase3_output = self.runner.run_with_checkpoint(
                lambda: self.phase3.run(phase2_output),
                checkpoint_name="phase3_complete"
            )
            # Phase3统计信息
            phase3_stats = phase3_output.get('statistics', {})
            logger.info(f"Phase 3 complete: Generated {phase3_stats.get('total_consensus', 0)} consensus sequences")
            logger.info(f"  - High quality: {phase3_stats.get('high_quality_consensus', 0)}")
            logger.info(f"  - Medium quality: {phase3_stats.get('medium_quality_consensus', 0)}")
            logger.info(f"  - Low quality: {phase3_stats.get('low_quality_consensus', 0)}")
            if 'deduplicated_sequences' in phase3_stats:
                logger.info(f"  - Removed duplicates: {phase3_stats['deduplicated_sequences']}")
            
            # Phase 4: 基因组屏蔽
            logger.info("\n>>> Phase 4: Genome Masking")
            phase4_output = self.runner.run_with_checkpoint(
                lambda: self.phase4.run(phase3_output),
                checkpoint_name="phase4_complete"
            )
            # Phase4统计信息
            phase4_stats = phase4_output.get('statistics', {})
            logger.info(f"Phase 4 complete: Masked genome ready")
            if 'masking_coverage' in phase4_stats:
                logger.info(f"  - Masking coverage: {phase4_stats['masking_coverage']:.2%}")
            
            # 生成最终结果
            final_results = self._compile_final_results(
                phase1_output, phase2_output, phase3_output, phase4_output
            )
            
            logger.info("\n" + "="*60)
            logger.info("TE Consensus Pipeline Completed Successfully!")
            logger.info(f"Final consensus library: {len(phase3_output['consensus_sequences'])} sequences")
            logger.info(f"Masked genome ready for RECON")
            logger.info("="*60)
            
            return final_results
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise
    
    def _load_and_convert_phase1_results(self) -> Dict[str, Any]:
        """加载并转换旧的Phase1结果"""
        try:
            logger.info(f"Loading legacy Phase1 results from {self.phase1_result_path}")
            
            # 检查是否需要转换
            phase1_path = Path(self.phase1_result_path)
            converted_path = phase1_path.parent / f"{phase1_path.stem}_converted{phase1_path.suffix}"
            
            if converted_path.exists():
                logger.info(f"Found converted file: {converted_path}")
                # 直接加载转换后的结果
                if converted_path.suffix == '.pkl':
                    with open(converted_path, 'rb') as f:
                        converted_data = pickle.load(f)
                elif converted_path.suffix == '.json':
                    with open(converted_path, 'r') as f:
                        converted_data = json.load(f)
                else:
                    raise ValueError(f"Unsupported file format: {converted_path}")
            else:
                logger.info("Converting legacy format to new pipeline format...")
                # 使用适配器转换
                adapter = LegacyPhase1Adapter(self.phase1_result_path)
                if not adapter.load_old_results():
                    raise RuntimeError("Failed to load legacy Phase1 results")
                
                converted_data = adapter.convert_to_new_format()
                
                # 保存转换后的结果以备下次使用
                adapter.save_converted_results(str(converted_path), converted_data)
                logger.info(f"Saved converted results for future use: {converted_path}")
            
            # 验证转换后的数据
            self._validate_phase1_data(converted_data)
            
            logger.info("Legacy Phase1 results loaded and converted successfully")
            return converted_data
            
        except Exception as e:
            logger.error(f"Failed to load/convert Phase1 results: {e}")
            raise
    
    def _validate_phase1_data(self, data: Dict[str, Any]):
        """验证Phase1数据格式"""
        required_keys = ['consensus_candidates', 'scores', 'rm_detailed_results']
        for key in required_keys:
            if key not in data:
                raise ValueError(f"Missing required key in Phase1 data: {key}")
        
        candidates = data['consensus_candidates']
        if not isinstance(candidates, list):
            raise ValueError("consensus_candidates must be a list")
        
        # 检查候选序列格式
        for candidate in candidates[:5]:  # 检查前5个
            if 'id' not in candidate or 'sequence' not in candidate:
                raise ValueError("Each candidate must have 'id' and 'sequence' fields")
        
        logger.info(f"Phase1 data validation passed: {len(candidates)} candidates")
    
    def _compile_final_results(self, phase1_output, phase2_output, phase3_output, phase4_output) -> Dict[str, Any]:
        """编译最终结果"""
        return {
            'phase1_summary': phase1_output.get('summary', ''),
            'phase2_statistics': phase2_output.get('statistics', {}),
            'phase3_statistics': phase3_output.get('statistics', {}),
            'phase4_statistics': phase4_output.get('statistics', {}),
            'final_consensus_sequences': phase3_output.get('consensus_sequences', []),
            'masked_genome_path': phase4_output.get('masked_genome_path', ''),
            'pipeline_complete': True
        }


def main():
    parser = argparse.ArgumentParser(description='TE Consensus Building Pipeline with Legacy Support')
    
    # 基本参数
    parser.add_argument('-r', '--repeatscout', help='RepeatScout output file')
    parser.add_argument('-g', '--genome', required=True, help='Genome file (FASTA)')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    
    # Legacy支持
    parser.add_argument('--load-phase1', help='Load legacy Phase1 results (.pkl or .json)')
    
    # 可选参数
    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of threads (default: 8)')
    parser.add_argument('-c', '--config', help='Configuration file (JSON)')
    parser.add_argument('--resume', action='store_true', help='Resume from checkpoints')
    parser.add_argument('--clear-cache', action='store_true', help='Clear cache before running')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # 设置日志级别
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # 验证参数
    if not args.load_phase1 and not args.repeatscout:
        parser.error("Either --repeatscout or --load-phase1 must be specified")
    
    if args.load_phase1 and not Path(args.load_phase1).exists():
        parser.error(f"Phase1 result file not found: {args.load_phase1}")
    
    try:
        # 创建配置
        if args.config and Path(args.config).exists():
            config = PipelineConfig.load(args.config)
            # 更新命令行参数
            config.genome_file = args.genome
            config.output_dir = args.output
            config.threads = args.threads
            if args.repeatscout:
                config.repeatscout_file = args.repeatscout
        else:
            # 如果加载旧结果，repeatscout文件可以为None
            repeatscout_file = args.repeatscout if args.repeatscout else "dummy.fa"
            config = PipelineConfig(
                repeatscout_file=repeatscout_file,
                genome_file=args.genome,
                output_dir=args.output,
                threads=args.threads
            )
        
        # 设置选项
        if args.clear_cache and hasattr(config, 'clear_cache'):
            config.clear_cache = True
        
        # 创建并运行流水线
        pipeline = TEConsensusPipelineWithLegacy(
            config=config,
            phase1_result_path=args.load_phase1
        )
        
        results = pipeline.run()
        
        print("\n✅ Pipeline completed successfully!")
        print(f"📊 Final consensus sequences: {len(results['final_consensus_sequences'])}")
        if results.get('masked_genome_path'):
            print(f"🧬 Masked genome: {results['masked_genome_path']}")
        
    except Exception as e:
        print(f"❌ Pipeline failed: {e}")
        logger.exception("Detailed error information:")
        exit(1)


if __name__ == "__main__":
    main()