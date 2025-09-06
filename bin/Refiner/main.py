#!/usr/bin/env python3
import logging
import argparse
import json
from pathlib import Path
from typing import Dict, Any

from config import PipelineConfig
from phase1_screening_optimized import SequenceScreenerOptimized
from phase2_consensus_expansion import ConsensusExpansionBuilder
from phase3_finalization_relaxed import ConsensusFilterRelaxed
from phase4_genome_masking import GenomeMasker
from utils.robust_runner import RobustRunner

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

class TEConsensusPipeline:
    """TE共识序列构建主流程"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.runner = RobustRunner(config)
        
        # 使用最新的优化版本（Expansion-focused pipeline）
        logger.info("Using expansion-focused pipeline for comprehensive genome annotation")
        self.phase1 = SequenceScreenerOptimized(config)
        self.phase2 = ConsensusExpansionBuilder(config)
        self.phase3 = ConsensusFilterRelaxed(config)
        self.phase4 = GenomeMasker(config)
        
        # 创建必要的目录
        for dir_path in [config.output_dir, config.temp_dir, 
                        config.cache_dir, config.checkpoint_dir]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    def run(self) -> Dict[str, Any]:
        """执行完整流程"""
        logger.info("="*60)
        logger.info("Starting TE Consensus Pipeline")
        logger.info(f"Input: {self.config.repeatscout_file}")
        logger.info(f"Genome: {self.config.genome_file}")
        logger.info(f"Output: {self.config.output_dir}")
        logger.info("="*60)
        
        try:
            # Phase 1: 筛选与评分
            logger.info("\n>>> Phase 1: Screening and Scoring")
            phase1_output = self.runner.run_with_checkpoint(
                self.phase1.run,
                checkpoint_name="phase1_complete"
            )
            logger.info(f"Phase 1 complete: {phase1_output['summary']}")
            
            # Phase 2: 扩展与共识构建
            logger.info("\n>>> Phase 2: Consensus Expansion and Building")
            consensus_list = self.runner.run_with_checkpoint(
                lambda: self.phase2.run(phase1_output),
                checkpoint_name="phase2_complete"
            )
            logger.info(f"Phase 2 complete: {len(consensus_list)} consensus sequences (expanded from {len(phase1_output['a_sequences'])+len([s for s in phase1_output['b_sequences'] if True])} input sequences)")
            
            # Phase 3: 质量控制与输出
            logger.info("\n>>> Phase 3: Quality Control and Output")
            phase3_output = self.runner.run_with_checkpoint(
                lambda: self.phase3.run(consensus_list),
                checkpoint_name="phase3_complete"
            )
            logger.info(f"Phase 3 complete: {phase3_output['summary']}")
            
            # Phase 4: 基因组Masking（可选）
            if self.config.__dict__.get('enable_masking', True):
                logger.info("\n>>> Phase 4: Genome Masking")
                phase4_output = self.runner.run_with_checkpoint(
                    lambda: self.phase4.run(phase1_output, phase3_output),
                    checkpoint_name="phase4_complete"
                )
                logger.info(f"Phase 4 complete: {phase4_output['summary']}")
                
                # 合并输出
                final_output = {
                    **phase3_output,
                    'masking_results': phase4_output
                }
            else:
                logger.info("Phase 4 skipped (masking disabled)")
                final_output = phase3_output
            
            logger.info("\n" + "="*60)
            logger.info("Pipeline completed successfully")
            logger.info("="*60)
            
            # 清理临时文件（可选）
            if not self.config.__dict__.get('keep_temp', False):
                self.cleanup_temp_files()
            
            return final_output
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise
    
    def cleanup_temp_files(self):
        """清理临时文件"""
        logger.info("Cleaning up temporary files...")
        # 清理临时目录
        temp_dir = Path(self.config.temp_dir)
        if temp_dir.exists():
            for file in temp_dir.glob("*"):
                file.unlink()
        
        # 可选：清理检查点
        if not self.config.__dict__.get('keep_checkpoints', False):
            self.runner.cleanup_checkpoints()
        
        logger.info("Cleanup complete")

def main():
    """命令行接口"""
    parser = argparse.ArgumentParser(
        description='TE Consensus Builder - Build high-quality TE consensus sequences'
    )
    
    # 必需参数
    parser.add_argument('-r', '--repeatscout', required=True, 
                       help='RepeatScout output file (filtered)')
    parser.add_argument('-g', '--genome', required=True,
                       help='Reference genome file (FASTA)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory')
    
    # 可选参数
    parser.add_argument('-t', '--threads', type=int, default=8,
                       help='Number of threads (default: 8)')
    parser.add_argument('-c', '--config', 
                       help='Configuration file (JSON)')
    parser.add_argument('--keep-temp', action='store_true',
                       help='Keep temporary files')
    parser.add_argument('--keep-checkpoints', action='store_true',
                       help='Keep checkpoint files')
    parser.add_argument('--resume', action='store_true',
                       help='Resume from checkpoints if available')
    parser.add_argument('--clear-cache', action='store_true',
                       help='Clear cache before running')
    parser.add_argument('--repeatmasker-quick', action='store_true',
                       help='Enable RepeatMasker quick mode (-q flag)')
    
    args = parser.parse_args()
    
    # 加载或创建配置
    if args.config:
        config = PipelineConfig.load(args.config)
        # 覆盖命令行参数
        config.repeatscout_file = args.repeatscout
        config.genome_file = args.genome
        config.output_dir = args.output
        config.threads = args.threads
    else:
        config = PipelineConfig(
            repeatscout_file=args.repeatscout,
            genome_file=args.genome,
            output_dir=args.output,
            threads=args.threads
        )
    
    # 添加额外的配置选项
    config.__dict__['keep_temp'] = args.keep_temp
    config.__dict__['keep_checkpoints'] = args.keep_checkpoints
    config.__dict__['repeatmasker_quick'] = args.repeatmasker_quick
    
    # 确保输出目录存在
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 保存配置供参考
    config_file = output_dir / "pipeline_config.json"
    config.save(str(config_file))
    
    # 清理缓存（如果需要）
    if args.clear_cache:
        from utils.cache_utils import clear_cache
        clear_cache()
    
    # 运行pipeline
    pipeline = TEConsensusPipeline(config)
    results = pipeline.run()
    
    # 打印结果摘要
    print("\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)
    print(f"Output directory: {args.output}")
    print(f"Masking library: {results['masking_library']}")
    print(f"Analysis library: {results['analysis_library']}")
    print(f"Statistics file: {results['statistics']}")
    print("="*60)

if __name__ == '__main__':
    main()
