#!/usr/bin/env python3
import logging
import argparse
import json
import multiprocessing as mp
from pathlib import Path
from typing import Dict, Any

# Set multiprocessing start method once at module load, before any forking.
# 'spawn' is safer for subprocess-heavy pipelines (avoids fork+exec issues).
try:
    mp.set_start_method('spawn', force=True)
except RuntimeError:
    pass  # Already set by another module or previous call

from config import PipelineConfig
from phase1_screening_optimized import SequenceScreenerOptimized
from phase2_chimera_splitting import ChimeraSplitter
from phase3_consensus_building import ConsensusBuilder
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
    """
    TE共识序列构建主流程
    
    重构后的清晰工作流：
    Phase 1: 识别共识构建候选序列
    Phase 2: 嵌合体检测和拆分
    Phase 3: 全面的共识序列构建
    Phase 4: 基因组屏蔽为RECON做准备
    """
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.runner = RobustRunner(config)

        logger.info("Initializing restructured TE Consensus Pipeline:")
        logger.info("  Phase 1: Candidate identification")
        logger.info("  Phase 2: Chimera detection & splitting")
        logger.info("  Phase 3: Comprehensive consensus building")
        logger.info("  Phase 4: Genome masking preparation")

        # Build genome BLAST database for Phase 3 BLASTN-based copy recruitment
        self._setup_genome_blast_db(config)

        self.phase1 = SequenceScreenerOptimized(config)
        self.phase2 = ChimeraSplitter(config)
        self.phase3 = ConsensusBuilder(config)
        self.phase4 = GenomeMasker(config)

        # 创建必要的目录
        for dir_path in [config.output_dir, config.temp_dir,
                        config.cache_dir, config.checkpoint_dir]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    def _setup_genome_blast_db(self, config: PipelineConfig):
        """Build genome BLAST database and detect RMBlast + scoring matrices
        for Phase 3 copy recruitment with RepeatMasker-level sensitivity."""
        import subprocess
        import os
        import shutil

        # --- Step 1: Detect rmblastn and scoring matrices ---
        self._detect_rmblast_environment(config)

        # --- Step 2: Build genome BLAST database ---
        genome_file = config.genome_file
        db_path = os.path.join(config.temp_dir or "temp_work", "genome_blastdb", "genome")
        db_dir = os.path.dirname(db_path)

        # Check if DB already exists
        if os.path.exists(db_path + ".nhr") or os.path.exists(db_path + ".nsq"):
            logger.info(f"Genome BLAST database already exists: {db_path}")
            config.genome_blast_db = db_path
            return

        try:
            os.makedirs(db_dir, exist_ok=True)
            cmd = [
                config.makeblastdb_exe,
                "-in", genome_file,
                "-dbtype", "nucl",
                "-out", db_path,
                "-parse_seqids"
            ]
            logger.info(f"Building genome BLAST database: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

            if result.returncode == 0:
                config.genome_blast_db = db_path
                logger.info(f"Genome BLAST database built successfully: {db_path}")
            else:
                logger.warning(f"makeblastdb failed (rc={result.returncode}): {result.stderr}")
                logger.warning("Copy recruitment will fall back to RepeatMasker")
                config.genome_blast_db = ""
        except Exception as e:
            logger.warning(f"Failed to build genome BLAST database: {e}")
            logger.warning("Copy recruitment will fall back to RepeatMasker")
            config.genome_blast_db = ""

    def _detect_rmblast_environment(self, config: PipelineConfig):
        """Detect rmblastn binary and RepeatMasker scoring matrices.

        RMBlast + RM scoring matrices provide RepeatMasker-level sensitivity
        for copy recruitment, unlike generic blastn which lacks:
        - Divergence-specific scoring matrices (14p/18p/20p/25p)
        - Complexity-adjusted scoring (-complexity_adjust)
        - Appropriate gap penalties for TE detection
        """
        import shutil
        import os

        # Detect rmblastn
        rmblastn_path = shutil.which(config.rmblastn_exe)
        if rmblastn_path:
            config.rmblastn_exe = rmblastn_path
            logger.info(f"Found rmblastn: {rmblastn_path}")
        else:
            logger.info("rmblastn not in PATH, copy recruitment will use standard blastn (lower sensitivity)")

        # Detect RepeatMasker scoring matrix directory
        # Search common locations relative to RepeatMasker installation
        rm_path = shutil.which(config.repeatmasker_exe)
        candidate_dirs = []
        if rm_path:
            rm_dir = os.path.dirname(os.path.realpath(rm_path))
            # RepeatMasker installed via conda: share/RepeatMasker/Matrices/ncbi/nt/
            for base in [rm_dir, os.path.join(rm_dir, '..', 'share', 'RepeatMasker')]:
                candidate_dirs.append(os.path.join(base, 'Matrices', 'ncbi', 'nt'))
                candidate_dirs.append(os.path.join(base, 'Libraries', 'Matrices', 'ncbi', 'nt'))

        # Also check conda env
        conda_prefix = os.environ.get('CONDA_PREFIX', '')
        if conda_prefix:
            candidate_dirs.append(os.path.join(conda_prefix, 'share', 'RepeatMasker', 'Matrices', 'ncbi', 'nt'))

        for d in candidate_dirs:
            d = os.path.normpath(d)
            test_matrix = os.path.join(d, '20p41g.matrix')
            if os.path.isfile(test_matrix):
                config.rmblast_matrix_dir = d
                logger.info(f"Found RepeatMasker scoring matrices: {d}")
                return

        logger.info("RepeatMasker scoring matrices not found; rmblastn will use default matrix")

    def run(self) -> Dict[str, Any]:
        """执行完整流程"""
        logger.info("="*60)
        logger.info("Starting TE Consensus Pipeline")
        logger.info(f"Input: {self.config.repeatscout_file}")
        logger.info(f"Genome: {self.config.genome_file}")
        logger.info(f"Output: {self.config.output_dir}")
        logger.info("="*60)
        
        try:
            # Phase 1: 候选序列识别
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
            phase3_stats = phase3_output.get('statistics', {})
            logger.info(f"Phase 3 complete: Generated {phase3_stats.get('total_consensus', 0)} consensus sequences")
            logger.info(f"  - High quality: {phase3_stats.get('high_quality_consensus', 0)}")
            logger.info(f"  - Medium quality: {phase3_stats.get('medium_quality_consensus', 0)}")
            logger.info(f"  - Low quality: {phase3_stats.get('low_quality_consensus', 0)}")
            
            # Phase 4: 基因组Masking（可选）
            # Changed default to False to avoid RepeatMasker issues
            if self.config.__dict__.get('enable_masking', False):
                logger.info("\n>>> Phase 4: Genome Masking")
                phase4_output = self.runner.run_with_checkpoint(
                    lambda: self.phase4.run(phase1_output, phase3_output),
                    checkpoint_name="phase4_complete"
                )
                logger.info(f"Phase 4 complete: Generated masked genome")
                
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
        try:
            temp_dir = Path(self.config.temp_dir)
            if temp_dir.exists():
                import shutil
                shutil.rmtree(temp_dir)
            logger.info("Cleanup complete")
        except Exception as e:
            logger.warning(f"Cleanup failed (non-fatal): {e}")

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
                       help='(Deprecated) Checkpoints are now always kept')
    parser.add_argument('--resume', action='store_true',
                       help='Resume from checkpoints if available')
    parser.add_argument('--clear-cache', action='store_true',
                       help='Clear cache before running')
    parser.add_argument('--repeatmasker-quick', action='store_true',
                       help='Enable RepeatMasker quick mode (-q flag)')
    parser.add_argument('--enable-masking', action='store_true',
                       help='Enable Phase 4 genome masking (disabled by default)')
    
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
    # keep_checkpoints is deprecated - checkpoints are always kept
    config.__dict__['repeatmasker_quick'] = args.repeatmasker_quick
    config.__dict__['enable_masking'] = args.enable_masking
    
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
