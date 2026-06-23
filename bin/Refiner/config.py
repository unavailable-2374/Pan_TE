from dataclasses import dataclass
from pathlib import Path
import json

@dataclass
class PipelineConfig:
    """Pipeline配置参数"""
    
    # 输入输出
    input_file: str
    genome_file: str
    output_dir: str
    temp_dir: str = "temp_work"
    cache_dir: str = "cache"
    checkpoint_dir: str = "checkpoints"

    # Phase 1 parameters
    min_length: int = 50
    max_length: int = 50000  # Maximum length limit
    max_n_percent: float = 0.2
    dust_window: int = 64
    dust_threshold: int = 7

    # Phase 1 sensitivity parameters
    # Trust the repeat finder's initial filtering
    trust_input: bool = True
    # Rescue C-class sequences that have biological evidence
    rescue_c_class: bool = True
    # Relaxed identity threshold for ancient TEs (from 50 to 40)
    low_identity_threshold: float = 40.0

    # Phase 2 parameters
    identity_threshold: float = 0.85
    coverage_threshold: float = 0.60
    max_recruits_per_family: int = 30
    msa_algorithm: str = "localpair"  # MAFFT L-INS-i
    column_coverage_threshold: float = 0.6

    # Phase 3 parameters
    redundancy_threshold_masking: float = 0.95
    redundancy_threshold_analysis: float = 0.90
    # Reduced from 5 to 2 - trust the repeat finder's detection
    min_copy_number: int = 2
    min_consensus_quality: float = 0.85
    skip_phase3_redundancy_removal: bool = True  # Skip Phase 3 CD-HIT deduplication
    
    # 性能参数
    threads: int = 8
    batch_size: int = 100
    use_parallel: bool = True
    max_retries: int = 3
    retry_backoff: int = 2
    
    # 外部工具路径
    repeatmasker_exe: str = "RepeatMasker"
    mafft_exe: str = "mafft"
    blastn_exe: str = "blastn"
    rmblastn_exe: str = "rmblastn"
    makeblastdb_exe: str = "makeblastdb"
    cdhit_exe: str = "cd-hit-est"
    rmblast_matrix_dir: str = ""  # Auto-detected from RepeatMasker installation
    
    # RepeatMasker优化参数
    repeatmasker_quick: bool = False  # 快速模式，自动基于基因组大小设置
    large_genome_threshold: int = 1073741824  # 1GB阈值，超过此大小启用快速模式

    # Stratified processing parameters (Phase 3 optimization)
    enable_stratified_processing: bool = True
    fast_path_min_copies: int = 10
    fast_path_min_identity: float = 75.0
    minimal_path_max_copies: int = 1

    # BLASTN-based copy recruitment (Phase 3 optimization)
    genome_blast_db: str = ""  # Set at runtime by main.py
    
    def save(self, filepath: str):
        """保存配置到文件"""
        with open(filepath, 'w') as f:
            json.dump(self.__dict__, f, indent=2)
    
    @classmethod
    def load(cls, filepath: str):
        """从文件加载配置"""
        with open(filepath, 'r') as f:
            return cls(**json.load(f))