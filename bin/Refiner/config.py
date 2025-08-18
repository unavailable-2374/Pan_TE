from dataclasses import dataclass
from pathlib import Path
import json

@dataclass
class PipelineConfig:
    """Pipeline配置参数"""
    
    # 输入输出
    repeatscout_file: str
    genome_file: str
    output_dir: str
    temp_dir: str = "temp_work"
    cache_dir: str = "cache"
    checkpoint_dir: str = "checkpoints"
    
    # Phase 1 参数
    min_length: int = 80
    max_length: int = 20000
    max_n_percent: float = 0.2
    dust_window: int = 64
    dust_threshold: int = 7
    
    # Phase 2 参数
    identity_threshold: float = 0.85
    coverage_threshold: float = 0.60
    max_recruits_per_family: int = 30
    msa_algorithm: str = "localpair"  # MAFFT L-INS-i
    column_coverage_threshold: float = 0.6
    
    # Phase 3 参数
    redundancy_threshold_masking: float = 0.95
    redundancy_threshold_analysis: float = 0.90
    min_copy_number: int = 5
    min_consensus_quality: float = 0.85
    
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
    makeblastdb_exe: str = "makeblastdb"
    cdhit_exe: str = "cd-hit-est"
    
    def save(self, filepath: str):
        """保存配置到文件"""
        with open(filepath, 'w') as f:
            json.dump(self.__dict__, f, indent=2)
    
    @classmethod
    def load(cls, filepath: str):
        """从文件加载配置"""
        with open(filepath, 'r') as f:
            return cls(**json.load(f))