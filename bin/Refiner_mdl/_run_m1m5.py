import sys, os
# Resolve bin/ and bin/Refiner relative to this file (dev/ablation harness; portable).
_BIN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _BIN); sys.path.insert(0, os.path.join(_BIN, 'Refiner'))
from Refiner_mdl.config import RefinerMdlConfig
from Refiner_mdl.main import RefinerMdlPipeline
D="/scratch/shuoc/TE/Arabidopsis_thaliana/demo"
cfg=RefinerMdlConfig(input_file=f"{D}/mdl-repeat/tmp/mdl_repeat.raw.fa", genome_file=f"{D}/genome/genome.fa",
    output_dir="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/m1m5_refine", temp_dir="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/m1m5_refine/temp_work", checkpoint_dir="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/m1m5_refine/checkpoints",
    bed_file=f"{D}/mdl-repeat/tmp/mdl_repeat.instances.bed", stats_file=f"{D}/mdl-repeat/tmp/mdl_repeat.stats.tsv",
    threads=16, enable_masking=False)
cfg.subfamily_divergence_cut=0.999   # single cluster -> single consensus (no subfamily) = M1-M5 ablation
cfg.keep_temp=False
RefinerMdlPipeline(cfg).run()
