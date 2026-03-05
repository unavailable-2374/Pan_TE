# Pan_TE Pipeline Architecture Memory Document

## Project Overview

Pan_TE is a **pan-genome transposable element (TE) annotation pipeline** that combines three independent TE detection strategies (LTR retrotransposon detection, RepeatScout, and RECON) into a unified workflow. It supports optional VCF-derived structural variant sequences and TE classification via ClassifyTE + RepeatClassifier.

**Entry point**: `bin/Pan_TE` (Python3)
**Version**: 1.0.0
**Dev repo**: `/home/shuoc/tool/tmp/Pan_TE/`
**Production repo (on PATH)**: `/home/shuoc/tool/Pan_TE/`

---

## Pipeline Execution Flow

```
Pan_TE (main orchestrator)
  |
  |-- Step 1: Genome Processing (process_genome)
  |     |-- [Optional] VCF Processing (vcf_processor.py)
  |     |-- clean_seq_fast: uppercase + degenerate base cleaning, rename to chr1..chrN
  |     |-- samtools faidx + index
  |     +-- Checkpoint: genome.ok
  |
  |-- Step 2+3: Parallel Detection (ThreadPoolExecutor, 2 workers)
  |     |
  |     |-- [Thread A] LTR Detection (process_ltr)
  |     |     |-- LTR_detect: Look4LTRs C++ tool (genome size-based splitting)
  |     |     |-- build_ltr_library.py: GraphGroup family assignment + consensus
  |     |     +-- Checkpoint: LTR.ok
  |     |
  |     +-- [Thread B] RepeatScout (process_repeatscout)
  |           |-- build_RS (Perl): multi-parameter RepeatScout
  |           |     |-- Genome sampling (BED-based, size-adaptive)
  |           |     |-- Dust + TRF soft masking
  |           |     |-- Parallel RepeatScout jobs (multiple l-mer x samples)
  |           |     |-- Merge results
  |           |     +-- Checkpoint: repeatscout_complete
  |           |-- Refiner Pipeline (Python3, 4 phases)
  |           |     |-- Phase 1: Candidate screening (complexity, length, N%)
  |           |     |-- Phase 2: Chimera detection & splitting
  |           |     |-- Phase 3: Consensus building (MAFFT MSA, quality tiers)
  |           |     +-- Phase 4: Genome masking (optional, RepeatMasker hard mask)
  |           +-- Checkpoint: RepeatScout.ok
  |
  |-- Step 4: RECON (process_recon)
  |     |-- Uses Refiner-masked genome from RepeatScout step
  |     |-- Uses original genome (genome/genome.fa) for flanking extension
  |     |-- run_RECON_advanced (Perl):
  |     |     |-- Small genome: single_track (masked track processing)
  |     |     |-- Large genome: dual_track
  |     |     |     |-- Track 1: Masked genome processing
  |     |     |     +-- Track 2: Sampling with progressive masking (90→270→540→1200 MB)
  |     |     |-- Unmasked region extraction with flanking extension (±500bp from original genome)
  |     |     |-- RMBlastN self-alignment (word_size=9, -complexity_adjust) -> MSP format
  |     |     |-- RECON (imagespread, eledef, eleredef, edgeredef, famdef)
  |     |     |-- Consensus: RepeatModeler Refiner (preferred) or Refiner.py (fallback)
  |     |     +-- Modules: lib/RECON/{Core,Config,Logger,Utils,MaskedTrack,SamplingTrack}.pm
  |     +-- Checkpoint: RECON.ok
  |
  |-- Step 5: Combine & Classify (combine_results)
  |     |-- Collect consensus from all sources:
  |     |     |-- LTR: Look4LTRs/consensus_lib/ltr_library.fa
  |     |     |-- RepeatScout: RepeatScout/refiner_output/phase3_analysis_library.fa
  |     |     +-- RECON: RECON/{consensi.fa, masked_track/consensi.fa, sampling_track/round_*/consensi.fa}
  |     |-- CD-HIT-EST deduplication (80% identity, local alignment, both strands)
  |     |-- [Optional] Classification via run_Classifier:
  |     |     |-- RepeatClassifier (parallel)
  |     |     +-- ClassifyTE (conda env: ClassifyTE_env, parallel)
  |     |     +-- Combine_for_Two: merge both classification results
  |     +-- Checkpoint: Combine.ok
  |
  +-- Output: Combine/TEs.fa (classified TE library)
```

---

## Directory Structure (Runtime)

```
<output_dir>/
  genome/
    genome.fa              # Cleaned, indexed genome (unmasked, used as original_genome by RECON)
    genome.fa.fai          # samtools index
  Look4LTRs/
    l4s.ok                 # LTR_detect checkpoint
    consensus_lib/
      ltr_library.fa       # Final LTR consensus library
  RepeatScout/
    tmp/                   # RepeatScout working files
      repeats.fa           # Raw RepeatScout output (merged)
    refiner_output/
      phase3_analysis_library.fa     # Refiner output for Combine
      consensus_masking.fa           # Refiner output for masking
      genome_final_masked.fa         # Masked genome for RECON
    consensi.fa            # Copy of Refiner masking output
  RECON/
    consensi.fa            # Single-track result (small genomes)
    masked_track/
      consensi.fa          # Dual-track: masked track result
    sampling_track/
      round_*/consensi.fa  # Dual-track: sampling track results
  Combine/
    raw_TEs_combined.fa    # All sources concatenated
    raw_TEs.fa             # Deduplicated TEs
    TEs.fa                 # Final classified TE library
  *.ok                     # Pipeline checkpoints
  pan_te.log               # Main pipeline log
```

---

## Key Components Detail

### 1. Genome Processor (`clean_seq_fast`)
- **Language**: Python3
- **Function**: Single-pass uppercase conversion + degenerate base (WRMKYSHBVDX) -> N replacement
- **Performance**: Streaming binary I/O with translation table, ~Mbp/s speed
- **Renaming**: All sequences renamed to chr1, chr2, ... chrN

### 2. LTR Detection (`LTR_detect` + `build_ltr_library.py`)
- **LTR_detect**: Python3 wrapper around Look4LTRs C++ tool
  - Genome <5GB: run on whole genome
  - 5-10GB: split into 2 parts
  - 10-15GB: 3 parts; 15-20GB: 4 parts
  - >20GB: sample 4 chunks of ~5GB each
  - Uses ProcessPoolExecutor for parallel part processing
- **build_ltr_library.py**: Builds consensus from .rtr files
  - Uses Look4LTRs GraphGroup family assignment
  - LTR identity-weighted consensus
  - Cross-family merging via FastGA
  - Multi-part: build per-part then CD-HIT-EST merge (80% identity)

### 3. RepeatScout (`build_RS`)
- **Language**: Perl (v2.0.0)
- **Strategy**: Size-adaptive multi-parameter approach
  - Small (<200MB): full genome, l-mer 14+16, 1 sample
  - Medium (200MB-1GB): l-mer 14+16+18, 1 sample (full or 800MB)
  - Large (1-2GB): 1x800MB sample, 3 l-mers
  - Very large (2-3GB): 2x800MB samples, 3 l-mers (6 jobs)
  - Huge (3-5GB): 3x800MB samples (9 jobs)
  - Massive (>5GB): 4x800MB samples (12 jobs)
- **Sampling**: BED-based progressive non-overlapping sampling
- **Masking**: Dust + TRF soft masking before RepeatScout
- **Parallelism**: Parallel::ForkManager for parallel l-mer x sample jobs
- **Dependencies**: RMBlast (via RepModelConfig or environment)

### 4. Refiner (`bin/Refiner/`)
- **Language**: Python3
- **Config**: `config.py` (PipelineConfig dataclass, JSON serializable)
- **Phase 1** (`phase1_screening_optimized.py`): Sequence screening
  - Complexity scoring (DUST + Shannon entropy)
  - Length filtering (50-50000bp)
  - N% filtering (max 20%)
  - Trusts RepeatScout's initial filtering
- **Phase 2** (`phase2_chimera_splitting.py`): Chimera detection
  - RepeatMasker-based chimera identification
  - Conservative splitting strategy
  - Parallel processing (I/O intensive)
- **Phase 3** (`phase3_consensus_building.py`): Consensus building
  - BLASTN recruitment of genomic copies
  - MAFFT L-INS-i multiple sequence alignment
  - Quality-tiered consensus (high/medium/low)
  - Hierarchical clustering (scipy)
  - Output: masking library + analysis library
- **Phase 4** (`phase4_genome_masking.py`): Genome masking (optional, default off)
  - RepeatMasker hard masking with consensus library
  - ID renaming for RepeatMasker compatibility
  - Output: masked genome for RECON

### 5. RECON (`run_RECON_advanced` + `lib/RECON/`)
- **Language**: Perl
- **Modules**: Config, Logger, Utils, Core, MaskedTrack, SamplingTrack
- **Processing modes**:
  - `single_track`: Small genomes - masked track only
  - `dual_track`: Large genomes - masked track + sampling track
- **RECON Optimizations** (dev repo):
  - **Flanking extension**: `extract_unmasked_regions_extended()` in MaskedTrack.pm extends unmasked regions ±500bp, merges within 200bp, extracts from original genome (not masked genome) to recover real DNA bases at boundaries
  - **Original genome threading**: `Config.pm:find_input_files()` discovers `../genome/genome.fa` as `$config->{original_genome_file}`; passed through `run_RECON_advanced` → `run_masked_track_independent()` / `run_sampling_track_independent()` → extraction functions
  - **Sensitive rmblastn**: `word_size=9` (was 11), `-complexity_adjust` added (Shannon entropy-weighted scoring)
  - **Lower thresholds**: `familySizeCutoff=3` (was 15), `minAlignmentSize=3` (was 8) — appropriate for masked genome with fewer remaining copies
  - **RepeatModeler Refiner**: `detect_rm_refiner()` in `build_for_RECON` checks `$CONDA_PREFIX/share/RepeatModeler/Refiner` then `which Refiner`; per-family fallback to `Refiner.py` on failure; handles both output header formats
  - **Unified MSP output**: All three rmblastn paths (single/parallel/sequential) produce MSP format via `MSPCollect.pl`; `run_single_rmblastn_process` outputs to temp file then converts (consistent with parallel path)
- **Core steps**: extract unmasked regions → gi|N rename → create BLAST DB → RMBlastN self-alignment → MSPCollect.pl → determine K → imagespread → eledef → eleredef → edgeredef → famdef → build_for_RECON (Refiner) → consensi.fa
- **Input**: Refiner-generated masked genome (`RepeatScout/refiner_output/genome_final_masked.fa`)
- **Exported functions per module** (after cleanup):
  - `Core.pm`: `run_recon_pipeline`, `run_rmblastn_self_alignment`, `run_parallel_rmblastn_chunked`, `create_blast_database`, `create_seq_name_list`, `determine_k_parameter`, `cleanup_intermediate_files`, `merge_msp_files`, `split_fasta_by_size`
  - `MaskedTrack.pm`: `run_masked_track_independent`, `extract_unmasked_regions`, `extract_unmasked_regions_extended`
  - `SamplingTrack.pm`: `run_sampling_track_independent`, `analyze_round_metrics`, `evaluate_stopping_criteria`, `update_accumulated_mask`, `log_debugging_info`, `extract_unmasked_regions_with_gi_format`

### 6. Classifier (`run_Classifier`)
- **Language**: Bash
- **Dual classification** (parallel):
  - RepeatClassifier (RepeatModeler tool)
  - ClassifyTE (Python, conda env: ClassifyTE_env, feature generation + ML evaluation)
- **Combination**: `Combine_for_Two` merges both results
- **Post-processing**: Standardize TE class names (LINE, LTR/Gypsy, etc.)
- **Timeout**: 2 hours per classification process

### 7. VCF Processor (`vcf_processor.py`)
- **Language**: Python3
- **Function**: Process structural variant VCF files to extract TE-related sequences
- **Thread allocation**: Proportional to file size
- **Parallel processing**: ThreadPoolExecutor

---

## Thread Allocation Strategy

The pipeline dynamically allocates threads based on genome size:
- **chunk_count**: Determined by genome file size (1-4 chunks for RepeatScout)
- **ltr_threads**: `max(1, total_threads - chunk_count * 3)`
- **rs_threads**: `max(1, total_threads)` (full allocation)
- LTR and RepeatScout run in **parallel** via ThreadPoolExecutor(max_workers=2)
- RECON runs **sequentially** after both complete

---

## Checkpoint System

All major steps have `.ok` checkpoint files in the output directory:
- `genome.ok`, `LTR.ok`, `RepeatScout.ok`, `RECON.ok`, `Combine.ok`
- Sub-checkpoints: `l4s.ok` (Look4LTRs), `repeatscout_complete` (build_RS)
- Refiner uses its own checkpoint system (`checkpoints/` directory)
- RECON sub-checkpoints per track: `seq_naming.ok`, `rmblastn.ok`, `msp_collection.ok`, `recon.ok`, `consensus.ok`, `round_completed.ok`
- Pipeline supports **resume**: skips completed steps on re-run
- **Important**: RECON checkpoints should be cleared for re-runs after code updates

---

## External Tool Dependencies

### Complete Tool Chain (verified available)

| Tool | Used By | Purpose | Location |
|------|---------|---------|----------|
| **Genome Processing** ||||
| clean_seq_fast | Pan_TE | Uppercase + degenerate base cleaning | Pan_TE/bin/ (Python3) |
| samtools | Pan_TE, build_RS, Refiner | FASTA indexing/extraction | hap.py-build/bin/ |
| index | Pan_TE | Additional FASTA indexing | Pan_TE/bin/ (Perl) |
| **LTR Detection** ||||
| look4ltrs | LTR_detect | LTR retrotransposon detection (C++) | Look4LTRs/bin/ |
| bedtools | build_ltr_library, build_RS, SamplingTrack | BED ops, sequence extraction | PGTA conda env |
| FastGA | build_ltr_library | Cross-family LTR alignment | /home/shuoc/tool/FASTGA/ |
| **RepeatScout** ||||
| build_lmer_table | build_RS | L-mer frequency table generation | PGTA conda env |
| RepeatScout | build_RS | De novo repeat discovery | PGTA conda env |
| filter-stage-1.prl | build_RS | RepeatScout output filtering | PGTA conda env |
| dustmasker | build_RS | Low-complexity masking | PGTA conda env |
| trf | build_RS | Tandem repeat masking | PGTA conda env |
| **Alignment/Masking** ||||
| rmblastn | build_RS, RECON, Refiner | Sensitive repeat alignment | PGTA conda env |
| blastn | Refiner Phase 3 | Genomic copy recruitment | PGTA conda env |
| makeblastdb | Refiner, RECON | BLAST database creation | PGTA conda env |
| RepeatMasker | Refiner Phase 2/4, SamplingTrack | TE annotation/masking | PGTA conda env |
| mafft | Refiner Phase 3, build_ltr_library | Multiple sequence alignment | PGTA conda env |
| **RECON** ||||
| imagespread | Core.pm | RECON step 1: image spreading | PGTA conda env |
| eledef | Core.pm | RECON step 2: element definition | PGTA conda env |
| eleredef | Core.pm | RECON step 3: element redefinition | PGTA conda env |
| edgeredef | Core.pm | RECON step 4: edge redefinition | PGTA conda env |
| famdef | Core.pm | RECON step 5: family definition | PGTA conda env |
| MSPCollect.pl | Core.pm, MaskedTrack | BLAST→MSP format conversion | Pan_TE/bin/ |
| Refiner (RM) | build_for_RECON | Iterative pairwise TE consensus | PGTA/share/RepeatModeler/ |
| Refiner.py | build_for_RECON (fallback) | K-mer Jaccard + MAFFT consensus | Pan_TE/bin/ |
| **Deduplication** ||||
| cd-hit-est | Combine, build_RS, build_ltr_library | Sequence clustering (80% identity) | PGTA conda env |
| **Classification** ||||
| renameTE | run_Classifier | TE sequence renaming | Pan_TE/bin/ (Perl) |
| RepeatClassifier | run_Classifier | RepeatModeler-based classification | PGTA conda env |
| ClassifyTE tools | run_Classifier | ML-based classification | ClassifyTE_env conda env |
| process_for_classify.py | run_Classifier | Classification result processing | Pan_TE/bin/ |
| Combine_for_Two | run_Classifier | Merge dual classification results | Pan_TE/bin/ |

### Conda Environments
- **PGTA**: Main environment with most bioinformatics tools (rmblastn, RepeatMasker, RECON, etc.)
- **ClassifyTE_env**: Isolated environment for ClassifyTE ML classification (Python 3.8)

---

## Input/Output

### Required Input
- `--genome`: FASTA file (.fa/.fasta/.fna/.fas), optionally gzip/bzip2 compressed

### Optional Input
- `--model-dir`: ClassifyTE model directory (for classification)
- `--vcf-dir`: VCF files directory (structural variants)
- `--threads`: Number of threads (default: 4, recommend multiples of 4)
- `-M`: Memory limit in MB
- `--fragment_size`: Fragment size (default: 40000)

### Final Output
- `Combine/TEs.fa`: Classified TE consensus library (if model-dir provided)
- `Combine/raw_TEs.fa`: Deduplicated TE library (without classification)

---

## Key Design Decisions

1. **Three-pronged detection**: LTR-specific (Look4LTRs) + de novo repeat (RepeatScout) + RECON for complementary coverage
2. **Sequential dependency**: RECON uses the masked genome from Refiner, reducing redundant work
3. **Adaptive scaling**: All major steps adjust strategy based on genome size
4. **Checkpoint-based resume**: Every major step checkpointed, full pipeline resumable
5. **Refiner as quality gate**: RepeatScout raw output refined through 4-phase pipeline before use
6. **Dual classification**: Two independent classifiers (RepeatClassifier + ClassifyTE) combined for robustness
7. **RECON flanking extension**: Compensates for TE fragmentation in masked genomes by extracting ±500bp flanking sequence from original genome
8. **RepeatModeler Refiner preferred**: Purpose-built iterative pairwise seed alignment produces better TE consensus than generic MAFFT; per-family fallback to Refiner.py ensures robustness

---

## Dev vs Production Note

The dev repo (`/home/shuoc/tool/tmp/Pan_TE/`) contains RECON optimizations not yet deployed to production (`/home/shuoc/tool/Pan_TE/`). PATH points to production. `run_RECON_advanced` uses `use lib "$Bin/../lib"` to resolve modules relative to its own location, so the production version loads production modules. Changes must be copied to production to take effect.
