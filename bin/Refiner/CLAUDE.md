# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a TE (Transposable Element) consensus sequence builder that constructs high-quality TE consensus sequences from RepeatScout output. It uses a three-phase "screen-represent-extend" strategy to build accurate consensus sequences through a sophisticated copy-based clustering approach.

**Input Requirements:**
- RepeatScout output file (filtered with filter-stage-1.prl)
- Reference genome file (≤400Mb recommended)

**Output:**
- Two consensus libraries: 95% redundancy removed (masking) and 90% redundancy removed (analysis)
- Statistics report with quality metrics

## Commands

### Basic Usage
```bash
# Standard run with correct logic (default)
python main.py -r repeatscout_output.fa -g genome.fa -o output_dir

# With custom threads
python main.py -r repeatscout.fa -g genome.fa -o output_dir -t 16

# Use custom configuration file
python main.py -r repeatscout.fa -g genome.fa -o output_dir -c config.json
```

### Recovery and Debugging
```bash
# Resume from checkpoint after interruption
python main.py -r repeatscout.fa -g genome.fa -o output_dir --resume --keep-checkpoints

# Clear cache and run fresh
python main.py -r repeatscout.fa -g genome.fa -o output_dir --clear-cache

# Keep temporary and checkpoint files for debugging
python main.py -r repeatscout.fa -g genome.fa -o output_dir --keep-temp --keep-checkpoints
```


### Configuration Management
```bash
# Generate default configuration file for customization
python -c "from config import PipelineConfig; PipelineConfig('dummy.fa', 'dummy.fa', 'dummy').save('config.json')"

# Use custom configuration (overrides command line params where specified)
python main.py -r repeatscout.fa -g genome.fa -o output_dir -c custom_config.json
```

### Installation
```bash
# Install Python dependencies
pip install -r requirements.txt

# Required external tools (must be in PATH or specified in config):
# - RepeatMasker (>= 4.1.0)
# - MAFFT (>= 7.450)  
# - BLAST+ (>= 2.10.0)
# - CD-HIT (>= 4.8.1)
```

### Testing
```bash
# Test current implementation
python test_optimizations.py
```

## Architecture

The pipeline has evolved to include multiple implementations with a focus on correct copy-based clustering:

### Current Implementation
The pipeline uses the latest optimized version with copy-based clustering for maximum accuracy.

### Three-Phase Architecture

#### Phase 1: Screening and Scoring (`phase1_screening_optimized.py`)
- Loads RepeatScout sequences and applies basic filters (length, N-content, GC-content)
- Calculates multi-dimensional complexity scores using DUST and Shannon entropy
- Runs RepeatMasker to evaluate genome coverage
- Categorizes sequences into A/B/C classes based on combined scores
  - A-class: High-quality seeds (score ≥ 0.75)
  - B-class: Medium quality (0.50 ≤ score < 0.75)  
  - C-class: Low quality (score < 0.50)

#### Phase 2: Copy-Based Consensus Building (`phase2_consensus_correct.py`)
**Key Innovation**: Clusters genome copies rather than RepeatScout sequences themselves

For each high-quality RepeatScout sequence:
1. **Find genome copies** using RepeatMasker
2. **Extract genome fragments** with location and orientation info
3. **Cluster copies** using hierarchical clustering with adaptive thresholds
4. **Build subfamily consensus** for each cluster using MAFFT MSA
5. **Detect TSDs** and refine boundaries
6. **Quality assessment** with boundary evaluation

#### Phase 3: Quality Control (`phase3_finalization_relaxed.py`)
- Relaxed chimera detection (preserves more sequences)
- Multi-track adaptive filtering (copy number, quality, boundary scores)
- Two-tier redundancy removal (95% and 90% thresholds)
- Comprehensive metadata generation

### Utility Modules (`utils/`)
- `robust_runner.py`: Checkpoint/resume capability, retry logic with exponential backoff
- `cache_utils.py`: Intelligent caching with MD5 hashing for expensive computations
- `sequence_utils.py`: Core sequence operations, TSD detection, orientation handling
- `blast_utils.py`: BLAST operations, similarity network construction
- `alignment_utils.py`: RepeatMasker and MAFFT wrappers with parameter tuning
- `complexity_utils.py`: DUST scoring, Shannon entropy, IUPAC ambiguity codes

## Configuration System

The `PipelineConfig` class manages 40+ parameters across multiple categories:

### Key Parameter Groups
- **Phase 1 Filtering**: `min_length` (80), `max_length` (20000), `max_n_percent` (0.2), `dust_threshold` (7)
- **Phase 2 Clustering**: `identity_threshold` (0.85), `coverage_threshold` (0.60), `max_recruits_per_family` (30)
- **Phase 3 Output**: `redundancy_threshold_masking` (0.95), `redundancy_threshold_analysis` (0.90), `min_copy_number` (5)
- **Performance**: `threads` (8), `max_retries` (3), `use_parallel` (True), `batch_size` (100)
- **External Tools**: Tool paths for RepeatMasker, MAFFT, BLAST+, CD-HIT

### Configuration Files
JSON-based configuration with inheritance and validation:
```python
# Save current config
config.save('my_config.json')

# Load and modify config
config = PipelineConfig.load('my_config.json')
config.threads = 32
```

## Error Recovery and Robustness

Multi-layered recovery system:

### Checkpoint System
- Phase-level checkpoints (`checkpoints/phase1_complete.pkl`)
- Automatic resume on restart with `--resume`
- Granular checkpoints for expensive operations

### Caching Strategy  
- Result caching with `@cache_result` decorator
- MD5-based cache keys for parameter sensitivity
- Cache invalidation and cleanup with `--clear-cache`

### Retry Logic
- Exponential backoff for transient failures
- Configurable retry attempts and delays
- Graceful degradation for non-critical failures

## Pipeline Variants and Performance

### Current Implementation Features
- **Copy-based clustering**: Clusters genome copies rather than RepeatScout sequences
- **Optimized performance**: Streamlined algorithms and parallel processing
- **Relaxed filtering**: Preserves more sequences while maintaining quality
- **High accuracy**: Sophisticated boundary detection and TSD analysis

### Memory and Performance Optimization
- Streaming file processing for large genomes
- Parallel processing with configurable thread pools
- Memory-efficient data structures for large sequence sets
- Batch processing to control memory usage

## Output Structure

```
output_dir/
├── consensus_masking.fa      # 95% redundancy removed, optimized for genome masking
├── consensus_analysis.fa     # 90% redundancy removed, with rich metadata headers
├── statistics.txt            # Comprehensive statistics and quality metrics
├── pipeline_config.json      # Complete configuration used for reproducibility
├── pipeline.log             # Detailed execution log
├── checkpoints/             # Phase completion markers (if --keep-checkpoints)
├── cache/                   # Cached computation results
└── temp_work/              # Temporary files (if --keep-temp)
```

## Key Algorithmic Innovations

### Copy-Based Clustering Logic
The critical difference from naive approaches:
- **Wrong**: Cluster RepeatScout sequences by similarity
- **Correct**: Find genome copies of each RepeatScout sequence, then cluster the copies

This produces more accurate subfamily detection and consensus building.

### Adaptive Thresholding
- Dynamic clustering thresholds based on sequence divergence patterns  
- P95-based filtering to handle varying repeat abundance
- Context-sensitive quality scoring

### TSD-Aware Boundary Refinement
- Systematic TSD (Target Site Duplication) detection
- Boundary adjustment based on TSD evidence
- Quality scoring incorporating structural features

## Development and Debugging

### Code Architecture Patterns
- Phase classes inherit common patterns: `__init__(config)`, `run()` method, `RobustRunner` integration
- Utility functions are stateless and well-tested
- Configuration drives behavior without hard-coded parameters

### Debugging Workflow
1. Run with `--keep-temp --keep-checkpoints` for full intermediate file retention
2. Check `pipeline.log` for detailed execution trace  
3. Examine `cache/` directory for cached intermediate results
4. Use `test_performance.py` for systematic performance comparison

### Testing
- `test_optimizations.py`: Validates current implementation correctness
- Test logging in `test_optimizations.log`

## External Tool Requirements

All tools must be in PATH or specified in configuration:
- **RepeatMasker** (≥ 4.1.0): Primary homology search engine
- **MAFFT** (≥ 7.450): Multiple sequence alignment with L-INS-i algorithm
- **BLAST+** (≥ 2.10.0): Sequence similarity search and network construction
- **CD-HIT** (≥ 4.8.1): Redundancy removal and clustering

## Parameter Tuning Guidelines

### For High-Repeat Genomes (>40% repetitive)
```python
config.min_copy_number = 5
config.identity_threshold = 0.80
config.redundancy_threshold_masking = 0.95
```

### For Low-Repeat Genomes (<20% repetitive)  
```python
config.min_copy_number = 2
config.identity_threshold = 0.70
config.redundancy_threshold_masking = 0.90
```

### For Ancient/Diverged TEs
```python
config.identity_threshold = 0.65
config.coverage_threshold = 0.40
config.dust_threshold = 10
```