# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ClassifyTE is a stacking-based machine learning framework for hierarchical classification of transposable elements (TEs). The system uses k-mer features extracted from DNA sequences and applies hierarchical classification algorithms to predict TE taxonomic labels.

## Environment Setup

This project requires conda and Java JDK. Use the provided environment file:

```bash
conda env create -f environment.yml python==3.7
conda activate ClassifyTE_env
```

## Core Workflow

### 1. Feature Generation
Extract k-mer features from FASTA sequences:

**Sequential (Original):**
```bash
python generate_feature_file.py -f input.fasta -d feature_directory -o output.csv
```

**Parallel (Optimized for multi-core systems):**
```bash
python generate_feature_file_parallel.py -f input.fasta -d feature_directory -o output.csv -j 200
# Or use the wrapper script:
./run_parallel_features.sh input.fasta output.csv feature_directory 200
```

### 2. Classification/Evaluation
Classify sequences using pre-trained models:
```bash
python evaluate.py -f feature_file.csv -n node.txt -d feature_directory -m model.pkl -a lcpnb
```

### 3. Training (Optional)
Train new models on datasets:
```bash
python train.py -f dataset.csv -n node_file.txt -m model_name -c cost_param -g gamma_param
```

## Key Architecture Components

### HierStack Package
- `hierarchy.py`: Builds taxonomic tree structure from node files using NetworkX
- `model.py`: Implements stacking classifiers with SVM base learners
- `classification.py`: Contains hierarchical classification algorithms (LCPNB, nLLCPN)
- `lcpnb.py` & `nllcpn.py`: Specific algorithm implementations
- `stackingClassifier.py`: Meta-classifier combining multiple base classifiers

### Data Flow
1. FASTA sequences → k-mer feature extraction (2-mer, 3-mer, 4-mer frequencies)
2. Features → hierarchical classification using stacking ensemble
3. Predictions follow taxonomic tree structure defined in node files

### File Structure
- `data/`: Input FASTA files and feature CSV files
- `models/`: Pre-trained pickle files (ClassifyTE_combined.pkl, ClassifyTE_repbase.pkl, ClassifyTE_pgsb.pkl)
- `nodes/`: Taxonomic hierarchy definitions (node.txt, node_pgsb.txt, node_repbase.txt)
- `features/`: Contains Java-based k-mer extraction tools (kanalyze-2.0.0/)
- `output/`: Classification results

### Model Types
- **ClassifyTE_combined.pkl**: Trained on combined dataset (C=512.0, gamma=0.0078125)
- **ClassifyTE_repbase.pkl**: Trained on RepBase dataset (C=128.0, gamma=0.0078125) 
- **ClassifyTE_pgsb.pkl**: Trained on PGSB dataset (C=32, gamma=0.03125)

### Classification Algorithms
- **LCPNB**: Local Classifier per Parent Node and Branch
- **nLLCPN**: non-Leaf Local Classifier per Parent Node

## Feature Generation Details

The system uses Java-based kanalyze tools for k-mer extraction. The feature generation process:
1. Splits FASTA into individual sequence files
2. Runs kanalyze to generate k-mer counts for 2-mer, 3-mer, and 4-mer
3. Combines features using Java collectors into single CSV file
4. Total features: 4² + 4³ + 4⁴ = 336 k-mer frequencies

## Parallel Processing Optimization

### Performance Features
- **Multi-core support**: Utilizes up to 256 CPU cores efficiently
- **Parallel k-mer extraction**: Each sequence processed independently using GNU parallel or xargs
- **Concurrent feature collection**: Java-based parallel processing with thread pools
- **Memory optimization**: Reduced memory footprint per process (1GB vs 3GB per kanalyze instance)

### Parallel Scripts
- `generate_feature_file_parallel.py`: Main parallel feature generation script
- `run_parallel_features.sh`: Convenient wrapper script with timing and validation
- `runKanalyzer_parallel`: Optimized shell script for parallel k-mer analysis
- `ParallelKmersFeaturesCollector.java`: Multi-threaded feature aggregation
- `benchmark_parallel.py`: Performance comparison tool

### Usage Recommendations
- Use parallel version for datasets with >10 sequences
- Adjust `-j` parameter based on available memory (default: 80% of CPU cores)
- Monitor system resources during large batch processing
- Expected speedup: 5-20x depending on sequence count and system specs

## Development Notes

- Models expect exactly 336 features (k-mer frequencies)
- Node files define hierarchical structure - each corresponds to specific trained model
- Output includes both CSV (sequence_id, predicted_label) and detailed text files
- Cross-validation training script available as CV_train.py
- Parallel processing maintains identical output format to sequential version