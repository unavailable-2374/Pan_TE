#!/bin/bash
# Wrapper script for parallel feature generation
# Usage: ./run_parallel_features.sh input.fasta [output.csv] [feature_dir] [num_jobs]

set -e  # Exit on any error

# Default values
FASTA_FILE="$1"
OUTPUT_FILE="${2:-feature_file.csv}"
FEATURE_DIR="${3:-features}"
NUM_JOBS="${4:-$(( $(nproc) * 4 / 5 ))}"

# Validate input
if [ -z "$FASTA_FILE" ]; then
    echo "Usage: $0 <fasta_file> [output_file] [feature_dir] [num_jobs]"
    echo "Example: $0 demo.fasta demo_features.csv demo_features 200"
    exit 1
fi

if [ ! -f "data/$FASTA_FILE" ]; then
    echo "Error: FASTA file 'data/$FASTA_FILE' not found!"
    exit 1
fi

echo "=================================================="
echo "ClassifyTE Parallel Feature Generation"
echo "=================================================="
echo "Input FASTA: data/$FASTA_FILE"
echo "Output CSV: data/$OUTPUT_FILE"
echo "Feature directory: $FEATURE_DIR"
echo "Parallel jobs: $NUM_JOBS"
echo "Available CPU cores: $(nproc)"
echo "=================================================="

# Check if required files exist
if [ ! -f "generate_feature_file_parallel.py" ]; then
    echo "Error: generate_feature_file_parallel.py not found!"
    exit 1
fi

if [ ! -f "features/ParallelKmersFeaturesCollector.java" ]; then
    echo "Error: ParallelKmersFeaturesCollector.java not found!"
    exit 1
fi

if [ ! -f "features/kanalyze-2.0.0/code/runKanalyzer_parallel" ]; then
    echo "Error: runKanalyzer_parallel script not found!"
    exit 1
fi

# Make scripts executable
chmod +x generate_feature_file_parallel.py
chmod +x features/kanalyze-2.0.0/code/runKanalyzer_parallel

# Record start time
START_TIME=$(date +%s)

# Run parallel feature generation
echo "Starting parallel feature generation..."
python3 generate_feature_file_parallel.py \
    -f "$FASTA_FILE" \
    -o "$OUTPUT_FILE" \
    -d "$FEATURE_DIR" \
    -j "$NUM_JOBS"

# Calculate elapsed time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "=================================================="
echo "Feature generation completed successfully!"
echo "Total time: ${ELAPSED} seconds"
echo "Output file: data/$OUTPUT_FILE"
echo "=================================================="

# Optional: Show file size and first few lines
if [ -f "data/$OUTPUT_FILE" ]; then
    echo "Output file info:"
    ls -lh "data/$OUTPUT_FILE"
    echo ""
    echo "First 3 lines of output:"
    head -n 3 "data/$OUTPUT_FILE"
    echo ""
    echo "Number of sequences processed: $(($(wc -l < "data/$OUTPUT_FILE") - 1))"
fi