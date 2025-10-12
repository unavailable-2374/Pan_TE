#!/usr/bin/env python3
"""
Parallel Feature Generation Script for ClassifyTE
Optimized for multi-core processing with up to 256 CPU cores

Author: Optimized version for parallel processing
License: MIT
"""

import os
import sys
import shutil
import subprocess
import multiprocessing as mp
from optparse import OptionParser
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time
from pathlib import Path

def setup_directories(feature_destpath, kanalyzer_input_destpath, kanalyzer_output_destpath):
    """Setup required directories for parallel processing"""
    
    # Clean and create input directory
    if os.path.isdir(kanalyzer_input_destpath):
        shutil.rmtree(kanalyzer_input_destpath)
    os.makedirs(kanalyzer_input_destpath)

    # Clean and create output directory
    if os.path.isdir(kanalyzer_output_destpath):
        shutil.rmtree(kanalyzer_output_destpath)
    os.makedirs(kanalyzer_output_destpath)
    
    # Create k-mer specific directories
    for kmer_size in ['2mer', '3mer', '4mer']:
        kmer_dir = os.path.join(kanalyzer_output_destpath, kmer_size)
        os.makedirs(kmer_dir, exist_ok=True)

def parallel_sequence_splitting(fasta_file, data_filepath, kanalyzer_input_destpath, max_workers=None):
    """Split FASTA file into individual sequences using parallel processing"""
    
    print(f"Reading and splitting FASTA file: {fasta_file}")
    
    with open(os.path.join(data_filepath, fasta_file), 'rt') as fp:
        content = fp.read()
        sequences = content.split(">")[1:]  # Skip empty first element
    
    print(f"Found {len(sequences)} sequences to process")
    
    def write_sequence(args):
        seq_data, seq_id, output_dir = args
        output_file = os.path.join(output_dir, f"seq{seq_id}.fasta")
        with open(output_file, 'w') as of:
            of.write(">" + seq_data)
    
    # Use optimal number of workers for I/O operations
    if max_workers is None:
        max_workers = min(32, mp.cpu_count())  # Limit for I/O operations
    
    # Prepare arguments for parallel processing
    args_list = [(seq_data, i+1, kanalyzer_input_destpath) 
                 for i, seq_data in enumerate(sequences) if seq_data.strip()]
    
    # Write sequences in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        executor.map(write_sequence, args_list)
    
    return len(args_list)

def create_file_list(kanalyzer_input_destpath, feature_destpath, kanalyzer_destpath):
    """Create list of sequence files for processing"""
    
    files = sorted([f for f in os.listdir(kanalyzer_input_destpath) 
                   if f.endswith('.fasta')])
    
    list_file_path = os.path.join(kanalyzer_input_destpath, 'list.txt')
    with open(list_file_path, 'w') as ff:
        for f in files:
            ff.write(f + '\n')
    
    # Copy list file to required locations
    shutil.copy2(list_file_path, os.path.join(feature_destpath, 'list.txt'))
    shutil.copy2(list_file_path, os.path.join(kanalyzer_destpath, 'list.txt'))
    
    return len(files)

def parallel_feature_generation(curr_dir1, output_file, feature_dir):
    """Run parallel k-mer feature generation"""
    
    feature_destpath = os.path.join(curr_dir1, feature_dir)
    kanalyzer_destpath = os.path.join(feature_destpath, "kanalyze-2.0.0", "code")
    
    print(f"Starting parallel k-mer analysis with {mp.cpu_count()} available cores")
    
    # Change to kanalyzer directory
    original_dir = os.getcwd()
    os.chdir(kanalyzer_destpath)
    
    try:
        # Make parallel script executable
        parallel_script = "./runKanalyzer_parallel"
        subprocess.run(['chmod', '775', parallel_script], check=True)
        
        # Run parallel k-mer analysis
        print("Running parallel k-mer feature extraction...")
        start_time = time.time()
        result = subprocess.run([parallel_script], check=True, capture_output=True, text=True)
        
        kmer_time = time.time() - start_time
        print(f"K-mer extraction completed in {kmer_time:.2f} seconds")
        
        if result.stdout:
            print("K-mer analysis output:", result.stdout)
            
    except subprocess.CalledProcessError as e:
        print(f"Error running parallel k-mer analysis: {e}")
        if e.stderr:
            print("Error output:", e.stderr)
        raise
    finally:
        os.chdir(original_dir)
    
    # Change to feature directory for Java compilation and execution
    os.chdir(feature_destpath)
    
    try:
        print("Compiling and running parallel feature collector...")
        start_time = time.time()
        
        # Compile Java files
        subprocess.run(['javac', 'ParallelKmersFeaturesCollector.java'], check=True)
        subprocess.run(['javac', 'BufferReaderAndWriter.java'], check=True)
        
        # Run parallel feature collector
        subprocess.run(['java', 'ParallelKmersFeaturesCollector'], check=True)
        
        collection_time = time.time() - start_time
        print(f"Feature collection completed in {collection_time:.2f} seconds")
        
        # Rename output file if needed
        if output_file != "feature_file.csv":
            subprocess.run(['mv', 'feature_file.csv', output_file], check=True)
        
        # Copy to data directory
        curr_dir2 = os.getcwd()
        data_dir = os.path.join(curr_dir1, 'data')
        shutil.copy2(os.path.join(curr_dir2, output_file), 
                     os.path.join(data_dir, output_file))
        
        # Clean up temporary file
        os.remove(output_file)
        
    except subprocess.CalledProcessError as e:
        print(f"Error in feature collection: {e}")
        raise
    finally:
        os.chdir(curr_dir1)

def parallel_get_data(fasta_file, feature_dir, max_workers=None):
    """Process FASTA data with parallel sequence splitting"""
    
    curr_dir1 = os.getcwd()
    data_filepath = os.path.join(curr_dir1, "data")
    feature_destpath = os.path.join(curr_dir1, feature_dir)
    kanalyzer_destpath = os.path.join(feature_destpath, "kanalyze-2.0.0", "code")
    kanalyzer_input_destpath = os.path.join(feature_destpath, "kanalyze-2.0.0", "input_data")
    kanalyzer_output_destpath = os.path.join(feature_destpath, "kanalyze-2.0.0", "output_data")
    
    # Setup directories
    setup_directories(feature_destpath, kanalyzer_input_destpath, kanalyzer_output_destpath)
    
    # Split sequences in parallel
    num_sequences = parallel_sequence_splitting(fasta_file, data_filepath, 
                                               kanalyzer_input_destpath, max_workers)
    
    # Create file list
    num_files = create_file_list(kanalyzer_input_destpath, feature_destpath, kanalyzer_destpath)
    
    print(f"Successfully prepared {num_files} sequence files for parallel processing")
    return num_files

def main():
    """Main function with optimized parallel processing"""
    
    parser = OptionParser(description="Parallel k-mer feature generation for ClassifyTE")
    parser.add_option("-f", "--filename", dest="filename", 
                     help="Name of the fasta file.")
    parser.add_option("-o", "--output", dest="output_filename", 
                     help="Name of feature file", default="feature_file.csv")
    parser.add_option("-d", "--featuredir", dest="feature_dir", 
                     help="Feature directory.", default="features")
    parser.add_option("-j", "--jobs", dest="max_jobs", type="int",
                     help=f"Maximum number of parallel jobs (default: {mp.cpu_count() * 4 // 5})",
                     default=mp.cpu_count() * 4 // 5)
    
    (options, args) = parser.parse_args()
    
    if not options.filename:
        parser.error("FASTA filename is required. Use -f option.")
    
    print(f"ClassifyTE Parallel Feature Generator")
    print(f"Available CPU cores: {mp.cpu_count()}")
    print(f"Using parallel jobs: {options.max_jobs}")
    print(f"Processing: {options.filename}")
    print(f"Output: {options.output_filename}")
    print("-" * 50)
    
    curr_dir1 = os.getcwd()
    feature_dir = os.path.join(curr_dir1, options.feature_dir)
    src = os.path.join(curr_dir1, "features")
    
    # Copy feature directory if using custom name
    if options.feature_dir != "features":
        if os.path.exists(feature_dir):
            shutil.rmtree(feature_dir)
        shutil.copytree(src, feature_dir)
    
    start_time = time.time()
    
    # Process data with parallel splitting
    num_sequences = parallel_get_data(options.filename, feature_dir, 
                                     max_workers=min(32, options.max_jobs))
    
    # Generate features in parallel
    parallel_feature_generation(curr_dir1, options.output_filename, feature_dir)
    
    total_time = time.time() - start_time
    
    print("-" * 50)
    print(f"Parallel processing completed successfully!")
    print(f"Processed {num_sequences} sequences")
    print(f"Total processing time: {total_time:.2f} seconds")
    print(f"Average time per sequence: {total_time/num_sequences:.3f} seconds")
    print(f"Output file: data/{options.output_filename}")

if __name__ == '__main__':
    main()