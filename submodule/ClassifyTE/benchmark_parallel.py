#!/usr/bin/env python3
"""
Benchmark script to compare original vs parallel feature generation
"""

import os
import time
import subprocess
import sys
from pathlib import Path

def run_command_with_timing(command, description):
    """Run a command and measure execution time"""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(command)}")
    print(f"{'='*60}")
    
    start_time = time.time()
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        end_time = time.time()
        elapsed = end_time - start_time
        
        print(f"✓ Completed successfully in {elapsed:.2f} seconds")
        if result.stdout:
            print("Output:", result.stdout[-500:])  # Show last 500 chars
        return elapsed, True
    except subprocess.CalledProcessError as e:
        end_time = time.time()
        elapsed = end_time - start_time
        print(f"✗ Failed after {elapsed:.2f} seconds")
        print(f"Error: {e}")
        if e.stderr:
            print(f"Stderr: {e.stderr}")
        return elapsed, False

def check_file_exists(filepath, description):
    """Check if a file exists and show its size"""
    if os.path.exists(filepath):
        size = os.path.getsize(filepath)
        print(f"✓ {description}: {filepath} ({size:,} bytes)")
        return True
    else:
        print(f"✗ {description}: {filepath} (not found)")
        return False

def main():
    """Main benchmark function"""
    
    # Check if demo file exists
    demo_file = "data/demo.fasta"
    if not os.path.exists(demo_file):
        print(f"Error: Demo file {demo_file} not found!")
        print("Please ensure demo.fasta exists in the data/ directory")
        return
    
    print("ClassifyTE Feature Generation Benchmark")
    print(f"Available CPU cores: {os.cpu_count()}")
    print(f"Demo file: {demo_file}")
    
    # Count sequences in demo file
    with open(demo_file, 'r') as f:
        seq_count = f.read().count('>')
    print(f"Number of sequences: {seq_count}")
    
    results = {}
    
    # Test 1: Original implementation
    print(f"\n{'#'*60}")
    print("TEST 1: Original Sequential Implementation")
    print(f"{'#'*60}")
    
    # Clean up any existing output
    for file in ["data/demo_original.csv", "data/demo_parallel.csv"]:
        if os.path.exists(file):
            os.remove(file)
    
    original_cmd = [
        "python3", "generate_feature_file.py",
        "-f", "demo.fasta",
        "-o", "demo_original.csv",
        "-d", "demo_features_original"
    ]
    
    original_time, original_success = run_command_with_timing(
        original_cmd, "Original Sequential Feature Generation"
    )
    results['original'] = {'time': original_time, 'success': original_success}
    
    # Verify original output
    original_output_exists = check_file_exists("data/demo_original.csv", "Original output")
    
    # Test 2: Parallel implementation
    print(f"\n{'#'*60}")
    print("TEST 2: Parallel Implementation")
    print(f"{'#'*60}")
    
    parallel_cmd = [
        "python3", "generate_feature_file_parallel.py",
        "-f", "demo.fasta", 
        "-o", "demo_parallel.csv",
        "-d", "demo_features_parallel",
        "-j", str(min(os.cpu_count() * 4 // 5, 50))  # Limit jobs for demo
    ]
    
    parallel_time, parallel_success = run_command_with_timing(
        parallel_cmd, "Parallel Feature Generation"
    )
    results['parallel'] = {'time': parallel_time, 'success': parallel_success}
    
    # Verify parallel output
    parallel_output_exists = check_file_exists("data/demo_parallel.csv", "Parallel output")
    
    # Compare outputs if both exist
    if original_output_exists and parallel_output_exists:
        print(f"\n{'='*60}")
        print("OUTPUT COMPARISON")
        print(f"{'='*60}")
        
        # Compare file sizes
        orig_size = os.path.getsize("data/demo_original.csv")
        par_size = os.path.getsize("data/demo_parallel.csv")
        print(f"Original file size: {orig_size:,} bytes")
        print(f"Parallel file size: {par_size:,} bytes")
        print(f"Size difference: {abs(orig_size - par_size):,} bytes")
        
        # Compare line counts
        with open("data/demo_original.csv", 'r') as f:
            orig_lines = len(f.readlines())
        with open("data/demo_parallel.csv", 'r') as f:
            par_lines = len(f.readlines())
        print(f"Original lines: {orig_lines}")
        print(f"Parallel lines: {par_lines}")
        
        # Quick content comparison (first few lines)
        print("\nFirst 3 lines comparison:")
        with open("data/demo_original.csv", 'r') as f:
            orig_head = f.readlines()[:3]
        with open("data/demo_parallel.csv", 'r') as f:
            par_head = f.readlines()[:3]
        
        for i, (orig_line, par_line) in enumerate(zip(orig_head, par_head)):
            match = "✓" if orig_line.strip() == par_line.strip() else "✗"
            print(f"Line {i+1}: {match}")
    
    # Performance summary
    print(f"\n{'#'*60}")
    print("PERFORMANCE SUMMARY")
    print(f"{'#'*60}")
    
    if results['original']['success'] and results['parallel']['success']:
        speedup = results['original']['time'] / results['parallel']['time']
        print(f"Original time:     {results['original']['time']:.2f} seconds")
        print(f"Parallel time:     {results['parallel']['time']:.2f} seconds")
        print(f"Speedup:           {speedup:.2f}x")
        print(f"Time saved:        {results['original']['time'] - results['parallel']['time']:.2f} seconds")
        print(f"Efficiency:        {(speedup / os.cpu_count()) * 100:.1f}% of theoretical maximum")
        
        if speedup > 1:
            print("✓ Parallel implementation is faster!")
        else:
            print("⚠ Parallel implementation is slower (possibly due to overhead)")
    else:
        print("⚠ Cannot compare performance due to failed executions")
        if not results['original']['success']:
            print(f"  Original failed: {results['original']['time']:.2f} seconds")
        if not results['parallel']['success']:
            print(f"  Parallel failed: {results['parallel']['time']:.2f} seconds")
    
    print(f"\nBenchmark completed!")

if __name__ == '__main__':
    main()