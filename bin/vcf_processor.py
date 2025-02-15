#!/usr/bin/env python3
"""
VCF Processor for Pan_TE

Handles parallel processing of multiple VCF files with dynamic thread allocation.
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
from typing import List, Dict, Tuple
from dataclasses import dataclass
import tempfile
import shutil

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    handlers=[
        logging.FileHandler('vcf_processor.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

@dataclass
class VCFFile:
    """Data class for VCF file information."""
    path: Path
    size: int
    allocated_threads: int = 1

class VCFProcessor:
    """Handles parallel processing of VCF files."""
    
    def __init__(self, vcf_dir: str, output_dir: str, total_threads: int):
        """Initialize the VCF processor.
        
        Args:
            vcf_dir: Directory containing VCF files
            output_dir: Directory for output files
            total_threads: Total available threads
        """
        self.vcf_dir = Path(vcf_dir)
        self.output_dir = Path(output_dir)
        self.total_threads = total_threads
        self.temp_dir = Path(output_dir) / 'tmp'
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def get_vcf_files(self) -> List[VCFFile]:
        """Get list of VCF files with their sizes."""
        vcf_files = []
        for file_path in self.vcf_dir.glob("*.vcf"):
            size = file_path.stat().st_size
            vcf_files.append(VCFFile(file_path, size))
        return sorted(vcf_files, key=lambda x: x.size, reverse=True)
    
    def allocate_threads(self, vcf_files: List[VCFFile]) -> List[VCFFile]:
        """Allocate threads to VCF files based on their sizes.
        
        Args:
            vcf_files: List of VCF files to process
            
        Returns:
            List of VCF files with allocated thread counts
        """
        if not vcf_files:
            return []
        
        total_size = sum(vcf.size for vcf in vcf_files)
        if total_size == 0:
            logger.warning("Total size of VCF files is 0, using equal thread distribution")
            threads_per_file = max(1, self.total_threads // len(vcf_files))
            for vcf in vcf_files:
                vcf.allocated_threads = threads_per_file
            return vcf_files

        remaining_threads = self.total_threads
        for i, vcf in enumerate(vcf_files):
            if i == len(vcf_files) - 1:
                # Last file gets all remaining threads
                vcf.allocated_threads = max(1, remaining_threads)
            else:
                # Calculate proportional thread allocation
                proportion = vcf.size / total_size
                allocated = max(1, int(proportion * self.total_threads))
                allocated = min(allocated, remaining_threads - (len(vcf_files) - i - 1))
                vcf.allocated_threads = allocated
                remaining_threads -= allocated

        return vcf_files
    
    def process_vcf(self, vcf_file: VCFFile) -> str:
        """Process a single VCF file."""
        vcf_temp_dir = os.path.join(self.temp_dir, vcf_file.path.stem)
        os.makedirs(vcf_temp_dir, exist_ok=True)
        
        output_file = os.path.join(vcf_temp_dir, f"{vcf_file.path.stem}_processed.fa")
        
        try:
            logger.info(f"Processing {vcf_file.path.name} with {vcf_file.allocated_threads} threads")
            
            # Ensure allocated_threads is an integer
            threads = str(max(1, int(vcf_file.allocated_threads)))
            
            cmd = [
                "decode_gfa.pl",
                "--vcf", str(vcf_file.path),
                "--threads", threads,
                "--out", vcf_temp_dir
            ]
            
            result = subprocess.run(cmd, 
                                  check=True, 
                                  capture_output=True, 
                                  text=True)
            
            if result.stderr:
                logger.warning(f"Warning messages from decode_gfa.pl: {result.stderr}")
            
            final_output = os.path.join(self.output_dir, f"{vcf_file.path.stem}_sequences.fa")
            if os.path.exists(output_file):
                with open(output_file, 'r') as src, open(final_output, 'w') as dst:
                    dst.write(src.read())
            else:
                raise FileNotFoundError(f"Expected output file not found: {output_file}")
                
            # Clean up
            import shutil
            shutil.rmtree(vcf_temp_dir)
            logger.info(f"Processed {vcf_file.path.name} successfully")
            
            return str(final_output)
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running decode_gfa.pl: {e.stderr}")
            if os.path.exists(vcf_temp_dir):
                shutil.rmtree(vcf_temp_dir)
            raise
        except Exception as e:
            logger.error(f"Error processing {vcf_file.path.name}: {str(e)}")
            if os.path.exists(vcf_temp_dir):
                shutil.rmtree(vcf_temp_dir)
            raise
    
    def process_all_vcfs(self) -> str:
        """Process all VCF files in parallel.
        
        Returns:
            Path to combined output file
        """
        vcf_files = self.get_vcf_files()
        if not vcf_files:
            raise ValueError(f"No VCF files found in {self.vcf_dir}")
            
        logger.info(f"Found {len(vcf_files)} VCF files")
        
        # Allocate threads
        vcf_files = self.allocate_threads(vcf_files)
        
        # Process VCF files in parallel
        output_files = []
        with ThreadPoolExecutor(max_workers=len(vcf_files)) as executor:
            future_to_vcf = {
                executor.submit(self.process_vcf, vcf): vcf 
                for vcf in vcf_files
            }
            
            for future in as_completed(future_to_vcf):
                vcf = future_to_vcf[future]
                try:
                    output_file = future.result()
                    output_files.append(output_file)
                    logger.info(f"Completed processing {vcf.path.name}")
                except Exception as e:
                    logger.error(f"Failed to process {vcf.path.name}: {str(e)}")
                    raise
        
        # Combine results
        combined_output = self.output_dir / "combined_vcf_sequences.fa"
        with open(combined_output, 'w') as outfile:
            for file_path in output_files:
                with open(file_path) as infile:
                    outfile.write(infile.read())
                os.remove(file_path)
        # Cleanup
        if self.temp_dir.exists():
            import shutil
            shutil.rmtree(self.temp_dir)
            logger.info("Cleaned up all temporary directories")
        
        return str(combined_output)

def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Process multiple VCF files in parallel"
    )
    
    parser.add_argument(
        '--vcf-dir',
        required=True,
        help='Directory containing VCF files'
    )
    
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of threads to use'
    )
    
    args = parser.parse_args()
    
    try:
        processor = VCFProcessor(
            args.vcf_dir,
            args.output_dir,
            args.threads
        )
        
        output_file = processor.process_all_vcfs()
        logger.info(f"Processing complete. Results written to: {output_file}")
        
    except Exception as e:
        logger.error(f"Processing failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()