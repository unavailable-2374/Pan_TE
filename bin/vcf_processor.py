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
    
    def process_all_vcfs(self) -> str:
        vcf_files = self.get_vcf_files()
        if not vcf_files:
            raise ValueError(f"No VCF files found in {self.vcf_dir}")
                
        vcf_files = self.allocate_threads(vcf_files)
        progress_dict = {vcf.path.name: 0 for vcf in vcf_files}
        complete_dict = {vcf.path.name: False for vcf in vcf_files}
        
        def display_progress():
            if not hasattr(display_progress, "initialized"):
                display_progress.initialized = False
                display_progress.last_lines = 0
                display_progress.header_printed = False

            if not display_progress.header_printed:
                print("\nProcessing VCF files:")
                print("-" * 60)
                display_progress.header_printed = True

            if display_progress.last_lines > 0:
                print(f"\33[{display_progress.last_lines}A", end="")

            lines_printed = 0
            for name, prog in sorted(progress_dict.items()):
                bar_length = 40
                filled = int(bar_length * prog / 100)
                bar = '=' * filled + '>' + ' ' * (bar_length - filled - 1)
                
                print(f"\r\033[K{name:<20} [{bar}] {prog}%")
                lines_printed += 1

            display_progress.last_lines = lines_printed
            
            sys.stdout.flush()

        def process_single_vcf(vcf_file: VCFFile) -> str:
            vcf_temp_dir = os.path.join(self.temp_dir, vcf_file.path.stem)
            os.makedirs(vcf_temp_dir, exist_ok=True)
            output_file = os.path.join(vcf_temp_dir, f"{vcf_file.path.stem}.fa")
            final_output = os.path.join(self.output_dir, f"{vcf_file.path.stem}.fa")

            cmd = [
                "decode_gfa.pl",
                "--vcf", str(vcf_file.path),
                "--threads", str(vcf_file.allocated_threads),
                "--out", vcf_temp_dir
            ]
            
            try:
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    universal_newlines=True
                )
                
                while True:
                    line = process.stderr.readline()
                    if not line and process.poll() is not None:
                        break
                    if line.startswith('PROGRESS:'):
                        try:
                            percent = int(line.strip().split(':')[1])
                            progress_dict[vcf_file.path.name] = percent
                            display_progress()
                        except (ValueError, IndexError):
                            continue
                            
                if process.returncode == 0:
                    complete_dict[vcf_file.path.name] = True
                    
                return output_file
                
            except Exception as e:
                logger.error(f"Error processing {vcf_file.path.name}: {str(e)}")
                raise

        output_files = []
        with ThreadPoolExecutor(max_workers=len(vcf_files)) as executor:
            future_to_vcf = {executor.submit(process_single_vcf, vcf): vcf
                            for vcf in vcf_files}
            
            for future in as_completed(future_to_vcf):
                vcf = future_to_vcf[future]
                try:
                    output_file = future.result()
                    output_files.append(output_file)
                except Exception as e:
                    logger.error(f"Failed to process {vcf.path.name}: {str(e)}")
                    raise

        logger.info("All VCF files processed successfully")

        combined_output = os.path.join(self.output_dir, "combined_vcf_sequences.fa")
        with open(combined_output, 'w') as outfile:
            for file_path in output_files:
                with open(file_path) as infile:
                    outfile.write(infile.read())
                os.remove(file_path)

        return combined_output

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
