#!/usr/bin/env python3
"""
VCF Processor for Pan_TE

Handles parallel processing of multiple VCF files with dynamic thread allocation.
"""

import os
import sys
import logging
import argparse
import math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
from collections import deque
from threading import Lock
from typing import List
from dataclasses import dataclass

logger = logging.getLogger(__name__)


def vcf_stem(path: Path) -> str:
    """Return a stable sample stem for .vcf and .vcf.gz inputs."""
    name = path.name
    lower = name.lower()
    if lower.endswith('.vcf.gz'):
        return name[:-7]
    if lower.endswith('.vcf'):
        return name[:-4]
    return path.stem

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

        if self.total_threads < 1:
            raise ValueError("total_threads must be positive")
        if not self.vcf_dir.is_dir():
            raise ValueError(f"VCF directory does not exist: {self.vcf_dir}")
        
    def get_vcf_files(self) -> List[VCFFile]:
        """Get list of VCF files with their sizes."""
        vcf_files = []
        paths = list(self.vcf_dir.glob("*.vcf"))
        paths.extend(self.vcf_dir.glob("*.vcf.gz"))
        seen_stems = set()
        for file_path in paths:
            stem = vcf_stem(file_path)
            if stem in seen_stems:
                raise ValueError(
                    f"Ambiguous VCF inputs share output stem {stem!r}; "
                    "keep only one of the .vcf/.vcf.gz files")
            seen_stems.add(stem)
            size = file_path.stat().st_size
            vcf_files.append(VCFFile(file_path, size))
        return sorted(vcf_files, key=lambda x: (-x.size, x.path.name))
    
    def allocate_threads(self, vcf_files: List[VCFFile]) -> List[VCFFile]:
        if not vcf_files:
            return []
        
        # If there are more files than threads, files are queued and each
        # active decoder receives one worker.  This keeps the total worker
        # count bounded without ever passing zero/negative --threads.
        if len(vcf_files) > self.total_threads:
            for vcf in vcf_files:
                vcf.allocated_threads = 1
            return vcf_files

        total_size = sum(vcf.size for vcf in vcf_files)
        weights = ([vcf.size / total_size for vcf in vcf_files]
                   if total_size > 0
                   else [1 / len(vcf_files)] * len(vcf_files))
        extra_threads = self.total_threads - len(vcf_files)
        raw_extra = [weight * extra_threads for weight in weights]
        assigned_extra = [math.floor(value) for value in raw_extra]
        remainder = extra_threads - sum(assigned_extra)
        order = sorted(
            range(len(vcf_files)),
            key=lambda i: (-(raw_extra[i] - assigned_extra[i]),
                           vcf_files[i].path.name))
        for i in order[:remainder]:
            assigned_extra[i] += 1
        for vcf, extra in zip(vcf_files, assigned_extra):
            vcf.allocated_threads = 1 + extra

        return vcf_files
    
    def process_all_vcfs(self) -> str:
        vcf_files = self.get_vcf_files()
        if not vcf_files:
            raise ValueError(f"No VCF files found in {self.vcf_dir}")
                
        vcf_files = self.allocate_threads(vcf_files)
        progress_dict = {vcf.path.name: 0 for vcf in vcf_files}
        progress_lock = Lock()
        
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
                filled = min(bar_length, int(bar_length * prog / 100))
                bar = '=' * filled
                if filled < bar_length:
                    bar += '>' + ' ' * (bar_length - filled - 1)
                
                print(f"\r\033[K{name:<20} [{bar}] {prog}%")
                lines_printed += 1

            display_progress.last_lines = lines_printed
            
            sys.stdout.flush()

        def process_single_vcf(vcf_file: VCFFile) -> str:
            stem = vcf_stem(vcf_file.path)
            vcf_temp_dir = os.path.join(self.temp_dir, stem)
            os.makedirs(vcf_temp_dir, exist_ok=True)
            output_file = os.path.join(vcf_temp_dir, f"{stem}.fa")

            cmd = [
                "decode_gfa.pl",
                "--vcf", str(vcf_file.path),
                "--threads", str(vcf_file.allocated_threads),
                "--out", vcf_temp_dir
            ]
            
            try:
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE,
                    text=True
                )
                stderr_tail = deque(maxlen=30)
                while True:
                    line = process.stderr.readline()
                    if not line and process.poll() is not None:
                        break
                    if line:
                        stderr_tail.append(line.rstrip())
                    if line.startswith('PROGRESS:'):
                        try:
                            percent = int(line.strip().split(':')[1])
                            progress_dict[vcf_file.path.name] = percent
                            with progress_lock:
                                display_progress()
                        except (ValueError, IndexError):
                            continue

                returncode = process.wait()
                if returncode != 0:
                    detail = '\n'.join(stderr_tail) or 'no stderr output'
                    raise RuntimeError(
                        f"decode_gfa.pl failed for {vcf_file.path.name} "
                        f"(exit {returncode}):\n{detail}")
                if not os.path.isfile(output_file):
                    raise RuntimeError(
                        f"decode_gfa.pl produced no output for "
                        f"{vcf_file.path.name}")
                return output_file
                
            except Exception as e:
                logger.error(f"Error processing {vcf_file.path.name}: {str(e)}")
                raise

        output_files = {}
        max_workers = min(len(vcf_files), self.total_threads)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_vcf = {executor.submit(process_single_vcf, vcf): vcf
                            for vcf in vcf_files}
            
            for future in as_completed(future_to_vcf):
                vcf = future_to_vcf[future]
                try:
                    output_file = future.result()
                    output_files[vcf.path] = output_file
                except Exception as e:
                    logger.error(f"Failed to process {vcf.path.name}: {str(e)}")
                    raise

        logger.info("All VCF files processed successfully")

        combined_output = os.path.join(self.output_dir, "combined_vcf_sequences.fa")
        with open(combined_output, 'w') as outfile:
            for vcf in vcf_files:
                file_path = output_files[vcf.path]
                with open(file_path) as infile:
                    outfile.write(infile.read())
                os.remove(file_path)

        return combined_output

def main():
    """Main execution function."""
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s')
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
