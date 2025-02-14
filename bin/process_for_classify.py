#!/usr/bin/env python3
"""
Process and classify transposable elements based on their characteristics.

This script processes TE sequences and assigns them to appropriate classification
hierarchies based on predefined rules and characteristics.
"""

import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class ClassificationRule:
    """Data class for TE classification rules."""
    source_name: str
    target_class: str

class TEClassifier:
    """Handles the classification of transposable elements."""
    
    CLASSIFICATION_HIERARCHY = {
        # LTR Elements
        'gypsy': 'LTR/Gypsy',
        'Copia': 'LTR/Copia',
        'Bel-Pao': 'LTR/Bel-Pao',
        'Retrovirus': 'LTR/Retrovirus',
        
        # Non-LTR Elements
        'DIRS': 'non-LTR/DIRS',
        'Ngaro': 'non-LTR/Ngaro',
        'VIPER': 'non-LTR/VIPER',
        'LINE': 'non-LTR/LINE',
        'SINE': 'non-LTR/SINE',
        'R2': 'non-LTR/R2',
        'RTE': 'non-LTR/RTE',
        'Jockey': 'non-LTR/Jockey',
        'L1': 'non-LTR/L1',
        
        # DNA Transposons
        'TIR': 'DNA/TIR',
        'Helitron': 'DNA/Helitron',
        'Maverick': 'DNA/Maverick',
        'Tc1-Mariner': 'DNA/Tc1-Mariner',
        'hAT': 'DNA/hAT',
        'Mutator': 'DNA/Mutator',
        'Merlin': 'DNA/Merlin',
        'PiggyBac': 'DNA/PiggyBac',
        'CACTA': 'DNA/CACTA',
        
        # General Categories
        'LTR': 'LTR',
        'non-LTR': 'non-LTR',
        'SubclassI': 'DNA',
        'SubclassII': 'DNA',
    }

    def __init__(self):
        """Initialize the classifier."""
        self.rules = [
            ClassificationRule(k, v) 
            for k, v in self.CLASSIFICATION_HIERARCHY.items()
        ]

    def get_classification(self, keyword: str) -> str:
        """
        Get the full classification for a given keyword.
        
        Args:
            keyword: The classification keyword to look up
            
        Returns:
            The full classification path
        """
        return self.CLASSIFICATION_HIERARCHY.get(keyword, keyword)

class TEProcessor:
    """Processes TE sequences and their classifications."""
    
    def __init__(self, classifier: TEClassifier):
        """
        Initialize the processor.
        
        Args:
            classifier: TEClassifier instance for handling classifications
        """
        self.classifier = classifier
        self.processed_lines: List[str] = []

    def process_file(self, input_file: Path) -> None:
        """
        Process the input file containing TE sequences.
        
        Args:
            input_file: Path to the input file
            
        Raises:
            FileNotFoundError: If input file doesn't exist
            ValueError: If file format is invalid
        """
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        logger.info(f"Processing file: {input_file}")
        
        with input_file.open('r') as file:
            # Skip header line
            next(file)
            
            # Process file two lines at a time
            while True:
                try:
                    id_line = next(file).strip()
                    class_line = next(file, '').strip()
                    
                    if not class_line:  # End of file
                        break
                        
                    if id_line.startswith('TE_'):
                        # Get classification keyword (remove comma)
                        classification_keyword = class_line[1:]
                        
                        # Map to full classification
                        full_classification = self.classifier.get_classification(
                            classification_keyword
                        )
                        
                        # Create processed line
                        processed_line = f"{id_line}#{full_classification}"
                        self.processed_lines.append(processed_line)
                        
                except StopIteration:
                    break
                except Exception as e:
                    logger.error(f"Error processing lines: {e}")
                    continue

        logger.info(f"Processed {len(self.processed_lines)} entries")

    def write_output(self, output_file: Path) -> None:
        """
        Write processed data to output file.
        
        Args:
            output_file: Path to the output file
            
        Raises:
            IOError: If writing to output file fails
        """
        logger.info(f"Writing output to: {output_file}")
        
        try:
            with output_file.open('w') as file:
                for line in self.processed_lines:
                    file.write(f"{line}\n")
                    
            logger.info("Output written successfully")
            
        except IOError as e:
            logger.error(f"Failed to write output: {e}")
            raise

def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Process and classify transposable elements"
    )
    
    parser.add_argument(
        'input_file',
        type=Path,
        help="Path to input file containing TE sequences"
    )
    
    parser.add_argument(
        'output_file',
        type=Path,
        help="Path to output file for processed sequences"
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Enable verbose logging"
    )
    
    return parser.parse_args()

def main() -> None:
    """Main execution function."""
    args = parse_arguments()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    try:
        # Initialize classifier and processor
        classifier = TEClassifier()
        processor = TEProcessor(classifier)
        
        # Process input file
        processor.process_file(args.input_file)
        
        # Write results
        processor.write_output(args.output_file)
        
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()