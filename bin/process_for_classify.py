#!/usr/bin/env python3
"""
Process and classify transposable elements based on their characteristics.

This script processes TE sequences and assigns them to appropriate classification
hierarchies based on predefined rules and characteristics.
"""

import csv
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

# Calibrated on A. thaliana PAF families that RepeatClassifier did classify,
# rerun with private KAnalyze temp dirs: at 0.9 precision against
# RepeatClassifier rises from 0.40 to 0.74 while about 40% fall back to Unknown.
DEFAULT_MIN_CONFIDENCE = 0.9

@dataclass
class ClassificationRule:
    """Data class for TE classification rules."""
    source_name: str
    target_class: str

class TEClassifier:
    """Handles the classification of transposable elements."""
    
    # ClassifyTE predicts Wicker nodes (nodes/tree.txt). Map them onto the
    # RepeatMasker Class/Family vocabulary so RepeatMasker's .tbl summary can
    # bucket them; bare node names are silently counted as Unclassified.
    # SubclassI is the DNA-transposon (TIR) subclass, not retrotransposon Class I.
    CLASSIFICATION_HIERARCHY = {
        # nodes/tree.txt
        'Retrotransposon': 'Unknown',
        'LTR': 'LTR',
        'Copia': 'LTR/Copia',
        'gypsy': 'LTR/Gypsy',
        'Bel-Pao': 'LTR/Pao',
        'DIRS': 'LTR/DIRS',
        'LINE': 'LINE',
        'R2': 'LINE/R2',
        'RTE': 'LINE/RTE',
        'Jockey': 'LINE/I-Jockey',
        'L1': 'LINE/L1',
        'I': 'LINE/I',
        'SINE': 'SINE',
        'tRNA': 'SINE/tRNA',
        '7SL': 'SINE/7SL',
        '5S': 'SINE/5S',
        'DNA transposon': 'DNA',
        'SubclassI': 'DNA',
        'TIRS': 'DNA',
        'Tc1-Mariner': 'DNA/TcMar-Tc1',
        'hAT': 'DNA/hAT',
        'Mutator': 'DNA/MULE-MuDR',
        'Merlin': 'DNA/Merlin',
        'Transib': 'DNA/CMC-Transib',
        'P': 'DNA/P',
        'PiggyBac': 'DNA/PiggyBac',
        'PIF-Harbinger': 'DNA/PIF-Harbinger',
        'CACTA': 'DNA/CMC-EnSpm',
        # Other Wicker nodes some ClassifyTE models emit
        'Retrovirus': 'LTR/ERV',
        'Ngaro': 'LTR/Ngaro',
        'TIR': 'DNA',
        'SubclassII': 'DNA',
        'Helitron': 'RC/Helitron',
        'Maverick': 'DNA/Maverick',
    }

    def __init__(self):
        """Initialize the classifier."""
        self.rules = [
            ClassificationRule(k, v) 
            for k, v in self.CLASSIFICATION_HIERARCHY.items()
        ]

    def get_classification(self, keyword: str) -> Optional[str]:
        """
        Get the RepeatMasker classification for a ClassifyTE label.
        
        Args:
            keyword: The ClassifyTE predicted label
            
        Returns:
            The RepeatMasker Class/Family, or None if the label is not known
        """
        return self.CLASSIFICATION_HIERARCHY.get(keyword)

class TEProcessor:
    """Processes TE sequences and their classifications."""
    
    def __init__(self, classifier: TEClassifier, min_confidence: float = 0.0):
        """
        Initialize the processor.
        
        Args:
            classifier: TEClassifier instance for handling classifications
            min_confidence: Minimum cumulative path probability for a node to be
                reported; deeper nodes below it back off to their ancestor
        """
        self.classifier = classifier
        self.min_confidence = min_confidence
        self.processed_lines: List[str] = []
        self.backoff_counts: Dict[str, int] = {}

    def choose_label(self, label: str, path: Optional[str], path_probs: Optional[str]) -> str:
        """
        Pick the deepest node on the predicted path whose cumulative probability
        reaches min_confidence; return 'Unknown' if even the first level does not.
        Without path probabilities (older ClassifyTE) the predicted label is kept.
        """
        if self.min_confidence <= 0 or not path or not path_probs:
            return label
        nodes = path.split('|')
        probs = [float(p) for p in path_probs.split('|')]
        if len(nodes) != len(probs):
            raise ValueError(f"Path and probabilities differ in length: {path} / {path_probs}")
        chosen = 'Unknown'
        cumulative = 1.0
        for node, prob in zip(nodes, probs):
            cumulative *= prob
            if cumulative < self.min_confidence:
                break
            chosen = node
        if chosen != label:
            key = f"{label}->{chosen}"
            self.backoff_counts[key] = self.backoff_counts.get(key, 0) + 1
        return chosen

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
        
        # ClassifyTE writes one prediction per line: "<sequence id>,<label>".
        unmapped: Dict[str, int] = {}
        with input_file.open('r') as file:
            reader = csv.reader(file)
            header = next(reader, None)
            if header is None:
                raise ValueError(f"Empty prediction file: {input_file}")
            has_probs = len(header) >= 4
            if self.min_confidence > 0 and not has_probs:
                logger.warning("Prediction file has no path probabilities; "
                               "--min-confidence is ignored (update ClassifyTE's evaluate.py)")
            for row in reader:
                if len(row) < 2 or not row[0].startswith('TE_'):
                    continue
                seq_id = row[0].split('#', 1)[0]
                label = row[1].strip()
                if has_probs and len(row) >= 4:
                    label = self.choose_label(label, row[2].strip(), row[3].strip())
                full_classification = self.classifier.get_classification(label)
                if full_classification is None:
                    # Keep the sequence visible as Unknown rather than guess.
                    unmapped[label] = unmapped.get(label, 0) + 1
                    full_classification = 'Unknown'
                # Combine_for_Two reads "<id>#Unknown,<class>"
                self.processed_lines.append(f"{seq_id}#Unknown,{full_classification}")

        for label, count in sorted(unmapped.items()):
            logger.warning(f"Unmapped ClassifyTE label '{label}' written as Unknown ({count} sequences)")
        for change, count in sorted(self.backoff_counts.items()):
            logger.info(f"Backed off below confidence {self.min_confidence}: {change} ({count} sequences)")

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
        '--min-confidence',
        type=float,
        default=DEFAULT_MIN_CONFIDENCE,
        help="Minimum cumulative ClassifyTE path probability for a node to be "
             "reported; less confident predictions back off to their ancestor "
             "(0 disables; default: %(default)s)"
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
        processor = TEProcessor(classifier, args.min_confidence)
        
        # Process input file
        processor.process_file(args.input_file)
        
        # Write results
        processor.write_output(args.output_file)
        
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()