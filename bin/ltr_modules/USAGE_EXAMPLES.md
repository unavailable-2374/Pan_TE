# LTR_QualityValidator Usage Examples

## Quick Start

```python
from ltr_modules.LTR_QualityValidator import LTRQualityValidator
from Bio import SeqIO

# Initialize validator
validator = LTRQualityValidator(
    min_ltr_len=100,      # Minimum LTR length
    max_ltr_len=10000,    # Maximum LTR length
    min_ltr_similarity=0.50  # Minimum similarity for ancient elements
)

# Load sequences
sequences = SeqIO.parse("ltr_candidates.fasta", "fasta")

# Identify LTR sequences
for seq_record in sequences:
    result = validator._score_ltr_structure(seq_record)

    print(f"Sequence: {seq_record.id}")
    print(f"  Type: {result['sequence_type']}")
    print(f"  Category: {result['category']}")
    print(f"  Terminal similarity: {result['terminal_similarity']:.3f}")
    print(f"  Detection method: {result['detection_method']}")
    print(f"  Score: {result['total_score']:.1f}")
    print()
```

## Understanding Results

### Result Dictionary Structure

```python
{
    'total_score': 75.5,              # Overall quality score (0-100)
    'category': 'confirmed',           # confirmed|probable|possible|not_ltr
    'sequence_type': 'complete_element', # single_ltr|complete_element
    'detection_method': 'similarity_based_medium',  # Detection strategy used
    'terminal_similarity': 0.850,     # LTR similarity (0-1)
    'ltr_length': 350,                # Detected LTR length
    'components': {                   # Individual feature scores
        'ltr_similarity': 85.0,
        'tsd': 100.0,
        'boundary_motifs': 100.0
    },
    'issues': []                      # List of detected issues
}
```

### Categories Explained

- **confirmed**: High-confidence LTR (score ≥70, similarity ≥80%)
- **probable**: Likely LTR (score ≥50, similarity ≥50%)
- **possible**: Low-confidence LTR (score ≥35)
- **not_ltr**: Does not meet LTR criteria (score <35)

### Detection Methods

- **length_based_too_short**: Sequence <300bp (too short for complete element)
- **similarity_based_short**: 300-800bp with 35% threshold
- **similarity_based_medium**: 800-3000bp with 40% threshold
- **similarity_based_long**: >3000bp with 40% threshold

## Example Scenarios

### Scenario 1: Recent Complete Element

```python
# 5kb element with 95% LTR similarity
result = {
    'sequence_type': 'complete_element',
    'category': 'confirmed',
    'terminal_similarity': 0.95,
    'detection_method': 'similarity_based_long',
    'total_score': 97.5
}
# Interpretation: High-quality recent LTR retrotransposon
```

### Scenario 2: Ancient Complete Element

```python
# 3kb element with 55% LTR similarity
result = {
    'sequence_type': 'complete_element',
    'category': 'probable',
    'terminal_similarity': 0.55,
    'detection_method': 'similarity_based_medium',
    'total_score': 52.0
}
# Interpretation: Ancient element with degraded LTRs
```

### Scenario 3: Single LTR Consensus

```python
# 800bp consensus sequence with TG...CA motifs but low terminal similarity
result = {
    'sequence_type': 'single_ltr',
    'category': 'confirmed',
    'terminal_similarity': 0.28,
    'detection_method': 'similarity_based_medium',
    'total_score': 78.0
}
# Interpretation: High-quality single LTR (not complete element)
```

### Scenario 4: Borderline Case

```python
# 1500bp sequence with 42% similarity (just above 40% threshold)
result = {
    'sequence_type': 'complete_element',  # ≥40% threshold
    'category': 'possible',
    'terminal_similarity': 0.42,
    'detection_method': 'similarity_based_medium',
    'total_score': 35.0
}
# Interpretation: Borderline complete element, may need manual review
```

## Filtering Recommendations

### Conservative Filter (High Precision)

```python
# Keep only high-confidence complete elements
def is_high_confidence(result):
    return (result['sequence_type'] == 'complete_element' and
            result['category'] in ['confirmed', 'probable'] and
            result['terminal_similarity'] >= 0.60)
```

### Balanced Filter (Moderate Precision/Recall)

```python
# Include ancient elements but filter obvious noise
def is_likely_ltr(result):
    return (result['category'] in ['confirmed', 'probable', 'possible'] and
            result['total_score'] >= 40)
```

### Permissive Filter (High Recall)

```python
# Include all potential LTR sequences for further analysis
def is_potential_ltr(result):
    return result['category'] != 'not_ltr'
```

## Batch Processing Example

```python
def classify_ltr_library(input_fasta, output_prefix):
    """
    Classify LTR sequences and write to separate files
    """
    validator = LTRQualityValidator()

    confirmed = []
    probable = []
    possible = []
    not_ltr = []

    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        result = validator._score_ltr_structure(seq_record)
        category = result['category']

        # Add classification info to sequence description
        seq_record.description = (
            f"{seq_record.description} "
            f"type={result['sequence_type']} "
            f"sim={result['terminal_similarity']:.3f} "
            f"score={result['total_score']:.1f}"
        )

        if category == 'confirmed':
            confirmed.append(seq_record)
        elif category == 'probable':
            probable.append(seq_record)
        elif category == 'possible':
            possible.append(seq_record)
        else:
            not_ltr.append(seq_record)

    # Write to separate files
    SeqIO.write(confirmed, f"{output_prefix}_confirmed.fasta", "fasta")
    SeqIO.write(probable, f"{output_prefix}_probable.fasta", "fasta")
    SeqIO.write(possible, f"{output_prefix}_possible.fasta", "fasta")
    SeqIO.write(not_ltr, f"{output_prefix}_filtered.fasta", "fasta")

    # Print statistics
    total = len(confirmed) + len(probable) + len(possible) + len(not_ltr)
    print(f"Confirmed: {len(confirmed)} ({len(confirmed)*100/total:.1f}%)")
    print(f"Probable: {len(probable)} ({len(probable)*100/total:.1f}%)")
    print(f"Possible: {len(possible)} ({len(possible)*100/total:.1f}%)")
    print(f"Filtered: {len(not_ltr)} ({len(not_ltr)*100/total:.1f}%)")

# Usage
classify_ltr_library("all_candidates.fasta", "classified")
```

## Diagnostic Logging

### Enable DEBUG Logging

```python
import logging

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Now validator will output detailed detection decisions
validator = LTRQualityValidator()
result = validator._score_ltr_structure(seq_record)

# Example DEBUG output:
# Sequence type detection for seq001: type=complete_element,
# method=similarity_based_medium, similarity=0.85, ltr_len=350bp, total_len=2000bp
# Using pre-detected LTR similarity: 0.850 (350bp)
# LTR identification: seq001 = confirmed (score=85.5)
```

## Advanced: Custom Thresholds

```python
# For specific use cases, you can adjust parameters
validator_strict = LTRQualityValidator(
    min_ltr_len=200,       # Longer minimum LTR
    max_ltr_len=5000,      # Shorter maximum LTR
    min_ltr_similarity=0.70  # Higher similarity requirement
)

validator_permissive = LTRQualityValidator(
    min_ltr_len=50,        # Very short LTRs allowed
    max_ltr_len=15000,     # Very long LTRs allowed
    min_ltr_similarity=0.40  # Lower similarity threshold
)
```

## Integration with Downstream Analysis

### Example: RepeatMasker Library Preparation

```python
def prepare_repeatmasker_library(input_fasta, output_fasta):
    """
    Filter and prepare high-quality LTR library for RepeatMasker
    """
    validator = LTRQualityValidator()
    high_quality = []

    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        result = validator._score_ltr_structure(seq_record)

        # Only include confirmed and probable complete elements
        if (result['sequence_type'] == 'complete_element' and
            result['category'] in ['confirmed', 'probable']):

            # Rename with classification info
            seq_record.id = f"{seq_record.id}#LTR/Unknown"
            high_quality.append(seq_record)

    SeqIO.write(high_quality, output_fasta, "fasta")
    print(f"Prepared {len(high_quality)} high-quality LTR sequences")

# Usage
prepare_repeatmasker_library("raw_ltrs.fasta", "repeatmasker_lib.fasta")
```

## Troubleshooting

### Issue: All sequences classified as single_ltr

**Possible causes:**
- Input contains consensus sequences, not complete elements
- Terminal similarity calculation fails (check sequence quality)
- Sequences are too short (<300bp)

**Solution:** Check `terminal_similarity` values in results

### Issue: Low scores despite high similarity

**Possible causes:**
- Missing TSD (reduces score by 35%)
- Missing boundary motifs (reduces score by 15%)

**Solution:** This is expected for solo LTRs or incomplete elements

### Issue: Inconsistent classification

**Possible causes:**
- Sequences near 40% similarity threshold
- High N content affecting similarity calculation

**Solution:** Check `detection_method` and `terminal_similarity` in results
