# Deployment Guide for Improved TE Consensus Pipeline

## Summary of Changes

The pipeline has been fundamentally redesigned with two major improvements:

### 1. Phase 1: Fixed Scoring System
- **Previous issue**: Extreme skewing with A=1, B=721, C=1703 due to dynamic percentile-based thresholds
- **Solution**: Fixed thresholds based on TE biological characteristics
  - A-class: score ≥ 0.65 (clear high-quality TEs)
  - B-class: 0.45 ≤ score < 0.65 (possible TEs needing verification)
  - C-class: score < 0.45 (likely not true TEs)
- **New weights**:
  - 40% copy number (most important for TEs)
  - 25% identity score (family conservation)
  - 20% complexity score (exclude simple repeats)
  - 10% length score (TE type appropriateness)
  - 5% boundary quality (structural features)

### 2. Phase 2: Expansion-Focused Design
- **Previous issue**: Reducing 960 inputs to 587 outputs (reduction approach)
- **Solution**: Complete redesign focusing on expansion and improvement
  - Direct use of Phase 1 RepeatMasker results (no re-running)
  - No artificial limits on copy numbers or chromosome distribution
  - Generates multiple consensus for different TE variants
  - Expected to produce 2000+ consensus sequences for comprehensive annotation

## Files to Upload

Upload the following modified files to your server:

### Core Files (Modified)
1. `/Users/caoshuo/Downloads/11-way/TCB/god/phase1_screening_optimized.py`
   - Fixed scoring system with TE-appropriate thresholds
   - Saves detailed RepeatMasker results for Phase 2

2. `/Users/caoshuo/Downloads/11-way/TCB/god/phase2_consensus_expansion.py` (NEW)
   - Completely new implementation focused on expansion
   - Uses Phase 1 RepeatMasker results directly
   - No artificial limits, generates more consensus

3. `/Users/caoshuo/Downloads/11-way/TCB/god/main.py`
   - Updated to use ConsensusExpansionBuilder instead of TEBiologyConsensusBuilder

### Optional Files (For Reference)
4. `/Users/caoshuo/Downloads/11-way/TCB/god/phase2_consensus_te_biology.py`
   - Previous Phase 2 implementation (keep for comparison)

5. `/Users/caoshuo/Downloads/11-way/TCB/god/phase2_consensus_correct.py`
   - Alternative Phase 2 implementation (keep for comparison)

## Upload Commands

```bash
# On your local machine, navigate to the directory
cd /Users/caoshuo/Downloads/11-way/TCB/god/

# Upload the modified files to your server
scp phase1_screening_optimized.py user@server:/path/to/TCB/god/
scp phase2_consensus_expansion.py user@server:/path/to/TCB/god/
scp main.py user@server:/path/to/TCB/god/
```

## Running the Improved Pipeline

### Step 1: Clear Previous Checkpoints
Since both Phase 1 and Phase 2 have been significantly modified, you must clear the checkpoints:

```bash
cd /path/to/TCB/god/
rm -rf checkpoints/phase1_complete.pkl
rm -rf checkpoints/phase2_complete.pkl
# Or clear all checkpoints:
rm -rf checkpoints/*.pkl
```

### Step 2: Run the Pipeline
```bash
python main.py \
  -r /path/to/repeatscout_output.fa \
  -g /path/to/genome.fa \
  -o output_dir \
  -t 32 \
  --keep-checkpoints
```

### Step 3: Monitor Progress
The pipeline will now show:
- Phase 1: More balanced A/B/C distribution (expected: A≈2000+, B≈50-200, C≈200-500)
- Phase 2: Expansion ratio >2x (expected: 3000-5000 consensus from ~2000 inputs)

## Expected Output

### Phase 1 Output
```
Classification thresholds (fixed) - A: >=0.65, B: >=0.45
Categorization: A=2128, B=48, C=249  # Example of balanced distribution
```

### Phase 2 Output
```
Phase 2 complete: Generated 4256 consensus sequences from 2176 inputs
Expansion ratio: 1.96x
```

### Final Output Files
```
output_dir/
├── consensus_masking.fa      # Comprehensive TE library for genome masking
├── consensus_analysis.fa     # Detailed TE library with metadata
├── statistics.txt            # Quality metrics and statistics
└── pipeline.log             # Detailed execution log
```

## Verification Steps

1. **Check Phase 1 distribution**:
   ```bash
   grep "Categorization:" pipeline.log
   ```
   Should show reasonable A/B/C distribution, not extreme skewing.

2. **Check Phase 2 expansion**:
   ```bash
   grep "Expansion ratio:" pipeline.log
   ```
   Should show ratio > 1.0, indicating expansion not reduction.

3. **Count final consensus sequences**:
   ```bash
   grep -c ">" output_dir/consensus_masking.fa
   ```
   Should be significantly more than the number of RepeatScout inputs.

## Troubleshooting

### If Phase 1 Still Shows Extreme Skewing
- Check that `phase1_screening_optimized.py` has the fixed thresholds:
  ```python
  A_THRESHOLD = 0.65
  B_THRESHOLD = 0.45
  ```

### If Phase 2 Generates Too Few Sequences
- Verify that `main.py` imports `ConsensusExpansionBuilder`
- Check that Phase 2 is using Phase 1 RepeatMasker results:
  ```bash
  grep "Using Phase 1 RepeatMasker results" pipeline.log
  ```

### If Pipeline Fails at Phase 2
- Ensure all utility modules are present and unchanged
- Check that `phase2_consensus_expansion.py` is in the correct location

## Performance Considerations

The new Phase 2 expansion approach will:
- Generate MORE consensus sequences (good for comprehensive annotation)
- Take slightly more time due to processing more variants
- Use more memory for storing multiple consensus per input
- Produce a more complete TE library for genome annotation

## Next Steps After Running

1. **Validate the expanded library**:
   - Check that consensus sequences cover expected TE families
   - Verify that subfamilies are properly represented

2. **Use for genome annotation**:
   - The expanded library is ideal for comprehensive genome masking
   - More variants mean better detection of diverged TE copies

3. **Optional post-processing**:
   - If the library is too large, can apply stricter filtering in Phase 3
   - Can cluster highly similar sequences post-hoc if needed

## Contact for Issues

If you encounter any issues:
1. Check the `pipeline.log` for detailed error messages
2. Ensure all modified files are correctly uploaded
3. Verify that checkpoints were cleared before running
4. Check that external tools (RepeatMasker, MAFFT) are accessible