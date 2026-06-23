# Pan_TE Pipeline Architecture Memory Document

## Project Overview

Pan_TE is a **pan-genome transposable element (TE) annotation pipeline** that combines three independent TE detection strategies (LTR retrotransposon detection, mdl-repeat, and TE-looker) into a unified workflow. It supports optional VCF-derived structural variant sequences and TE classification via ClassifyTE + RepeatClassifier.

**Step 4 has been migrated from RepGraph (Perl, BLAST-self-alignment + Leiden) to TE-looker (Rust, `dtr` CLI, tri-track aggregation: nhmmer + sourmash + FastGA → window-graph → Leiden CPM → consensus). TE-looker is targeted at D_cc 45–75% (the "twilight zone" / older divergent TEs); active/young TEs (D_cc > 75%) are caught by the LTR and mdl-repeat tracks.**

**Entry point**: `bin/Pan_TE` (Python3)
**Version**: 1.0.0
**Dev repo**: `/home/shuoc/tool/tmp/Pan_TE/`
**Production repo (on PATH)**: `/home/shuoc/tool/Pan_TE/`

---

## Pipeline Execution Flow

```
Pan_TE (main orchestrator)
  |
  |-- Step 1: Genome Processing (process_genome)
  |     |-- [Optional] VCF Processing (vcf_processor.py)
  |     |-- clean_seq_fast: uppercase + degenerate base cleaning, rename to chr1..chrN
  |     |-- samtools faidx + index
  |     +-- Checkpoint: genome.ok
  |
  |-- Step 2+3: Parallel Detection (ThreadPoolExecutor, 2 workers)
  |     |
  |     |-- [Thread A] LTR Detection (process_ltr)
  |     |     |-- LTR_detect: Look4LTRs C++ tool (genome size-based splitting)
  |     |     |-- build_ltr_library.py: GraphGroup family assignment + consensus
  |     |     +-- Checkpoint: LTR.ok
  |     |
  |     +-- [Thread B] mdl-repeat (process_mdl_repeat)
  |           |-- build_mdl (Perl): mdl-repeat wrapper
  |           |     |-- Dust + TRF masking
  |           |     |-- mdl-repeat (handles large genomes internally)
  |           |     +-- Checkpoint: mdl_repeat_complete
  |           |-- Refiner_mdl Pipeline (Python3, 3 phases)
  |           |     |-- Phase 0: Metadata triage + tiering
  |           |     |-- Phase 1: Tiered consensus polishing
  |           |     +-- Phase 2: Library split (masking + analysis)
  |           +-- Checkpoint: mdl_repeat.ok
  |
  |-- Step 4: TE-looker (process_te_looker)
  |     |-- Uses processed genome (genome/genome.fa)
  |     |-- dtr run (Rust binary at ~/tool/te-looker/target/release/dtr):
  |     |     |-- Stage 0: dustmasker ∪ TRF ∪ NUMT → soft-mask BED
  |     |     |-- Stage 1: tri-track evidence aggregation
  |     |     |     - Track 1: nhmmer with HMM library auto-built from mdl-repeat consensus (--dfam-hmm)
  |     |     |     - Track 2: sourmash gather (de novo by default; user can opt-in via --te-looker-extra-args "--track2-library …")
  |     |     |     - Track 3: FastGA self-alignment
  |     |     |-- Stage 2: window-graph (k-NN, default k=4, stride 100bp)
  |     |     |-- Stage 3: Leiden CPM dual-γ community detection
  |     |     |-- Stage 4: per-cluster consensus (abPOA + MAFFT + Refiner-HMM)
  |     |     |-- Stage 4.5: family-merger (sourmash pre-filter → blastn → IQ-TREE phylo gate)
  |     |     |-- Stage 5: boundary refinement + structural features (LTR/TIR/Helitron)
  |     |     +-- Stage 6: G1–G6 hard gates → families.fasta + provenance.json
  |     +-- Checkpoint: te-looker.ok (legacy RepGraph.ok still recognized on resume)
  |
  |-- Step 5: Combine & Classify (combine_results)
  |     |-- Collect consensus from all sources:
  |     |     |-- LTR: Look4LTRs/consensus_lib/ltr_library.fa
  |     |     |-- mdl-repeat: mdl_repeat/refiner_output/phase3_analysis_library.fa
  |     |     +-- TE-looker: te-looker/run/families.fasta (consensi.fa kept as legacy alias)
  |     |-- CD-HIT-EST deduplication (90% identity, local alignment, both strands; tight because each upstream source is already internally merged)
  |     |-- [Optional] Classification via run_Classifier:
  |     |     |-- RepeatClassifier (parallel)
  |     |     +-- ClassifyTE (conda env: ClassifyTE_env, parallel)
  |     |     +-- Combine_for_Two: merge both classification results
  |     +-- Checkpoint: Combine.ok
  |
  +-- Output: Combine/TEs.fa (classified TE library)
```

---

## Directory Structure (Runtime)

```
<output_dir>/
  genome/
    genome.fa              # Cleaned, indexed genome (unmasked, also passed to TE-looker)
    genome.fa.fai          # samtools index
  Look4LTRs/
    l4s.ok                 # LTR_detect checkpoint
    consensus_lib/
      ltr_library.fa       # Final LTR consensus library
  mdl_repeat/
    tmp/                   # mdl-repeat working files
      repeats.fa           # mdl-repeat output
      repeats.bed          # Instance BED file
      repeats.stats        # Family statistics
    refiner_output/
      phase3_analysis_library.fa     # Refiner_mdl output for Combine
      consensus_masking.fa           # Refiner_mdl output (fed to TE-looker as --track2-library)
    consensi.fa            # Copy of Refiner_mdl masking output (legacy alias)
  te-looker/                 # (or RepGraph/ if resuming from a legacy run)
    run/
      families.fasta       # TE-looker primary output (Class A/D/E gated families)
      provenance.json      # RFC 8785 canonicalized run provenance (BLAKE3 hashed)
      stage{1,3,5,6}_*.json # Per-stage diagnostic JSON
      summary.json         # Run summary
    consensi.fa            # Copy of families.fasta for Combine-stage discovery (legacy filename)
  Combine/
    raw_TEs_combined.fa    # All sources concatenated
    raw_TEs.fa             # Deduplicated TEs
    TEs.fa                 # Final classified TE library
  *.ok                     # Pipeline checkpoints
  pan_te.log               # Main pipeline log
```

---

## Key Components Detail

### 1. Genome Processor (`clean_seq_fast`)
- **Language**: Python3
- **Function**: Single-pass uppercase conversion + degenerate base (WRMKYSHBVDX) -> N replacement
- **Performance**: Streaming binary I/O with translation table, ~Mbp/s speed
- **Renaming**: All sequences renamed to chr1, chr2, ... chrN

### 2. LTR Detection (`LTR_detect` + `build_ltr_library.py`)
- **LTR_detect**: Python3 wrapper around Look4LTRs C++ tool
  - Genome <5GB: run on whole genome
  - 5-10GB: split into 2 parts
  - 10-15GB: 3 parts; 15-20GB: 4 parts
  - >20GB: sample 4 chunks of ~5GB each
  - Uses ProcessPoolExecutor for parallel part processing
- **build_ltr_library.py**: Builds consensus from .rtr files
  - Uses Look4LTRs GraphGroup family assignment
  - LTR identity-weighted consensus
  - Cross-family merging via FastGA
  - Multi-part: build per-part then CD-HIT-EST merge (80% identity)

### 3. mdl-repeat (`build_mdl`)
- **Language**: Perl (v3.0.0) wrapper around mdl-repeat C binary
- **Strategy**: Direct execution — mdl-repeat handles large genomes internally
  - No manual sampling or parallel splitting required
  - Dust + TRF masking applied to full genome as preprocessing
  - mdl-repeat uses MDL (Minimum Description Length) principle for family detection
- **Masking**: Dust + TRF hard masking before mdl-repeat (parallel TRF via ForkManager)
- **Output format**: `>R=N length=L copies=C mdl=M` FASTA headers
- **Auxiliary output**: Instance BED file, family statistics TSV
- **Binary location**: `/home/shuoc/tool/mdl-repeat/bin/mdl-repeat`
- **Dependencies**: dustmasker, trf, RMBlast (via RepModelConfig or environment)

### 4. Refiner (`bin/Refiner/`)
- **Language**: Python3
- **Config**: `config.py` (PipelineConfig dataclass, JSON serializable)
- **Phase 1** (`phase1_screening_optimized.py`): Sequence screening
  - Complexity scoring (DUST + Shannon entropy)
  - Length filtering (50-50000bp)
  - N% filtering (max 20%)
  - Trusts the repeat finder's initial filtering
- **Phase 2** (`phase2_chimera_splitting.py`): Chimera detection
  - RepeatMasker-based chimera identification
  - Conservative splitting strategy
  - Parallel processing (I/O intensive)
- **Phase 3** (`phase3_consensus_building.py`): Consensus building
  - BLASTN recruitment of genomic copies
  - MAFFT L-INS-i multiple sequence alignment
  - Quality-tiered consensus (high/medium/low)
  - Hierarchical clustering (scipy)
  - Output: masking library + analysis library
- **Phase 4** (`phase4_genome_masking.py`): Genome masking (optional, default off, legacy artifact)
  - RepeatMasker hard masking with consensus library
  - ID renaming for RepeatMasker compatibility
  - Was used for the old RECON/RepGraph path; TE-looker takes the unmasked genome and uses `--track2-library` instead

### 4b. Refiner_mdl (`bin/Refiner_mdl/`)
- **Language**: Python3
- **Purpose**: Refinement pipeline for mdl-repeat output (replaces Refiner for mdl-repeat data)
- **Purpose**: Refinement pipeline for mdl-repeat output, called directly by `build_mdl`
- **Config**: `config.py` (RefinerMdlConfig dataclass, 43 parameters)
- **Phase 0** (`phase0_triage.py`): Metadata parsing + quality triage
  - Parse mdl-repeat FASTA headers (R, length, copies, mdl)
  - Hard filter: length<50, single-copy, low entropy, high N%, mdl<=0
  - Adaptive tiering: P25/P75 percentile on mdl_per_copy → T1/T2/T3
  - Fragment assembly: BED co-occurrence (primary) or BLASTN coordinate co-occurrence (fallback)
  - Two-round CD-HIT-EST dedup (95% then 90%)
- **Phase 1** (`phase1_consensus.py`): Tiered consensus polishing
  - Batch BLASTN copy recruitment → samtools faidx extraction
  - T1 (fast): top-20 hits + MAFFT --auto + chimera ≥300bp
  - T2 (standard): top-50 hits + MAFFT L-INS-i + chimera ≥150bp
  - T3 (lightweight): top-10 hits + MAFFT --auto, no chimera
  - Chimera detection: bin coverage CV > 1.5 → split
  - Quality protection: low-occupancy alignment → keep original consensus
  - ProcessPoolExecutor parallel, sorted by length*copies descending
- **Phase 2** (`phase2_library_split.py`): Final QC + library split
  - QC: length, entropy, N% filters
  - **mdl-repeat native quality gate (`enable_mdl_quality_gate`, default on)**: mdl-repeat scores every family in its FASTA header (`tier=core/warn/reject`, `accept=exclusive/standalone/…`); `phase0_triage.parse_mdl_fasta` now parses these into `mdl_tier` / `mdl_accept`. On a large repeat-rich genome the bulk of written families are low-confidence (warn / standalone / few copies) and inflate the library, so this gate keeps a family iff `mdl_tier ∈ mdl_quality_keep_tiers ("core")` OR it carries TE structural signal (protein/LTR/TIR/known-homology) OR `copies ≥ mdl_quality_highcopy_keep (50)`; it drops only the no-signal ∩ low-tier ∩ low-copy intersection (audited to `mdl_quality_dropped.tsv`). SAFETY: skipped entirely if no structural signal could be computed (missing RepeatPeps/RepeatMasker libs) and never gates headers lacking `tier=` (older mdl-repeat builds). The structural-signal rescue protects divergent/dark-matter TEs that mdl flagged warn (pure `tier=core` alone over-drops). Validated: maize 20,159 → ~7.2 k / 15 Mb (all 5,382 TE-signal families kept); Arabidopsis end-to-end gate 17,775 → 4,720 (only warn dropped, 0 TE-signal families lost). This is the mdl analogue of the te-looker host-side gate — both follow the rule "a de novo tool over-produces; curate with the tool's own confidence + structural/copy evidence."
  - `consensus_masking.fa` = T1 + T2 (consumed by TE-looker as `--track2-library` for Stage 1 Track 2 sourmash evidence)
  - `phase3_analysis_library.fa` = T1 + T2 + T3 (for Combine classification)
  - `refiner_mdl_stats.json`: per-phase statistics (incl. `mdl_quality_gate` block)
- **Output**: `consensus_masking.fa`, `phase3_analysis_library.fa`
- **Reused from Refiner**: `utils/complexity_utils.py`, `utils/alignment_utils.py`

### 5. TE-looker (`dtr` Rust binary)
- **Language**: Rust (cargo workspace, 12 crates), CLI binary `dtr` (Dark TE Rebuilder)
- **Location**: `~/tool/te-looker/target/release/dtr` (resolved via `--te-looker-bin`, `TE_LOOKER_BIN` env var, or `dtr` on PATH)
- **Targeted operating range**: D_cc 45–75% (older divergent TEs). Active/young TEs are caught upstream by Look4LTRs and mdl-repeat.
- **Pipeline**: Stage 0 (dustmasker / TRF / NUMT) → Stage 1 (tri-track: nhmmer + sourmash + FastGA) → Stage 2 (window-graph k-NN) → Stage 3 (Leiden CPM dual-γ) → Stage 4 (abPOA + MAFFT + Refiner-HMM consensus) → Stage 4.5 (family merger: sourmash pre-filter + blastn + IQ-TREE phylo gate) → Stage 5 (boundary refine + LTR/TIR/Helitron structural features) → Stage 6 (Class A/D/E validation gates)
- **Stage 1 Track 1 seeding**: Pan_TE auto-builds an HMM library from mdl-repeat's `consensus_masking.fa` (via `hmmbuild --dna`, one HMM per consensus, packaged as a Stockholm-blocks file) and passes it as `--dfam-hmm`. This makes Track 1 (nhmmer) fast on mdl-repeat-known families and leaves Tracks 2 (sourmash) + 3 (FastGA) free for de novo discovery of the divergent / fragmented TEs that TE-looker is designed for. Suppress this auto-seeding by passing `--dfam-hmm` or `--seed-from-repeatscout` via `--te-looker-extra-args`.
- **Window-stride auto-adaptation**: Pan_TE picks `--window-stride` by genome-size tier (50 / 100 / 200 / 500 / 1000 bp for XS / S-M / L / XL / XXL) before invoking `dtr`. Stage 2 graph memory scales as |V| ≈ genome_bp / stride, so this is the dominant scaling knob (see `design_v2/SCALABILITY_STRATEGY.md`).
- **Min-count auto-scaling**: Pan_TE picks `--min-count` (seed = canonical 16-mer occurring ≥C times) by genome-size tier (50 / 100 / 200 / 300 / 500 for <300 Mb / <1 Gb / <3 Gb / <10 Gb / ≥10 Gb) via `choose_te_looker_min_count`. dtr's built-in default (20) over-seeds large genomes (validated floor on 135 Mb TAIR10 was ≥200). Override via `--te-looker-extra-args "--min-count N"`.
- **⚠ Installed `dtr` is DISCOVERY-ONLY (reality vs the spec above)**: the binary actually built (`~/tool/te-looker/core/target/release/dtr`, te-core 0.3.0) is the discovery core only — `te-discover` (A1 ungapped extend-to-consensus) → `te-refine` (spoa). The tri-track / Leiden / Stage 4.5 merger / Stage 6 gates / `provenance.json` described above are `docs/V4_DESIGN.md` design intent, **not implemented**. It accepts but ignores `--dfam-hmm` / `--window-stride`; it emits `disc.consensi.fa` + `disc.members.bed` + `families.fasta` with NO cross-seed merge and NO copy/TE gate, so on a large repeat-rich genome it MASSIVELY over-produces (maize: 58,711 families / 92 Mb, ~10× inflated — one element fragmented into ~10 adjacent sub-families because ungapped extension halts at the first divergence/indel breakpoint).
- **Host-side discovery gate (Pan_TE `_gate_and_merge_te_looker`, default on)**: because the installed dtr does not gate, `process_te_looker` does it. If `provenance.json` is present (a future gated dtr) the families.fasta is trusted as-is; otherwise the output is treated as discovery-grade and gated: drop `members < TE_LOOKER_MIN_COPIES (5)` or `len < TE_LOOKER_MIN_LEN (100)`, then collapse cross-seed redundancy with `cd-hit-est -c 0.80 -aS 0.6`. A genome-scaled family-count sanity bound (`> genome_Mb × 15`) hard-aborts rather than ship a bloated library. Validated: maize 58,711 → ~4.5 k / 8.7 Mb; Arabidopsis 522 → 143 / 0.30 Mb. (NOTE: this is a redundancy-collapse stage — total bp, not family count, is the false-positive metric; fragment stitching that conserves bp does not help. The A1 fragmentation root cause is not fixed at source; gapped-extension was tried and failed — fragmentation is divergence- not indel-driven.)
- **Resume**: if `te-looker/run/families.fasta` already exists and is non-empty, `process_te_looker` skips the multi-hour dtr discovery (and the mdl→HMM seed build) and only re-applies the gate. Delete that file to force full re-discovery.
- **Output (in `te-looker/run/`)**:
  - `families.fasta`: final TE consensus library (multi-class evidence gated)
  - `provenance.json`: RFC 8785 canonicalized provenance, BLAKE3 hashed (G10 byte-level reproducibility)
  - `stage1_*.json`, `stage3_*.json`, `stage5_*.json`, `stage6_*.json`: per-stage diagnostic JSON
  - `summary.json`: run summary
- **Key default parameters**: `--min-copies 5`, `--evidence-min-classes 2`, `--cluster-gamma 0.005` (macro), `--graph-k 4`, `--window-stride 100`, `--merger-method three-tier`, `--enable-structural true`. Override via `--te-looker-extra-args`.
- **Determinism**: All PRNG via `ChaCha20Rng` with seed derivation; `--seed` flag exposed.
- **Skip flags**: `--skip-sourmash`, `--skip-fastga`, `--skip-stage0`, `--skip-merger`, `--skip-stage5` for environments missing optional tools.
- **Design spec (authoritative)**: `/home/shuoc/tool/repgraph/design_v2/FINAL_PLAN.md`, `RUST_ARCHITECTURE.md`, `BR_round2.md`, `ED_round2.md`.

### 6. Classifier (`run_Classifier`)
- **Language**: Bash
- **Dual classification** (parallel):
  - RepeatClassifier (RepeatModeler tool)
  - ClassifyTE (Python, conda env: ClassifyTE_env, feature generation + ML evaluation)
- **Combination**: `Combine_for_Two` merges both results
- **Post-processing**: Standardize TE class names (LINE, LTR/Gypsy, etc.)
- **Timeout**: 2 hours per classification process

### 7. VCF Processor (`vcf_processor.py`)
- **Language**: Python3
- **Function**: Process structural variant VCF files to extract TE-related sequences
- **Thread allocation**: Proportional to file size
- **Parallel processing**: ThreadPoolExecutor

---

## Thread Allocation Strategy

The pipeline allocates threads as follows:
- **ltr_threads**: `max(1, total_threads // 2)`
- **mdl_threads**: `max(1, total_threads)` (full allocation, mdl-repeat is multi-threaded internally)
- LTR and mdl-repeat run in **parallel** via ThreadPoolExecutor(max_workers=2)
- TE-looker runs **sequentially** after both complete (uses full `--threads`; Stage 1 Track 2 needs mdl-repeat's consensus_masking.fa)

---

## Checkpoint System

All major steps have `.ok` checkpoint files in the output directory:
- Canonical: `genome.ok`, `LTR.ok`, `mdl-repeat.ok`, `te-looker.ok`, `combine.ok`
- Legacy names still recognized on resume: `mdl_repeat.ok`, `RepGraph.ok`, `Combine.ok`
- Sub-checkpoints: `l4s.ok` (Look4LTRs), `mdl_repeat_complete` (build_mdl), `dust_trf_masking` (build_mdl)
- Refiner / Refiner_mdl use their own checkpoint system (`checkpoints/` directory)
- TE-looker has its own provenance/state internally (no Pan_TE sub-checkpoints; re-run by deleting `te-looker.ok` and the `te-looker/run/` directory)
- Pipeline supports **resume**: skips completed steps on re-run

---

## External Tool Dependencies

### Complete Tool Chain (verified available)

| Tool | Used By | Purpose | Location |
|------|---------|---------|----------|
| **Genome Processing** ||||
| clean_seq_fast | Pan_TE | Uppercase + degenerate base cleaning | Pan_TE/bin/ (Python3) |
| samtools | Pan_TE, build_mdl, Refiner | FASTA indexing/extraction | hap.py-build/bin/ |
| index | Pan_TE | Additional FASTA indexing | Pan_TE/bin/ (Perl) |
| **LTR Detection** ||||
| look4ltrs | LTR_detect | LTR retrotransposon detection (C++) | Look4LTRs/bin/ |
| bedtools | build_ltr_library, SamplingTrack | BED ops, sequence extraction | PGTA conda env |
| FastGA | build_ltr_library | Cross-family LTR alignment | /home/shuoc/tool/FASTGA/ |
| **mdl-repeat** ||||
| mdl-repeat | build_mdl | De novo repeat discovery (MDL principle) | /home/shuoc/tool/mdl-repeat/bin/ |
| dustmasker | build_mdl | Low-complexity masking | PGTA conda env |
| trf | build_mdl | Tandem repeat masking | PGTA conda env |
| **Alignment/Masking** ||||
| rmblastn | Refiner | Sensitive repeat alignment | PGTA conda env |
| blastn | Refiner Phase 3, te-looker Stage 4.5 | Genomic copy recruitment / merger | PGTA conda env |
| makeblastdb | Refiner, te-looker Stage 4.5 | BLAST database creation | PGTA conda env |
| RepeatMasker | Refiner Phase 2/4 | TE annotation/masking | PGTA conda env |
| mafft | Refiner Phase 3, build_ltr_library, te-looker Stage 4 | Multiple sequence alignment | PGTA conda env |
| **TE-looker (`dtr`)** ||||
| dtr | Pan_TE Step 4 | Tri-track de novo TE discovery (Rust) | ~/tool/te-looker/target/release/ |
| hmmbuild | Pan_TE Step 4 (host-side, mdl-repeat → HMM seed) + te-looker Stage 4 | HMM model construction | PGTA conda env |
| nhmmer | te-looker Stage 1 Track 1 + Stage 5 boundary refine | HMM profile sequence search | PGTA conda env |
| sourmash | te-looker Stage 1 Track 2, Stage 4.5 | k-mer sketch gather | PGTA conda env |
| FastGA | te-looker Stage 1 Track 3 | Self-alignment (Track 3) | ~/.cargo/lib/sweepga/ |
| abpoa | te-looker Stage 4 | Partial-order MSA consensus | ~/tool/abPOA/bin/ |
| iqtree | te-looker Stage 4.5 phylo gate | SH-aLRT + UFBoot phylogenetic gate | PGTA conda env |
| **Deduplication** ||||
| cd-hit-est | Combine, build_ltr_library | Sequence clustering (80% identity) | PGTA conda env |
| **Classification** ||||
| renameTE | run_Classifier | TE sequence renaming | Pan_TE/bin/ (Perl) |
| RepeatClassifier | run_Classifier | RepeatModeler-based classification | PGTA conda env |
| ClassifyTE tools | run_Classifier | ML-based classification | ClassifyTE_env conda env |
| process_for_classify.py | run_Classifier | Classification result processing | Pan_TE/bin/ |
| Combine_for_Two | run_Classifier | Merge dual classification results | Pan_TE/bin/ |

### Conda Environments
- **PGTA**: Main environment with most bioinformatics tools (rmblastn, RepeatMasker, nhmmer, sourmash, iqtree, etc.)
- **ClassifyTE_env**: Isolated environment for ClassifyTE ML classification (Python 3.8)
- TE-looker's `dtr` binary expects FastGA + abpoa on `PATH`; Pan_TE's `te_looker_env()` prepends the `dtr` binary directory plus optional `FASTGA_BIN` / `ABPOA_BIN` / `TE_LOOKER_TOOL_PATH` env vars to PATH before invocation.

---

## Input/Output

### Required Input
- `--genome`: FASTA file (.fa/.fasta/.fna/.fas), optionally gzip/bzip2 compressed

### Optional Input
- `--model-dir`: ClassifyTE model directory (for classification)
- `--vcf-dir`: VCF files directory (structural variants)
- `--threads`: Number of threads (default: 4, recommend multiples of 4)
- `-M`: Memory limit in MB
- `--fragment_size`: Fragment size (default: 40000)

### Final Output
- `Combine/TEs.fa`: Classified TE consensus library (if model-dir provided)
- `Combine/raw_TEs.fa`: Deduplicated TE library (without classification)

---

## Key Design Decisions

1. **Three-pronged detection**: LTR-specific (Look4LTRs, active/young) + de novo repeat (mdl-repeat, broad-spectrum MDL) + TE-looker (D_cc 45–75% twilight zone) for complementary coverage across the divergence spectrum
2. **Tri-track aggregation in TE-looker**: Replaces the old single-pass BLAST self-alignment (RepGraph) with three orthogonal evidence tracks (nhmmer profile / sourmash sketch / FastGA alignment) aggregated on a window grid before community detection. This is the algorithmic upgrade that closes the "dark TE matter" gap.
3. **mdl-repeat → TE-looker seeding via Track 1 HMM**: mdl-repeat's `consensus_masking.fa` is converted to a multi-HMM library by `hmmbuild` and passed to `dtr run --dfam-hmm`, so the nhmmer track recognizes known families efficiently. Tracks 2 (sourmash) and 3 (FastGA) run de novo on the same genome. This replaces both the old "soft-mask the genome before RepGraph" pattern AND an earlier intermediate design that fed the same FASTA to Track 2 as a sketch query (which biased Track 2 and left Track 1 silent).
4. **mdl-repeat handles scaling**: mdl-repeat internally manages large genomes — no manual sampling or parallel splitting needed
5. **Checkpoint-based resume**: Every major step checkpointed, full pipeline resumable (legacy `RepGraph.ok` checkpoint name still recognized for backwards compat)
6. **Refiner_mdl as quality gate**: mdl-repeat raw output refined through 3-phase pipeline (triage, consensus polishing, library split)
7. **Dual classification**: Two independent classifiers (RepeatClassifier + ClassifyTE) combined for robustness
8. **TE-looker byte-level reproducibility (G10)**: All PRNG via ChaCha20 seed derivation, `BTreeMap` for deterministic iteration, RFC 8785 canonicalized JSON + BLAKE3 provenance hashing. Same input + same seed → identical bytes.
9. **TE-looker scope discipline**: `dtr` aborts below D_cc 45% (Eddy 2008 random-match floor) and abstains above 75% (defer to Look4LTRs / mdl-repeat). The CLI banner enforces this; do not silence it.
10. **De novo confidence gating (both de novo tracks)**: every de novo repeat finder over-produces (one element fragmented into many sub-families; divergent copies escaping recruitment re-seeded as redundant families). The fix is NOT downstream consensus rebuilding/stitching (which conserves total bp) and NOT discovery-parameter tuning alone — it is a confidence/redundancy gate that drops or collapses the low-confidence bulk while protecting structurally-real TEs. **The false-positive metric is total library bp (vs RepeatModeler2/EDTA ≤15 MB for maize), not family count.** mdl-repeat: gate on its native `tier`/`accept` ∪ TE-protein/structural signal ∪ copy support (Refiner_mdl Phase 2). TE-looker: copy/length gate + cd-hit-0.80 redundancy collapse (`_gate_and_merge_te_looker`). Use TE-protein homology (not EDTA membership, not bp) as the real-TE discriminator (see memory note [[edta-not-reference]]).

---

## Dev vs Production Note

The dev repo (`/home/shuoc/tool/tmp/Pan_TE/`) contains the TE-looker migration (Step 4: RepGraph → `dtr`). Production at `/home/shuoc/tool/Pan_TE/` may still expect the legacy RepGraph paths until this work is propagated. PATH points to production; copy `bin/Pan_TE` and remove the same legacy Perl files (`bin/run_repgraph`, `bin/MSPCollect.pl`, `bin/build_consensus`, `bin/Refiner.py`, `lib/RepGraph/`) before running in production. The `dtr` binary itself lives in `~/tool/te-looker/` and is shared between dev and production.
