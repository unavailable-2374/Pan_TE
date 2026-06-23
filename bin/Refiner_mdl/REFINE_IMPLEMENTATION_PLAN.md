# mdl-repeat Refinement — Implementation Plan (code-level spec)

## Conclusion first

This document turns `REFINE_STRATEGY_DESIGN.md` into an implementation-level
specification: a full Phase 1 rewrite from **whole-genome BLAST recruitment** to
**BED-seeded extend–align–trim**, plus the minimal Phase 0 / `main.py` / `config.py`
changes that support it. The Phase 2 contract (`consensus_masking.fa` +
`phase3_analysis_library.fa`) and the `build_mdl` ↔ `main.py` CLI are unchanged.

The rewrite splits the current monolithic `phase1_consensus.py` into six modules
(orchestrator + extract + align + boundary + chimera + fallback). Copies come from
`mdl_repeat.instances.bed` via `samtools faidx` (linear in repeat bp, no
genome-scale alignment); boundaries are corrected by occupancy-trimming a padded
per-family MSA (extension *and* trim from one mechanism); chimeras are split only
when an MSA occupancy valley and copy-span bimodality agree; and the family loop is
sharded across processes with per-shard resume.

**The single most important invariant, preserved verbatim from the current code:
refinement may only replace a family's consensus when it can demonstrably improve
it; on any failure or low-confidence path it keeps the original mdl-repeat seed
(`consensus_source='original'`).**

Environment facts this plan is pinned to (verified on this system, not assumed):
`samtools 1.4.1` (no `-r` region-file flag → positional regions, batched);
`abpoa 1.5.4` at `/home/shuoc/tool/abPOA/bin/abpoa` (MSA via `-r 1`); genome cleaned
to `chr1..chrN` so BED chrom names already match `genome.fa.fai` (no `ref|…|`
prefix on the BED path; the prefix-strip rule is retained only for the BLAST
fallback path). All performance statements are qualitative or marked
"待 benchmark"; no wall-clock or speedup numbers are asserted.

---

## 1. Module / file layout

The current `phase1_consensus.py` (single file, BLAST-centric) is **replaced** by an
orchestrator plus five focused modules. Each is independently testable (milestones,
§12). Nothing else under `bin/Refiner_mdl/` changes structurally except the named
hooks in §10.

```
bin/Refiner_mdl/
  phase1_consensus.py      # run_phase1 orchestrator: shard plan, BED routing, pool, concat, post-merge
  phase1_extract.py        # coordinate math, samtools extraction, reverse-complement, subsampling
  phase1_align.py          # build_msa (MAFFT/abPOA dispatch), Alignment row utilities, occupancy profile
  phase1_boundary.py       # extend-align-trim, occupancy-aware consensus, alignment QC, original-seed fallback
  phase1_chimera.py        # MSA occupancy-valley + span-bimodality chimera detection (recursive)
  phase1_fallback.py       # bounded per-family blastn recruitment for under-instanced families
```

Shared lightweight type alias used across modules (define in `phase1_align.py`,
import elsewhere):

```python
# An MSA as a list of equal-length uppercase strings, one per copy, '-' = gap.
# We deliberately do NOT thread Bio.Align objects through the pipeline: MAFFT and
# abPOA produce different containers, so both are normalized to List[str] at the
# build_msa boundary and every downstream routine (occupancy, trim, consensus,
# chimera) operates on List[str]. This decouples boundary/chimera logic from the
# MSA engine.
Alignment = List[str]                  # rows[i] is copy i, all len == aln_len
AlignedFamily = Tuple[Alignment, List["CopyMeta"]]   # rows + parallel per-copy meta
```

### 1.1 `phase1_extract.py`

```python
def bed_to_samtools_region(chrom: str, b_start: int, b_end: int,
                           pad: int, chrom_len: int) -> Tuple[str, int, int]:
    """Convert a BED interval (+pad) to a 1-based-inclusive samtools region.

    BED is 0-based half-open [b_start, b_end). samtools faidx 'chrom:X-Y' is
    1-based inclusive. Core (no pad) -> chrom:(b_start+1)-(b_end). With pad and
    clamping to [1, chrom_len]:

        s1 = max(1, b_start + 1 - pad)
        e1 = min(chrom_len, b_end + pad)

    Returns (region_str, s1, e1). region_str == f"{chrom}:{s1}-{e1}".
    Raises ValueError if e1 < s1 (only possible with corrupt coords; caller drops).
    Pure function — unit-tested for off-by-one in both directions (§11, §12.M1).
    """

def strip_seqid_prefix(sid: str) -> str:
    """Strip an NCBI-style 'ref|NC_x.y|' wrapper -> 'NC_x.y'.

    No-op for clean chr names (BED path). Retained for the blastn fallback whose
    subject ids inherit makeblastdb -parse_seqids pipes. Identical logic to the
    current _extract_hits prefix handling (MEMORY: BLAST DB seqid format mismatch).
        if '|' in sid: return first non-empty token after the first '|', else sid.
    """

def load_fai_lengths(fai_path: str) -> Dict[str, int]:
    """Parse genome.fa.fai -> {chrom: length}. Column 0 = name, column 1 = length.
    Loaded once per worker process (small: one row per chromosome)."""

def reverse_complement(seq: str) -> str:
    """ACGTN-safe RC with full IUPAC complement table; unknown chars -> 'N'."""

def extract_padded_copies(instances: List["Instance"], genome_file: str,
                          fai_lengths: Dict[str, int], pad_fn,
                          config) -> List["Copy"]:
    """Batch-extract co-oriented, padded genomic copies for one family.

    instances : list of Instance(chrom, start, end, strand, divergence)
    pad_fn    : callable(family_len) -> int pad (see compute_pad)

    Algorithm:
      1. For each instance compute region via bed_to_samtools_region; skip + count
         any that raise (corrupt/out-of-range coords) — never fabricate a copy.
      2. Build region_to_meta : Dict[str, deque[Instance]] (a region string may
         repeat if two instances share coords; use a deque and pop in order).
      3. Sort unique regions by (chrom, start) so faidx seeks are monotonic
         (page-cache friendly; §7).
      4. Ensure genome_file + '.fai' exists (samtools faidx genome_file once if
         missing); fail loudly if indexing errors (do not silently continue).
      5. Batch positional samtools calls (config.extract_batch_size, default 200):
             cmd = [samtools, 'faidx', genome_file, *batch_regions]
         NOTE: samtools 1.4.1 has NO -r region-file flag (verified). Positional
         args are mandatory here. If config.samtools_supports_region_file is later
         auto-detected True (samtools >= 1.9 pinned in PGTA), collapse to a single
         '-r region_file' call per family — keep positional batching as the
         detected fallback (§7).
      6. Parse FASTA output; for each '>chrom:s-e' header pop one meta from
         region_to_meta[key]; uppercase the sequence; reverse_complement if
         meta.strand == '-'.
      7. Return List[Copy(id, sequence, divergence, strand)].

    Resident memory is O(sum of extracted bp for this family) only.
    """

def compute_pad(family_len: int, config) -> int:
    """pad = int(clamp(config.pad_fraction * family_len, config.pad_min,
                       config.pad_cap)).  Default 0.2*L clamped to [50, 500] bp.
    Extension is intrinsically capped at `pad` per side (only pad bp of flank is
    ever extracted), so no separate extension cap is needed."""

def stratified_subsample(instances: List["Instance"], cap: int,
                         seed: int) -> List["Instance"]:
    """Down-sample to <= cap, stratified by divergence then length.

    if len(instances) <= cap: return list(instances)
    rng = random.Random(seed)                      # deterministic (G10 posture)
    sort instances by divergence ascending
    n_strata = min(10, max(1, len(instances)//cap + 1))
    split into n_strata contiguous equal-count divergence bins
    within each bin: sort by length descending, then rng.shuffle ties
    round-robin pull one per bin until cap reached
    return selection
    Rationale: uniform random over-represents the youngest/most-numerous subfamily;
    stratifying by divergence keeps the consensus representative of the full age
    spectrum (strategy §4.1)."""
```

`Instance`, `Copy`, `CopyMeta` are plain dataclasses (or `SimpleNamespace`) — see §2.

### 1.2 `phase1_align.py`

```python
def select_regime(n_copies: int, mean_div: float, config) -> str:
    """Return 'mafft_linsi' | 'mafft_auto' | 'abpoa'.
       n <= small_family_threshold and mean_div <= small_family_max_div -> 'mafft_linsi'
       n <= large_family_threshold                                       -> 'mafft_auto'
       else                                                              -> 'abpoa'
    Defaults: small=30, small_max_div=0.1, large=200 (§9)."""

def build_msa(copies: List["Copy"], regime: str, config) -> Optional[AlignedFamily]:
    """Align co-oriented copies; return (rows, metas) or None on failure.

    Copies are ALREADY strand-co-oriented (§1.1 step 6), so MAFFT is called with
    adjustdirection=False — this is the key change from the current code, which
    relied on MAFFT --adjustdirection to guess orientation and could mis-flip short
    or divergent copies.

    mafft_*:
        subsample to config.msa_subsample_cap_mafft (default 100) if needed
        aln = run_mafft(seqs, algorithm='localpair'|'auto', adjustdirection=False,
                        thread=1, config=mcfg)        # one thread per family (§7)
        if aln is None: return None
        rows, metas = _bio_alignment_to_rows(aln, copies)   # reorder to input order by id
    abpoa:
        subsample to config.msa_subsample_cap_abpoa (default 300) if needed
        write copies to tmp.fa; run:
            [config.abpoa_exe, '-r', '1', '-o', out.msa, tmp.fa]
        parse out.msa (FASTA-with-gaps; abPOA labels it PIR but emits '>name' +
        gapped rows — use a tolerant FASTA reader, NOT Bio PIR parser); order is
        input order. On non-zero exit or empty output: return None.
    All rows uppercased; gap char normalized to '-'. Validate all rows equal length
    (assert; on mismatch log + return None — never silently truncate).
    """

def occupancy_profile(rows: Alignment) -> "np.ndarray":
    """occ[col] = (# rows with base in {A,C,G,T}) / n_rows, for each column.
    N and '-' and '.' count as absent. Returned as float ndarray length aln_len."""

def column_majority(rows: Alignment, col: int) -> Tuple[Optional[str], int, int]:
    """Return (majority_base, majority_count, n_present) over ACGT at a column."""
```

`_bio_alignment_to_rows` converts a Bio `MultipleSeqAlignment` (from `run_mafft`) to
`List[str]`, reordering rows to match the input `copies` order by sequence id (MAFFT
may reorder; abPOA preserves order). The parallel `metas` list carries
`(id, divergence, strand)` per row for chimera span analysis.

### 1.3 `phase1_boundary.py`

```python
def trim_termini(rows: Alignment, config) -> Tuple[int, int]:
    """Two-pointer occupancy trim. Returns (left, right) inclusive column bounds.

    occ = occupancy_profile(rows); n = len(rows)
    floor = max(config.trim_min_copies, config.trim_occupancy_floor * n)  # present-count
    def solid(col):
        base, cnt, npres = column_majority(rows, col)
        return npres >= floor and base is not None and cnt / max(npres,1) >= 0.5
    left = 0
    while left < aln_len and not solid(left): left += 1
    right = aln_len - 1
    while right >= left and not solid(right): right -= 1
    return left, right
    # ONLY termini move inward. Interior low-occupancy columns between left..right
    # are never trimmed (guard against truncating AT-rich / indel-rich interiors,
    # strategy §3.1). Extension beyond the original BED edge happens automatically
    # when shared pad columns stay solid; per-copy-unique flank collapses below the
    # floor and is trimmed."""

def build_consensus(rows: Alignment, left: int, right: int,
                    metas: List["CopyMeta"], config) -> str:
    """Occupancy-aware majority consensus over columns [left, right].

    Reuses the proven logic of the current build_majority_consensus, restricted to
    the kept span and operating on List[str]:
      for col in [left..right]:
        counts over ACGT (ignore -, N, .)
        total = sum(counts)
        if total == 0: continue                       # pure insertion column -> drop
        if total < max(2, 0.1*n): continue            # near-empty insertion -> drop
        base, cnt = most common
        emit base if cnt/total >= config.consensus_min_occupancy else 'N'
    Optional (config.consensus_divergence_weighted, default False): weight each
    vote by (1 - copy.divergence) so cleaner copies dominate. IUPAC codes are NOT
    emitted (downstream RepeatClassifier / HMM build prefer ACGTN; strategy §4.2)."""

def alignment_qc_ok(rows: Alignment, config) -> bool:
    """Identical intent to current check_alignment_quality: fraction of columns with
    occupancy < config.consensus_min_occupancy must be <=
    config.consensus_occupancy_fail_ratio. Operates on List[str]."""

def refine_family(rec: Dict, copies: List["Copy"], config) -> Dict:
    """Produce the refined record for ONE family. Never regresses below mdl-repeat.

    original_seq = rec['sequence']
    out = dict(rec); out['sequence'] = original_seq; out['consensus_source']='original'
    if len(copies) < config.min_copies_for_msa: return out          # fallback
    mean_div = mean(c.divergence for c in copies)
    regime = select_regime(len(copies), mean_div, config)
    fam = build_msa(copies, regime, config)
    if fam is None: return out                                       # MSA failed
    rows, metas = fam
    if not alignment_qc_ok(rows, config): return out                # poor alignment
    left, right = trim_termini(rows, config)
    if right - left + 1 < config.min_chimera_fragment: return out    # collapsed
    cons = build_consensus(rows, left, right, metas, config)
    if len(cons) < config.min_chimera_fragment: return out
    out['sequence'] = cons
    out['actual_length'] = len(cons)
    out['consensus_source'] = regime                                 # 'mafft_linsi'|'mafft_auto'|'abpoa'
    out['_aln'] = (rows, metas, left, right)   # transient handle for chimera step; stripped before return
    return out
    """
```

### 1.4 `phase1_chimera.py`

```python
def detect_chimera(rec: Dict, config, depth: int = 0) -> List[Dict]:
    """MSA-based, BLAST-free chimera splitter. AND-gate of two signals.

    Requires rec['_aln'] = (rows, metas, left, right). If absent (original-seed
    fallback families have no alignment) -> return [rec] (cannot split safely).

    if depth >= config.chimera_max_depth: return [strip_aln(rec)]
    rows_full, metas, left, right = rec['_aln']
    rows = [r[left:right+1] for r in rows_full]      # work on the kept span only
    W = right - left + 1; n = len(rows)
    if n < config.chimera_min_hits or W < 2*config.min_chimera_fragment:
        return [strip_aln(rec)]

    occ = occupancy_profile(rows); med = median(occ)
    # --- Signal 1: occupancy valley ---
    valley = occ < config.chimera_occupancy_ratio * med      # boolean per column
    suppress termini: ignore first/last config.chimera_min_valley_cols columns
    (rows, metas, left, right) = rec['_aln']                  # interior-only run
    find longest interior run of valley==True; need len >= config.chimera_min_valley_cols
    if no qualifying run: return [strip_aln(rec)]
    v = midpoint column of the run

    # --- Signal 2: copy-span bimodality coinciding with the valley ---
    for each row compute first_col, last_col of non-gap within span
    margin = config.chimera_min_valley_cols
    left_only  = fraction of rows with last_col  <  v + margin and first_col < v
    right_only = fraction of rows with first_col >  v - margin and last_col  > v
    spanning   = fraction of rows with first_col <= v-margin and last_col >= v+margin
    bimodal = (left_only  >= config.chimera_span_group_min_frac and
               right_only >= config.chimera_span_group_min_frac and
               spanning   <= config.chimera_max_span_frac)
    if not bimodal: return [strip_aln(rec)]      # valley without two copy groups
                                                 # = AT-rich/indel dip, NOT a chimera

    # --- Optional Signal 3 (config.enable_bed_independence_chimera, default off) ---
    # Corroborate that left and right halves occur as independent genomic insertions
    # (inverse of Phase 0 co-occurrence). Raises confidence; not required for split.

    # --- Split: map column v to a consensus base coordinate, recurse on columns ---
    split_base = count of emitted consensus positions in build_consensus over
                 columns [left .. left+v]   (re-walk the same emit rules)
    frag1 = make_fragment(rec, seq[:split_base], rows[:, left:left+v], '_chimfrag1')
    frag2 = make_fragment(rec, seq[split_base:], rows[:, left+v:right+1], '_chimfrag2')
    drop any fragment with len < config.min_chimera_fragment
    return detect_chimera(frag1, config, depth+1) + detect_chimera(frag2, config, depth+1)
    """
```

`make_fragment` rebuilds a child record whose `_aln` is the column sub-slice (so the
recursion re-trims and re-checks on the sub-alignment without re-extraction —
strictly cheaper and consistent). `strip_aln(rec)` deletes the transient `_aln`
handle before a record leaves Phase 1 (it must not reach Phase 2 / be pickled).

### 1.5 `phase1_fallback.py`

```python
def recruit_by_blastn(rec: Dict, config) -> List["Copy"]:
    """Bounded per-family copy recruitment for under-instanced families (§7 strategy).

    Triggered only when a family has < config.min_copies_for_msa BED instances.
    Requires config.genome_blast_db (built lazily by main.ensure_genome_blast_db).
    Runs ONE blastn of rec['sequence'] vs the genome DB:
        -max_target_seqs config.blastn_max_targets, -max_hsps 1,
        -evalue config.blastn_evalue, -dust no, pident filter:
            config.min_recruit_pident (70) by default;
            if config.enable_divergent_blast_recruitment: lower to recruit 45-75%
            copies (off by default — that band is TE-looker's remit, strategy §2.2).
    Extract subject spans WITH pad (strip_seqid_prefix on sseqid; strand from
    sign of sstart vs send -> RC accordingly), returning List[Copy].
    Cost is O(few families x genome), never O(all families x genome).
    """
```

### 1.6 `phase1_consensus.py` (orchestrator)

```python
def run_phase1(phase0_output: Dict, config) -> Dict:
    """Shard plan -> BED route -> parallel per-shard refine -> concat -> post-merge.

    records = phase0_output['records']
    index_path = phase0_output['instance_index']          # TSV from Phase 0 (§2)
    shard_dir = os.path.join(config.temp_dir, 'phase1_shards'); mkdir

    # 1. Shard plan: partition records into shards of config.shard_size, sorted by
    #    len*copies descending so heavy families spread across shards (load balance).
    shards = make_shards(records, config.shard_size)
    write each shard's records to shard_{i}.records.jsonl

    # 2. Route instances: stream index_path ONCE, append each row to the
    #    shard file of its owner record. index is pre-sorted by record_id (Phase 0),
    #    so this is a streaming group-by, O(total instances) IO, no RAM index.
    route_instances(index_path, record_to_shard, shard_dir)

    # 3. Parallel refine, one worker per shard (ProcessPoolExecutor,
    #    max_workers=config.threads). Resume: skip shards whose shard_{i}.done exists.
    with ProcessPoolExecutor(max_workers=config.threads) as ex:
        futs = [ex.submit(_process_shard, i, shard_dir, config_dict)
                for i in range(len(shards)) if not done(i)]
        await all; on worker exception: log, mark shard failed (records fall back to
        original seeds for that shard via _process_shard's own try/except per family).

    # 4. Concatenate shard_{i}.records.jsonl(refined) -> result_records
    # 5. Post-refine merge: single CD-HIT-EST at config.post_refine_merge_identity
    #    (0.95) to collapse boundary-convergent duplicates newly exposed by refine
    #    (strategy §6). Reuses phase0._run_cdhit. This is the looser merge DEFERRED
    #    from Phase 0.
    result_records = merge_post_refine(result_records, config)
    return {'records': result_records, 'stats': {...}}
    """

def _process_shard(shard_idx, shard_dir, config_dict) -> None:
    """Worker: load shard records + routed instances, refine each family, write
    shard_{i}.records.jsonl + shard_{i}.consensus.fa, touch shard_{i}.done.

    cfg = SimpleNamespace(**config_dict); fai = load_fai_lengths(cfg.genome_file+'.fai')
    for rec in shard_records:
        try:
            insts = instances_for(rec)                 # routed rows; scaffolds: §1.7
            if len(insts) < cfg.min_copies_for_msa and cfg.genome_blast_db \
               and cfg.enable_fallback_recruitment:
                copies = recruit_by_blastn(rec, cfg)
            else:
                pad_fn = lambda L: compute_pad(L, cfg)
                copies = extract_padded_copies(insts, cfg.genome_file, fai, pad_fn, cfg)
            refined = refine_family(rec, copies, cfg)
            do_chimera = chimera_eligible(rec, cfg)    # tier/length/topology gate (§1.7)
            out = detect_chimera(refined, cfg) if do_chimera else [strip_aln(refined)]
        except Exception as e:
            log warning; out = [dict(rec)]             # NEVER fake success: keep seed
        write out
    write shard_{i}.done
    """
```

### 1.7 Two carried-over behaviors

- **Chimera eligibility gate** (`chimera_eligible`): preserve the current tier /
  length / topology policy — T1 if `len>=300`, T2 if `len>=150`, T3 only if
  `topology=='complex'`; `topology=='complex'` overrides the length guard on any
  tier. This gate decides *whether* to call `detect_chimera`; the detection
  algorithm itself is the new MSA one.
- **Scaffold instances** (`build_scaffold_instances`): a scaffold record (from Phase
  0 chain assembly) has no single BED family. Reconstruct spanning loci from member
  instances: for each chromosome+strand, where the member families occur adjacent in
  chain order within `config.fragment_gap`, emit one spanning Instance
  `(chrom, first_member_start, last_member_end, strand, mean_member_div)`. If
  `>= min_copies_for_msa` spanning loci are found, refine like a normal family; else
  keep the concatenated seed (`consensus_source='original'`). This reuses the
  adjacency scan already in Phase 0.

---

## 2. Data structures and contracts

### 2.1 Per-family instance index (Phase 0 → Phase 1)

**Disk format** (`<temp_dir>/instances_index.tsv`, written by Phase 0, sorted by
`record_id`):

| col | field | source | notes |
|----|-------|--------|-------|
| 0 | `record_id` | Phase 0 record | owner; scaffold members carry the scaffold's id |
| 1 | `member_R` | BED col 3 `R=N` | original family of the row (== record's R for non-scaffold) |
| 2 | `chrom` | BED col 0 | already `chr1..chrN` (clean genome) |
| 3 | `start` | BED col 1 | 0-based half-open |
| 4 | `end` | BED col 2 | 0-based half-open |
| 5 | `strand` | BED col 5 | `+` / `-` |
| 6 | `divergence` | `1 - int(BED col 4)/1000` | from `score=int(1000·(1−div))` |

Built by a **single streaming pass** over `mdl_repeat.instances.bed`, routed through
`R= → record_id` (with scaffold-member remap), then `sort -k1,1` on disk. No
full in-RAM index — this is what keeps 10–30 Gb tractable. `phase0_output` gains one
key: `instance_index: str` (the TSV path). This is the *only* new Phase 0 output.

**In-worker structure** after routing (`Instance` dataclass):

```python
@dataclass
class Instance:
    chrom: str; start: int; end: int; strand: str; divergence: float
```

### 2.2 Record dict — fields across phases

| field | P0 | P1 | P2 | status | notes |
|-------|----|----|----|--------|-------|
| `id` | ✓ | ✓ | ✓ | keep | `mdl_R<N>` or `scaffold_<k>` |
| `R` | ✓ | ✓ | ✓ | keep | int family number |
| `length` | ✓ | ✓ | ✓ | keep | original mdl-repeat length |
| `copies` | ✓ | ✓ | ✓ | keep | from header; scaffolds: min over members |
| `mdl` | ✓ | ✓ | ✓ | keep | |
| `sequence` | ✓ | ✓→refined | ✓ | keep | P1 replaces with consensus (or original) |
| `actual_length` | ✓ | ✓ | ✓ | keep | recomputed by P1 |
| `divergence` | opt | opt | opt | keep | family-level (header), distinct from per-instance |
| `topology` | opt | opt | – | keep | drives chimera eligibility |
| `entropy`,`dust_score`,`n_frac` | ✓ | ✓ | ✓ | keep | P2 recomputes for QC anyway |
| `mdl_per_copy`,`tier` | ✓ | ✓ | ✓ | keep | |
| `is_scaffold`,`scaffold_members` | opt | opt | – | keep | scaffold handling |
| `consensus_source` | – | **new** | ✓ | **new** | `original`\|`mafft_linsi`\|`mafft_auto`\|`abpoa` |
| `is_chimera_fragment` | – | opt | – | keep | set on split fragments |
| `_aln` | – | transient | – | **internal** | MUST be stripped before record leaves P1 |

**Deprecated / removed from the data flow:** BLAST hit dicts (`sseqid/sstart/send/
qstart/qend/pident/length`) no longer flow between functions on the common path —
they exist only transiently inside `recruit_by_blastn`. The `hits` argument that the
old `_process_single`/`detect_chimera` threaded everywhere is gone.

`Copy` / `CopyMeta`:

```python
@dataclass
class Copy:      id: str; sequence: str; divergence: float; strand: str
@dataclass
class CopyMeta:  id: str; divergence: float; strand: str        # parallel to MSA rows
```

### 2.3 Phase 1 → Phase 2 contract (unchanged)

`run_phase1` returns `{'records': [...], 'stats': {...}}`; every record has
`id, sequence, tier, copies, mdl, length` populated. `run_phase2` is called exactly
as today and produces `consensus_masking.fa` + `phase3_analysis_library.fa`. No
Phase 2 edit is required.

---

## 3. Copy extraction (samtools faidx, this system)

1. **Coordinate conversion** — `bed_to_samtools_region` (§1.1). BED 0-based
   half-open → samtools 1-based inclusive: `s1=max(1,start+1-pad)`,
   `e1=min(chrom_len,end+pad)`. Unit-tested both directions (§12.M1).
2. **No `-r` flag** — samtools 1.4.1 verified to lack region-file support. Regions
   are passed positionally, batched at `extract_batch_size=200`
   (`samtools faidx genome.fa chr1:101-300 chr2:51-180 …`). This is the MEMORY-
   recorded constraint, honored exactly. `samtools_supports_region_file` defaults
   False; if a newer samtools is pinned and auto-detection flips it True, a single
   `-r` call per family replaces the batching (positional path kept as fallback).
3. **Seqid prefix** — BED chrom names are clean (`chr1…`), so `strip_seqid_prefix`
   is a no-op on this path; it is applied defensively and is load-bearing only on
   the blastn fallback (subject ids `ref|…|` from `-parse_seqids`).
4. **Strand** — minus-strand instances are reverse-complemented after extraction so
   all copies are co-oriented to the family's `+` frame **before** MSA. MAFFT is
   then called with `adjustdirection=False`.
5. **Pad** — `compute_pad` = `clamp(0.2·L, 50, 500)` bp. Large enough to expose a
   5′-truncated terminus, capped so it does not pull in the next (often nested) TE.
   Extension is bounded by `pad` because only `pad` bp of flank is extracted.

---

## 4. extend–align–trim (precise algorithm)

- **Pad default & basis**: `0.2·L ∈ [50,500]` bp (§3.5). LTR-scale (kb) termini are
  out of scope by design — Look4LTRs owns LTRs.
- **Trim floor (`solid`)**: a column survives if `n_present ≥ max(3, 0.3·n_copies)`
  **and** its majority base has `count/n_present ≥ 0.5`. The absolute floor of 3
  stops a single chance-similar flank from "extending"; the 0.3 fraction scales with
  family size (mirrors CIAlign/trimAl occupancy trimming).
- **Two-pointer scan**: from each end inward to the first `solid` column; keep
  `[left, right]`. Shared pad columns beyond the BED edge stay solid → **extension**;
  per-copy-unique flank collapses → **trim**. Both behaviors from one mechanism.
- **Interior protection**: columns strictly between `left` and `right` are never
  removed regardless of occupancy; low-occupancy interiors are emitted as the
  majority base or `N` (insertion-only columns dropped). This is the explicit guard
  against truncating AT-rich / indel-rich interiors.
- **Degenerate fallbacks** (→ keep original mdl-repeat seed, `consensus_source=
  'original'`): `< min_copies_for_msa` copies; `build_msa` returns None;
  `alignment_qc_ok` False; trimmed span `< min_chimera_fragment`; consensus
  `< min_chimera_fragment`. This realizes the never-regress invariant.

---

## 5. Consensus construction

- **Engine by regime** (`select_regime`): `mafft_linsi` (`--localpair`) for small
  low-divergence families (`≤30` copies, mean div `≤0.1`); `mafft_auto` for medium
  (`≤200`); `abpoa` for large (`>200`). Thresholds in config.
- **Subsample caps**: MAFFT 100, abPOA 300, both via `stratified_subsample`
  (divergence×length strata, fixed `subsample_seed=42` for determinism).
- **abPOA invocation**: `abpoa -r 1 -o out.msa in.fa` (MSA output mode 1). Output is
  FASTA-with-gaps; parse tolerantly (do not use a strict PIR parser). On non-zero
  exit / empty output → `build_msa` returns None → original-seed fallback.
- **Consensus call**: occupancy-aware majority over the kept span (`build_consensus`,
  §1.3); emit `N` below `consensus_min_occupancy` (0.5 effective via the existing
  rule), drop pure-insertion columns; ACGTN only (no IUPAC). Divergence-weighted
  voting is an opt-in (`consensus_divergence_weighted`, default False).
- **Quality protection**: the `refine_family` fallback chain (§4) guarantees the
  output is never worse than mdl-repeat's own all-copies weighted-DP consensus.

---

## 6. Chimera detection (new algorithm)

Two-signal **AND-gate** on the trimmed MSA (`phase1_chimera.detect_chimera`, §1.4):

1. **Occupancy valley** — longest interior run of columns with
   `occ < chimera_occupancy_ratio · median_occ` (ratio 0.5), run length
   `≥ chimera_min_valley_cols` (3), termini excluded.
2. **Copy-span bimodality at the valley** — with `v` = valley midpoint and
   `margin = chimera_min_valley_cols`: require `left_only ≥ 0.2` **and**
   `right_only ≥ 0.2` **and** `spanning ≤ 0.5`. A valley *without* two disjoint copy
   groups (e.g. an AT-rich interior dip) fails this and is **not** split.
3. **(Optional) BED independence** — `enable_bed_independence_chimera` (default
   False): corroborate that each half occurs as an independent genomic insertion
   (inverse of Phase 0 co-occurrence). Raises confidence; not required.

**Split**: map column `v` to a consensus base coordinate by re-walking the emit
rules up to `v`; recurse on the two **column sub-ranges** (no re-extraction) up to
`chimera_max_depth=2`; discard fragments `< min_chimera_fragment`. Conservative by
construction — a false split permanently fractures a real family, so precision is
favored over recall.

---

## 7. Bounded blastn fallback

- **Trigger**: family has `< min_copies_for_msa` BED instances **and**
  `enable_fallback_recruitment` **and** a genome DB is available. Per-family, never
  an all-records sweep.
- **DB scope**: the lazily-built **whole-genome** DB (built once, reused) — building
  a per-family subject DB would be pointless since we are searching the genome. The
  *cost* bound comes from running blastn for only the small minority of
  under-instanced families, not from the DB size.
- **Params**: `-max_target_seqs blastn_max_targets (50)`, `-max_hsps 1`,
  `-evalue blastn_evalue`, `-dust no`, pident `≥ min_recruit_pident (70)`;
  `enable_divergent_blast_recruitment` (default False) lowers pident into the
  45–75% band only for TE-looker-less runs. Extracted subject spans are padded and
  RC'd by HSP orientation, then fed to the same extend-align-trim.

---

## 8. Scalability, parallelism, checkpointing

- **Parallel model**: `ProcessPoolExecutor(max_workers=config.threads)`, one task
  per **shard** (`shard_size=500` families). Each MSA runs single-threaded
  (`thread=1`); parallelism is across families/shards, avoiding the thread explosion
  the shared `run_mafft` already guards against.
- **Memory**: never load the genome into Python — `samtools faidx` random access
  bounds resident memory to the extracted bp of the family in flight; MSA inputs
  bounded by subsample caps. No pickle-of-all-records is ever held.
- **I/O**: regions sorted by `(chrom,start)` for monotonic seeks; single streaming
  BED pass for routing; optional per-chromosome extraction for instance-dense
  genomes (future, not required for v1).
- **Checkpointing (per-shard, replaces the three coarse pickles)**:
  - Phase 0: one checkpoint — write `phase0.records.jsonl` + `instances_index.tsv`
    + a `phase0.done` marker (records jsonl is ~input size, acceptable; no giant
    pickle).
  - Phase 1: `shard_{i}.records.jsonl` (refined) + `shard_{i}.consensus.fa` +
    `shard_{i}.done` per shard. **Resume = skip shards with a `.done` marker** and
    concatenate the rest. A killed 30 Gb run resumes at shard granularity.
  - Phase 2: single `phase2.done` after writing the two libraries.
- **Temp management**: all shard files under `<temp_dir>/phase1_shards/`; removed by
  the existing `_cleanup_temp` unless `--keep-temp`. Per-family MSA tmp files are
  `tempfile`-created and unlinked in `finally` (as `run_mafft` already does).

---

## 9. `config.py` changes

**Add** (defaults + one-line basis):

| param | default | basis |
|-------|---------|-------|
| `pad_fraction` | `0.2` | pad ∝ element length |
| `pad_min` | `50` | floor exposes short truncated termini |
| `pad_cap` | `500` | ceiling avoids pulling in the next/nested TE |
| `trim_occupancy_floor` | `0.3` | fractional present-count floor (CIAlign-style) |
| `trim_min_copies` | `3` | absolute floor; blocks single-flank "extension" |
| `min_copies_for_msa` | `5` | below this, keep mdl-repeat seed (matches current gate) |
| `msa_subsample_cap_mafft` | `100` | progressive MSA tractable to ~100 seqs |
| `msa_subsample_cap_abpoa` | `300` | POA scales further; bound runtime |
| `small_family_threshold` | `30` | L-INS-i affordable below this |
| `small_family_max_div` | `0.1` | L-INS-i reserved for low-divergence families |
| `large_family_threshold` | `200` | MAFFT→abPOA switch |
| `abpoa_exe` | `"/home/shuoc/tool/abPOA/bin/abpoa"` | verified path (PATH fallback `abpoa`) |
| `chimera_occupancy_ratio` | `0.5` | valley = occ < 0.5·median |
| `chimera_min_valley_cols` | `3` | min interior valley run / span margin |
| `chimera_span_group_min_frac` | `0.2` | min support for each half (bimodality) |
| `chimera_max_span_frac` | `0.5` | max full-span copies allowed when splitting |
| `enable_bed_independence_chimera` | `False` | optional 3rd corroborating signal |
| `consensus_divergence_weighted` | `False` | optional cleaner-copy weighting |
| `post_refine_merge_identity` | `0.95` | collapse boundary-convergent dups post-refine |
| `shard_size` | `500` | families per parallel task |
| `samtools_supports_region_file` | `False` | auto-detected; 1.4.1 → False |
| `enable_fallback_recruitment` | `True` | allow per-family blastn for under-instanced |
| `enable_divergent_blast_recruitment` | `False` | 45–75% band is TE-looker's remit |
| `subsample_seed` | `42` | deterministic subsampling (G10 posture) |
| `extract_batch_size` | `200` | positional faidx batch (arg-length safe) |

**Change**: `dedup_round1_identity` `0.95 → 0.98` (pre-refine = near-exact only);
**stop using** `dedup_round2_identity` in Phase 0 (the 0.90 round is removed from
pre-refine; the looser merge moves to Phase 1's post-refine pass at 0.95). Keep the
field for back-compat but unreferenced in Phase 0.

**Deprecate (kept for back-compat, unused on common path)**: `blastn_batch_size`,
`chimera_cv_threshold`, `chimera_low_coverage_ratio`. `t1_max_hits`/`t2_max_hits`/
`t3_max_hits` and `blastn_max_targets` now apply **only** to the bounded fallback.

**Keep unchanged**: tiering params, Phase 2 thresholds, `min_recruit_pident` (now a
fallback/weighting threshold), external-tool exe names.

---

## 10. Interface impact (per file)

- **`phase0_triage.py`** — *one logic change + one new output*:
  (a) `deduplicate` → single CD-HIT-EST round at `dedup_round1_identity` (0.98);
  drop the 0.90 second round. (b) Add `write_instance_index(records, bed_path,
  out_path)` (single BED pass, R=→record_id with scaffold remap, score→divergence,
  `sort -k1,1`); `run_phase0` returns the new `instance_index` key. Everything else
  (parse, hard filter, tiering, `assemble_fragments_bed`, cyclic drop) unchanged.
- **`phase1_consensus.py` + 5 new modules** — full rewrite per §1.
- **`phase2_library_split.py` / `te_structure_filter.py`** — **no change**. They
  consume `phase1_output['records']` and the contract holds. (`te_structure_filter`
  is consensus-set-scale self-BLAST, not genome-scale — reused as-is.)
- **`main.py`** — (a) Remove the unconditional `makeblastdb` from `__init__`; replace
  `_setup_genome_blast_db` with **`ensure_genome_blast_db()`** called lazily by
  Phase 1 only when a fallback family needs it. (b) **Remove the 600 s timeout
  silent-failure trap**: build with no timeout (or a generous one) and on real
  failure set `genome_blast_db=""` **and log an error** — under-instanced families
  then keep their seeds; the common path is unaffected because it never touches the
  DB. (c) Replace the three `pickle` checkpoints with the per-shard markers of §8
  (Phase 1 resume lives inside `run_phase1`; `main.run` just checks `phase1.done`).
- **`build_mdl` (Perl)** — **no change**. It already passes `--bed`/`--stats` and
  consumes `consensus_masking.fa`; the BED becomes the *primary* input with no new
  CLI surface. The lazy fallback DB lands under `config.temp_dir/genome_blastdb`
  (same path as today), so cleanup is unchanged.

---

## 11. Edge cases and failure handling

Every path below either falls back to the original seed or raises loudly — **never
fabricates a copy, alignment, or success**.

- **BED missing / family absent from BED** → 0 instances → under-instanced path
  (fallback blastn if enabled, else original seed). Logged per family count.
- **Single-copy / `< min_copies_for_msa`** → original seed (no MSA attempted).
- **Cross-scaffold members** → handled by `build_scaffold_instances`; if spanning
  loci `< min_copies_for_msa`, scaffold keeps its concatenated seed.
- **Coordinate out of range / negative / `end ≤ start`** → `bed_to_samtools_region`
  raises `ValueError`; that instance is skipped and counted (`skipped_bad_coords`),
  remaining copies still used.
- **Clamping at chromosome ends** → `s1`/`e1` clamped to `[1, chrom_len]` from the
  `.fai`; reduces effective pad on that side (documented, not an error).
- **`.fai` missing** → build once with `samtools faidx genome_file`; if that fails,
  raise (cannot proceed without random access) — do not silently skip extraction.
- **samtools non-zero exit / empty output for a batch** → log stderr, skip that
  batch's copies, continue; if a family ends with `< min_copies_for_msa` usable
  copies → original seed.
- **MAFFT failure / timeout** (`run_mafft` returns None) → original seed.
- **abPOA non-zero exit / empty / unequal row lengths** → `build_msa` None →
  original seed.
- **CD-HIT-EST (post-refine merge) failure** → log, return the un-merged set (no
  data loss; downstream Combine still dedups at 90%).
- **Worker process crash** → per-family `try/except` keeps seeds; a whole-shard
  crash leaves no `.done`, so resume re-runs only that shard.
- **Genome DB build failure (fallback)** → `genome_blast_db=""`, error logged,
  under-instanced families keep seeds; common path unaffected.
- **Duplicate identical instance coordinates** → `region_to_meta` deque pop handles
  the collision deterministically.

---

## 12. Implementation order (milestones) + validation

Each milestone is independently testable; validation uses the **real** Arabidopsis
mdl-repeat output (MEMORY: ~23 K sequences from a ~120 Mb genome) — no synthetic
numbers, and correctness is judged by *observable invariants*, not fabricated
targets.

- **M1 — `phase1_extract.py` (pure coordinate + extraction).** Unit-test
  `bed_to_samtools_region` against hand-worked BED↔samtools cases incl. clamping and
  `end≤start`; test `reverse_complement` round-trips; test `strip_seqid_prefix`.
  Integration: extract a handful of known instances and verify the extracted
  sequence equals an independent `samtools faidx` by hand (and that minus-strand
  copies are RC'd). *Pass = byte-identical extraction; off-by-one caught here.*
- **M2 — `phase1_align.py` + `phase1_boundary.py`.** On a single medium Arabidopsis
  family: build MSA (MAFFT), compute occupancy, trim, build consensus. *Pass =
  consensus length is finite and within `[trimmed_span]`; on a family with shared
  flank the consensus extends past the mdl-repeat length; on a ragged family it
  trims; original-seed fallback fires when copies `< min_copies_for_msa`.* Compare
  refined vs original lengths as a distribution (report the distribution, do not
  assert a number).
- **M3 — `phase1_chimera.py`.** Feed a deliberately concatenated two-family
  alignment (constructed from two real Arabidopsis families) and confirm a split at
  the junction; feed a single AT-rich family and confirm **no** split (AND-gate
  precision). *Pass = split only on the true composite.*
- **M4 — orchestrator + sharding + Phase 0 index.** Run end-to-end on Arabidopsis
  with small `shard_size` to force multiple shards; kill mid-run and resume.
  *Pass = resumed run skips `.done` shards and the final record count equals a
  single-shot run; no family silently lost (count in ≥ count out only via explicit
  chimera splits / post-merge, both logged).*
- **M5 — fallback + `main.py` lazy DB + abPOA path.** Force a few under-instanced
  families (or a large family `>200` to exercise abPOA). *Pass = fallback blastn runs
  only for under-instanced families; abPOA produces a parseable MSA; lazy DB builds
  once and only when triggered; removing TE-looker does not break the common path.*
- **M6 — full pipeline + Phase 2 contract.** Run `build_mdl`-style invocation to
  Phase 2. *Pass = `consensus_masking.fa` + `phase3_analysis_library.fa` are
  produced with the same header schema as today; downstream Combine/TE-looker
  consumers see no interface change.* Spot-check a well-characterized Arabidopsis
  family (e.g. an ATHILA/Copia element) boundary against its RepBase consensus as a
  qualitative sanity check (report agreement, do not invent a percentage).

**Benchmarking** (separate, after M6): measure wall-clock and peak RSS on
Arabidopsis and one larger genome, comparing BED-seeded vs the old genome-BLAST
path. All such numbers are produced *by running*, never asserted here (待 benchmark).
