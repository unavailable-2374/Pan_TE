# mdl-repeat Refinement Strategy — Ground-Up Design

## Executive summary (conclusion first)

The refinement of mdl-repeat family seeds should be redesigned around **one
authoritative fact the current implementation ignores: the
`mdl_repeat.instances.bed` file already gives the exact genomic coordinates and
strand of every copy of every family.** The current Phase 1 throws this away and
re-discovers copies by BLASTing each seed against a whole-genome BLAST database —
an operation whose cost scales as roughly *O(number of families × genome size)*,
which becomes intractable at 10–30 Gb genomes with 10⁵–10⁶ families, and which is
also algorithmically *unable to extend* a consensus past its own length.

The new strategy is **BED-seeded extend–align–trim**:

1. **Copy source = BED instances, not genome BLAST.** Read every copy's
   coordinates directly. Extraction cost is linear in the total bp of repeat
   instances (random-access `samtools faidx`), with no genome-scale alignment and
   no whole-genome `makeblastdb`. BLAST is retained only as a *bounded, optional
   per-family fallback* for families with too few BED instances to polish.
2. **Boundary refinement = extend–align–trim.** Extract each instance with a
   flank pad, co-orient copies using the BED strand, build a per-family MSA, then
   trim each terminus inward to the first column meeting an occupancy floor. Pad
   columns that are shared across copies survive the trim and **extend** the
   consensus beyond mdl-repeat's ragged boundary; pad columns that are
   per-copy-unique genomic flank collapse in occupancy and are **trimmed**. This
   gives both behaviors (extend *and* trim) that the HSP-bounded BLAST path could
   never give.
3. **Chimera detection = MSA occupancy valley + copy-span bimodality + BED
   independence**, replacing the BLAST `qstart/qend` coverage heuristic.
4. **Family-parallel, sharded, resumable execution** — the family loop is
   embarrassingly parallel and never holds the genome or all alignments in RAM.

Expected gains: **the dominant near-quadratic genome×families coupling is removed
entirely** (replaced by linear instance extraction + capped, parallel per-family
MSA), and **boundary quality gains a true extension capability plus principled
occupancy trimming and strand-correct co-orientation.** Concrete wall-clock
numbers must be benchmarked on a real genome — they are deliberately not asserted
here.

Migration: Phase 0 is mostly reusable (one change to dedup aggressiveness, plus
emitting a per-family instance index). Phase 1 is a full rewrite. Phase 2 is
reusable as-is. `main.py` drops the unconditional whole-genome `makeblastdb` from
the critical path and replaces coarse pickle checkpoints with shard markers.

---

## 1. Problem statement and what mdl-repeat already guarantees

mdl-repeat is not RepeatScout-with-extra-output; it already runs an internal
refinement pipeline (merge under an 80-80-80 rule, Otsu subfamily split,
spatial-co-occurrence fragment assembly, exclusive-coverage prune) and only emits
a family when `mdl_score > 0`. Its three outputs are:

| Output | Content | Refine-relevant fact |
|---|---|---|
| consensus FASTA | `>R=N length=L copies=C mdl=M [div=D] [topo=…]` | one weighted-DP consensus per family, already built from *all* copies |
| `instances.bed` | `chr  start  end  R=N  score  strand` | **exact coordinates + strand of every copy**; `score = int(1000·(1−divergence))` ⇒ per-instance divergence is recoverable as `1 − score/1000` |
| `stats.tsv` | per-family summary | family-level QC metadata |

Two consequences drive the whole design:

- **The copy set is already known.** Re-deriving it by BLAST is redundant work
  that duplicates information mdl-repeat handed us for free, *and* it is the
  single most expensive step in the pipeline.
- **What mdl-repeat does *not* do well is what refinement must add:** (a) MSA-grade
  consensus polishing (its consensus is column-wise weighted-DP voting, not a true
  progressive/POA MSA); (b) **boundary extension/trim** — mdl-repeat's seed-extend
  boundaries are ragged and frequently 5′-truncated relative to the true
  full-length element; (c) coverage-valley chimera splitting (its Otsu split only
  handles bimodal *divergence*, not coverage discontinuity); (d) deeper redundancy
  collapse after boundaries converge.

So refinement is precisely a **consensus-polish + boundary + chimera + redundancy**
layer, exactly analogous to the RepeatModeler post-RECON consensus stage — but
seeded from known coordinates instead of a fresh search.

### 1.1 The two confirmed defects in the current implementation

**Defect A — performance.** `phase1_consensus.run_blastn_batched` BLASTs every
record (batched, but all of them) against `config.genome_blast_db`, the whole
genome. Both factors scale with genome size: the number of families grows and the
DB to scan grows. `main._setup_genome_blast_db` builds that DB with a hard
`timeout=600`, which a 10–30 Gb genome will silently blow past, leaving
`genome_blast_db=""` and *all* recruitment silently skipped (every family falls
back to its raw mdl-repeat seed). This is both slow and a silent-failure trap.

**Defect B — boundaries.** Even when BLAST succeeds, `_extract_hits` extracts the
*subject* (genomic) span of each HSP, but consensus is then built only from
`long_hits` requiring `h['length'] >= slen*0.8` — i.e. alignments that cover ≥80%
of the *query consensus*. The consensus can therefore never grow past the query
length: the boundary is **clamped to the existing (ragged, possibly truncated)
mdl-repeat seed**. Trimming is implicit and extension is impossible. Conversely,
if one were to "just extract the BED spans," the consensus would inherit
mdl-repeat's ragged ends with no correction at all. Neither extreme is right; the
correct answer is extend-then-trim.

---

## 2. Copy-source strategy

**Decision: BED-first, with a bounded optional BLAST fallback. Never whole-genome
BLAST of the full family set.**

### 2.1 Primary path — BED instance extraction

For each family `R=N`, take its instance list `[(chrom, start, end, strand,
divergence), …]` parsed once in Phase 0 and threaded through to Phase 1 (no
re-parsing). Extract each instance's sequence by `samtools faidx` random access.

- **Complexity.** Parsing BED: one pass, *O(total instances)*. Extraction:
  *O(total bp of instances)* = *O(genome × repeat fraction)*, linear, with random
  access amortized by the `.fai` index. There is **no genome-scale alignment and
  no whole-genome DB build** on this path.
- **Scalability.** This is the property that makes 10–30 Gb tractable: cost grows
  with *repeat content*, not with *families × genome*.
- **Boundary/divergence quality.** BED instances are mdl-repeat's accepted copies
  (instances within `max-divergence = 0.30`, i.e. ≥70% identity). They are
  exactly the copies that legitimately belong to the family. Their per-instance
  divergence (from the score column) lets us (a) stratify subsampling and (b)
  optionally down-weight the most divergent copies in consensus.
- **Strand.** The BED carries strand. We **reverse-complement minus-strand
  instances before MSA** so all copies are co-oriented. This is strictly better
  than the old path, which normalized `sstart/send` and *lost* strand, then leaned
  on `MAFFT --adjustdirection` to guess orientation — unreliable for short or
  divergent copies.

### 2.2 When BLAST is still warranted (bounded fallback only)

Two legitimate cases, both handled *per-family* and *capped*, never as an
all-records genome sweep:

1. **BED missing or family has too few instances to polish** (`< min_copies_for_msa`,
   e.g. <5). A family can be under-instanced if the BED was not emitted (older
   mdl-repeat) or if scaffolding/dedup in Phase 0 left a record without a clean
   coordinate set. For *only these families*, run a targeted `blastn` of that one
   consensus against the genome DB with a capped `max_target_seqs`, to recruit
   enough copies for an MSA. Cost is *O(few families × genome)*, not
   *O(all families × genome)* — typically a small minority of families.
2. **Deliberate divergent recruitment (off by default).** If one wants copies
   *below* mdl-repeat's 70% identity floor, only BLAST can find them. This is
   intentionally **disabled by default**: copies in the 45–75% twilight zone are
   the remit of TE-looker (Step 4), not mdl-repeat refine. Recovering them here
   would duplicate TE-looker and risk over-extending families with noisy copies.
   Expose it as an opt-in flag for users who run Pan_TE without the TE-looker step.

The whole-genome `makeblastdb` therefore moves **off the critical path**: build it
lazily only if the fallback is actually triggered (and once, reused), and remove
the 600 s timeout trap (build with no timeout, or stream-validate, and fail loudly
if it errors rather than silently emptying the DB path).

---

## 3. Boundary refinement — extend–align–trim (core design point)

**Decision: extract with a flank pad, co-orient, MSA, then occupancy-trim each
terminus inward. Extension and trimming are two outcomes of the *same* mechanism,
not separate steps.**

### 3.1 Mechanism

For a family with co-oriented padded copies:

1. **Pad extraction.** For instance `(chrom, b_start, b_end, strand)` (BED is
   0-based half-open), extract `samtools` region
   `chrom:(b_start+1−pad) … (b_end+pad)`, clamped to `[1, chrom_len]`
   (chrom_len from the `.fai`). Convert BED↔samtools coordinates explicitly — this
   is a classic off-by-one/coordinate-system bug surface and must be unit-tested.
   Reverse-complement if `strand == '−'`.
2. **MSA** of the padded, co-oriented copies (anchored conceptually on the BED
   core; see §4 for tool/tiering).
3. **Column occupancy profile.** For each alignment column, occupancy =
   (#copies with a non-gap, non-N base) / (#copies). The padded flanks of *true*
   element termini are shared by many copies → high occupancy. Genuine unique
   genomic flank differs per copy → occupancy collapses within a few bp of the
   real boundary.
4. **Two-pointer terminal trim.** Walk inward from the left end to the first
   column with `occupancy ≥ T_occ` *and* a supported majority base; symmetric from
   the right. Keep everything between the two pointers. **Crucially, only the
   termini are trimmed; interior low-occupancy columns are *kept*** (emitted as the
   majority base or `N`), so an internal indel-rich or AT-rich dip never truncates
   the element. This is the explicit guard against over-trimming internal regions.

Extension and trim fall out automatically:

- If the consensus core's true boundary lies **inside the original BED span**
  (mdl-repeat over-extended into flank), occupancy collapses *before* the BED edge
  and the pointer stops short → **trim**.
- If the true boundary lies **outside the BED span** (mdl-repeat truncated the
  element), the pad columns just beyond the BED edge are shared across copies, stay
  above `T_occ`, and are retained → **extension** beyond mdl-repeat's call.

### 3.2 Parameter choices and rationale

| Parameter | Recommended | Rationale / tradeoff |
|---|---|---|
| `pad` | `clamp(0.2 · L, 50, 500)` bp | Must be large enough to expose a truncated terminus but small enough not to pull in the *next* (often nested) TE insertion. Element-length-proportional with a hard cap; LTR-scale termini (kb) are deliberately out of scope — Look4LTRs owns LTRs. |
| `T_occ` (trim occupancy floor) | `max(3 copies, 0.3 · N_copies)` present | Absolute floor of 3 prevents a single chance-similar flank from "extending"; the fractional term scales with family size. Mirrors CIAlign/trimAl occupancy trimming. |
| consensus min-occupancy (interior `N` call) | ~0.5 | Below this an interior column is emitted as `N` rather than a base; keeps the element contiguous without asserting a base that a minority supports. |
| extension cap | ≤ `pad` per side | Extension is bounded by available pad material; document that elements whose termini exceed `pad` cannot be fully recovered. |

### 3.3 Quality protection (never regress below mdl-repeat)

If a family has `< min_copies_for_msa` co-oriented copies, or the MSA fails the
occupancy QC (too large a fraction of low-occupancy columns), **keep the original
mdl-repeat consensus unchanged.** mdl-repeat already built that consensus from all
copies via weighted DP; refinement must only *replace* it when it can demonstrably
*improve* it. This is the single most important correctness invariant and is
already honored by the current code's `consensus_source='original'` fallback —
preserve it.

---

## 4. Consensus construction

**Decision: MAFFT for small/medium families, abPOA for large families, with
stratified subsampling above a copy cap; majority-rule occupancy-aware base
calling; original-seed fallback on failure.**

### 4.1 MSA engine selection by family size and divergence

| Family regime | Engine / mode | Rationale |
|---|---|---|
| small, low-div (`≤ ~30` copies, mean div `≤ ~0.1`) | MAFFT **L-INS-i** (`--localpair`) | most accurate; affordable at this size |
| medium (`~30–200` copies) | MAFFT `--auto` / FFT-NS-2 | progressive is adequate; L-INS-i too slow |
| large (`> cap`, e.g. `> 100–200` copies) | **abPOA** (already a Pan_TE dependency at `~/tool/abPOA/bin/abpoa`) on a stratified subsample | partial-order alignment scales near-linearly in sequence count where progressive MSA degrades; consensus is POA-native |

- **Stratified subsampling above the cap.** Families can have 10⁴–10⁶ copies
  (abundant SINE/LINE). Aligning all of them is infeasible and unnecessary.
  Subsample to a cap (e.g. 50–100 for MAFFT, a few hundred for abPOA) **stratified
  by per-instance divergence (BED score) and length**, so the consensus reflects
  the full divergence spectrum and is not biased toward the dominant young
  subfamily. Uniform-random subsampling would over-represent the most numerous
  (usually youngest) copies.
- **Determinism.** Subsampling must use a fixed seed so reruns are reproducible
  (consistent with TE-looker's G10 reproducibility posture elsewhere in Pan_TE).

### 4.2 Consensus calling

Majority-rule with occupancy gating (the existing `build_majority_consensus` is a
reasonable base): emit the majority base when its support among non-gap bases meets
the threshold, else `N`; skip columns that are essentially all-gap. Optionally
weight votes by `(1 − instance_divergence)` so cleaner copies dominate — a modest,
defensible refinement over unweighted majority. IUPAC ambiguity codes are
*not* recommended for the masking/analysis libraries (downstream
RepeatClassifier/RepeatMasker and the TE-looker HMM build prefer ACGTN).

---

## 5. Chimera detection (BLAST-free)

**Decision: split only when *two independent signals* agree — an MSA occupancy
valley and bimodal copy-span endpoints — optionally corroborated by BED
independence. Recursion depth ≤2, hard minimum fragment length.**

The old detector binned the consensus and counted BLAST `qstart/qend` overlaps per
bin (CV threshold). With copies now coming from the MSA we have strictly richer
information:

1. **Occupancy valley.** From the same column-occupancy profile used for trimming,
   find the longest interior run of columns with `occupancy < ratio · median_occupancy`.
   A true composite/chimera shows a sharp interior occupancy trough where the two
   halves' copy sets stop overlapping.
2. **Copy-span bimodality.** For each copy record its `[first_col, last_col]` span
   in the MSA. A genuine single family: most spans cover the whole consensus. A
   chimera: spans cluster into two groups (left-supporting vs right-supporting)
   whose endpoints pile up at the valley. Require this bimodality to **coincide**
   with the occupancy valley before splitting — this is what prevents splitting on
   a merely AT-rich / indel-rich interior dip (which produces a valley but *not*
   two disjoint span groups).
3. **BED independence (strong corroboration, reuses Phase 0 machinery).** If the
   left segment occurs as an *independent* genomic insertion (the left half present
   at loci where the right half is absent), the split is biologically real. This is
   the inverse of Phase 0's spatial-co-occurrence fragment assembly and can share
   its BED-scanning code.

Split at the valley midpoint; recurse on each fragment up to
`chimera_max_depth = 2`; discard fragments below `min_chimera_fragment`. Splitting
is conservative by construction (two-signal AND-gate), because a false split
permanently fractures a real family.

---

## 6. Redundancy and family merging

**Decision: light, exactness-oriented dedup *before* refine; a moderate
identity-based merge *after* refine; leave cross-source and twilight-zone merging
to Combine (CD-HIT-EST 90%) and TE-looker (Stage 4.5).**

- **Before refine (Phase 0).** Only collapse near-identical duplicate seeds (e.g.
  CD-HIT-EST at ~98% identity, high coverage) to avoid polishing the same family
  twice. **Recommendation: do *not* run the current 90%-identity second round
  pre-refine** — merging at 90% *before* boundaries are corrected can fuse distinct
  subfamilies that extend–align–trim would have kept separate (over-merge destroys
  families, which is worse than transient redundancy). Defer the looser merge to
  after refinement.
- **After refine.** Boundary refinement makes two formerly-ragged seeds converge to
  the *same* polished consensus, newly exposing redundancy that only now is
  visible. Run a single moderate CD-HIT-EST pass (e.g. 95% identity, `-aS 0.8`) to
  collapse these boundary-convergent duplicates. Keep it at 95%, not 90%, to avoid
  collapsing genuine subfamilies.
- **Don't over-merge.** Combine already deduplicates all three sources together at
  90%, and TE-looker runs its own family-merger with a phylogenetic gate. Refiner_mdl
  should hand off a clean, internally-non-redundant set, not an aggressively merged
  one — aggressive merging here would pre-empt and degrade those downstream stages.

---

## 7. Scalability, parallelism, and checkpointing

**Decision: family-level parallelism over a sharded work plan; never materialize
the genome or all alignments in memory; per-shard checkpoints for resume.**

- **Parallel model.** The family loop is embarrassingly parallel. Partition
  families into shards; a `ProcessPoolExecutor` worker owns a shard and, per family,
  does: select/subsample instances → batch `samtools faidx` extract (padded) →
  reverse-complement minus strand → MSA (single-threaded per family; parallelism is
  across families, not within) → consensus → extend-trim → chimera. One MAFFT/abPOA
  thread per family avoids the thread-explosion the shared `run_mafft` already
  guards against.
- **Memory.** Never load the genome into Python — `samtools faidx` random access
  keeps resident memory at *O(extracted bp for the family in flight)*. MSA inputs
  are bounded by the subsample cap. The pipeline never holds all 10⁵–10⁶ records or
  their alignments simultaneously (the current pickle-of-everything checkpoint does,
  and must be replaced).
- **I/O efficiency of extraction.**
  - **Sort each worker's regions by `(chrom, start)`** so faidx seeks are
    monotonic and friendly to the OS page cache.
  - **Prefer `samtools faidx -r region_file`** (one call per shard) over
    hundreds of positional-argument calls. *Caveat (must resolve, not work around):*
    the system samtools (1.4.1, per project memory) does **not** support `-r`; the
    current code uses positional args in batches of 200 as a compatibility
    workaround. Per the project's degradation discipline this is a *resolvable*
    blocker — **pin a samtools ≥1.9** in the PGTA env so `-r` works, collapsing
    per-shard extraction to a single call. Keep the positional-batch path only as a
    detected fallback.
  - Optionally extract per chromosome to turn random access into a streamed pass for
    extremely instance-dense genomes.
- **Big-family cap.** Enforce the subsample cap (§4.1) so a handful of ultra-abundant
  families cannot dominate runtime or blow up an MSA.
- **Checkpoint granularity.** Replace the three coarse `pickle` blobs
  (`phase0/1/2_complete.pkl`, which serialize *all* records and are both
  memory-heavy and all-or-nothing) with **per-shard output files + a done-marker
  per shard**. Phase 1 then resumes by skipping shards whose markers exist and
  concatenating shard FASTAs at the end. This makes a multi-day 30 Gb run safely
  resumable at fine granularity.

---

## 8. Correctness and reviewer acceptability

### 8.1 Relationship to established paradigms

| Established method | What it does | This design's relationship |
|---|---|---|
| **RepeatScout** | seed-and-extend family seeds, no internal refine | mdl-repeat replaces it (MDL-arbitrated, with internal refine); we refine mdl-repeat's *output* |
| **RECON** | clustering of genomic copies into families | mdl-repeat already does family formation; we do the post-family consensus layer |
| **RepeatModeler2 `Refiner` / `alignAndCallConsensus`** | recruit copies (RMBlast) → MSA → **extend/trim against flanks** → iterate | our extend–align–trim mirrors this boundary logic, but **seeds copies from known BED coordinates instead of re-searching**, which is the performance win; the boundary mathematics are the standard, accepted ones |
| **CIAlign / trimAl** | occupancy-based alignment column trimming | our terminal occupancy trim is the same well-accepted principle, applied only at termini to preserve interior structure |
| **abPOA (POA)** | scalable partial-order MSA + consensus | adopted for large families where progressive MSA degrades |

The design uses no novel, hard-to-defend heuristic: it is the standard
recruit→align→consensus→boundary→dedup pipeline, with the single principled change
that recruitment reads coordinates mdl-repeat already produced.

### 8.2 Stated limitations (be explicit in any write-up)

1. **Divergence ceiling inheritance.** BED-seeded recruitment recovers only
   mdl-repeat's accepted copies (≥~70% identity). Copies in the 45–75% twilight
   zone are *deliberately* left to TE-looker; document this division of labor so a
   reviewer doesn't read it as a sensitivity gap. (The opt-in BLAST recruitment of
   §2.2 exists for Pan_TE runs without TE-looker.)
2. **Extension bounded by pad and by genomic reality.** If *no* instance of a
   family is full-length in the genome (e.g. all copies 5′-truncated), no MSA can
   reconstruct the missing terminus — a fundamental limit, not a design flaw.
   Extension is also capped at `pad` per side.
3. **Coordinate-system discipline.** BED (0-based, half-open) ↔ samtools (1-based,
   inclusive) conversion and strand-aware reverse-complement are correctness-critical
   and must be unit-tested; an off-by-one here silently shifts every boundary.
4. **Conservative chimera splitting** can leave some true composites unsplit (the
   two-signal AND-gate favors precision over recall) — the right tradeoff, since a
   false split permanently fractures a real family and is harder to recover
   downstream than an unsplit composite.

---

## 9. Migration impact relative to the current three-phase implementation

| Component | Disposition | Detail |
|---|---|---|
| **Phase 0** (`phase0_triage.py`) | **mostly reusable** | Keep parsing, hard filter, tiering, BED fragment assembly (`assemble_fragments_bed`), cyclic-topology drop. **Change dedup:** pre-refine round only at ~98% (exact/near-exact); drop the 90% second round to after refine. **Add:** build and return a per-family instance index `{R=N → [(chrom,start,end,strand,div), …]}` parsed once, so Phase 1 never re-reads the BED; ensure scaffolds (from chain assembly) carry their members' coordinate sets through to Phase 1. |
| **Phase 1** (`phase1_consensus.py`) | **full rewrite** | Remove `run_blastn_batched` + whole-genome recruitment from the common path. New worker: BED instances → subsample → padded `samtools faidx` extract → strand RC → MSA (MAFFT/abPOA by size) → occupancy-aware consensus → **extend–align–trim** → occupancy-valley + span-bimodality chimera. Keep the `consensus_source='original'` fallback invariant. Retain a *bounded per-family* `blastn` only for under-instanced families (§2.2). |
| **Phase 2** (`phase2_library_split.py`, `te_structure_filter.py`) | **reusable as-is** | QC filters, TE-structure verification (self-BLAST of the *consensus set* only — small, not genome-scale), and masking/analysis library split all operate on the refined consensus set and are unaffected. The `consensus_masking.fa` / `phase3_analysis_library.fa` contract to `build_mdl` and downstream is preserved. |
| **`main.py`** | **modify** | Move whole-genome `makeblastdb` off the critical path (lazy, only if §2.2 fallback fires; build once; **remove the 600 s timeout silent-failure trap** — fail loudly). Replace the three `pickle` checkpoints with per-shard markers + final concat. |
| **`build_mdl` (Perl)** | **no interface change** | It already passes `--bed`/`--stats` when present and consumes `consensus_masking.fa`. The BED it forwards (`mdl_repeat.instances.bed`) becomes the *primary* input rather than an optional aid — no new CLI surface needed. |
| **`config.py`** | **add / relax params** | *Add:* `pad_fraction`, `pad_min`, `pad_cap`, `trim_occupancy_floor`, `min_copies_for_msa`, `msa_subsample_cap`, `large_family_threshold` (MAFFT↔abPOA), `abpoa_exe`, `chimera_occupancy_ratio`, `post_refine_merge_identity`, `shard_size`, `samtools_supports_region_file` (auto-detected), `enable_divergent_blast_recruitment` (default false), `subsample_seed`. *Relax/deprecate:* `blastn_max_targets`, `blastn_batch_size`, `t1/t2/t3_max_hits` (now apply only to the bounded fallback). *Keep:* tiering, dedup identities (repurposed per §6), Phase 2 thresholds. |

### 9.1 Net effect on cost and quality

- **Cost.** The pipeline's dominant term changes from
  *≈ O(N_families × genome) BLAST + a whole-genome `makeblastdb`* to
  *≈ O(total instance bp) extraction + O(N_families × capped MSA), fully
  family-parallel and sharded.* The near-quadratic genome×families coupling — the
  thing that makes the current path intractable at 10–30 Gb — is removed. (Absolute
  speedups are genome-dependent and must be measured, not asserted.)
- **Quality.** Boundaries gain a real **extension** capability (impossible in the
  HSP-clamped BLAST path) plus principled occupancy trimming; co-orientation
  becomes strand-correct via BED strand + reverse-complement; chimera detection
  gains a precision-oriented two-signal gate; and the never-regress-below-mdl-repeat
  invariant guarantees refinement can only help or hold, never harm.
