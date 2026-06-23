# mdl-repeat Refinement Strategy — Representativeness-First Redesign (v2)

## Executive summary (conclusion first)

**The refine stage is being re-centered on two user-anchored goals that the
current single-consensus pipeline does not serve: (1) representativeness — find
every real family member and, whenever a family is divergent, emit one consensus
per subfamily so the library can recover all members on downstream genome search;
and (2) conservative false-positive reduction — actively trim mdl-repeat's own
noise (pseudo-families, low-complexity, over-fragmentation, chimeras) but only
where support is clearly absent, favoring sensitivity over aggression.** The
governing principle under both is that **copy-completeness is the shared
foundation**: only with each family's full, real member set in hand can the
pipeline (a) see the divergence structure and cluster subfamilies correctly, and
(b) judge a family real-vs-noise by whether it recruits a coherent copy set.

The single most important architectural change relative to the implemented M1–M5
pipeline (the `REFINE_STRATEGY_DESIGN.md` / `REFINE_IMPLEMENTATION_PLAN.md` design
that those modules realize): the implemented pipeline assumes **one consensus per
family** and treats copy recruitment as a *failsafe* (blastn fallback fires only
when BED instances `< min_copies_for_msa`). v2 inverts both assumptions:

1. **Always cluster.** After copies are assembled, the recruited copy set is
   always partitioned by pairwise divergence; every qualifying cluster
   (member-count and homogeneity thresholds met) emits its own subfamily
   consensus via the existing extend-align-trim (M2). A homogeneous family yields
   one consensus (one cluster); a divergent family yields several. This is the
   direct realization of goal 1 and the structural inversion of "single consensus
   unless a structural chimera is detected."
2. **Selective recall, never tolerated incompleteness.** Recall is promoted from
   a failsafe to a first-class, scalable mechanism whose target is always *all
   real members* of every family. "Selective" governs only *where compute is
   spent* (families whose BED copy set is already complete and homogeneous skip
   full-genome search; only families that *signal* missing/divergent members
   trigger bounded targeted recall) — it never means tolerating an incomplete
   copy set. The mechanism is a cheap k-mer/minimizer pre-screen followed by
   per-family-capped targeted alignment, so cost stays well below
   O(families × genome).
3. **Conservative FP reduction** stays sensitivity-weighted: only obviously
   unsupported families are removed, and the recall result feeds the real-vs-noise
   decision (a family that cannot recruit a coherent copy set is the pseudo-family
   signal).

**Two orthogonal split axes are both retained and explicitly coordinated.** The
existing M3 structural-chimera splitter ("two *different* elements concatenated"
— A+B junction) and the new subfamily-divergence splitter ("*one* element
diverged into sub-lineages") answer different biological questions on different
signals; v2 specifies their order and mutual non-interference rather than merging
them.

This document **supersedes** the v1 strategy/implementation pair on three specific
points (§9): the "fallback default-off / recall is optional" posture, the
"one-consensus-per-family" output contract, and the framing of the 45–75% band as
purely TE-looker's remit. It **reuses** M1–M5 as a near-complete substrate
(§7): M1 extraction, M2 extend-align-trim, M3 structural chimera, M4 sharding, and
M5 bounded recall all survive; the net new work is a **selective-recall trigger
layer** and a **subfamily clustering module** inserted between M1 and M2, plus a
representativeness metric harness.

---

## 1. Pipeline restructure (主轴: 找齐拷贝 → 总是聚类拆 subfamily → 每簇共识 → 保守降假阳)

### 1.1 Stage map and relation to the implemented Phase 0/1/2

The implemented code is Phase 0 (triage/tiering/BED index) → Phase 1 (per-family
extract → refine → chimera, sharded) → Phase 2 (QC + library split). v2 keeps
Phase 0 and Phase 2 essentially intact and **restructures the per-family body of
Phase 1** from a linear "extract → single consensus → chimera" into a
copy-completeness-first, cluster-always flow:

```
Phase 0  (unchanged in spirit; one added per-family signal column — §2.4)
  parse → hard filter → tiering → BED fragment assembly → near-exact dedup
  → write_instance_index  (+ emit per-family completeness pre-signals)

Phase 1  (restructured per-family body; same sharded orchestrator M4)
  for each family (shard-parallel):
    A. ASSEMBLE COPIES  (find-all-copies foundation — §2)
         A1. BED instances  (M1 extract, primary, cheap)
         A2. completeness triage: is this family's copy set already
             complete+homogeneous?  (signals §2.3)
         A3. if NOT → selective targeted recall  (§2.5; cheap pre-screen →
             per-family-capped alignment; upgraded M5)
         → full real member copy set for this family
    B. CLUSTER ALWAYS  (subfamily axis — §3)
         B1. pairwise divergence/identity matrix over assembled copies
         B2. hierarchical/graph clustering
         B3. qualifying clusters (min members, intra-cluster homogeneity)
             → 1 cluster (homogeneous) … k clusters (divergent)
    C. PER-CLUSTER CONSENSUS  (M2 extend-align-trim, once per cluster — §3.4)
         each qualifying cluster → one subfamily consensus
         (boundary extend/trim is the SAME M2 mechanism, per cluster)
    D. STRUCTURAL CHIMERA  (M3, per subfamily consensus — §5)
         orthogonal A+B split, applied to each emitted subfamily MSA
    E. CONSERVATIVE FP PRUNE  (§4)
         drop only clearly unsupported families/clusters; recall-informed
  → subfamily-resolved, noise-reduced consensus records

Phase 2  (unchanged contract; subfamily records flow through as normal records)
  QC + masking/analysis library split → consensus_masking.fa + analysis lib
```

The reorder is the design's spine: **A (copy-completeness) precedes B
(clustering) precedes C (consensus)**, because clustering on an incomplete copy
set would miss exactly the divergent distal members that motivate subfamily
splitting in the first place, and a consensus built before clustering is the
"one-consensus-sits-mid-spectrum" failure the user named.

### 1.2 Why the implemented order cannot be patched in place

The implemented `refine_family` (phase1_boundary.py) takes the *whole* copy set,
builds *one* MSA, trims it, and emits *one* consensus; M3 then optionally splits
on an A+B junction. There is no point in that flow where divergence-structured
sub-lineages of a *single* element are separated. Adding clustering as a
post-hoc split of the single consensus is wrong (the single MSA has already
averaged the sub-lineages together and the distal members may have failed
`alignment_qc_ok` or been outvoted). Clustering must happen on the **copies**,
before the MSA — hence a Structural restructure of the Phase 1 body, not a Local
patch. This is the v2 diagnosis: the root cause is the one-consensus mental model
baked into `refine_family`, so the fix is a scoped rewrite of the per-family body,
preserving M1/M2/M3/M4/M5 as called components.

---

## 2. Find-all-copies: selective, scalable recall (goal-1 foundation)

**Decision: BED-first copy assembly (cheap, complete for the young core), plus a
signal-triggered, genome-scalable targeted recall whose objective is the full
real member set — selective in *where* it spends compute, never in *whether* it
aims for completeness.**

### 2.1 What the real data says the recall must find

The chr4 validation (real, already run) is the design anchor: mdl-repeat families
are 84% low-copy (`copies < 5`, median 3); the BED `copies` column equals the BED
instance count exactly (2914/2914), so the BED is a faithful record of the copies
mdl-repeat *accepted*. Independent RepeatMasker/blastn checks show the BED is
**near-complete for full-length copies of low-copy families** (ratio ≈ 1) **but
systematically misses divergent/fragmented homologs**: for R=150, BED lists 2
instances, yet blastn at ≥70% identity finds 34 hits — only 4 full-length
(cov ≥ 80%), the remaining ~30 being 12–29%-divergent partials. So the members
that "finding all copies" must recover are **the divergent ones**, and those
divergent members are precisely the input that goal-1 subfamily clustering needs.
This couples the two goals at the data level: the recall is not chasing more
young copies (BED already has them), it is chasing the divergent tail that both
reveals subfamily structure and full-length-extends boundaries.

### 2.2 Primary path — BED assembly (M1, retained verbatim)

`phase1_extract.extract_padded_copies` is the cheap, correct base: linear in
repeat bp via `samtools faidx`, strand-co-oriented, padded for extend-align-trim.
For a family whose BED copy set is already complete and homogeneous (§2.3), this
is the *entire* copy-assembly step — no genome search, exactly the "don't spend
compute where it isn't needed" half of selectivity.

### 2.3 The completeness signal: "does this family still have un-recruited members?"

A family is flagged **incomplete → recall-eligible** when any of the following
fire (cheap, computed from data already in hand — no genome alignment to decide
*whether* to recall):

| Signal | Source | What it indicates | Cost |
|---|---|---|---|
| **BED internal divergence high** | per-instance divergence (BED score col, already in the index) — spread/variance across the family's instances | the family already spans a wide divergence band ⇒ a divergent tail likely exists below mdl-repeat's acceptance | O(instances), free |
| **Copy count vs expectation mismatch** | family `copies` vs length-class prior (a 300 bp element with 2 copies is suspicious; SINE-length with 2 is implausibly low) | mdl-repeat under-instanced this family | O(1), free |
| **Cheap k-mer pre-screen hit beyond BED** | sourmash/minimizer sketch of the seed vs a genome sketch index (built once) → estimated homolog count ≫ BED count | divergent/distal homologs exist that the BED's identity floor excluded | O(genome sketch), amortized once |
| **Boundary-truncation signal** | M2 occupancy at the BED-span edge stays solid (extension wanted) but few copies reach it | full-length members are under-represented | computed during a cheap first-pass MSA |
| **Under-instanced** | BED instances `< min_copies_for_msa` | the implemented M5 trigger — kept as one signal among several, no longer the *only* trigger | O(1), free |

The first three are computable in Phase 0 / early Phase 1 without per-family genome
alignment, so the *decision* to recall is cheap even though recall itself is not.
This is the crux of "selective ≠ tolerating incompleteness": every family still
*targets* its full member set, but the expensive recall is spent only on the
families that signal a missing tail. Families that pass all signals are declared
complete and proceed to clustering on their BED copies alone.

### 2.4 Phase 0 emits the pre-signals (the one new Phase 0 output)

`write_instance_index` already streams the BED once per family; v2 extends that
single pass to also emit a tiny per-family **completeness pre-signal** record
(BED divergence spread, instance count, length-class expectation flag) alongside
the index. No new BED pass, no genome touch. The k-mer pre-screen index (genome
sketch) is built once at Phase 1 start (shared, like the lazy blastn DB) and
queried per flagged family.

### 2.5 The recall mechanism: cheap pre-screen → per-family-capped targeted alignment

**Decision: a two-tier recall that never runs O(families × genome) full blastn.**

- **Tier 1 — cheap pre-screen (genome-scalable gate).** A single genome k-mer
  sketch index (sourmash, already a Pan_TE dependency for TE-looker) is built once.
  Each flagged family's seed is queried against it to (a) confirm a divergent tail
  exists and (b) localize candidate regions. Cost is sketch-comparison, sub-linear
  in genome bp, shared index. This replaces "blind whole-genome blastn per family."
- **Tier 2 — bounded targeted alignment.** Only for families the pre-screen
  confirms, run the existing **M5 `recruit_by_blastn`** but scoped: capped at
  `blastn_max_targets` recruited copies per family (M5 already enforces this via
  `_blast_hits_to_instances(max_keep=...)`), and — where the pre-screen localized
  candidate regions — restricted to those subjects rather than the whole DB. The
  recruited subject spans flow into the same M1 padded extractor, so recruited
  copies and BED copies are identical `Copy` objects downstream.
- **Per-family cap and global budget.** Each family is capped (copy count + a
  wall-budget); a global cap on the number of recall-triggered families bounds
  total cost. Over-budget families keep their BED copy set and are flagged in the
  metric output as "completeness not fully verified" (honest, not silently
  degraded — §6).

**Cost envelope.** BED assembly is O(repeat bp). The recall decision is O(instances
+ one sketch query per flagged family). Tier-2 alignment runs for the minority of
flagged families, each capped. The dominant term stays linear-ish in repeat
content, never the O(families × genome) blastn sweep the v1 design correctly
fought — but unlike v1, recall is *on by default for flagged families* rather than
an off-by-default failsafe.

### 2.6 Division of labor with TE-looker (D_cc 45–75% twilight)

The boundary is **refined**, not abandoned. The recall here aims at members that
are *the same family as a known mdl-repeat seed* but fell below mdl-repeat's ~70%
acceptance — i.e. divergent members reachable from the seed by targeted alignment
(the R=150 ~70%-identity tail). It does **not** attempt de novo discovery of
families that have *no* mdl-repeat seed at all in the 45–75% band — that remains
TE-looker's job (Stage 1 tri-track, de novo). Concretely:

- **Refine recalls**: members homologous to an existing mdl-repeat family seed,
  down to the pre-screen's detectable identity floor (configurable; default keeps
  the productive ≥~60–70% band where alignment is reliable and the consensus stays
  coherent).
- **TE-looker owns**: families with no mdl-repeat anchor, and the deep-twilight
  (D_cc 45–60%) where sketch homology to a seed becomes unreliable.

This is a productive overlap, not a duplication: recall makes each *known* family
complete (improving the masking-track HMM seed Pan_TE passes to TE-looker), while
TE-looker finds the *unknown* families. The metric output (§6) reports the recall's
identity-floor so the handoff boundary is auditable.

---

## 3. Subfamily clustering and multi-consensus (goal 1, the "always cluster" decision)

**Decision: always cluster the assembled copy set by pairwise divergence; every
qualifying cluster emits one subfamily consensus through the existing M2
extend-align-trim. Homogeneous family → 1 cluster → 1 consensus; divergent family
→ k clusters → k consensus. Conservative cluster-qualification thresholds prevent
over-splitting (goal 2 echo).**

### 3.1 Why single-consensus fails divergent families (the design motivation)

A single consensus sits at the centroid of the divergence spectrum. For a family
whose members range over, say, 70–95% mutual identity, that centroid is ~15%
diverged from both extremes; a downstream RepeatMasker/blastn search seeded by it
recovers the central members but misses both tails — exactly the "找不回远端拷贝"
failure. Splitting into subfamily consensi, each the centroid of a *tight*
cluster, gives each tail a near-by query. This is standard practice (RepeatModeler2
emits subfamily models; RECON/RepeatScout lineages are routinely sub-clustered),
so it is reviewer-defensible.

### 3.2 Clustering substrate: pairwise divergence over assembled copies

- **Distance.** Pairwise identity/divergence among the assembled copies. Two
  tractable sources: (a) a k-mer/minimizer distance (cheap, scales to large copy
  sets, no all-vs-all alignment) for the first cut, refined by (b) alignment-derived
  identity within candidate clusters. The BED per-instance divergence is a *1-D*
  proxy (distance to the seed) and is insufficient alone — two members equidistant
  from the seed can be far from each other — so the clustering uses pairwise copy×copy
  distance, not the 1-D BED score.
- **Algorithm.** Hierarchical clustering (scipy linkage + a divergence-cut) is the
  natural fit and is already precedented in the *old* `bin/Refiner/phase3_consensus_building.py`
  (BLASTN recruitment → scipy hierarchical clustering → tiered consensus); v2
  adapts that proven pattern. Graph community detection (the same Leiden family
  TE-looker uses) is the alternative for very large copy sets; the choice is by
  copy-set size, mirroring M2's MAFFT-vs-abPOA regime split.
- **Cut threshold.** The divergence cut that defines subfamily boundaries is the
  central tunable. A conservative default (a relatively *coarse* cut, e.g. split
  only when sub-lineages differ by more than a clear margin) errs toward *fewer*
  subfamilies, honoring goal 2's "don't over-split." The cut is a config parameter
  with a documented biological basis (subfamily-level divergence, not member-level
  noise), reviewable and tunable.

### 3.3 "Qualifying cluster" thresholds (conservative gate against over-splitting)

A cluster emits its own consensus only if:

| Criterion | Rationale |
|---|---|
| **member count ≥ `subfamily_min_members`** (e.g. ≥ `min_copies_for_msa`) | a 1–2-member cluster cannot support a trustworthy consensus; its members fold into the nearest qualifying cluster instead of spawning a singleton subfamily |
| **intra-cluster homogeneity ≥ floor** (mean pairwise identity above a bound, or divergence spread below a bound) | guarantees the cluster is actually tight enough that its centroid represents its members; a loose cluster is re-cut, not emitted |
| **cluster distinctness from siblings** (inter-cluster divergence ≥ the cut) | prevents emitting two near-identical subfamilies that post-refine merge (§7) would collapse anyway — cheaper to not split them |

Members in sub-threshold clusters are **not discarded** — they are reassigned to
the nearest qualifying cluster (so no real member is lost), or, if a whole family
is too sparse to form even one qualifying cluster, the family falls back to a
single all-copy consensus via the existing M2 path (the never-regress invariant
still holds). This makes "always cluster" safe: the worst case degrades gracefully
to today's single-consensus behavior, never to data loss.

### 3.4 Per-cluster consensus = M2 extend-align-trim, called once per cluster

Each qualifying cluster's copies are handed to the **existing**
`phase1_boundary.refine_family` / `build_msa` / `trim_termini` / `build_consensus`
pipeline unchanged — M2 is simply invoked per cluster instead of per family. Every
M2 invariant carries over per cluster: never-regress fallback, occupancy trim with
interior protection, MAFFT/abPOA regime by cluster size, ACGTN-only consensus. The
boundary extend/trim now operates on a *homogeneous* cluster, which is strictly
better for occupancy signal (the distal-member dilution that could blur a mixed
family's boundary is gone).

### 3.5 Subfamily naming and lineage record (provenance)

Each emitted consensus carries a stable id and a lineage record: parent family R,
subfamily index, member count, intra-cluster mean identity, and the cut threshold
used. Naming convention extends the existing `mdl_R<N>` scheme, e.g.
`mdl_R<N>_sf<k>` (and the existing `_chimfrag<j>` suffix composes orthogonally for
a subfamily that is *also* structurally split — §5). The lineage is written to the
record dict (like the existing `consensus_source` / `is_chimera_fragment` fields)
so Phase 2 and downstream Combine see normal records with provenance, and the
metric harness (§6) can trace a recovered member back to the subfamily that
recovered it. **No process metadata** (version strings, decision narration) enters
the FASTA headers — only biological provenance.

### 3.6 Single-vs-multi decision is emergent, not a branch

There is no "decide whether to split" branch. Clustering always runs; a homogeneous
family simply produces one qualifying cluster and therefore one consensus. This is
the literal realization of the locked decision "总是聚类，达标簇各出共识，不是仅当
单条不足才拆" — the number of consensi is an output of the divergence structure, not
a gated special case.

---

## 4. Conservative false-positive reduction (goal 2)

**Decision: prune only families/clusters with clearly absent support; weight every
threshold toward sensitivity; let the recall result (§2) inform the real-vs-noise
call. When in doubt, keep.**

### 4.1 mdl-repeat noise types and their conservative discriminants

| Noise type | Conservative discriminant | "Keep when in doubt" boundary |
|---|---|---|
| **Low-complexity / SSR** | DUST + Shannon entropy (already in Phase 0 hard_filter and Phase 2 QC) | thresholds already tuned permissive; v2 changes nothing here — this noise is reliably caught and rarely ambiguous |
| **Over-fragmented pieces** | Phase 0 BED co-occurrence fragment assembly already re-joins adjacent fragments into scaffolds; a residual short fragment that *also* recruits a coherent copy set is kept (it may be a real short element, e.g. a solo LTR or MITE) | only drop a short fragment if it recruits **no** coherent copy set AND fails TE-structure signal — never drop on length alone |
| **Pseudo-family (the key one)** | **a family is pseudo iff recall cannot assemble a coherent, connected copy set** — i.e. its purported copies do not form a clusterable group with mutual homology. This is the §1 governing principle operationalized: real families recruit coherent copies, noise does not | require *positive* evidence of incoherence (copies mutually unrelated / un-clusterable), not mere low copy count — an 84%-of-families-are-low-copy genome (chr4) means low copy count is the norm, not a noise signal |
| **Chimera (A+B)** | the existing M3 structural splitter — *splits* rather than drops, so no real sequence is lost; the two halves each survive if each has support | M3's two-signal AND-gate is already conservative (favors under-splitting); retained as-is |

### 4.2 The recall ↔ real/false coupling (governing principle, made concrete)

The pseudo-family discriminant is *defined* through recall: "find all copies"
yields, per family, the assembled copy set; a real family's copies form one or
more coherent clusters (§3), a pseudo-family's "copies" are an incoherent grab-bag
that no cut produces a qualifying cluster from. So the same clustering pass that
serves goal 1 (subfamily resolution) serves goal 2 (a family that produces *zero*
qualifying clusters and whose copies are mutually incoherent is the pseudo-family
flag). This is why §1 insists copy-completeness is the shared foundation:
**without the recalled copy set there is no principled way to tell a real low-copy
family from a pseudo-family** — copy count alone cannot, given the chr4 reality.

### 4.3 Conservative boundary, stated explicitly

- A family/cluster is dropped **only** when it has *positive* evidence of being
  unsupported (incoherent copy set, low-complexity, no TE structural signal) —
  never on a single weak signal, and never on low copy count alone.
- A family that *cannot be recall-verified* (over budget, §2.5) is **kept**, not
  dropped, and flagged — absence of verification is not evidence of absence.
- Splitting (M3 chimera, §3 subfamily) is preferred over dropping wherever the
  sequence might carry real signal, because a split is recoverable downstream and a
  drop is not.
- Downstream Combine (CD-HIT-EST 90% across all three sources) and TE-looker's own
  gates provide additional FP control, so refine deliberately under-prunes and
  hands off a *sensitive* set rather than a *clean* one.

---

## 5. Coordinating the two orthogonal split axes (structural chimera vs subfamily divergence)

These are independent axes and v2 keeps both:

- **M3 structural chimera** answers "are two *different* elements concatenated into
  one consensus?" Signal: an interior occupancy *valley* with copy-span
  *bimodality* (copies split into left-reaching and right-reaching groups). Splits
  *along the sequence* (A | B).
- **Subfamily divergence (new)** answers "is this *one* element diverged into
  sub-lineages?" Signal: clustered structure in the *pairwise copy×copy distance
  matrix* (copies fall into divergence groups, each spanning the *whole* element).
  Splits *across the copy set* (subfamily 1 vs subfamily 2), each spanning the full
  length.

They cannot be confused because their signals are orthogonal: a chimera shows
*positional* bimodality (copies cover different *parts*); a multi-subfamily family
shows *whole-length divergence* groups (copies cover the *same* part at different
identities). **Order: cluster first (§3), then run M3 per emitted subfamily
consensus.** Rationale: clustering on the full copy set first means each
subfamily's MSA is homogeneous, so M3's occupancy-valley signal is cleaner (a mixed
family's divergence could otherwise masquerade as a shallow valley and confuse the
chimera gate). Running M3 inside each subfamily also means an A+B element that
*also* has subfamilies is handled correctly: subfamily split first, then each
subfamily consensus is checked for an internal A+B junction, composing as
`mdl_R<N>_sf<k>_chimfrag<j>`. The two axes never operate on the same signal, so
there is no double-counting and no interference.

---

## 6. Representativeness / quality metrics (verifiable design targets)

The design goal "the library can recover all family members" must be *measured*,
not asserted. The metric harness (computed on real data, never fabricated):

- **Member recovery rate (primary).** Take the emitted consensus library, search
  it back against the genome (RepeatMasker or blastn), and compute, per family, the
  fraction of that family's *known true members* (the union of BED instances + the
  recall-confirmed divergent members from §2) that are hit at a coverage/identity
  threshold. The headline number is the distribution of per-family recovery rate
  **before vs after subfamily splitting** — the design predicts the divergent-family
  tail of that distribution lifts when multi-consensus replaces single-consensus.
  Reported as a distribution on real data; **no target number is invented here** —
  the harness computes it from the run.
- **Divergent-member recovery (the R=150 case).** Specifically track recovery of the
  divergent partials that single-consensus misses (the 12–29%-diverged tail), since
  that is the exact population goal 1 targets.
- **Subfamily count vs family divergence.** Sanity curve: families with wider copy
  divergence should emit more subfamilies; a flat curve would indicate the cut is
  too coarse, a runaway curve too fine. This is a tuning diagnostic, not a result.
- **FP-side metric.** Fraction of families dropped, with the per-family reason
  (low-complexity / incoherent-copies / no-TE-signal), and a manual spot-check
  audit on a sample of dropped families to confirm none had recoverable support —
  the conservative-pruning verification.
- **Boundary metric (inherited from M2).** Refined-vs-original length distribution,
  and qualitative agreement of a well-characterized Arabidopsis family
  (ATHILA/Copia) boundary against RepBase. Reported, not asserted.

All metrics are observable invariants on real runs (chr4 Arabidopsis first, then a
larger genome), consistent with the project's no-fabricated-numbers discipline.

---

## 7. M1–M5 reuse map (reuse / extend / new)

| Module | v2 disposition | Detail |
|---|---|---|
| **Phase 0** (`phase0_triage.py`) | **reuse + small extend** | Parse, hard filter, tiering, BED fragment assembly, near-exact dedup, `write_instance_index` all unchanged. **Extend** the existing single BED pass in `write_instance_index` to also emit per-family completeness pre-signals (§2.4). No new genome pass. |
| **M1 `phase1_extract.py`** | **reuse verbatim** | `bed_to_samtools_region`, `extract_padded_copies`, `reverse_complement`, `stratified_subsample`, `compute_pad` are the copy-assembly base for both BED and recalled copies. Unchanged. |
| **M2 `phase1_align.py` / `phase1_boundary.py`** | **reuse, called per cluster** | `build_msa`, `select_regime`, `trim_termini`, `build_consensus`, `refine_family`, `alignment_qc_ok` invoked **once per qualifying cluster** instead of once per family. The never-regress fallback and interior-protected trim carry over per cluster. No edit to the functions; the *caller* (the new per-family body) changes. |
| **M3 `phase1_chimera.py`** | **reuse, called per subfamily** | `detect_chimera` (two-signal AND-gate) runs on each emitted subfamily consensus's `_aln` (§5). Unchanged; just relocated after clustering. |
| **M4 `phase1_consensus.py`** | **reuse + extend orchestration** | Sharding, instance routing, resume, post-refine merge, data-integrity reconciliation all reused. **Extend** the per-family worker body (`_process_shard`'s inner loop) to call: completeness-triage → selective recall → cluster → per-cluster M2 → per-subfamily M3. The `n_out == n_in + extra` invariant generalizes (extra now also counts subfamily splits, tracked independently). The lazy-shared-DB pattern extends to the shared genome sketch index. |
| **M5 `phase1_fallback.py`** | **reuse + promote** | `recruit_by_blastn` + `ensure_genome_blast_db` reused as Tier-2 of selective recall (§2.5). **Promoted** from "fires only when BED < min_copies" to "fires whenever a completeness signal flags the family," and **fronted** by the cheap sketch pre-screen so it is not run blindly. `enable_divergent_blast_recruitment` semantics tighten to the seed-homologous divergent band (§2.6). |
| **NEW — `phase1_completeness.py`** | **new** | Completeness signals (§2.3), pre-signal consumption, sketch pre-screen gate (§2.5 Tier 1), recall trigger decision + per-family/global budget. |
| **NEW — `phase1_cluster.py`** | **new** | Pairwise copy distance, hierarchical/graph clustering, qualifying-cluster gate, member reassignment, subfamily lineage records (§3). Borrows the scipy-hierarchical pattern from `bin/Refiner/phase3_consensus_building.py`. |
| **NEW — `phase1_metrics.py`** (or a `tools/` harness) | **new** | Member-recovery-rate computation (§6), run on the emitted library; reporting only, off the critical path. |
| **Phase 2** (`phase2_library_split.py`) | **reuse verbatim** | Consumes records; subfamily records are normal records with `tier`/`copies`/`sequence`. `consensus_masking.fa` / `phase3_analysis_library.fa` contract unchanged. (One consideration: per-subfamily `copies` should be the cluster's member count, not the parent family's, so Phase 2's copy-threshold gates act per subfamily — a field-population detail, not a Phase 2 edit.) |
| **`config.py`** | **add params** | `subfamily_min_members`, `subfamily_divergence_cut`, `subfamily_intra_homogeneity_floor`, `cluster_method` (hierarchical/graph by size), completeness-signal thresholds (`bed_divergence_spread_trigger`, copy-count-vs-length-class table), `recall_global_family_budget`, `recall_per_family_wall_s`, `sketch_index_*`, recall identity floor. Keep all M1–M5 params. |

The reuse ratio is high: the net-new code is two modules (completeness trigger,
clustering) plus a metric harness; M1–M5 are called components, not rewrites.

---

## 8. Scalability and correctness

### 8.1 Large-genome scaling

- **Copy assembly** stays O(repeat bp) for the BED path (M1). The recall decision
  is O(instances + one sketch query per flagged family) — the sketch index is built
  once and shared like the lazy blastn DB. Tier-2 alignment runs for a capped,
  flagged minority. No O(families × genome) sweep.
- **Clustering** is the new cost. Per-family pairwise distance is O(c²) in copy
  count c; bounded by the existing subsample caps (`msa_subsample_cap_mafft/abpoa`),
  and the k-mer first-cut avoids all-vs-all alignment for large c. For
  ultra-abundant families the graph-clustering regime (Leiden-style) scales where
  O(c²) would not — chosen by copy-set size, mirroring M2's MAFFT↔abPOA split.
- **Parallelism / memory** inherit M4 unchanged: shard-parallel families,
  single-threaded MSA per family, `samtools faidx` random access bounds resident
  memory to the family in flight. Per-shard `.done` resume is unchanged. The new
  per-family work (cluster + per-cluster MSA) stays inside one worker, so the
  parallel model and memory envelope are preserved.

### 8.2 Conservative pruning does not lose real data

- Clustering never discards members (sub-threshold clusters reassign to nearest
  qualifying cluster; a too-sparse family degrades to single-consensus — §3.3).
- Pruning requires *positive* evidence of non-support; unverifiable families are
  kept and flagged (§4.3).
- The M4 data-integrity reconciliation (`output == input + fragments_added −
  merged`) generalizes: `fragments_added` now also counts subfamily splits, each
  tracked from an independent source so the cross-check stays a real audit, not a
  tautology. Any silent family loss still trips the existing `DATA INTEGRITY` error.

### 8.3 Reviewer acceptability

- Subfamily multi-consensus is **standard** (RepeatModeler2 subfamily models;
  RECON/RepeatScout sub-clustering). The divergence-cut clustering is the accepted
  approach; v2's only stance is a *conservative* cut.
- The recall is bounded, seed-homologous, and identity-floor-documented — not a
  fishing expedition; its boundary with TE-looker is auditable via the metric
  harness.
- The two-axis split (structural chimera vs subfamily divergence) is justified by
  orthogonal signals, and each axis individually mirrors accepted tools (M3 ≈
  composite-element splitting; subfamily clustering ≈ RepeatModeler2).
- **Limitations to state explicitly**: (a) recall cannot recover a family with *no*
  full-length genomic copy (a fundamental limit, not a flaw); (b) the divergence cut
  is a tunable with no universal optimum — the conservative default trades some
  subfamily resolution for fewer false splits, and the metric harness exposes the
  trade so a user can retune; (c) the pre-screen identity floor sets where refine's
  responsibility ends and TE-looker's begins — families in the deep twilight are out
  of scope by design.

---

## 9. What v2 supersedes from the v1 design

v1 (`REFINE_STRATEGY_DESIGN.md` / `REFINE_IMPLEMENTATION_PLAN.md`) was a sound
*performance* redesign (BED-seeded extend-align-trim replacing whole-genome BLAST),
and v2 keeps its entire engineering substrate. v2 overturns three v1 *assumptions*
that conflict with the now-explicit goals:

1. **"Recall is an optional, default-off failsafe."** v1 §2.2 / config
   `enable_divergent_blast_recruitment=False` framed recruitment of divergent
   copies as off-by-default and triggered only by under-instancing. **Superseded:**
   recall is goal-1 core and *cannot* be off; it is made *selective* (signal-driven,
   cheap-pre-screened, budgeted) instead of optional. Completeness is the target for
   every family, not a fallback for sparse ones.
2. **"One consensus per family."** v1's whole Phase 1 emits a single consensus per
   family (extend-align-trim), with M3 chimera as the only multiplier.
   **Superseded:** always-cluster multi-consensus is the default output shape; a
   single consensus is now the *special case* of a one-cluster family, not the norm.
3. **"The 45–75% band is purely TE-looker's remit; refine must not touch it."**
   v1 §2.2/§8.1 drew a hard line at mdl-repeat's ~70% acceptance. **Refined (not
   fully overturned):** refine *does* recall the seed-homologous divergent tail
   below ~70% (the R=150 population), because those members are needed for subfamily
   resolution and boundary completion; TE-looker still owns *de novo* discovery of
   anchorless families and the deep twilight. The line moves from "70% identity" to
   "homologous-to-a-known-seed vs de-novo," with an auditable identity floor.

The deeper reframe: v1 optimized *cost of building one consensus*; v2 reprioritizes
around *copy-completeness as the foundation* — once you commit to finding every
real member, single-consensus and default-off recall both become untenable, and the
pipeline reorders to A(find-all) → B(cluster) → C(per-cluster consensus) → D/E
(orthogonal chimera split + conservative prune).

---

## 10. Implementation milestone re-ordering (built on M1–M5)

Since M1–M5 exist and pass their suites, v2's new work is additive and ordered by
dependency, each step validated on real chr4 Arabidopsis data with observable
invariants (no fabricated targets):

- **N1 — completeness signals (Phase 0 extend + `phase1_completeness.py` signals,
  no recall yet).** Compute and emit the per-family completeness pre-signals; on
  chr4, verify the families flagged "incomplete" are enriched for the known
  divergent-tail families (e.g. R=150 flags; a clean low-copy homogeneous family
  does not). *Pass = flag set matches the independent blastn-tail evidence on a
  spot-check sample; decision cost is O(instances), no genome alignment.*
- **N2 — subfamily clustering (`phase1_cluster.py`) on BED copies only.** Wire
  clustering between M1 extract and M2, calling M2 per cluster; recall still off.
  *Pass = a homogeneous family yields exactly 1 cluster → 1 consensus (matches
  today's output); a known divergent family yields >1 subfamily; no member lost
  (reassignment invariant holds); M4 data-integrity reconciliation passes with
  subfamily splits counted in `fragments_added`.*
- **N3 — per-subfamily M3 ordering (§5).** Run M3 after clustering, per subfamily.
  *Pass = a constructed subfamily-of-a-chimera case yields `_sf<k>_chimfrag<j>`
  correctly; a pure multi-subfamily family is not additionally chimera-split; a pure
  A+B chimera with one subfamily behaves as today.*
- **N4 — selective recall (sketch pre-screen Tier 1 + promoted M5 Tier 2 +
  budgets).** Turn on signal-driven recall for flagged families. *Pass = recall runs
  only for flagged families (count matches N1 flags), each capped; on R=150 the
  recalled divergent partials appear and feed an additional/wider subfamily; the
  shared sketch index builds once; over-budget families keep BED copies and are
  flagged, never silently dropped.*
- **N5 — conservative FP prune coupling (§4).** Wire the pseudo-family discriminant
  (zero qualifying clusters + incoherent copies) into a drop decision, sensitivity-
  weighted. *Pass = a manufactured incoherent grab-bag family is dropped; every real
  low-copy family is kept; the dropped-family audit sample shows none had recoverable
  support.*
- **N6 — representativeness metric harness (`phase1_metrics.py`) + full run.** Run
  end-to-end through Phase 2; compute member-recovery-rate before/after subfamily
  splitting. *Pass = Phase 2 contract unchanged (`consensus_masking.fa` +
  `phase3_analysis_library.fa` same schema); the recovery-rate distribution is
  produced on real data and its divergent-family tail is observably higher with
  multi-consensus than with single-consensus (reported, not asserted).*

Each milestone is independently testable and degrades safely to the current
(passing) behavior if its new component is disabled, so the path from the working
M1–M5 pipeline to the v2 goals is incremental and reversible at every step.
