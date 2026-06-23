from dataclasses import dataclass
import json


@dataclass
class RefinerMdlConfig:
    """Configuration for Refiner_mdl pipeline."""

    # Input/Output
    input_file: str = ""         # mdl-repeat output FASTA
    genome_file: str = ""        # Reference genome FASTA
    output_dir: str = ""
    temp_dir: str = "temp_work"
    checkpoint_dir: str = "checkpoints"

    # Optional mdl-repeat auxiliary files
    bed_file: str = ""           # mdl-repeat -instances output
    stats_file: str = ""         # mdl-repeat -stats output

    # Phase 0: Hard filtering
    min_length: int = 50
    max_n_percent: float = 0.2
    entropy_threshold_long: float = 1.0   # Shannon entropy (k=1) for seq >= 100bp
    entropy_threshold_short: float = 0.8  # Shannon entropy (k=1) for seq < 100bp
    # k=2 (dinucleotide) Shannon entropy — catches dinucleotide-repeat low-complexity
    # (e.g. (CA)n) that the k=1 gate passes. Max possible = 4 bits; clean TEs ~3.1-3.9.
    entropy2_threshold_long: float = 2.0   # for seq >= 100bp
    entropy2_threshold_short: float = 1.6  # for seq < 100bp (fewer dinucleotides sampled)
    # Max fraction of the sequence in low-complexity windows (per-window DUST). Guards
    # against mostly-simple consensi that slip the average DUST gate.
    lowcomplexity_fraction_max: float = 0.6
    min_mdl_score: float = 0.0            # MDL > 0
    dust_score_threshold: float = 0.7     # DUST > 0.7 = SSR/simple repeat

    # Phase 2: Masking library — length-dependent copies threshold
    # Longer sequences need fewer copies (high specificity);
    # shorter sequences need more copies evidence.
    masking_copies_long: int = 2          # length >= 500bp
    masking_copies_medium: int = 3        # 200bp <= length < 500bp
    masking_copies_short: int = 5         # 100bp <= length < 200bp
    masking_copies_vshort: int = 10       # length < 100bp

    # Phase 2: Analysis library — T3 inclusion threshold
    analysis_t3_min_copies: int = 3       # T3 needs copies >= this for analysis lib

    # Phase 2: TE-structure handling. Structure/homology cannot decide true-vs-false
    # repeat (divergent / dark-matter TEs lack detectable structure), so by default
    # the four checks only TAG records (sidecar TSV) and drop nothing; real-vs-noise
    # is handled upstream by N5 (recurrence/coherence). Set True to restore the legacy
    # hard gate (drops sequences with no TE signal) for explicit comparison only.
    enable_te_structure_gate: bool = False

    # Phase 2: mdl-repeat NATIVE quality gate ─────────────────────────────
    # mdl-repeat emits a per-family confidence in the FASTA header (tier=core/warn/reject,
    # accept=exclusive/standalone/…). On large, repeat-rich genomes the bulk of written
    # families are low-confidence (warn-tier / standalone / few copies) and inflate the
    # library with redundant fragments — yet Refiner_mdl otherwise discards this signal.
    # This gate keeps a family iff it is mdl high-confidence (tier in keep_tiers) OR carries
    # independent TE structural signal (protein/LTR/TIR/known-homology) OR has high copy
    # support; it drops only the no-signal ∩ low-tier ∩ low-copy intersection. Records whose
    # header lacks tier= (older mdl-repeat builds) are never gated (kept), preserving old
    # behaviour. The structural-signal rescue protects divergent/dark-matter TEs that mdl
    # flagged warn — so it does not repeat the over-aggressive pure-tier=core mistake.
    enable_mdl_quality_gate: bool = True
    mdl_quality_keep_tiers: str = "core"   # comma-separated mdl tiers kept unconditionally
    mdl_quality_highcopy_keep: int = 50    # copies >= this kept even at a low mdl tier

    # Phase 2: Protein-homology gene exclusion (negative filter) ──────────
    # Multi-copy host gene families (F-box, NBS-LRR, kinases, transporters …) are
    # mistaken for "repeats" by mdl-repeat and over-mask coding sequence. This
    # negative filter drops a consensus iff it (a) strongly matches a CELLULAR
    # protein (external DB, TE proteins removed — gene-annotation-free, available
    # before the genome's own gene annotation) AND (b) shows ZERO TE evidence from
    # annotate_te_structure (no protein_domain / tir / ltr / known_te_homology).
    # Any TE signal keeps the sequence, so divergent / no-ORF / structural TEs are
    # never dropped. DB absent -> step is skipped (tag-only); never drops without DB.
    enable_gene_exclusion: bool = False
    cellular_protein_db: str = ""         # blast DB of cellular proteins, TE removed
    gene_excl_min_cov: float = 0.5        # min cellular-protein coverage to drop (needs no-TE-evidence)
    gene_excl_drop_cov: float = 0.7       # cellular coverage at/above which we drop UNCONDITIONALLY
                                          # (even with a TE domain) — a near-full host-gene hit to the
                                          # TE-removed DB = domesticated/host gene, not a real TE
    gene_excl_min_pident: float = 35.0    # min % identity (cross-species tolerant floor)
    gene_excl_evalue: float = 1e-5
    gene_excl_max_target_seqs: int = 3
    gene_excl_jobs: int = 0               # concurrent blastx query-chunks (0 -> threads)

    # Phase 2: protein-homology flank trim (cut terminal host-gene over-extension) ─
    enable_flank_trim: bool = False
    flank_trim_min_cov: float = 0.3       # min cellular cov to consider a flank (band [min, drop_cov))
    flank_trim_margin: int = 30           # gene span must touch a terminus within this many bp
    flank_trim_retain_min: int = 100      # retained TE core must stay >= this many bp
    flank_trim_jobs: int = 0              # concurrent blastx query-chunks (0 -> threads)

    # Phase 2: Non-TE repeat separation (classification / routing) ─────────
    # Route organelle (NUMT/NUPT), rDNA (45S/5S), and satellite/tandem consensi OUT of
    # the TE library into dedicated tracks — they are recurrent but not transposons, and
    # leaving them tagged "TE" overstates genome TE content (one 45S rDNA consensus masks
    # ~12.6 Mb). Each detector skips cleanly if its tool/reference is absent.
    enable_nonte_separation: bool = False
    organelle_ref: str = ""               # FASTA of chloroplast + mito genomes
    barrnap_exe: str = "barrnap"          # rRNA detection (eukaryote kingdom)
    trf_exe: str = "trf"                  # satellite/tandem detection
    organelle_min_pid: float = 80.0       # min % identity for organelle blastn
    organelle_min_cov: float = 0.5        # min summed coverage for organelle call
    satellite_tandem_frac: float = 0.5    # min fraction of length under tandem repeat

    # Phase 2: Short low-copy noise filter (genome-size adaptive) ──────────
    # Drop the short-AND-low-copy corner of mdl-repeat output (marginal MDL fragments,
    # ~50% of entries but <10% of masking). JOINT gate so long-low-copy full elements
    # and short-high-copy MITEs are both kept. min_copies scales with genome size.
    enable_lowcopy_filter: bool = False
    genome_size_bp: int = 0               # set by orchestrator; 0 -> conservative base
                                          # In large-genome SAMPLED mode this is the FULL
                                          # genome size (from --genome-full), NOT the sample,
                                          # so the genome-size-adaptive copy tiers stay correct.
    lowcopy_max_len: int = 200            # below this, a low-copy entry is fragment noise
    lowcopy_min_copies_override: int = 0  # >0 forces min_copies (else genome-size tier)
    hard_min_copies: int = 0              # >0: drop ALL entries with copies < this,
                                          # regardless of length (hard recurrence floor;
                                          # aligns mdl with TE-looker's --min-copies handoff)

    # Phase 2: genomic copy recruitment via rmblastn (reliable copy count) ─
    # mdl `copies` is an unreliable MDL estimate; recruit the real genomic copy number by
    # rmblastn vs the genome and write rec['genomic_copies'], which the low-copy / hard-
    # floor filter then reads in preference to `copies`. Skips -> falls back to mdl copies.
    enable_copy_recruit: bool = False
    rmblastn_exe: str = "rmblastn"
    rm_matrix_dir: str = ""               # ncbi/nt matrix dir; auto from CONDA_PREFIX
    copy_recruit_divergence: int = 25     # rmblastn matrix 'p' (divergence tolerance)
    copy_recruit_cov: float = 0.5         # min consensus coverage for an instance = 1 copy
    copy_recruit_gapopen: int = 30
    copy_recruit_gapextend: int = 6
    copy_recruit_evalue: float = 1e-5
    copy_recruit_timeout_s: int = 21600
    copy_recruit_jobs: int = 0            # concurrent rmblastn query-chunks (0 -> threads)

    # Phase 1: rmblastn copy source (#1) + nested-seed de-nesting (N7) ─────
    # When enabled, Phase 1 recruits each family's copies by rmblastn (instead of mdl BED)
    # — the diverged members the BED misses make a more representative consensus. The same
    # per-family rmblastn drives N7 seed-chimera resolution. Both default OFF (ship dark);
    # OFF => byte-identical to the BED path.
    enable_rmblastn_copies: bool = False
    copy_extract_min_cov: float = 0.5     # min consensus coverage for an rmblastn instance
                                          # to be used as a COPY (drops partial gene-paralog
                                          # matches from MSA; N7 still sees all loci)
    enable_seed_chimera: bool = False
    seed_chimera_blast_timeout_s: int = 600
    chimera_seed_bin_bp: int = 20
    chimera_seed_margin_bp: int = 30
    chimera_seed_max_span_frac: float = 0.15
    chimera_seed_min_side_copies: int = 4
    chimera_seed_min_disjoint: float = 0.8
    chimera_seed_min_component_bp: int = 80
    chimera_seed_max_breakpoints: int = 3
    chimera_seed_require_standalone: bool = True   # post-split improvement gate

    # Phase 1: per-copy structural (TSD-anchored) boundary trim ────────────
    # Trim each copy to its TSD-defined TE boundary before MSA, so the consensus does not
    # over-extend into a host-gene flank shared across copies (the dominant CDS over-mask).
    enable_structural_trim: bool = False
    tsd_min: int = 4                      # min TSD length (bp)
    tsd_max: int = 20                     # max TSD length (bp)
    tsd_slack: int = 15                   # junction search slack each side (bp)
    struct_trim_min_te: int = 50          # min TE length kept after a TSD trim
    struct_trim_min_copies: int = 4       # need >= this many copies to vote a consensus TSD
    struct_trim_min_frac: float = 0.3     # consensus TSD must be in >= this fraction of copies
    struct_trim_offset_bin: int = 3       # bin TSD boundary offsets (bp) to tolerate jitter

    # Copy-gate rescue: recover dropped-but-TE-structured consensi (TEsorter HMM +
    # LTR/TIR/known-TE). Rescue-only (never drops), so safe for divergent TEs.
    enable_rescue: bool = False
    tesorter_exe: str = "TEsorter"
    tesorter_db: str = "rexdb"            # REXdb viridiplantae+metazoa HMM
    tesorter_timeout_s: int = 14400

    # Phase 0: Fragment assembly — directed greedy chain algorithm
    fragment_gap: int = 50                # Max gap between adjacent instances (bp)
    cooccurrence_threshold: int = 3       # Min independent co-occurrence sites
    cooccurrence_ratio: float = 0.3       # Min co-occurrence / min(copies) ratio
    min_direction_consistency: float = 0.8  # Min fraction of instances in dominant order
    min_copy_ratio: float = 0.3           # Min min(copies)/max(copies) between pair
    max_overlap_fraction: float = 0.5     # Max overlap as fraction of shorter seq
    max_scaffold_length: int = 0          # 0 = adaptive (P95 input length × 2, cap 50kb)
    # Fragment-chain stitching (shared with the te-looker track via fragment_stitcher).
    stitch_gap: int = 300                 # Max gap between co-linear fragments (was fragment_gap=50;
                                          # recovers fragments split by a short divergent linker)
    stitch_window_k: int = 4              # Look-ahead window: skip an interleaved 3rd family
    stitch_max_len: int = 25000           # Full-length ceiling (let fragmented LTR-RTs reassemble)

    # Phase 0: Deduplication
    # Single CD-HIT-EST round at near-exact identity (pre-refine = collapse only
    # duplicate descriptions of the same family).  The looser merge that used to be
    # the 0.90 second round is DEFERRED to Phase 1's post-refine pass (boundary-
    # convergent duplicates are only exposed after extend-align-trim).
    dedup_round1_identity: float = 0.98
    dedup_round2_identity: float = 0.90  # back-compat field; UNUSED in Phase 0

    # Phase 1 (rewrite M1): BED-seeded copy extraction + coordinate handling
    pad_fraction: float = 0.2             # pad per side as fraction of element length
    pad_min: int = 50                     # floor: exposes short truncated termini
    pad_cap: int = 500                    # ceiling: avoids pulling in next/nested TE
    extract_batch_size: int = 200         # positional faidx batch (arg-length safe)
    samtools_supports_region_file: bool = False  # samtools 1.4.1 has no -r flag
    subsample_seed: int = 42              # deterministic subsampling (G10 posture)
    min_copies_for_msa: int = 5           # below this, keep mdl-repeat seed

    # Phase 1 (rewrite M2): MSA regime selection + extend-align-trim boundary refine
    msa_subsample_cap_mafft: int = 100    # progressive MAFFT tractable to ~100 seqs
    msa_subsample_cap_abpoa: int = 300    # POA scales further; bound large-family runtime
    abpoa_exe: str = "abpoa"   # resolved from PATH (conda provides it); override via config
    # MSA subprocess wall-clock cap (seconds).  HARD timeout on BOTH the abPOA and the
    # MAFFT subprocess (the previous code left abPOA at a 1800 s timeout and MAFFT
    # un-bounded).  A pathological large family — e.g. chr4 R=2 (a 20 kb LTR × 65 copies)
    # — drives abPOA into a >4 TB band-matrix allocation that OOM-kills or hangs the
    # whole shard for minutes; this was the exact cause of the chr4 full-run stall.  On
    # timeout build_msa returns None and M2's never-regress fallback keeps the original
    # mdl-repeat seed (never faked, never hung).  300 s comfortably covers every real
    # chr4 family while bounding the pathological tail.
    msa_wall_s: int = 60       # per-MSA subprocess wall; backstop for any slow family
                               # that slips the sweet-spot guard (then falls back to seed)
    # Copy-set guards for BOTH MSA engines (abPOA and MAFFT) — refuse to even LAUNCH the
    # aligner when the subsampled work set is large enough that abPOA's banded DP would
    # blow memory before the wall clock could fire (an OOM-kill is not catchable by a
    # subprocess timeout) OR MAFFT --auto would grind for its full 600 s wall before
    # falling back (the dominant per-family cost on a full run).  The danger case for both
    # is a small number of VERY LONG copies (giant LTR/retrotransposon elements, e.g. the
    # chr4 R=2 family: 65 copies up to ~400 kb with padding, ~3.6 Mbp total), not copy
    # count alone.  When EITHER guard trips, build_msa returns None up front and the caller
    # keeps the original mdl-repeat seed (never-regress, no wasted minutes).  0 disables a
    # guard.  Defaults sized from the chr4 R=2 blow-up measured on the real genome: the
    # per-copy length cap (50 kb) catches the giant element, the total-bp cap (8 Mbp)
    # catches a many-medium-copy aggregate while leaving ordinary large families (the
    # largest ordinary chr4 family is well under 1 Mbp total) untouched.
    abpoa_max_total_bp: int = 8_000_000   # total residues across the subsampled work set
    abpoa_max_seq_bp: int = 50_000        # any single copy longer than this -> skip MSA
    small_family_threshold: int = 30      # L-INS-i affordable below this
    small_family_max_div: float = 0.1     # L-INS-i reserved for low-divergence families
    large_family_threshold: int = 200     # MAFFT -> abPOA switch (by copy count)
    # MAFFT -> abPOA switch by TOTAL MSA work (n_copies × avg seq length): MAFFT on many
    # long sequences is minutes (R=1567: 100 × 8.5 kb refine >300 s), so route big-bp
    # clusters to the far faster POA regardless of count. 150 kb ≈ 100 copies × 1.5 kb.
    msa_abpoa_bp_threshold: int = 150_000
    # Hard ceiling on MSA work (n_copies × avg seq length): clusters above this are
    # down-sampled so abPOA/MAFFT stay bounded. Large divergent long-element families
    # (e.g. R=1567: 114 × 8.5 kb) cost ~80 s/cluster and fall back to the seed anyway;
    # a representative subset gives the same consensus far cheaper.
    msa_max_total_bp: int = 300_000
    # Sweet-spot guard: skip subfamily MSA refinement (keep the mdl-repeat seed) when the
    # element is longer than this, or n_copies × consensus_len exceeds the work budget.
    # Such families are minutes-scale to MSA and fall back to the seed anyway (measured),
    # and long elements are Look4LTRs' remit. Value lives in shorter/moderate families.
    subfamily_max_consensus_len: int = 2500   # element-length ceiling for subfamily MSA;
                                              # longer -> instant seed (value lives ≲2.5 kb,
                                              # e.g. R=8/R=17 ~1 kb; longer elements are
                                              # Look4LTRs' remit and MSA-slow / seed-bound)
    subfamily_work_budget: int = 500_000      # n_copies × consensus_len ceiling
    trim_occupancy_floor: float = 0.3     # fractional present-count floor (CIAlign-style)
    trim_min_copies: int = 3              # absolute floor; blocks single-flank "extension"
    consensus_divergence_weighted: bool = False  # optional cleaner-copy vote weighting

    # Phase 1 (rewrite M3): MSA-based BLAST-free chimera detection (AND-gate)
    chimera_occupancy_ratio: float = 0.5  # valley = occupancy < ratio * median(occupancy)
    chimera_min_valley_cols: int = 3      # min interior valley run length / span margin
    chimera_span_group_min_frac: float = 0.2  # min copy support for each half (bimodality)
    chimera_max_span_frac: float = 0.5    # max full-span copies tolerated when splitting
    enable_bed_independence_chimera: bool = False  # optional 3rd corroborating signal (off)

    # N2 (v2): subfamily clustering — "always cluster, each qualifying cluster -> one
    # consensus" (REFINE_STRATEGY_DESIGN_v2.md §3, §7).  Clustering ALWAYS runs between
    # M1 extract and M2 consensus; a homogeneous family yields one cluster (one
    # consensus, original id, byte-identical to today), a divergent family yields
    # several (`mdl_R<N>_sf<k>`).  These thresholds govern the conservative gate that
    # keeps "always cluster" from over-splitting (§3.3) and degrade safely to today's
    # single-consensus behavior.
    # ────────────────────────────────────────────────────────────────────
    # PRIMARY distance substrate is all-vs-all BLASTN LOCAL-alignment identity among a
    # family's own copies (phase1_cluster.pairwise_distance, §3.2 "Algorithm" — the
    # proven Refiner/phase3 BLASTN-recruit→hierarchical pattern).  k-mer Jaccard is
    # RETAINED only as the size-routed coarse fallback for copy sets above
    # `cluster_graph_size_threshold` (too large for an O(c²) all-vs-all BLAST).  This
    # k-mer size governs that fallback only; k=13 is the standard repeat/TE scale.
    subfamily_kmer_size: int = 13
    # Minimum HSP coverage of the SHORTER copy for a BLASTN HSP to count as evidence in
    # the primary distance (phase1_cluster._blastn_identity_distance).  A coverage GATE,
    # not a weight: it rejects short spurious high-identity hits, then the LOCAL %identity
    # of the surviving HSP is the distance — flank-length-invariant, so M1's per-copy
    # padding does not deflate it (a coverage-WEIGHTED score does, and was MEASURED to
    # corrupt the split on real chr4: R=8 padded → intra-identity 0.25 vs 0.86 with the
    # gate).  0.4 keeps the validated chr4 splits (R=8 → 4 subfamilies, R=17 → 3) while
    # rejecting fragmentary chance hits.
    subfamily_min_hsp_coverage: float = 0.4
    # The divergence cut that defines subfamily boundaries — the central tunable
    # (§3.2 "Cut threshold").  A BLASTN local-identity DISTANCE cut (distance =
    # 1 − pident): members whose pairwise distance stays <= this merge into one
    # subfamily; sub-lineages beyond it separate.  Default 0.30 == "split only when two
    # sub-lineages are below ~70% mutual local identity" — the conservative few-splits
    # side (a false split permanently fractures a real family; a missed split only
    # reproduces today's single-consensus output, recoverable).
    #
    # CALIBRATED ON REAL chr4 (measured, see phase1_cluster module docstring +
    # tests/integration_n2_cluster.py):  at 0.30 the genuinely divergent multi-copy
    # families RESOLVE — R=8 (40 copies) → 4 subfamilies (intra-identity 0.82–0.92,
    # inter-subfamily 0.18); R=17 → 3 (intra 0.80–0.92, inter 0.08); R=40 → 2 — while
    # homogeneous / low-copy families STAY one cluster (R=113/89/94/149/216 → 1), so the
    # degrade-to-today invariant holds with no over-splitting.  An EARLIER k-mer-Jaccard
    # substrate (and a global-MSA core substrate) both saturated at ~0.5–1.0 for every
    # family and could not separate; BLASTN local identity is the substrate that works
    # because mdl-stage copies carry rearrangement/internal-repeat that defeats colinear
    # MSA.  Retune via the N6 metric harness once recall (N4) adds the divergent tail.
    subfamily_divergence_cut: float = 0.30
    # Minimum members for a cluster to emit its own subfamily consensus (§3.3).
    # Reuses min_copies_for_msa's rationale (a 1-2-member cluster cannot support a
    # trustworthy consensus); kept as an independent field so it can be tuned without
    # moving the MSA floor.  Sub-min clusters are reassigned to the nearest qualifying
    # cluster, never dropped.
    subfamily_min_members: int = 5
    # Intra-cluster homogeneity floor: mean pairwise BLASTN local IDENTITY (1 - distance)
    # within a cluster must clear this for the cluster to qualify (§3.3).  Now a REAL
    # identity scale (the BLASTN substrate resolves true identity, unlike the saturating
    # k-mer one), so the floor is meaningful: 0.6 == "an emitted subfamily's members are
    # at least ~60% mutually identical on average".  Calibrated to the real chr4 splits
    # (emitted subfamilies measured at intra-identity 0.79–0.92, well clear of 0.6) while
    # staying conservative enough not to reject a real but internally-diverse subfamily.
    # A cluster below this floor is not emitted; its members reassign to a tighter one.
    subfamily_intra_homogeneity_floor: float = 0.6
    # Copy-set size at/below which the all-vs-all BLASTN + scipy hierarchical path is
    # used; above it, the coarse k-mer connected-component graph regime (avoids the
    # O(c²) all-vs-all BLAST + linkage memory on ultra-abundant families).  Mirrors M2's
    # MAFFT↔abPOA size route.  Both cut at subfamily_divergence_cut.
    cluster_graph_size_threshold: int = 500
    # Total-bp cap for the all-vs-all distance BLAST (the dominant scaling knob — BLAST
    # cost grows with total residues, not just copy count).  A few-copy GIANT-element
    # family (e.g. a 65-copy 20 kb LTR ≈ 1.3 Mbp) would otherwise run a multi-minute
    # all-vs-all BLAST; above this cap the copies are divergence-stratified subsampled for
    # the distance (M2-style), with unsampled copies folded into the main cluster
    # (never-regress).  500 kbp keeps ordinary families whole (the largest chr4 family,
    # 428×136 bp ≈ 58 kbp, runs in ~25 s all-vs-all) while bounding the giant-element tail.
    distance_blast_max_total_bp: int = 500_000
    # Hard COUNT cap for the O(copies²) all-vs-all clustering distance (the bp cap above
    # does not bind when copies are individually small). Ultra-abundant families are
    # down-sampled to a divergence-stratified representative subset before clustering.
    cluster_distance_max_copies: int = 150
    # Drop copies longer than this multiple of the family's MEDIAN copy length before
    # clustering — a mis-annotated giant instance (measured: a 671 kb 'copy' of an 8.5 kb
    # family) is noise that bloats every pairwise comparison.
    cluster_outlier_len_mult: float = 5.0

    # N1 (v2): family-completeness pre-signals (cheap, zero genome alignment)
    # ────────────────────────────────────────────────────────────────────
    # All thresholds below decide ONLY whether a family is "recall-eligible"
    # (worth a later targeted recall in N4); nothing is dropped here.  Costs are
    # O(instances)/O(1) per family — no blastn, no sketch — so the decision is
    # always cheap even though recall itself (N4) is not.  See
    # REFINE_STRATEGY_DESIGN_v2.md §2.3, §2.4, §7.
    #
    # BED divergence-spread trigger: a family whose per-instance divergence
    # (1 - BED_score/1000) is spread wider than this already straddles a broad
    # divergence band, which signals a divergent tail exists below mdl-repeat's
    # ~70% acceptance floor (the R=150 case: instance divs {0.02, 0.676},
    # spread 0.328).  Default 0.31 is calibrated to the REAL chr4 spread
    # distribution (measured, not guessed): per-family spread has median 0.151
    # and p90 0.309, so a 0.31 floor fires only for the extreme-spread tail.
    # This is deliberately conservative — a lower floor (e.g. the median 0.15)
    # over-fires badly: on chr4 it flags ~50% of all families and ~78% of the
    # *well-instanced* families, including R=2 (a 20 kb LTR element, 65 BED
    # copies, spread 0.263 from internal mosaic divergence) whose BED is already
    # complete (blastn finds 33 hits <= 65 BED, zero divergent partials).  High
    # *internal* spread in a well-sampled family is not the same as a *missing*
    # tail, so the trigger sits at the distribution's extreme to stay selective
    # (REFINE_STRATEGY_DESIGN_v2.md §2.1 anchors thresholds to the real data;
    # §2.3 stresses selective ≠ tolerating incompleteness — N4 still targets the
    # full member set, this only governs where recall compute is spent).
    bed_divergence_spread_trigger: float = 0.31
    # Minimum instances required before the spread signal is trusted (std of a
    # 2-point sample is noisy but still meaningful for the extreme R=150 case;
    # we keep it at 2 so two-instance wide-band families are not missed, but the
    # spread metric itself uses population std which is defined for n>=2).
    bed_spread_min_instances: int = 2
    # Copy-count-vs-length-class expectation table.  mdl-repeat tends to
    # under-instance longer elements (a full-length element leaves fewer
    # near-identical full copies than a short high-copy SINE/MITE), so a long
    # element with very few BED instances is a "missing members" signal.  Each
    # entry: (min_length_inclusive, expected_min_copies).  A family whose
    # BED instance count is below the expectation for its length class is
    # flagged.  Evaluated as: walk the table descending by length floor, take
    # the first class whose floor <= family length, compare n_instances to its
    # expected_min_copies.  Defaults are conservative (favor NOT flagging) and
    # consistent with Phase 2's length-dependent masking_copies_* gates:
    #   >=1000 bp  : expect >=4 copies  (long elements; 2-3 is suspicious)
    #   500-999 bp : expect >=3 copies
    #   200-499 bp : expect >=3 copies
    #   <200 bp    : expect >=2 copies  (short SINE/MITE class; low bar)
    length_class_copy_expectation: tuple = (
        (1000, 4),
        (500, 3),
        (200, 3),
        (0, 2),
    )

    # N4 (v2): selective recall — find the divergent tail of incomplete families
    # ────────────────────────────────────────────────────────────────────
    # REFINE_STRATEGY_DESIGN_v2.md §2.5 (two-tier recall), §2.6 (TE-looker
    # boundary / identity floor), §7 (M5 promote row), §10 N4.  Recall fires
    # ONLY for families a cheap, SELECTIVE completeness gate flags (see
    # phase1_completeness.recall_eligible): on real chr4 the gate
    # (under_instanced AND a specific signal) fires for 20.9% of families vs the
    # 85.2% that `under_instanced` alone would — the bottleneck the gate exists to
    # avoid.  The recalled copies (R=150's 12-29%-divergent partials) flow back
    # into N2 clustering so true subfamily structure emerges.
    enable_selective_recall: bool = True   # master switch; False degrades to N3 behaviour
    # Tier-2 recall pident floor.  PLACEHOLDER 0.65 (== 65% identity); the v2
    # productive band is ~60-70% (§2.6: "keeps the productive ≥~60-70% band where
    # alignment is reliable and the consensus stays coherent").  This is the
    # refine↔TE-looker handoff identity floor — recall recruits members homologous
    # to a KNOWN mdl-repeat seed down to this floor; below it is TE-looker's
    # de-novo twilight (D_cc 45-60%).  TUNABLE: retune against a TE-looker
    # divergence benchmark in N6 (the metric harness reports the floor so the
    # handoff is auditable).
    recall_identity_floor: float = 0.65
    # Global cap on the NUMBER of families that may trigger a Tier-2 recall in one
    # run.  A hard ceiling on total recall compute independent of how many families
    # the gate flags: once this many families have recalled, further eligible
    # families KEEP their BED copies and are flagged completeness_verified=False
    # (honest, never silently dropped — §2.5, §6).  Default 5000 comfortably covers
    # chr4's ~610 co-fire-eligible families; lower it to bound a whole-genome run.
    recall_global_family_budget: int = 5000
    # Per-family wall-clock budget (seconds) for the Tier-2 blastn recruit.  A
    # family whose recruit blastn exceeds this is treated as a recall failure:
    # the BED copies are kept and the family flagged completeness_verified=False.
    # Bounds the tail cost of a pathological giant-element family without faking a
    # recall result.  Passed through to recruit_by_blastn's subprocess timeout.
    recall_per_family_wall_s: int = 300
    # Genome k-mer sketch pre-screen (§2.5 Tier 1).  DISABLED by default on the
    # basis of a real-chr4 empirical test (run during N4, recorded in the N4
    # report): a sourmash MinHash sketch of the 472 bp R=150 seed against the
    # 18.5 Mbp chr4 genome sketch returns ZERO matches under `sourmash gather` and
    # `sourmash search --containment` even at scaled=1 (all k-mers kept) — the
    # surviving shared 21-mers per divergent partial are too few (median ~1-13% of
    # the 446 seed k-mers, many at 1-3%) for a genome-scale containment estimate to
    # clear any usable threshold.  This is the N2 k-mer-saturation result
    # generalized: at 75-90% identity the sketch is blind to the divergent tail it
    # would need to gate on, so it cannot serve as the cheap Tier-1 confirm.  The
    # SELECTIVE gate is instead the N1 specific-signal co-fire
    # (recall_eligible).  Flip to True only if a future benchmark shows a sketch
    # configuration that detects the tail.
    enable_genome_sketch_prescreen: bool = False

    # N5 (v2): conservative false-positive pruning — the pseudo-family discriminant
    # ────────────────────────────────────────────────────────────────────
    # REFINE_STRATEGY_DESIGN_v2.md §4 (all), §4.1 pseudo-family row, §4.2, §4.3, §8.2,
    # §10 N5.  A family is dropped ONLY when, after recall (N4) and clustering (N2/N3),
    # it forms NO qualifying cluster AND its copies are POSITIVELY incoherent
    # (mutually un-homologous).  The drop NEVER rests on low copy count — chr4 is 84%
    # low-copy REAL families (median 3 copies), so low copy count is the NORM, not a
    # noise signal (§4.1/§4.2).  Disabling the master switch degrades exactly to N4.
    enable_conservative_fp_prune: bool = True   # master switch; False degrades to N4
    # Minimum assembled copies before an incoherence call is TRUSTED (§4.3 "when in
    # doubt, keep").  A family with fewer copies than this is KEPT unconditionally —
    # there is too little pairwise evidence to declare a grab-bag, and the chr4 reality
    # (84% of REAL families are low-copy) means a low count must never by itself
    # condemn a family.  Set to min_copies_for_msa (5): below the MSA floor a family
    # already cannot form a qualifying cluster for benign reasons, so we refuse to also
    # call it pseudo.  This is the single most important guard against误删 real
    # low-copy families.
    pseudo_family_min_copies_for_call: int = 5
    # Median pairwise BLASTN local IDENTITY (1 - distance) over the copy set, BELOW
    # which the copies count as "mutually unrelated" (§4.1/§4.2 incoherence).  0.45 is
    # deliberately FAR below the subfamily homogeneity floor (0.6): a copy set whose
    # MEDIAN pair is below 45% identity is not a loose-but-real subfamily, it is a
    # grab-bag.  Conservative by design — a real divergent family's copies still share
    # well above this (the chr4 measured emitted subfamilies sit at 0.79–0.92, and even
    # the divergent recall tail targets the ≥~60-70% band, §2.6), so a true family is
    # never close to this floor.  Median (not mean) so one chance high-identity pair in
    # a grab-bag cannot rescue it from the drop, and one low-identity outlier in a real
    # family cannot condemn it.
    pseudo_family_max_intra_identity: float = 0.45
    # Fraction of copy PAIRS that must share NO significant local alignment (no
    # gate-passing HSP, distance == 1.0) for the set to count as "unconnected".  The
    # drop requires the homologous-pair fraction to be BELOW this, i.e. the copies are
    # mostly mutually un-alignable.  0.25 == "fewer than a quarter of copy pairs share
    # any homology" — a real family (even a divergent one) has nearly every pair
    # connected through the all-vs-all BLAST; only a genuine grab-bag drops this low.
    # This is the "is the copy set even CONNECTED?" axis of §4.2; it must hold TOGETHER
    # with the low-median-identity axis (both conjuncts) for a drop, so a family that is
    # connected-but-low-identity (a fast-evolving but real family) is KEPT.
    pseudo_family_min_homologous_pair_frac: float = 0.25
    # ── Seed-coherence protection (§4.3 "存疑即留 / when in doubt, keep") ──────
    # A family's purported copies can be mutually incoherent (each pairwise-unrelated)
    # while the family STILL has a real consensus that ≥1 copy clearly belongs to — the
    # canonical case is a star-topology divergent family (every copy matches the seed but
    # not each other) or a real family carrying one anomalous instance (e.g. chr4 R=286:
    # a 2579 bp consensus whose copies are 2578/2578/2578/41659/7307 bp, mutually
    # un-homologous, yet one copy matches the consensus at 99.8% identity).  The copy×copy
    # incoherence conjuncts ALONE judge such a family pseudo and would误删 it.  This guard
    # adds a NECESSARY pre-condition for the drop: the family's own seed/consensus
    # (rec['sequence']) is blastn'd against its copies, and the family is only allowed to
    # be called pseudo when the seed itself recruits FEWER than seed_coherence_min copies
    # at >= seed_coherence_min_pident identity over >= seed_coherence_min_coverage of the
    # seed.  seed_coherence_min = 1: a single coherent seed→copy match is positive evidence
    # the family has real signal, so it is KEPT (§4.3).  A truly spurious family — whose
    # seed matches NONE of its copies AND whose copies are mutually incoherent — still
    # drops.  The seed→copies blastn is per-family on a tiny copy set (no genome-scale
    # cost) and is computed ONLY when the family has already passed every incoherence
    # conjunct and is about to be dropped (short-circuit), so it adds no cost to the
    # overwhelming KEEP majority.
    enable_seed_coherence_protection: bool = True   # part of the N5 master switch family
    seed_coherence_min: int = 1            # >=1 seed→copy match ⇒ keep (real signal exists)
    seed_coherence_min_pident: float = 80.0   # %identity floor for a seed→copy match to count
    seed_coherence_min_coverage: float = 0.50  # min fraction of the SEED covered by the HSP

    # Phase 1 (rewrite M4): orchestrator sharding + post-refine merge + fallback
    shard_size: int = 500                 # families per parallel shard task
    post_refine_merge_identity: float = 0.95  # collapse boundary-convergent dups post-refine
    enable_fallback_recruitment: bool = True   # per-family blastn for under-instanced families
    enable_divergent_blast_recruitment: bool = False  # 45-75% band is TE-looker's remit

    # Phase 1: Consensus polishing
    blastn_max_targets: int = 50
    blastn_evalue: float = 1e-5
    blastn_batch_size: int = 500
    min_recruit_pident: float = 70.0    # Min %identity for copy recruitment
                                        # Aligned with mdl-repeat max-divergence=0.30
                                        # (instances accepted at identity >= 70%);
                                        # stricter values would drop the most
                                        # divergent — and most informative — copies
                                        # from MAFFT consensus polish.
    t1_max_hits: int = 20
    t2_max_hits: int = 50
    t3_max_hits: int = 10
    chimera_cv_threshold: float = 1.5
    chimera_low_coverage_ratio: float = 0.3
    min_chimera_fragment: int = 50
    chimera_max_depth: int = 2            # Max recursive chimera splitting depth
    chimera_min_hits: int = 5             # Min BLASTN hits for chimera detection
    consensus_min_occupancy: float = 0.3
    consensus_occupancy_fail_ratio: float = 0.5

    # Phase 1: Tier thresholds (fallback for small inputs)
    small_input_threshold: int = 100
    tier_absolute_high: float = 20.0   # mdl_per_copy for T1 (small input)
    tier_absolute_low: float = 5.0     # mdl_per_copy for T3 (small input)

    # Phase 1: T2 downgrade threshold
    t2_downgrade_count: int = 10000     # If T2 > this, downgrade short seqs
    t2_downgrade_length: int = 200      # Length threshold for downgrade

    # Performance
    threads: int = 8
    max_retries: int = 3
    retry_backoff: int = 2

    # External tools
    blastn_exe: str = "blastn"
    makeblastdb_exe: str = "makeblastdb"
    mafft_exe: str = "mafft"
    cdhit_exe: str = "cd-hit-est"
    blastx_exe: str = "blastx"
    samtools_exe: str = "samtools"
    repeatmasker_exe: str = "RepeatMasker"

    # Large-genome SAMPLED mode ───────────────────────────────────────────
    # When the genome exceeds the large-genome threshold, build_mdl runs mdl-repeat and
    # all of Refiner_mdl's rmblastn against a SAMPLE of the genome (passed as --genome /
    # config.genome_file), and passes the COMPLETE genome here as --genome-full. The only
    # consumer of the full genome is phase2b (genome-wide copy count + final hard floor);
    # every other rmblastn surface (copy gate, Phase-1 per-family, recall, fallback)
    # follows config.genome_file = the sample automatically. Empty -> not sampled, the
    # ≤2Gb path is byte-identical to before.
    genome_full_file: str = ""
    # Defer the copy-count hard floor / short-low-copy drop in the copy gate (set with
    # genome_full_file). Sample-based copy counts undercount real families, so the copy
    # gate only recruits (for Phase-1 ordering) + rescues nothing; the REAL floor is
    # applied by phase2b on genome-wide counts. False -> copy gate drops as today.
    defer_copy_floor: bool = False

    # Runtime (set by main.py)
    genome_blast_db: str = ""
    enable_masking: bool = False
    keep_temp: bool = False

    def save(self, filepath: str):
        with open(filepath, 'w') as f:
            json.dump(self.__dict__, f, indent=2)

    @classmethod
    def load(cls, filepath: str):
        with open(filepath, 'r') as f:
            data = json.load(f)
        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in data.items() if k in valid_fields}
        return cls(**filtered)
