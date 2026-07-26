# Pan_TE — Zhou Lab @ AGIS

Pan_TE is a **pan-genome transposable-element (TE) annotation pipeline**. It combines three
complementary de novo detection tracks across the divergence spectrum and unifies them into a
single classified TE library:

- **Look4LTRs** — LTR retrotransposon detection (active / young elements)
- **mdl-repeat** — broad-spectrum de novo repeat discovery (Minimum Description Length)
- **TE-looker (`dtr`)** — divergent, "twilight-zone" TEs (older elements other tools miss)

It runs on a single genome, on a pan-genome (reference-free **PAF** or a **GFA** graph), or in
enrichment mode (a reference genome augmented with pan-genome insertions), and optionally
classifies the result with RepeatClassifier + ClassifyTE.

> **Platform:** Linux (x86-64 / aarch64). look4ltrs requires GCC + OpenMP and is not built for macOS.

## Installation

### Option A - conda / bioconda

Once the published `pan_te` package is available from the configured channels, installation is
self-contained: Conda supplies `look4ltrs`, `mdl-repeat`, `te-looker` and their runtime helper
scripts. No Git submodules or local compilation are required.

    conda create -n pan_te -c conda-forge -c bioconda pan_te
    conda activate pan_te

For source development, use Option B.

### Option B — from source

    git clone --recursive https://github.com/unavailable-2374/Pan_TE.git
    cd Pan_TE
    conda env create -f environment.yml        # creates the `pan_te` env (or: mamba env create)
    conda activate pan_te
    ./build_tools.sh                            # development fallback: builds custom tools locally
    echo 'export PATH=$PWD/bin:$PATH' >> ~/.bashrc

For development only, `build_tools.sh` resolves each custom tool from a `PAN_TE_*_SRC`
environment override, an optional repo-local submodule, or a clone. The te-looker source lives at
<https://github.com/unavailable-2374/te-looker> — set `PAN_TE_TE_LOOKER_SRC=/path/to/te-looker`
if it is not already checked out.

The LTR consensus helper is resolved in this order: an explicit
`PAN_TE_LTR_LIBRARY_BUILDER` path, the Conda `look4ltrs` package, then the legacy
`submodule/Look4LTRs` development checkout.

**Classification (optional)** runs in its own isolated environment because of pinned ML deps:

    conda env create -f environment-classifyte.yml   # creates `ClassifyTE_env`

`run_Classifier` invokes that environment's Python directly to avoid nested Conda
activation conflicts. If the environment is installed outside Conda's standard env
directories, set `CLASSIFYTE_ENV_PREFIX=/absolute/path/to/the/environment`.

### Data assets (required for real runs)

The Dfam libraries (for RepeatMasker) and the ClassifyTE models cannot be redistributed via
conda. Fetch them after install with the bundled helper:

    pan_te-setup-data --dfam-list                    # how to pick a Dfam partition
    pan_te-setup-data \
        --classifyte-dir ~/Pan_TE_data/ClassifyTE \
        --dfam-url https://www.dfam.org/releases/Dfam_3.9/families/dfam39_full.0.h5.gz

(From a source checkout, call `bin/pan_te-setup-data`. Verify the current Dfam release/filenames
at <https://www.dfam.org/releases/>.)

## Usage

    # Single genome
    Pan_TE --genome genome.fasta --threads 80 --model-dir /path/to/ClassifyTE

    # Pan-genome, reference-free (all-vs-all PAF)
    Pan_TE --paf all_vs_all.paf --paf-fasta asm1.fna asm2.fna ... --threads 80

    # PGGB graph (embedded P/W paths)
    Pan_TE --gfa pggb.gfa --threads 80

    # rGFA graph: primary minigraph GAF + the query assemblies are required
    Pan_TE --gfa pangenome.rgfa --gfa-gaf assemblies_to_graph.gaf \
        --gfa-gaf-fasta assemblies.fa --threads 80

    # Also annotate every GFA node and produce a per-sample presence matrix
    Pan_TE --gfa pangenome.gfa --anno --threads 80 --model-dir /path/to/ClassifyTE

    # Enrichment (reference genome + pan-genome insertions)
    Pan_TE --genome ref.fasta --paf all_vs_all.paf --paf-fasta assemblies.fa --threads 80
    Pan_TE --genome ref.fasta --gfa pangenome.gfa --threads 80

    # Add explicit insertion alleles from plain or gzip-compressed SV VCFs
    Pan_TE --genome ref.fasta --vcf-dir sv_vcfs --threads 80

### Key options

| Option | Description |
|--------|-------------|
| `--genome <fa>` | Genome FASTA. Required unless `--paf` or `--gfa` is given. |
| `--model-dir <dir>` | Path to ClassifyTE (enables classification). |
| `--paf <file>` / `--paf-fasta <fa...>` | PAF from all-vs-all alignment + its assembly FASTA(s). FASTA identifiers must be unique across assemblies. |
| `--paf-mapping <tsv>` | Optional chromosome mapping (seqname, assembly_id, chr_id). |
| `--gfa <file>` / `--gfa-ref <name>` | Pangenome graph (rGFA / PGGB smoothxg); reference sample/path auto-detected if omitted. rGFA default allele mode requires `--gfa-gaf`. |
| `--gfa-mode {alleles,legacy}` | `alleles` (default) uses PGGB paths or rGFA+GAF to extract locus-aware alternates; `legacy` retains the historical topology walker and is not sample-aware. |
| `--gfa-gaf <file>` / `--gfa-gaf-fasta <fa>` | Primary minigraph GAF and its query FASTA. Required for default rGFA allele extraction; FASTA names/lengths are validated against GAF. |
| `--gfa-gaf-min-mapq <n>` / `--gfa-gaf-min-identity <f>` | Quality thresholds for rGFA candidate and presence calls (defaults: 20 and 0.90). |
| `--gfa-window-overlap <bp>` | Overlap for alleles longer than 30 kb (default 1000). Terminal windows are shifted so no sub-100-bp tail is lost. |
| `--gfa-spacer-length <bp>` | Artificial N spacer between packed candidates (default 0, so every allele is an independent FASTA record). |
| `--gfa-clean-intermediates` | Remove candidate FASTA intermediates after writing the provenance manifest and summary. |
| `--anno` | After building the TE library, run sensitive RepeatMasker (`-s`) on GFA P/W path contexts and write node-level TE coverage plus per-sample presence. rGFA without P/W uses `--gfa-gaf`: node columns are conservative `1/NA`, while locus calls are `0/1/NA`. Default simple/low-complexity hits absent from `combine/TEs.fa` are excluded. Requires `--gfa`. |
| `--vcf-dir <dir>` | Directory of `.vcf`/`.vcf.gz` structural variants to fold into `--genome`. Explicit ALT insertions and symbolic insertions with `SVINSSEQ`, `INSSEQ`, `SEQ`, or `CONSENSUS` are supported. |
| `--te-looker-bin <path>` / `--te-looker-extra-args "<args>"` | Override the `dtr` binary / pass extra args to `dtr run`. |
| `--te-looker-w <int>` | te-discover minimizer window, default 1 (exact: every k-mer). Values >1 emit only ~`2/(w+1)` of k-mer occurrences, cutting memory roughly proportionally — but they **change the seed set**, so this is opt-in only and `--min-count` is re-derived to match. |
| `--check-env` | Run the environment preflight (tools, Perl/Python modules, data assets, memory) and exit without running the pipeline. |
| `--out-dir <dir>` | Work directory (default: current directory). |
| `--threads <int>` | Threads (default 4; multiples of 4 recommended). |
| `-M <int>` | Memory limit in MB (default 0 = unlimited). |
| `--fragment-size <int>` | Fragment length (default 40000). |

## Output

- `combine/raw_te_library.fa` — deduplicated de novo TE library (all three tracks merged).
- `combine/TEs.fa` — the classified TE library (when `--model-dir` is provided).
- `genome/gfa_enrichment/gfa_allele_manifest.tsv` — one row per locus-aware allele occurrence, with locus, coordinates, type, exact representative, sample support, QC status, and rGFA+GAF provenance when applicable.
- `genome/gfa_enrichment/gfa_rgfa_locus_presence.tsv` — rGFA candidate-locus `0/1/NA` matrix. The SR=0 reference sample is fixed to `0`; other `0` calls require one quality GAF record spanning both candidate anchors; unresolved/overlapping calls remain `NA`.
- `genome/gfa_enrichment/gfa_unique_locus_alleles.fa` — exact locus-aware representative sequences before QC (kept by default).
- `genome/gfa_enrichment/gfa_processing_stats.json` — graph, allele, deduplication, QC, and artificial-spacer statistics.
- `anno/gfa_node_annotations.tsv` — one row per GFA node, including TE-covered/uncovered bp, coverage fraction, annotation status, and sample columns. P/W graphs use 0/1; rGFA+GAF uses conservative 1/NA node evidence.
- `anno/gfa_rgfa_node_presence.tsv` — long-form complete-node GAF traversal evidence for rGFA.
- `anno/gfa_node_te_hits.tsv` — every RepeatMasker hit projected to 0-based, half-open node coordinates.
- `anno/gfa_path_node_map.tsv` — auditable mapping from RepeatMasker path contexts back to GFA nodes.

Reference-path validation against an independent RepeatMasker run can be performed with:

    python bin/compare_gfa_reference_annotation.py \
        --gfa pangenome.gfa --reference-fasta reference.fa \
        --node-hits anno/gfa_node_te_hits.tsv \
        --direct-out reference.fa.out --repeat-library combine/TEs.fa \
        --out-dir anno/reference_comparison

The comparison verifies every `SR:i:0` node against its `SN`/`SO` reference
substring (including orientation), filters the direct run to custom-library
families, and reports coverage differences, class/family agreement, and direct-hit fragmentation.

## PGGB scaling model

In the default `alleles` mode, reference occurrences are used as ordered anchors and each
alternate P path is consumed one at a time. Exact alleles are aggregated online by
`(graph locus, sequence SHA-256)`; occurrence metadata is kept in a temporary SQLite store
and emitted to the manifest at the end. Memory therefore follows the embedded graph
segments, selected reference paths, one alternate path, and unique locus alleles—not all
cohort path steps and duplicate allele sequences. Identical sequences at different loci
remain separate, while copies of the same allele at one locus contribute sample and
occurrence support. Rearrangements use a maximum collinear occurrence chain, and missing
contigs, missing segment sequence, or unsupported non-zero P overlaps stop the run.

## Environment preflight

Every run verifies its environment before doing any work, and `--check-env` runs that check
on its own:

    Pan_TE --check-env --genome genome.fasta

It resolves every external tool the requested run will invoke, the Perl and Python modules the
stages import, the ClassifyTE assets when `--model-dir` is given, and whether the te-looker
step can fit in the memory this job actually has. Problems are reported together, each naming
its fix, and nothing is created until the checks pass.

The memory check takes the smallest ceiling it can find — `/proc/meminfo`, cgroup v1/v2 limits,
`SLURM_MEM_PER_NODE` / `SLURM_MEM_PER_CPU`, and `-M` — because a scheduler allocation is
usually far below the node's physical RAM. te-discover holds one 16-byte entry per genomic
position at the default `w=1`, so a 3.1 Gb genome needs roughly 56 GB (about 110 GB against
te-looker ≤ 0.3.0, which held two copies of that buffer). If that does not fit, Pan_TE stops
and tells you the requirement, the detected ceiling and its source, rather than quietly
subsampling: `w` is never chosen from available memory, because the resulting library would
then depend on which machine happened to run it.

## Notes

- The pipeline is checkpointed and resumable. GFA-backed `genome.ok` checkpoints include an input/configuration fingerprint; stale or missing fingerprints are rejected instead of being silently reused.
- `PartialTE` means that only part of a node is TE; `te_covered_bp`, `te_fraction`, and `te_uncovered_bp` preserve the boundary instead of labeling the whole node as TE.
- PAF sequence names and lengths must match `--paf-fasta`; malformed or failed PAF/VCF enrichment stops the run instead of silently omitting data.
- Symbolic VCF alleles without an embedded sequence are skipped because no nucleotide sequence can be added to TE discovery.
- Packaging recipes and the bioconda submission notes live under [`recipes/`](recipes/).
