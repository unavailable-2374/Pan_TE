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

### Option A — conda / bioconda

> **Status (in progress):** the tool dependencies (`look4ltrs`, `mdl-repeat`, `te-looker`) are
> submitted to bioconda and passing CI — [bioconda-recipes #66613](https://github.com/bioconda/bioconda-recipes/pull/66613),
> awaiting maintainer review. The `pan_te` recipe follows once those merge. The command below
> works once all four are merged into the bioconda channel:

    conda create -n pan_te -c conda-forge -c bioconda pan_te
    conda activate pan_te

Until then, use Option B.

### Option B — from source

    git clone --recursive https://github.com/unavailable-2374/Pan_TE.git
    cd Pan_TE
    conda env create -f environment.yml        # creates the `pan_te` env (or: mamba env create)
    conda activate pan_te
    ./build_tools.sh                            # builds look4ltrs, FastGA, mdl-repeat, dtr into the env
    echo 'export PATH=$PWD/bin:$PATH' >> ~/.bashrc

`build_tools.sh` resolves each custom tool from its repo-local submodule, a `PAN_TE_*_SRC`
environment override, or a clone. The te-looker source lives at
<https://github.com/unavailable-2374/te-looker> — set `PAN_TE_TE_LOOKER_SRC=/path/to/te-looker`
if it is not already checked out.

**Classification (optional)** runs in its own isolated environment because of pinned ML deps:

    conda env create -f environment-classifyte.yml   # creates `ClassifyTE_env`

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

    # Pan-genome, graph (GFA from minigraph-cactus or PGGB)
    Pan_TE --gfa pangenome.gfa --threads 80

    # Enrichment (reference genome + pan-genome insertions)
    Pan_TE --genome ref.fasta --paf all_vs_all.paf --paf-fasta assemblies.fa --threads 80
    Pan_TE --genome ref.fasta --gfa pangenome.gfa --threads 80

### Key options

| Option | Description |
|--------|-------------|
| `--genome <fa>` | Genome FASTA. Required unless `--paf` or `--gfa` is given. |
| `--model-dir <dir>` | Path to ClassifyTE (enables classification). |
| `--paf <file>` / `--paf-fasta <fa...>` | PAF from all-vs-all alignment + its assembly FASTA(s). |
| `--paf-mapping <tsv>` | Optional chromosome mapping (seqname, assembly_id, chr_id). |
| `--gfa <file>` / `--gfa-ref <name>` | Pangenome graph (rGFA / smoothxg); reference path auto-detected if omitted. |
| `--vcf-dir <dir>` | Directory of structural-variant VCFs to fold in. |
| `--te-looker-bin <path>` / `--te-looker-extra-args "<args>"` | Override the `dtr` binary / pass extra args to `dtr run`. |
| `--out-dir <dir>` | Work directory (default: current directory). |
| `--threads <int>` | Threads (default 4; multiples of 4 recommended). |
| `-M <int>` | Memory limit in MB (default 0 = unlimited). |
| `--fragment-size <int>` | Fragment length (default 40000). |

## Output

- `combine/raw_te_library.fa` — deduplicated de novo TE library (all three tracks merged).
- `combine/TEs.fa` — the classified TE library (when `--model-dir` is provided).

## Notes

- The pipeline is checkpointed and resumable: re-running skips completed steps.
- Packaging recipes and the bioconda submission notes live under [`recipes/`](recipes/).
