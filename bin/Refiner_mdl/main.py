#!/usr/bin/env python3
"""
Refiner_mdl: Refinement pipeline for mdl-repeat output.

Architecture: 3 phases
  Phase 0: Metadata parsing + hard filtering + tiering + fragment assembly + dedup
  Phase 1: Tiered consensus polishing (T1/T2/T3) + chimera detection
  Phase 2: Final QC + library split (masking + analysis)
"""

import logging
import argparse
import os
import sys
import shutil
import multiprocessing as mp
from pathlib import Path
from typing import Dict, Any

try:
    mp.set_start_method('spawn', force=True)
except RuntimeError:
    pass

# Add parent dirs to path for imports
_bin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_refiner_dir = os.path.join(_bin_dir, 'Refiner')
for _d in [_bin_dir, _refiner_dir]:
    if _d not in sys.path:
        sys.path.insert(0, _d)

from Refiner_mdl.config import RefinerMdlConfig
from Refiner_mdl.phase0_triage import run_phase0
from Refiner_mdl.phase1_consensus import run_phase1
from Refiner_mdl.phase2_library_split import run_phase2

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class RefinerMdlPipeline:
    """Orchestrator for the Refiner_mdl pipeline."""

    def __init__(self, config: RefinerMdlConfig):
        self.config = config

        # Create output dirs
        for d in [config.output_dir, config.temp_dir, config.checkpoint_dir]:
            Path(d).mkdir(parents=True, exist_ok=True)

        # Add file handler for logging
        log_file = os.path.join(config.output_dir, 'refiner_mdl.log')
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.INFO)
        fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logging.getLogger().addHandler(fh)

        # Genome BLAST DB is NOT built here. It is built lazily, once, inside Phase 1
        # (run_phase1 → ensure_genome_blast_db) and ONLY when under-instanced families
        # actually need the bounded blastn fallback. The old unconditional build (with
        # its 600 s silent-failure trap) sat on the critical path of every run — even
        # ones where the common BED-seeded path never touches the DB — so it is gone.

        logger.info("Refiner_mdl pipeline initialized")
        logger.info(f"  Input:   {config.input_file}")
        logger.info(f"  Genome:  {config.genome_file}")
        logger.info(f"  Output:  {config.output_dir}")
        logger.info(f"  Threads: {config.threads}")

    def ensure_genome_blast_db(self):
        """Lazy genome BLAST DB entry point (delegates to the canonical builder).

        Not called on the common path: Phase 1 triggers the build itself in the parent
        process before spawning shard workers, only when an under-instanced family
        needs the fallback. Exposed here as an explicit hook; on failure it sets
        config.genome_blast_db = "" and logs an error (under-instanced families then
        keep their mdl-repeat seeds — never a faked success)."""
        from phase1_fallback import ensure_genome_blast_db as _ensure
        return _ensure(self.config)

    def _load_checkpoint(self, name: str):
        """Load checkpoint if exists."""
        import pickle
        cp_file = os.path.join(self.config.checkpoint_dir, f'{name}.pkl')
        if os.path.exists(cp_file):
            try:
                with open(cp_file, 'rb') as f:
                    data = pickle.load(f)
                logger.info(f"Restored checkpoint: {name}")
                return data
            except Exception as e:
                logger.warning(f"Failed to load checkpoint {name}: {e}")
        return None

    def _save_checkpoint(self, name: str, data):
        """Save checkpoint."""
        import pickle
        cp_file = os.path.join(self.config.checkpoint_dir, f'{name}.pkl')
        try:
            with open(cp_file, 'wb') as f:
                pickle.dump(data, f)
            logger.info(f"Saved checkpoint: {name}")
        except Exception as e:
            logger.warning(f"Failed to save checkpoint {name}: {e}")

    def run(self) -> Dict[str, Any]:
        """Execute the full pipeline."""
        logger.info("=" * 60)
        logger.info("Starting Refiner_mdl Pipeline")
        logger.info("=" * 60)

        # Phase 0
        phase0_output = self._load_checkpoint('phase0_complete')
        if phase0_output is None:
            phase0_output = run_phase0(self.config)
            self._save_checkpoint('phase0_complete', phase0_output)

        # Copy gate (recruit genomic copies + recurrence hard-floor + TE-structure
        # rescue) — runs BEFORE Phase 1 so the floor applies to whole-family copy
        # counts and only survivors get the expensive representativeness work.
        gated_output = self._load_checkpoint('copy_gate_complete')
        if gated_output is None:
            from Refiner_mdl.phase05_copy_gate import run_copy_gate
            gated_output = run_copy_gate(phase0_output, self.config)
            self._save_checkpoint('copy_gate_complete', gated_output)

        # Phase 1
        phase1_output = self._load_checkpoint('phase1_complete')
        if phase1_output is None:
            phase1_output = run_phase1(gated_output, self.config)
            phase1_output['phase0_stats'] = phase0_output.get('stats', {})
            phase1_output['copy_gate_stats'] = gated_output.get('copy_gate_stats', {})
            self._save_checkpoint('phase1_complete', phase1_output)

        # Phase 2
        phase2_output = self._load_checkpoint('phase2_complete')
        if phase2_output is None:
            phase2_output = run_phase2(phase1_output, self.config)
            self._save_checkpoint('phase2_complete', phase2_output)

        # Phase 2b (large-genome SAMPLED mode only): genome-wide copy count + the
        # deferred recurrence floor, applied once on the final library against the
        # COMPLETE genome. No-op when --genome-full is unset (≤2Gb path unchanged).
        final_output = phase2_output
        if getattr(self.config, 'genome_full_file', ''):
            phase2b_output = self._load_checkpoint('phase2b_genome_count_complete')
            if phase2b_output is None:
                from Refiner_mdl.phase2b_genome_count import run_genome_count
                phase2b_output = run_genome_count(phase2_output, self.config)
                self._save_checkpoint('phase2b_genome_count_complete', phase2b_output)
            final_output = phase2b_output

        # Cleanup temp
        if not self.config.keep_temp:
            self._cleanup_temp()

        logger.info("=" * 60)
        logger.info("Refiner_mdl Pipeline Complete")
        logger.info(f"  Masking library:  {final_output.get('masking_library', 'N/A')}")
        logger.info(f"  Analysis library: {final_output.get('analysis_library', 'N/A')}")
        logger.info("=" * 60)

        return final_output

    def _cleanup_temp(self):
        """Remove temporary files."""
        temp_dir = Path(self.config.temp_dir)
        if temp_dir.exists() and str(temp_dir) != self.config.output_dir:
            try:
                shutil.rmtree(temp_dir)
                logger.info("Temporary files cleaned up")
            except Exception as e:
                logger.warning(f"Cleanup failed: {e}")


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Refiner_mdl: Refinement pipeline for mdl-repeat output'
    )
    parser.add_argument('-r', '--input', required=True,
                        help='mdl-repeat output FASTA file')
    parser.add_argument('-g', '--genome', required=True,
                        help='Reference genome FASTA. In large-genome SAMPLED mode this is '
                             'a SAMPLE of the genome (all rmblastn aligns against it); pass '
                             'the complete genome via --genome-full for the final count.')
    parser.add_argument('--genome-full', default='',
                        help='Complete genome FASTA for the large-genome SAMPLED mode. When '
                             'set, --genome is treated as a sample: the copy-count hard floor '
                             'is DEFERRED past the per-family work and applied once at the end '
                             '(phase2b) against this full genome. Genome-size copy tiers use '
                             'this size. Leave empty (default) for the standard full-genome run.')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=8,
                        help='Number of threads (default: 8)')
    parser.add_argument('--bed', default='',
                        help='mdl-repeat instances BED file')
    parser.add_argument('--stats', default='',
                        help='mdl-repeat stats file')
    parser.add_argument('--keep-temp', action='store_true',
                        help='Keep temporary files')
    parser.add_argument('--enable-copy-recruit', action='store_true',
                        help='Recruit real genomic copy number via rmblastn (recommended; '
                             'replaces mdl copies for the recurrence filter)')
    parser.add_argument('--hard-min-copies', type=int, default=0,
                        help='Drop families below this genomic copy count (0=off; '
                             '5 aligns with the TE-looker handoff)')
    parser.add_argument('--copy-recruit-cov', type=float, default=0.5,
                        help='Min consensus coverage for an rmblastn instance to count '
                             'as one genomic copy (default 0.5)')
    parser.add_argument('--copy-recruit-divergence', type=int, default=25,
                        help='rmblastn matrix divergence level p (14/18/20/25); lower = '
                             'stricter recruitment, fewer gene-paralog false copies')
    parser.add_argument('--enable-rescue', action='store_true',
                        help='Rescue dropped-but-TE-structured consensi (TEsorter HMM + '
                             'LTR/TIR/known-TE); rescue-only, never drops')
    parser.add_argument('--cellular-protein-db', default='',
                        help='BLAST protein DB of cellular proteins (TE removed); enables '
                             'gene-exclusion of host-gene-derived consensi')
    parser.add_argument('--organelle-ref', default='',
                        help='FASTA of chloroplast+mito genomes for NUMT/NUPT removal')
    parser.add_argument('--enable-nonte', action='store_true',
                        help='Separate organelle / rRNA / satellite into non-TE tracks')
    parser.add_argument('--barrnap-exe', default='barrnap',
                        help='barrnap executable for rRNA detection (non-TE separation)')
    parser.add_argument('--enable-rmblastn-copies', action='store_true',
                        help='Recruit Phase-1 copies via rmblastn (instead of mdl BED) — '
                             'diverged members make a more representative consensus')
    parser.add_argument('--enable-seed-chimera', action='store_true',
                        help='De-nest chimeric SEED consensi (N7) via rmblastn query-'
                             'coverage breakpoints (needs --enable-rmblastn-copies)')
    parser.add_argument('--enable-flank-trim', action='store_true',
                        help='Trim terminal host-gene flanks off TE consensi (protein '
                             'homology; needs --cellular-protein-db)')
    parser.add_argument('--enable-structural-trim', action='store_true',
                        help='Per-copy TSD-anchored boundary trim (family-consensus) so '
                             'the consensus does not over-extend into shared host flanks')

    args = parser.parse_args()

    # Genome size for the genome-size-adaptive low-copy gate. In sampled mode the COPY
    # TIERS must reflect the FULL genome (a 200 Mb sample of a 16 Gb genome must still use
    # the XXL tier), so prefer --genome-full when set.
    genome_full = args.genome_full or ''
    size_src = genome_full if (genome_full and os.path.exists(genome_full)) else args.genome
    try:
        genome_size_bp = os.path.getsize(size_src)
    except OSError:
        genome_size_bp = 0
    sampled_mode = bool(genome_full and os.path.exists(genome_full)
                        and os.path.abspath(genome_full) != os.path.abspath(args.genome))

    config = RefinerMdlConfig(
        input_file=args.input,
        genome_file=args.genome,
        output_dir=args.output,
        temp_dir=os.path.join(args.output, 'temp_work'),
        checkpoint_dir=os.path.join(args.output, 'checkpoints'),
        bed_file=args.bed,
        stats_file=args.stats,
        threads=args.threads,
        enable_masking=False,
        genome_size_bp=genome_size_bp,
        enable_copy_recruit=args.enable_copy_recruit,
        copy_recruit_cov=args.copy_recruit_cov,
        copy_recruit_divergence=args.copy_recruit_divergence,
        hard_min_copies=args.hard_min_copies,
        enable_lowcopy_filter=(args.enable_copy_recruit or args.hard_min_copies > 0),
        enable_rescue=args.enable_rescue,
        cellular_protein_db=args.cellular_protein_db,
        enable_gene_exclusion=bool(args.cellular_protein_db),
        organelle_ref=args.organelle_ref,
        enable_nonte_separation=args.enable_nonte,
        barrnap_exe=args.barrnap_exe,
        enable_rmblastn_copies=args.enable_rmblastn_copies,
        enable_seed_chimera=args.enable_seed_chimera,
        enable_flank_trim=args.enable_flank_trim,
        enable_structural_trim=args.enable_structural_trim,
        genome_full_file=(genome_full if sampled_mode else ''),
        defer_copy_floor=sampled_mode,
    )
    config.keep_temp = args.keep_temp

    # Create output directory before saving config
    os.makedirs(args.output, exist_ok=True)

    # Save config
    config.save(os.path.join(args.output, 'refiner_mdl_config.json'))

    # Ensure the genome(s) Phase 1 extracts family copies from are samtools-indexed.
    # A missing .fai makes every Phase-1 shard worker crash (No such file: *.fai) and
    # silently fall back to raw seeds — Phase 1 no-ops and ships unrefined fragments.
    import subprocess as _sp
    for _g in (config.genome_file, getattr(config, 'genome_full_file', '')):
        if _g and os.path.isfile(_g) and not os.path.isfile(_g + '.fai'):
            logger.info(f"Indexing genome for Phase-1 copy extraction: {_g}")
            try:
                _sp.run([getattr(config, 'samtools_exe', 'samtools'), 'faidx', _g],
                        check=True)
            except Exception as e:
                logger.warning(f"samtools faidx failed for {_g}: {e}")

    pipeline = RefinerMdlPipeline(config)
    results = pipeline.run()

    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print(f"Output directory:  {args.output}")
    print(f"Masking library:   {results.get('masking_library', 'N/A')}")
    print(f"Analysis library:  {results.get('analysis_library', 'N/A')}")
    print(f"Statistics:        {results.get('statistics', 'N/A')}")
    print("=" * 60)


if __name__ == '__main__':
    main()
