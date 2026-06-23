#!/usr/bin/env python3
"""
N6 end-to-end driver: Phase0 -> Phase1 -> Phase2 on real chr4, then the §6 metric
harness (member recovery multi-vs-single, divergent-member recovery, subfamily-vs-
divergence, FP-side, boundary).  build_mdl-style invoke (constructs the config directly).

Emits a single JSON (n6_e2e_report.json) plus stdout logging of every phase count, the
conservation equation, and the two library `grep -c '>'`-equivalent counts.  Nothing is
fabricated — every number comes from the real run / real blastn.
"""
import json
import logging
import os
import sys
import time

# Refiner_mdl uses intra-package relative imports (phase2 -> .te_structure_filter), so
# the modules must load as the Refiner_mdl package.  Put bin/ on the path and import via
# the package; ALSO put the package dir on the path so the package's own top-level
# imports (phase1_align etc.) resolve.
_pkg_dir = os.path.dirname(os.path.abspath(__file__))
_bin_dir = os.path.dirname(_pkg_dir)
for _d in (_pkg_dir, _bin_dir):
    if _d not in sys.path:
        sys.path.insert(0, _d)

from Refiner_mdl.config import RefinerMdlConfig
from Refiner_mdl.phase0_triage import run_phase0
from Refiner_mdl.phase1_consensus import run_phase1
from Refiner_mdl.phase2_library_split import run_phase2
import Refiner_mdl.phase1_metrics as metrics

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
log = logging.getLogger('n6_e2e')

GENOME = os.environ.get('N6_GENOME', '/tmp/m6_e2e/genome/chr4.fa')
FAMILIES = os.environ.get(
    'N6_FAMILIES', '/home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_families.fa')
BED = os.environ.get(
    'N6_BED', '/home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_instances.bed')
OUTDIR = os.environ.get('N6_OUTDIR',
                        sys.argv[1] if len(sys.argv) > 1 else '/tmp/m6_e2e/n6_run')


def main():
    t0 = time.time()
    os.makedirs(OUTDIR, exist_ok=True)
    cfg = RefinerMdlConfig(
        input_file=FAMILIES,
        genome_file=GENOME,
        output_dir=os.path.join(OUTDIR, 'refiner_output'),
        temp_dir=os.path.join(OUTDIR, 'refiner_output', 'temp_work'),
        checkpoint_dir=os.path.join(OUTDIR, 'refiner_output', 'checkpoints'),
        bed_file=BED,
        stats_file='',
        threads=8,
        enable_masking=False,
    )
    for d in (cfg.output_dir, cfg.temp_dir, cfg.checkpoint_dir):
        os.makedirs(d, exist_ok=True)

    # ── Phase 0 ───────────────────────────────────────────────────────────
    log.info("=== PHASE 0 ===")
    p0 = run_phase0(cfg)
    p0_records = p0['records']
    # Original seed lengths for the boundary metric (parent-family id -> length).
    original_lengths = {r['id']: len(r['sequence']) for r in p0_records}
    p0_n = len(p0_records)
    log.info("Phase 0: emitted %d records (post hard-filter + tiering + dedup)", p0_n)

    # ── Phase 1 ───────────────────────────────────────────────────────────
    log.info("=== PHASE 1 ===")
    p1 = run_phase1(p0, cfg)
    p1['phase0_stats'] = p0.get('stats', {})
    p1_records = p1['records']
    p1_stats = p1['stats']
    pruned = p1.get('pruned', [])
    log.info("Phase 1 stats: %s", json.dumps(p1_stats))

    # ── Phase 2 ───────────────────────────────────────────────────────────
    log.info("=== PHASE 2 ===")
    p2 = run_phase2(p1, cfg)
    p2_stats = p2['stats']

    masking_lib = p2['masking_library']
    analysis_lib = p2['analysis_library']

    def count_fa(path):
        if not os.path.exists(path):
            return 0
        with open(path) as fh:
            return sum(1 for ln in fh if ln.startswith('>'))

    masking_seqs = count_fa(masking_lib)
    analysis_seqs = count_fa(analysis_lib)

    # ── consensus_source distribution (incl. subfamily multi-consensus) ─────
    src_dist = {}
    sf_records = 0
    for r in p1_records:
        s = r.get('consensus_source', 'original')
        src_dist[s] = src_dist.get(s, 0) + 1
        if '_sf' in r['id']:
            sf_records += 1

    # ── §6 metric harness ───────────────────────────────────────────────────
    log.info("=== METRIC HARNESS (§6) ===")
    report = metrics.run_metrics(
        p1_records, cfg,
        bed_path=BED,
        pruned_records=pruned,
        original_lengths=original_lengths,
        n_input_families=p0_n,
        min_cov=0.5, pident_min=70.0)

    # ── assemble the e2e summary ───────────────────────────────────────────
    summary = {
        'elapsed_s': round(time.time() - t0, 1),
        'phase_counts': {
            'phase0_records': p0_n,
            'phase1_output': p1_stats['output_count'],
            'phase1_fragments_added': p1_stats['fragments_added'],
            'phase1_split_families': p1_stats['split_families'],
            'phase1_merged': p1_stats['merged'],
            'phase1_n_pruned': p1_stats['n_pruned'],
            'phase1_n_recalled': p1_stats['n_recalled'],
            'phase1_n_over_budget': p1_stats['n_over_budget'],
            'phase1_conservation_ok': p1_stats['conservation_ok'],
            'phase2_masking_seqs': masking_seqs,
            'phase2_analysis_seqs': analysis_seqs,
            'phase2_qc_failed': p2_stats.get('qc_failed', 0),
            'phase2_no_te_signal': p2_stats.get('no_te_signal', 0),
        },
        'conservation_equation': {
            'input': p1_stats['input_count'],
            'fragments_added': p1_stats['fragments_added'],
            'pruned': p1_stats['n_pruned'],
            'merged': p1_stats['merged'],
            'output': p1_stats['output_count'],
            'holds': (p1_stats['output_count']
                      == p1_stats['input_count'] + p1_stats['fragments_added']
                      - p1_stats['n_pruned'] - p1_stats['merged']),
        },
        'consensus_source_distribution': src_dist,
        'subfamily_records': sf_records,
        'masking_library': masking_lib,
        'analysis_library': analysis_lib,
        'metrics': report,
        'phase2_header_sample_note': 'see stdout grep below',
    }

    out_json = os.path.join(OUTDIR, 'n6_e2e_report.json')
    with open(out_json, 'w') as fh:
        json.dump(summary, fh, indent=2, default=str)
    log.info("E2E report written: %s", out_json)

    # ── stdout digest ───────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("N6 END-TO-END DIGEST (real chr4)")
    print("=" * 70)
    print(json.dumps(summary['phase_counts'], indent=2))
    print("--- conservation ---")
    print(json.dumps(summary['conservation_equation'], indent=2))
    print("--- consensus_source ---")
    print(json.dumps(src_dist, indent=2), "subfamily_records=", sf_records)
    mr = report['member_recovery']
    if 'error' in mr:
        print("MEMBER RECOVERY ERROR:", mr['error'])
    else:
        print("--- member recovery (all families) ---")
        print("multi :", json.dumps(mr['recovery_multi_dist']))
        print("single:", json.dumps(mr['recovery_single_dist']))
        print("--- member recovery (split families only, n=%d) ---"
              % mr['n_split_families'])
        print("multi :", json.dumps(mr['split_family_recovery_multi_dist']))
        print("single:", json.dumps(mr['split_family_recovery_single_dist']))
        print("--- divergent-member recovery ---")
        print(json.dumps(mr['divergent_member_recovery'], indent=2))
    if 'fp_side' in report:
        print("--- FP-side (N5 prune) ---")
        fp = report['fp_side']
        print("dropped=%d rate=%.4f reasons=%s"
              % (fp['n_dropped'], fp['drop_rate'], json.dumps(fp['by_reason'])))
    if 'boundary' in report:
        print("--- boundary (refined vs original length) ---")
        print(json.dumps(report['boundary'], indent=2))
    print("=" * 70)
    print("DONE in %.1fs" % summary['elapsed_s'])


if __name__ == '__main__':
    main()
