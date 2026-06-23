"""
M4 integration harness — real testA mdl-repeat output through Phase 0 -> Phase 1.

NOT a unit test: invokes MAFFT / CD-HIT-EST / samtools on the real testA fixture
(/home/shuoc/tool/mdl-repeat/tests). Validates the data-integrity invariant
(output == input + chimera_fragments_added - post_merge_merged, no silent loss) and
the per-shard resume (a half-completed run skips .done shards and reproduces the
single-shot record count).

Run (PGTA env), from a clean working dir:
    python3 bin/Refiner_mdl/tests/integration_m4_orchestrator.py
"""

import logging
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig          # noqa: E402
from phase0_triage import run_phase0          # noqa: E402
import phase1_consensus as orch               # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s %(name)s: %(message)s')
log = logging.getLogger('m4_harness')

DATA = "/home/shuoc/tool/mdl-repeat/tests"
GENOME = f"{DATA}/data/testA.fa"
FAMILIES = f"{DATA}/results/testA_families.fa"
BED = f"{DATA}/results/testA_instances.bed"


def _mk_config(work_dir, shard_size, threads=3):
    return RefinerMdlConfig(
        input_file=FAMILIES,
        genome_file=GENOME,
        output_dir=os.path.join(work_dir, 'out'),
        temp_dir=os.path.join(work_dir, 'temp'),
        checkpoint_dir=os.path.join(work_dir, 'ckpt'),
        bed_file=BED,
        threads=threads,
        shard_size=shard_size,
        genome_blast_db="",          # no fallback DB -> under-instanced keep seed
    )


def _shard_state(shard_dir, n):
    return {i: os.path.exists(orch._shard_paths(shard_dir, i)['done'])
            for i in range(n)}


def main():
    for f in (GENOME, FAMILIES, BED):
        if not os.path.isfile(f):
            log.error("missing fixture: %s", f)
            sys.exit(1)
    if not os.path.exists(GENOME + '.fai'):
        log.error("missing %s.fai (run: samtools faidx %s)", GENOME, GENOME)
        sys.exit(1)

    work = tempfile.mkdtemp(prefix='m4_harness_')
    log.info("work dir: %s", work)
    try:
        # ── Phase 0 (shared by both runs) ──────────────────────────────
        cfg = _mk_config(work, shard_size=1, threads=3)
        os.makedirs(cfg.temp_dir, exist_ok=True)
        p0 = run_phase0(cfg)
        input_count = len(p0['records'])
        log.info("PHASE0: %d input families -> %d records; index=%s",
                 p0['stats']['input_count'], input_count, p0['instance_index'])
        n_index_rows = sum(1 for _ in open(p0['instance_index'])) \
            if p0['instance_index'] else 0
        log.info("PHASE0: instance index rows = %d", n_index_rows)

        # ── Run A: single shot, shard_size=1 (force one shard per family) ──
        runA = os.path.join(work, 'runA')
        cfgA = _mk_config(runA, shard_size=1, threads=3)
        os.makedirs(cfgA.temp_dir, exist_ok=True)
        # re-point the index into runA temp (write_instance_index lives in p0 temp;
        # re-run phase0 per run dir to keep temp self-contained)
        p0A = run_phase0(cfgA)
        outA = orch.run_phase1(p0A, cfgA)
        sA = outA['stats']
        shard_dir_A = os.path.join(cfgA.temp_dir, 'phase1_shards')
        log.info("RUN A stats: %s", sA)
        log.info("RUN A shard .done state: %s",
                 _shard_state(shard_dir_A, sA['n_shards']))

        expected = sA['input_count'] + sA['fragments_added'] - sA['merged']
        assert sA['output_count'] == expected, (
            f"conservation broken: {sA['output_count']} != {expected}")
        assert sA['conservation_ok'], "conservation_ok flag false"
        log.info("RUN A CONSERVATION OK: output(%d) == input(%d) + "
                 "fragments_added(%d) - merged(%d)",
                 sA['output_count'], sA['input_count'],
                 sA['fragments_added'], sA['merged'])
        countA = sA['output_count']
        srcs = {}
        for r in outA['records']:
            srcs[r.get('consensus_source', '?')] = \
                srcs.get(r.get('consensus_source', '?'), 0) + 1
        log.info("RUN A consensus_source breakdown: %s", srcs)

        # ── Run B: resume — complete fully, wipe some shards' .done, rerun ──
        runB = os.path.join(work, 'runB')
        cfgB = _mk_config(runB, shard_size=1, threads=3)
        os.makedirs(cfgB.temp_dir, exist_ok=True)
        p0B = run_phase0(cfgB)
        outB1 = orch.run_phase1(p0B, cfgB)
        shard_dir_B = os.path.join(cfgB.temp_dir, 'phase1_shards')
        n_shards = outB1['stats']['n_shards']
        log.info("RUN B first pass done; shards=%d state=%s",
                 n_shards, _shard_state(shard_dir_B, n_shards))

        # Simulate a kill that left shard 0 done but shards 1..n-1 incomplete:
        # remove their .done + refined + stats (their inputs/instances remain).
        wiped = []
        for i in range(1, n_shards):
            p = orch._shard_paths(shard_dir_B, i)
            for k in ('done', 'refined', 'stats', 'consensus'):
                if os.path.exists(p[k]):
                    os.remove(p[k])
            wiped.append(i)
        log.info("RUN B simulated kill: wiped .done for shards %s "
                 "(shard 0 left complete)", wiped)
        log.info("RUN B pre-resume .done state: %s",
                 _shard_state(shard_dir_B, n_shards))

        outB2 = orch.run_phase1(p0B, cfgB)
        sB2 = outB2['stats']
        log.info("RUN B resume stats: %s", sB2)
        log.info("RUN B post-resume .done state: %s",
                 _shard_state(shard_dir_B, n_shards))

        assert sB2['output_count'] == countA, (
            f"resume count {sB2['output_count']} != single-shot {countA}")
        assert sB2['conservation_ok'], "resume conservation_ok false"
        log.info("RUN B RESUME OK: resumed output(%d) == single-shot output(%d)",
                 sB2['output_count'], countA)

        print("\n" + "=" * 64)
        print("M4 INTEGRATION HARNESS — REAL testA")
        print("=" * 64)
        print(f"input families            : {input_count}")
        print(f"instance index rows       : {n_index_rows}")
        print(f"shards (shard_size=1)     : {sA['n_shards']}")
        print(f"chimera fragments added   : {sA['fragments_added']} "
              f"({sA['split_families']} families split)")
        print(f"post-refine merged        : {sA['merged']}")
        print(f"final output records      : {sA['output_count']}")
        print(f"conservation (in+add-merge): "
              f"{sA['input_count']}+{sA['fragments_added']}-{sA['merged']} "
              f"= {sA['output_count']}  [OK]")
        print(f"consensus_source breakdown: {srcs}")
        print(f"resume: wiped shards {wiped}, shard 0 kept; "
              f"resumed count {sB2['output_count']} == single-shot {countA}  [OK]")
        print("=" * 64)
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
