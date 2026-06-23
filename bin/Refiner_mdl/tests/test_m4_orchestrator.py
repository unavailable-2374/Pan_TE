"""
M4 unit tests — phase1_consensus orchestrator (sharding, routing, scaffold
reconstruction, chimera eligibility, post-refine merge, resume).

Inputs are SYNTHETIC small records / index rows, constructed to exercise the
scheduler's bookkeeping and the data-integrity accounting. They are unit-test
fixtures, NOT real biological data — real-family end-to-end behavior (counts
conservation + resume) is validated separately by the M4 integration harness on
testA.fa (tests/integration_m4_orchestrator.py).

merge_post_refine and the resume test invoke CD-HIT-EST, so this module must run in
an environment where `cd-hit-est` is on PATH (PGTA conda env).

Run:  cd bin/Refiner_mdl && python3 -m unittest tests.test_m4_orchestrator
"""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig  # noqa: E402
from phase1_extract import Instance  # noqa: E402
import phase1_consensus as orch  # noqa: E402


def _cfg(**over):
    c = RefinerMdlConfig()
    for k, v in over.items():
        setattr(c, k, v)
    return c


def _rec(rid, seq, copies=5, **extra):
    r = {'id': rid, 'R': int(rid.split('R')[-1]) if 'R' in rid else 0,
         'sequence': seq, 'copies': copies, 'mdl': 100.0, 'length': len(seq),
         'tier': 'T2'}
    r.update(extra)
    return r


# ═══════════════════════════════════════════════════════════════════════
# make_shards — count conservation + heavy-family spread
# ═══════════════════════════════════════════════════════════════════════

class TestMakeShards(unittest.TestCase):

    def test_count_conservation(self):
        recs = [_rec(f'mdl_R{i}', 'A' * (10 + i), copies=i + 1) for i in range(7)]
        shards = orch.make_shards(recs, shard_size=3)
        self.assertEqual(len(shards), 3)                  # ceil(7/3)
        total = sum(len(s) for s in shards)
        self.assertEqual(total, 7)                        # nothing lost or duplicated
        ids = {r['id'] for s in shards for r in s}
        self.assertEqual(ids, {f'mdl_R{i}' for i in range(7)})

    def test_heavy_families_spread(self):
        # Two clearly heaviest families (len*copies) must not both land in shard 0.
        recs = [_rec(f'mdl_R{i}', 'A' * 10, copies=1) for i in range(6)]
        recs[0] = _rec('mdl_R0', 'A' * 1000, copies=100)   # heaviest
        recs[1] = _rec('mdl_R1', 'A' * 900, copies=100)    # second heaviest
        shards = orch.make_shards(recs, shard_size=2)       # 3 shards
        shard_of = {r['id']: i for i, s in enumerate(shards) for r in s}
        self.assertNotEqual(shard_of['mdl_R0'], shard_of['mdl_R1'])

    def test_empty(self):
        self.assertEqual(orch.make_shards([], 3), [])

    def test_single_shard_when_size_exceeds_n(self):
        recs = [_rec(f'mdl_R{i}', 'A' * 10) for i in range(3)]
        shards = orch.make_shards(recs, shard_size=500)
        self.assertEqual(len(shards), 1)
        self.assertEqual(len(shards[0]), 3)


# ═══════════════════════════════════════════════════════════════════════
# route_instances — each row to the owner shard, streaming group-by
# ═══════════════════════════════════════════════════════════════════════

class TestRouteInstances(unittest.TestCase):

    def test_rows_land_in_owner_shard(self):
        with tempfile.TemporaryDirectory() as d:
            index = os.path.join(d, 'idx.tsv')
            # Sorted by record_id (col 0). Two families, owners in different shards.
            rows = [
                "mdl_R0\t0\tchr1\t100\t200\t+\t0.0",
                "mdl_R0\t0\tchr1\t300\t400\t-\t0.1",
                "mdl_R1\t1\tchr2\t500\t600\t+\t0.2",
            ]
            with open(index, 'w') as fh:
                fh.write('\n'.join(rows) + '\n')
            record_to_shard = {'mdl_R0': 0, 'mdl_R1': 1}
            counts = orch.route_instances(index, record_to_shard, d)
            self.assertEqual(counts.get(0), 2)
            self.assertEqual(counts.get(1), 1)
            by0 = orch._load_shard_instances(orch._shard_paths(d, 0)['instances'])
            by1 = orch._load_shard_instances(orch._shard_paths(d, 1)['instances'])
            self.assertEqual(len(by0['mdl_R0']), 2)
            self.assertEqual(len(by1['mdl_R1']), 1)
            self.assertEqual(by0['mdl_R0'][0].chrom, 'chr1')
            self.assertEqual(by0['mdl_R0'][1].strand, '-')

    def test_orphan_rows_skipped(self):
        with tempfile.TemporaryDirectory() as d:
            index = os.path.join(d, 'idx.tsv')
            with open(index, 'w') as fh:
                fh.write("mdl_R9\t9\tchr1\t1\t2\t+\t0.0\n")   # no owner
            counts = orch.route_instances(index, {'mdl_R0': 0}, d)
            self.assertEqual(sum(counts.values()), 0)


# ═══════════════════════════════════════════════════════════════════════
# chimera_eligible — tier / length / topology gate
# ═══════════════════════════════════════════════════════════════════════

class TestChimeraEligible(unittest.TestCase):

    def test_t1_length_gate(self):
        self.assertTrue(orch.chimera_eligible(
            {'tier': 'T1', 'sequence': 'A' * 300}, _cfg()))
        self.assertFalse(orch.chimera_eligible(
            {'tier': 'T1', 'sequence': 'A' * 299}, _cfg()))

    def test_t2_length_gate(self):
        self.assertTrue(orch.chimera_eligible(
            {'tier': 'T2', 'sequence': 'A' * 150}, _cfg()))
        self.assertFalse(orch.chimera_eligible(
            {'tier': 'T2', 'sequence': 'A' * 149}, _cfg()))

    def test_t3_not_eligible_without_complex(self):
        self.assertFalse(orch.chimera_eligible(
            {'tier': 'T3', 'sequence': 'A' * 5000}, _cfg()))

    def test_complex_topology_overrides_length(self):
        self.assertTrue(orch.chimera_eligible(
            {'tier': 'T3', 'sequence': 'A' * 10, 'topology': 'complex'}, _cfg()))
        self.assertTrue(orch.chimera_eligible(
            {'tier': 'T2', 'sequence': 'A' * 10, 'topology': 'complex'}, _cfg()))


# ═══════════════════════════════════════════════════════════════════════
# build_scaffold_instances — spanning reconstruction + fallback
# ═══════════════════════════════════════════════════════════════════════

class TestBuildScaffoldInstances(unittest.TestCase):

    def test_reconstructs_spanning_loci(self):
        # Two assembled element copies: members within fragment_gap merge into one
        # spanning locus each; the two copies sit far apart -> two spanning loci.
        cfg = _cfg(fragment_gap=50, min_copies_for_msa=2)
        gap = cfg.fragment_gap
        routed = [
            Instance('chr1', 100, 200, '+', 0.1),
            Instance('chr1', 200 + gap, 300, '+', 0.1),   # adjacent -> merges
            Instance('chr1', 100000, 100100, '+', 0.2),
            Instance('chr1', 100100 + gap, 100200, '+', 0.2),  # adjacent -> merges
        ]
        spanning = orch.build_scaffold_instances(
            {'id': 'scaffold_1', 'is_scaffold': True}, routed, cfg)
        self.assertEqual(len(spanning), 2)
        self.assertEqual(spanning[0].start, 100)
        self.assertEqual(spanning[0].end, 300)
        self.assertEqual(spanning[1].start, 100000)

    def test_fallback_when_too_few_loci(self):
        # Only one spanning locus < min_copies_for_msa -> empty (keeps seed).
        cfg = _cfg(fragment_gap=50, min_copies_for_msa=5)
        routed = [Instance('chr1', 100, 200, '+', 0.1),
                  Instance('chr1', 240, 300, '+', 0.1)]
        spanning = orch.build_scaffold_instances(
            {'id': 'scaffold_1', 'is_scaffold': True}, routed, cfg)
        self.assertEqual(spanning, [])

    def test_strand_separates_loci(self):
        cfg = _cfg(fragment_gap=50, min_copies_for_msa=2)
        routed = [Instance('chr1', 100, 200, '+', 0.1),
                  Instance('chr1', 210, 300, '-', 0.1)]   # different strand
        spanning = orch.build_scaffold_instances(
            {'id': 'scaffold_1', 'is_scaffold': True}, routed, cfg)
        self.assertEqual(len(spanning), 2)


# ═══════════════════════════════════════════════════════════════════════
# merge_post_refine — no data loss (requires cd-hit-est)
# ═══════════════════════════════════════════════════════════════════════

class TestMergePostRefine(unittest.TestCase):

    def test_distinct_sequences_all_retained(self):
        import random
        rng = random.Random(3)
        recs = [_rec(f'f{i}', ''.join(rng.choice('ACGT') for _ in range(400)))
                for i in range(6)]
        out = orch.merge_post_refine(recs, _cfg())
        self.assertEqual(out['merged'], 0)
        self.assertEqual(len(out['records']), 6)

    def test_identical_collapse_counts_reconcile(self):
        # Three identical + three distinct: identicals collapse, total reconciles.
        import random
        rng = random.Random(5)
        dup = 'ACGT' * 100
        recs = [_rec('d0', dup), _rec('d1', dup), _rec('d2', dup)]
        recs += [_rec(f'u{i}', ''.join(rng.choice('ACGT') for _ in range(400)))
                 for i in range(3)]
        out = orch.merge_post_refine(recs, _cfg())
        self.assertEqual(len(out['records']) + out['merged'], 6)
        self.assertGreaterEqual(out['merged'], 2)   # at least 2 dups collapsed

    def test_under_two_noop(self):
        recs = [_rec('only', 'A' * 400)]
        out = orch.merge_post_refine(recs, _cfg())
        self.assertEqual(out['merged'], 0)
        self.assertEqual(len(out['records']), 1)


# ═══════════════════════════════════════════════════════════════════════
# Resume — fully pre-completed shards are skipped and concatenated
# ═══════════════════════════════════════════════════════════════════════

class TestResumeSkipsDone(unittest.TestCase):
    """All shards pre-marked .done with refined.jsonl + stats.json: run_phase1 must
    spawn NO workers, just concatenate. Validates the resume short-circuit and the
    counts-conservation accounting without invoking any MSA."""

    def test_all_done_concatenates_without_workers(self):
        import random
        rng = random.Random(11)
        recs = [_rec(f'mdl_R{i}',
                     ''.join(rng.choice('ACGT') for _ in range(400)), copies=5)
                for i in range(5)]
        with tempfile.TemporaryDirectory() as d:
            cfg = _cfg(temp_dir=d, threads=2, shard_size=2)
            shards = orch.make_shards(recs, cfg.shard_size)
            shard_dir = os.path.join(d, 'phase1_shards')
            os.makedirs(shard_dir, exist_ok=True)
            # Pre-complete every shard exactly as a prior run would have.
            for i, shard in enumerate(shards):
                p = orch._shard_paths(shard_dir, i)
                orch._write_jsonl(p['records'], shard)
                refined = [dict(r, consensus_source='original') for r in shard]
                orch._write_jsonl(p['refined'], refined)
                with open(p['stats'], 'w') as fh:
                    json.dump({'shard': i, 'n_in': len(shard),
                               'n_out': len(refined), 'extra': 0,
                               'split_families': 0}, fh)
                with open(p['done'], 'w') as fh:
                    fh.write('done\n')

            # Sentinel: if a worker were spawned it would call this and fail loudly.
            orig = orch._process_shard

            def _boom(*a, **k):
                raise AssertionError("worker spawned despite .done marker")

            orch._process_shard = _boom
            try:
                out = orch.run_phase1(
                    {'records': recs, 'instance_index': None}, cfg)
            finally:
                orch._process_shard = orig

            self.assertEqual(out['stats']['input_count'], 5)
            self.assertEqual(out['stats']['fragments_added'], 0)
            self.assertEqual(out['stats']['output_count'] + out['stats']['merged'], 5)
            self.assertTrue(out['stats']['conservation_ok'])


if __name__ == "__main__":
    unittest.main(verbosity=2)
