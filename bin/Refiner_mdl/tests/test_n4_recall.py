"""
N4 unit tests — selective recall: co-fire gate, Tier-2 recruit floor, global
budget keep-and-flag, recruited-Copy type identity, no-DB safety, safe-disable.

All inputs here are SYNTHETIC (clearly labeled). The pieces that would otherwise
hit blastn are exercised either through the no-DB / disabled paths (no subprocess)
or by monkeypatching `recruit_by_blastn` to a deterministic stub. The ONE place a
real blastn matters — that the Tier-2 floor is genuinely lower than the M5 floor —
is asserted at the pident-floor level (`_recall_pident_floor`) without a subprocess.
Real-genome integration (R=150 recall → subfamily emergence) lives in the separate
N4 integration script.

Run:  cd bin/Refiner_mdl && python3 -m unittest tests.test_n4_recall
"""

import os
import sys
import unittest
from types import SimpleNamespace

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig                       # noqa: E402  (import config
# before the phase modules: phase0_triage prepends ../Refiner to sys.path at import,
# which would otherwise shadow this local config.py.)
from phase1_extract import Copy, Instance                 # noqa: E402
from phase1_completeness import recall_eligible, assess_completeness  # noqa: E402
import phase1_fallback as fb                              # noqa: E402
import phase1_consensus as pc                             # noqa: E402


def _cfg(**over):
    c = RefinerMdlConfig(genome_file="/nonexistent/genome.fa")
    for k, v in over.items():
        setattr(c, k, v)
    return c


def _inst(div, chrom="chr1", start=0, end=400, strand="+"):
    """SYNTHETIC instance helper. div = 1 - BED_score/1000."""
    return Instance(chrom=chrom, start=start, end=end, strand=strand,
                    divergence=div)


# ═══════════════════════════════════════════════════════════════════════
# 1. recall_eligible CO-FIRE selectivity (the core selectivity contract)
# ═══════════════════════════════════════════════════════════════════════

class TestRecallEligibleCoFire(unittest.TestCase):
    """under_instanced ALONE must NOT light recall; co-fire with a specific signal must."""

    def test_under_instanced_alone_not_eligible(self):
        # 2 instances (< min_copies_for_msa=5) but TIGHT divergence (no spread) and a
        # SHORT element whose length-class expectation is met -> only under_instanced.
        cfg = _cfg()
        rec = {'id': 'R=under', 'length': 120, 'copies': 2}
        # length 120 (<200 class, expect 2) with n=2 -> length_class does NOT fire;
        # divergences identical -> spread 0.0 -> does NOT fire.
        insts = [_inst(0.05, end=120), _inst(0.05, start=500, end=620)]
        v = recall_eligible(rec, insts, cfg)
        self.assertTrue(v.flags['under_instanced'])
        self.assertFalse(v.flags['bed_divergence_spread'])
        self.assertFalse(v.flags['copy_count_vs_length_class'])
        self.assertFalse(v.incomplete_recall_eligible,
                         "under_instanced alone must NOT trigger expensive recall")

    def test_cofire_divergence_spread_eligible(self):
        # The R=150 archetype: 2 instances (under_instanced) with a WIDE divergence
        # band (0.02 and 0.676 -> spread 0.328 >= 0.31) -> co-fire -> eligible.
        cfg = _cfg()
        rec = {'id': 'R=150like', 'length': 472, 'copies': 2}
        insts = [_inst(0.02), _inst(0.676)]
        v = recall_eligible(rec, insts, cfg)
        self.assertTrue(v.flags['under_instanced'])
        self.assertTrue(v.flags['bed_divergence_spread'])
        self.assertTrue(v.incomplete_recall_eligible,
                        "under_instanced + wide divergence spread must be recall-eligible")

    def test_cofire_length_class_eligible(self):
        # A LONG element (1500 bp, expect >=4 copies) with only 2 instances of TIGHT
        # divergence -> under_instanced + length_class co-fire (no spread needed).
        cfg = _cfg()
        rec = {'id': 'R=long', 'length': 1500, 'copies': 2}
        insts = [_inst(0.05, end=1500), _inst(0.06, start=5000, end=6500)]
        v = recall_eligible(rec, insts, cfg)
        self.assertTrue(v.flags['under_instanced'])
        self.assertTrue(v.flags['copy_count_vs_length_class'])
        self.assertFalse(v.flags['bed_divergence_spread'])
        self.assertTrue(v.incomplete_recall_eligible)

    def test_well_instanced_specific_signal_not_eligible(self):
        # A WELL-instanced family (>= min_copies) with a genuinely WIDE divergence
        # band is NOT recall-eligible: co-fire REQUIRES under_instanced (a well-
        # sampled family that merely has internal divergence already has its members
        # — §2.3 R=2 case, the chr4 20 kb LTR with 65 BED copies and high internal
        # spread that blastn confirms is already complete).
        cfg = _cfg()
        rec = {'id': 'R=2like', 'length': 20000, 'copies': 8}
        # Half near 0.0, half near 0.7 -> population spread well above 0.31.
        insts = [_inst(d, start=i * 30000, end=i * 30000 + 20000)
                 for i, d in enumerate([0.0, 0.02, 0.05, 0.05, 0.68, 0.70, 0.72, 0.70])]
        v = recall_eligible(rec, insts, cfg)
        self.assertFalse(v.flags['under_instanced'])  # 8 >= 5
        self.assertTrue(v.flags['bed_divergence_spread'])  # wide band
        self.assertFalse(v.incomplete_recall_eligible,
                         "well-instanced family must not pay for recall on spread alone")

    def test_disabled_recall_never_eligible(self):
        # enable_selective_recall=False -> nothing eligible (degrades to N3).
        cfg = _cfg(enable_selective_recall=False)
        rec = {'id': 'R=150like', 'length': 472, 'copies': 2}
        insts = [_inst(0.02), _inst(0.676)]
        v = recall_eligible(rec, insts, cfg)
        self.assertFalse(v.incomplete_recall_eligible)
        # The underlying signal flags are still computed (audit), just not acted on.
        self.assertTrue(v.flags['bed_divergence_spread'])


# ═══════════════════════════════════════════════════════════════════════
# 2. Tier-2 recruit floor: recall floor is genuinely lower than the M5 floor
# ═══════════════════════════════════════════════════════════════════════

class TestRecallPidentFloor(unittest.TestCase):

    def test_recall_floor_below_m5_floor(self):
        # recall_identity_floor=0.65 -> 65% pident; the M5 floor is min_recruit_pident
        # =70%. The recall floor MUST be lower so it captures the divergent tail the
        # M5 floor would exclude (the R=150 75-90% partials sit above 65 and the
        # most-divergent below 70).
        cfg = _cfg(recall_identity_floor=0.65, min_recruit_pident=70.0)
        self.assertEqual(fb._recall_pident_floor(cfg), 65.0)
        self.assertLess(fb._recall_pident_floor(cfg), cfg.min_recruit_pident)

    def test_recall_floor_fraction_to_percent(self):
        cfg = _cfg(recall_identity_floor=0.60)
        self.assertEqual(fb._recall_pident_floor(cfg), 60.0)

    def test_recall_mode_uses_recall_floor_in_filter(self):
        # _blast_hits_to_instances with the recall floor keeps a 68%-identity hit that
        # the 70% M5 floor would drop (root cause: lower floor recovers more tail).
        line = "chr4\t100\t300\t68.0\t200\t150.0"
        kept_recall = fb._blast_hits_to_instances(line, pident_min=65.0)
        kept_m5 = fb._blast_hits_to_instances(line, pident_min=70.0)
        self.assertEqual(len(kept_recall), 1)
        self.assertEqual(len(kept_m5), 0)


# ═══════════════════════════════════════════════════════════════════════
# 3. No-DB safety: recall never crashes / never fabricates when DB absent
# ═══════════════════════════════════════════════════════════════════════

class TestRecallNoDBSafety(unittest.TestCase):

    def test_recall_mode_no_db_returns_empty(self):
        cfg = _cfg(genome_blast_db="")
        rec = {'id': 'x', 'sequence': 'ACGT' * 100}
        self.assertEqual(fb.recruit_by_blastn(rec, cfg, recall_mode=True), [])

    def test_recall_mode_empty_sequence_returns_empty(self):
        cfg = _cfg(genome_blast_db="/tmp/somedb")
        self.assertEqual(
            fb.recruit_by_blastn({'id': 'x', 'sequence': ''}, cfg, recall_mode=True),
            [])


# ═══════════════════════════════════════════════════════════════════════
# 4. Worker-level recall: merge, over-budget keep+flag, Copy-type identity
#    (recruit_by_blastn monkeypatched to a deterministic stub — no subprocess)
# ═══════════════════════════════════════════════════════════════════════

class _StubExtract:
    """Deterministic stand-ins so _process_shard runs without a genome/DB."""


def _make_eligible_rec(rid='R=150like'):
    # Header record as Phase 0 would emit; 472 bp seed, 2 wide-divergence BED insts.
    return {'id': rid, 'sequence': 'ACGTACGTAC' * 47, 'length': 472, 'copies': 2,
            'tier': 'T2'}


def _bed_copy(i):
    return Copy(id=f"chr1:{i}-{i+400}(+)", sequence='ACGT' * 100,
                divergence=0.02, strand='+')


def _recalled_copy(i):
    # A DIVERGENT partial as Tier-2 recall would return — SAME Copy type as BED.
    return Copy(id=f"chr1:{900+i}-{900+i+200}(+)", sequence='ACGTACGT' * 25,
                divergence=0.25, strand='+')


class TestWorkerRecallIntegration(unittest.TestCase):
    """Drive _process_shard's recall branch with stubbed extract + recruit."""

    def setUp(self):
        # Patch the heavy I/O functions in phase1_consensus to deterministic stubs.
        self._orig_extract = pc.extract_padded_copies
        self._orig_recruit = pc.recruit_by_blastn
        self._orig_fai = pc.load_fai_lengths
        self._orig_refine = pc._refine_one_family

        pc.load_fai_lengths = lambda p: {'chr1': 10_000_000}
        # BED extract returns 2 copies for the eligible family (under-instanced).
        pc.extract_padded_copies = lambda insts, *a, **k: [_bed_copy(0), _bed_copy(1)]
        # Capture the copies handed to clustering so we can assert recall merged in.
        self.captured = {}

        def _capture_refine(rec, copies, cfg):
            self.captured[rec['id']] = list(copies)
            # Emit one record per family (degenerate single-consensus) so the worker's
            # n_out == n_in invariant holds and we isolate the recall behaviour.
            out = dict(rec)
            out.pop('_aln', None)
            return [out]

        pc._refine_one_family = _capture_refine

    def tearDown(self):
        pc.extract_padded_copies = self._orig_extract
        pc.recruit_by_blastn = self._orig_recruit
        pc.load_fai_lengths = self._orig_fai
        pc._refine_one_family = self._orig_refine

    def _run_worker(self, cfg, recall_stub):
        import json
        import tempfile
        pc.recruit_by_blastn = recall_stub
        shard_dir = tempfile.mkdtemp(prefix='n4_worker_')
        paths = pc._shard_paths(shard_dir, 0)
        rec = _make_eligible_rec()
        with open(paths['records'], 'w') as fh:
            fh.write(json.dumps(rec) + '\n')
        # 2 BED instances with wide divergence -> recall_eligible co-fires.
        with open(paths['instances'], 'w') as fh:
            fh.write("R=150like\t.\tchr1\t0\t400\t+\t0.02\n")
            fh.write("R=150like\t.\tchr1\t500\t900\t+\t0.676\n")
        cfg_dict = dict(cfg.__dict__)
        stats = pc._process_shard(0, shard_dir, cfg_dict)
        refined = pc._read_jsonl(paths['refined'])
        return stats, refined

    def test_eligible_within_budget_recalls_and_merges(self):
        cfg = _cfg(genome_blast_db='/tmp/fakedb', recall_worker_budget=10,
                   enable_selective_recall=True)
        stub = lambda rec, cfg, recall_mode=False: [_recalled_copy(0), _recalled_copy(1)]
        stats, refined = self._run_worker(cfg, stub)
        # Recall fired once; 2 BED + 2 recalled copies reached clustering.
        self.assertEqual(stats['n_recalled'], 1)
        self.assertEqual(stats['n_over_budget'], 0)
        self.assertEqual(len(self.captured['R=150like']), 4)
        # The record is flagged verified with the recalled count.
        self.assertTrue(refined[0]['completeness_verified'])
        self.assertEqual(refined[0]['recalled_copies'], 2)

    def test_recalled_copies_same_type_as_bed(self):
        cfg = _cfg(genome_blast_db='/tmp/fakedb', recall_worker_budget=10)
        stub = lambda rec, cfg, recall_mode=False: [_recalled_copy(0)]
        self._run_worker(cfg, stub)
        merged = self.captured['R=150like']
        # Every object handed to clustering is a phase1_extract.Copy — BED and recalled
        # alike — so N2 cannot distinguish provenance (the §2.5 invariant).
        self.assertTrue(all(isinstance(c, Copy) for c in merged))
        self.assertEqual(len(merged), 3)  # 2 BED + 1 recalled

    def test_over_budget_keeps_bed_and_flags(self):
        # Budget 0 -> eligible family cannot recall -> keep BED copies, flag False.
        cfg = _cfg(genome_blast_db='/tmp/fakedb', recall_worker_budget=0)
        called = {'n': 0}

        def _stub(rec, cfg, recall_mode=False):
            called['n'] += 1
            return [_recalled_copy(0)]

        stats, refined = self._run_worker(cfg, _stub)
        self.assertEqual(stats['n_recalled'], 0)
        self.assertEqual(stats['n_over_budget'], 1)
        self.assertEqual(called['n'], 0, "over-budget family must NOT call blastn")
        # BED copies kept (2), record flagged not-verified — never dropped.
        self.assertEqual(len(self.captured['R=150like']), 2)
        self.assertFalse(refined[0]['completeness_verified'])
        # Conservation: the family still produces exactly one output record.
        self.assertEqual(stats['n_out'], stats['n_in'])

    def test_disabled_recall_degrades_to_bed_only(self):
        # enable_selective_recall=False -> no recall, BED copies only, no recall flag,
        # and the recruit stub is never called (degrade-to-N3 safety at the worker
        # level). enable_fallback_recruitment is also off here to isolate the recall
        # path from the legacy M5 fallback (which legitimately recruits when BED is
        # too sparse and is NOT what N4's disable switch governs).
        cfg = _cfg(genome_blast_db='/tmp/fakedb', recall_worker_budget=10,
                   enable_selective_recall=False, enable_fallback_recruitment=False)
        called = {'n': 0}

        def _stub(rec, cfg, recall_mode=False):
            called['n'] += 1
            return [_recalled_copy(0)]

        stats, refined = self._run_worker(cfg, _stub)
        self.assertEqual(stats['n_recalled'], 0)
        self.assertEqual(stats['n_over_budget'], 0)
        self.assertEqual(called['n'], 0)
        self.assertEqual(len(self.captured['R=150like']), 2)
        self.assertNotIn('completeness_verified', refined[0])


# ═══════════════════════════════════════════════════════════════════════
# 5. Budget split arithmetic (parent splits the global cap across workers)
# ═══════════════════════════════════════════════════════════════════════

class TestBudgetSplit(unittest.TestCase):

    def test_ceil_division_covers_global_budget(self):
        # The parent computes ceil(budget / n_workers) per worker; n_workers copies
        # of that per-worker cap must be >= the global budget (no under-provision).
        for budget, nw in [(5000, 8), (610, 4), (3, 1), (100, 7)]:
            per = (budget + nw - 1) // nw
            self.assertGreaterEqual(per * nw, budget)


if __name__ == '__main__':
    unittest.main(verbosity=2)
