"""
M2 unit tests — phase1_align + phase1_boundary.

Inputs are SYNTHETIC List[str] alignments and synthetic Copy objects, constructed
to exercise specific geometric behaviors (extend / trim / interior-protection /
fallback). They are unit-test fixtures, NOT real biological data — the real-family
behavior is validated separately by the M2 integration run on testA.fa.

The one MSA-runnable refine_family case actually invokes MAFFT, so this module must
be run inside an environment where `mafft` is on PATH (PGTA conda env).

Run (from repo root, PGTA env):
    python3 -m unittest bin.Refiner_mdl.tests.test_m2_boundary
or: cd bin/Refiner_mdl && python3 -m unittest tests.test_m2_boundary
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig  # noqa: E402
from phase1_align import (  # noqa: E402
    column_majority,
    occupancy_profile,
    select_regime,
)
from phase1_boundary import (  # noqa: E402
    alignment_qc_ok,
    build_consensus,
    refine_family,
    trim_termini,
)
from phase1_extract import Copy  # noqa: E402


def _cfg(**over):
    c = RefinerMdlConfig()
    for k, v in over.items():
        setattr(c, k, v)
    return c


class TestOccupancyProfile(unittest.TestCase):

    def test_basic(self):
        occ = occupancy_profile(["A-CG", "ATCG"])
        self.assertEqual(list(occ), [1.0, 0.5, 1.0, 1.0])

    def test_n_counts_as_absent(self):
        # N is not a definite base -> absent.
        occ = occupancy_profile(["ANCG", "ATCG"])
        self.assertEqual(list(occ), [1.0, 0.5, 1.0, 1.0])

    def test_empty(self):
        self.assertEqual(len(occupancy_profile([])), 0)


class TestColumnMajority(unittest.TestCase):

    def test_unanimous(self):
        self.assertEqual(column_majority(["AAC", "A-C", "ATC"], 0), ("A", 3, 3))

    def test_majority_with_gap(self):
        # col1: A, -, T -> present A=1,T=1 ... use a clearer majority instead.
        self.assertEqual(column_majority(["AAC", "AAC", "ATC"], 1), ("A", 2, 3))

    def test_all_absent(self):
        self.assertEqual(column_majority(["-N.", "-N."], 0), (None, 0, 0))


class TestSelectRegime(unittest.TestCase):

    def test_small_lowdiv_linsi(self):
        self.assertEqual(select_regime(10, 0.05, _cfg()), "mafft_linsi")

    def test_small_highdiv_auto(self):
        # n small but divergence above the L-INS-i ceiling -> auto.
        self.assertEqual(select_regime(10, 0.2, _cfg()), "mafft_auto")

    def test_medium_auto(self):
        self.assertEqual(select_regime(50, 0.05, _cfg()), "mafft_auto")

    def test_large_abpoa(self):
        self.assertEqual(select_regime(250, 0.05, _cfg()), "abpoa")

    def test_boundary_at_thresholds(self):
        c = _cfg()
        # n == small_family_threshold and div == small_family_max_div -> linsi.
        self.assertEqual(select_regime(30, 0.1, c), "mafft_linsi")
        # n == large_family_threshold -> auto.
        self.assertEqual(select_regime(200, 0.5, c), "mafft_auto")
        # one past large -> abpoa.
        self.assertEqual(select_regime(201, 0.5, c), "abpoa")


class TestTrimTermini(unittest.TestCase):
    """Shared flank survives (extend), unique flank collapses (trim), interior
    low-occupancy column is protected."""

    def setUp(self):
        # 5 rows x 10 cols. Present('A')/absent('-') pattern per column:
        #   col0: row0           (1 present)  -> unique flank, trimmed
        #   col1: row0,1         (2 present)  -> below floor(3), trimmed
        #   col2,3: all 5        (solid)      -> left boundary at col2
        #   col4: row0 only      (1 present)  -> INTERIOR low-occ, must be PROTECTED
        #   col5,6,7: all 5      (solid)      -> right boundary at col7
        #   col8: row0,1         (2 present)  -> trimmed
        #   col9: row0           (1 present)  -> trimmed
        self.rows = [
            "AAAAAAAAAA",  # row0 present everywhere
            "-AAA-AAAA-",  # row1: cols 1,2,3,5,6,7,8
            "--AA-AAA--",  # row2: cols 2,3,5,6,7
            "--AA-AAA--",  # row3
            "--AA-AAA--",  # row4
        ]

    def test_trim_bounds(self):
        left, right = trim_termini(self.rows, _cfg())
        self.assertEqual((left, right), (2, 7))

    def test_interior_lowocc_not_trimmed(self):
        # col4 (interior, present=1) lies strictly inside [left, right] -> kept.
        left, right = trim_termini(self.rows, _cfg())
        self.assertLess(left, 4)
        self.assertGreater(right, 4)

    def test_all_gap_collapses(self):
        rows = ["----", "----"]
        left, right = trim_termini(rows, _cfg())
        self.assertGreater(left, right)  # nothing solid -> empty span


class TestBuildConsensus(unittest.TestCase):

    def test_pure_insertion_column_dropped(self):
        # col0 is all gaps (pure insertion) -> dropped; only col1 'A' emitted.
        rows = ["-A", "-A", "-A", "-A", "-A"]
        cons = build_consensus(rows, 0, 1, [], _cfg())
        self.assertEqual(cons, "A")

    def test_kept_span_solid(self):
        # The trim fixture's [2,7] span: cols 2,3,5,6,7 solid 'A', col4 dropped.
        rows = [
            "AAAAAAAAAA",
            "-AAA-AAAA-",
            "--AA-AAA--",
            "--AA-AAA--",
            "--AA-AAA--",
        ]
        cons = build_consensus(rows, 2, 7, [], _cfg())
        self.assertEqual(cons, "AAAAA")  # col4 (present=1<2) dropped

    def test_low_occupancy_emits_N(self):
        # 8 present, evenly split 2A/2C/2G/2T -> max 0.25 < min_occupancy(0.3) -> N.
        rows = ["A", "A", "C", "C", "G", "G", "T", "T"]
        cons = build_consensus(rows, 0, 0, [], _cfg(consensus_min_occupancy=0.3))
        self.assertEqual(cons, "N")

    def test_no_iupac_emitted(self):
        rows = ["A", "A", "A", "C", "C"]  # 3A/2C, max 0.6 -> 'A', never an IUPAC code
        cons = build_consensus(rows, 0, 0, [], _cfg())
        self.assertEqual(cons, "A")
        self.assertTrue(set(cons) <= set("ACGTN"))


class TestAlignmentQc(unittest.TestCase):

    def test_good_alignment_passes(self):
        rows = ["ACGT", "ACGT", "ACGT", "ACGT", "ACGT"]
        self.assertTrue(alignment_qc_ok(rows, _cfg()))

    def test_mostly_gap_fails(self):
        # Every column has occupancy 1/5 = 0.2 < 0.3 -> all low -> fail.
        rows = ["ACGT", "----", "----", "----", "----"]
        self.assertFalse(alignment_qc_ok(rows, _cfg()))

    def test_empty_fails(self):
        self.assertFalse(alignment_qc_ok([], _cfg()))


def _mk_copies(seqs, divergence=0.02):
    return [Copy(id=f"c{i}", sequence=s, divergence=divergence, strand="+")
            for i, s in enumerate(seqs)]


class TestRefineFamilyFallback(unittest.TestCase):

    def test_under_min_copies_keeps_seed(self):
        rec = {"id": "mdl_R1", "sequence": "ACGT" * 50, "length": 200}
        copies = _mk_copies(["ACGT" * 50] * 3)  # 3 < min_copies_for_msa(5)
        out = refine_family(rec, copies, _cfg())
        self.assertEqual(out["consensus_source"], "original")
        self.assertEqual(out["sequence"], rec["sequence"])
        self.assertNotIn("_aln", out)

    def test_zero_copies_keeps_seed(self):
        rec = {"id": "mdl_R1", "sequence": "ACGT" * 50}
        out = refine_family(rec, [], _cfg())
        self.assertEqual(out["consensus_source"], "original")

    def test_runnable_msa_refines(self):
        """Real MAFFT invocation: 6 low-divergence ~180bp copies should produce a
        non-original consensus via the mafft_linsi regime. Requires mafft on PATH."""
        import random
        rng = random.Random(7)
        base = "".join(rng.choice("ACGT") for _ in range(180))
        seqs = []
        for _ in range(6):
            b = list(base)
            for _ in range(4):  # a handful of point substitutions per copy
                p = rng.randrange(len(b))
                b[p] = rng.choice("ACGT")
            seqs.append("".join(b))
        rec = {"id": "mdl_RX", "sequence": base, "length": len(base)}
        out = refine_family(rec, _mk_copies(seqs, divergence=0.03), _cfg())
        self.assertEqual(out["consensus_source"], "mafft_linsi")
        self.assertGreaterEqual(len(out["sequence"]), _cfg().min_chimera_fragment)
        self.assertIn("_aln", out)
        self.assertTrue(set(out["sequence"]) <= set("ACGTN"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
