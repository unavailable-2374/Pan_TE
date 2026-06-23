"""
M3 unit tests — phase1_chimera.detect_chimera (MSA-based, BLAST-free).

Inputs are SYNTHETIC List[str] alignments, hand-constructed to exercise the
AND-gate geometry (true composite vs AT-rich single family vs degenerate splits).
They are unit-test fixtures, NOT real biological data — real-family behavior is
validated separately by the M3 integration run on testA.fa (see
tests/integration_m3_chimera.py).

Each fixture's geometry is documented inline (occupancy per column block, the two
copy-span groups). No external tools are invoked here; detect_chimera operates on
the supplied alignment rows directly.

Run (from repo root):
    python3 -m unittest bin.Refiner_mdl.tests.test_m3_chimera
or: cd bin/Refiner_mdl && python3 -m unittest tests.test_m3_chimera
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig  # noqa: E402
from phase1_chimera import detect_chimera, make_fragment, strip_aln  # noqa: E402


def _cfg(**over):
    c = RefinerMdlConfig()
    for k, v in over.items():
        setattr(c, k, v)
    return c


def _rec(rows, sequence, left=None, right=None, **extra):
    """Build a refined-style record with a transient _aln handle.

    metas are None placeholders: detect_chimera's span logic reads only rows, and
    build_consensus only consults metas when consensus_divergence_weighted is on
    (default off), so None is safe here and keeps fixtures terse."""
    if left is None:
        left = 0
    if right is None:
        right = len(rows[0]) - 1
    rec = {'id': 'mdl_RX', 'sequence': sequence,
           '_aln': (rows, None, left, right)}
    rec.update(extra)
    return rec


# ═══════════════════════════════════════════════════════════════════════
# True composite -> should split
# ═══════════════════════════════════════════════════════════════════════

class TestTrueComposite(unittest.TestCase):
    """A 120-col MSA where copies fall into a left-reaching 'A' group and a
    right-reaching 'C' group, joined only by a few spanning copies. The junction
    (cols 55-64) is covered by the spanning copies alone -> occupancy valley AND
    span bimodality both fire -> split."""

    def _build_rows(self):
        # left-only (6): A in cols 0-54, gap 55-119
        left_only = ["A" * 55 + "-" * 65 for _ in range(6)]
        # right-only (6): gap 0-64, C in cols 65-119
        right_only = ["-" * 65 + "C" * 55 for _ in range(6)]
        # spanning (2): full-length G across the whole element
        spanning = ["G" * 120 for _ in range(2)]
        return left_only + right_only + spanning  # n = 14, W = 120

    def test_splits_into_two(self):
        rows = self._build_rows()
        # Consensus the fixture imitates: A*55 (junction drops) C*55. The exact seq
        # is irrelevant to the split decision; supply a plausible full consensus.
        seq = "A" * 55 + "C" * 55
        out = detect_chimera(_rec(rows, seq), _cfg(), depth=0)
        self.assertEqual(len(out), 2, "true composite must split into 2 fragments")
        # Fragments are returned _aln-free (must not leak the MSA past Phase 1).
        for frag in out:
            self.assertNotIn('_aln', frag)
            self.assertTrue(frag.get('is_chimera_fragment'))

    def test_split_point_near_junction(self):
        rows = self._build_rows()
        seq = "A" * 55 + "C" * 55
        out = detect_chimera(_rec(rows, seq), _cfg(), depth=0)
        # frag1 sequence length == split_base; the A/C junction sits at base 55.
        # The valley spans cols 55-64 so the midpoint maps a few bases past 55;
        # assert it lands in the junction neighborhood, not deep in either arm.
        frag1_len = len(out[0]['sequence'])
        self.assertTrue(45 <= frag1_len <= 75,
                        f"split base {frag1_len} should sit near the junction (~55)")
        # Concatenation is conserved (no bases invented or lost in the split).
        self.assertEqual(out[0]['sequence'] + out[1]['sequence'], seq)


# ═══════════════════════════════════════════════════════════════════════
# AT-rich single family -> must NOT split (AND-gate precision)
# ═══════════════════════════════════════════════════════════════════════

class TestATRichSingleFamily(unittest.TestCase):
    """A 120-col MSA with a low-occupancy interior dip (cols 55-64), but EVERY copy
    spans the full element (A flank on the left, T flank on the right). Signal 1
    (valley) fires, Signal 2 (bimodality) does not -> no split. This is the AND-gate
    guarding a real AT-rich / indel-rich interior from being fractured."""

    def test_no_split(self):
        rows = []
        for i in range(10):
            junction = ['-'] * 10
            junction[i % 10] = 'A'           # one base per row -> occ 0.1 in junction
            rows.append("A" * 55 + "".join(junction) + "T" * 55)
        seq = "A" * 55 + "T" * 55
        out = detect_chimera(_rec(rows, seq), _cfg(), depth=0)
        self.assertEqual(len(out), 1, "single AT-rich family must NOT split")
        self.assertNotIn('_aln', out[0])

    def test_uniform_full_occupancy_no_split(self):
        # No valley at all: 8 identical full-length copies -> Signal 1 never fires.
        rows = ["ACGT" * 30 for _ in range(8)]   # W = 120
        out = detect_chimera(_rec(rows, "ACGT" * 30), _cfg(), depth=0)
        self.assertEqual(len(out), 1)


# ═══════════════════════════════════════════════════════════════════════
# Prerequisite / guard behavior
# ═══════════════════════════════════════════════════════════════════════

class TestGuards(unittest.TestCase):

    def test_missing_aln_returns_self(self):
        rec = {'id': 'mdl_R1', 'sequence': 'ACGT' * 50}   # no _aln
        out = detect_chimera(rec, _cfg())
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]['sequence'], rec['sequence'])
        self.assertNotIn('_aln', out[0])

    def test_depth_limit_returns_self(self):
        rows = ["A" * 55 + "-" * 65 for _ in range(6)] + \
               ["-" * 65 + "C" * 55 for _ in range(6)] + \
               ["G" * 120 for _ in range(2)]
        rec = _rec(rows, "A" * 55 + "C" * 55)
        out = detect_chimera(rec, _cfg(), depth=_cfg().chimera_max_depth)
        self.assertEqual(len(out), 1)
        self.assertNotIn('_aln', out[0])

    def test_too_few_hits_returns_self(self):
        # n = 4 < chimera_min_hits (5): geometry is composite but evidence too thin.
        rows = ["A" * 55 + "-" * 65, "A" * 55 + "-" * 65,
                "-" * 65 + "C" * 55, "-" * 65 + "C" * 55]
        out = detect_chimera(_rec(rows, "A" * 55 + "C" * 55), _cfg())
        self.assertEqual(len(out), 1)

    def test_narrow_span_returns_self(self):
        # W = 80 < 2 * min_chimera_fragment (100): too short to host two fragments.
        rows = ["A" * 38 + "-" * 42 for _ in range(6)] + \
               ["-" * 42 + "C" * 38 for _ in range(6)]
        out = detect_chimera(_rec(rows, "A" * 40 + "C" * 40), _cfg())
        self.assertEqual(len(out), 1)


# ═══════════════════════════════════════════════════════════════════════
# Short-fragment discard
# ═══════════════════════════════════════════════════════════════════════

class TestShortFragmentDropped(unittest.TestCase):
    """An asymmetric composite whose LEFT arm is shorter than min_chimera_fragment.
    The split is real (valley + bimodality), but the sub-min left fragment is
    discarded; only the right fragment survives."""

    def test_left_fragment_discarded(self):
        # left-only 'A' over cols 0-24 (25 bp), right-only 'C' over cols 35-99,
        # spanning 'G' across 100 cols. Junction cols 25-34 -> valley, v ~ 30.
        left_only = ["A" * 25 + "-" * 75 for _ in range(6)]
        right_only = ["-" * 35 + "C" * 65 for _ in range(6)]
        spanning = ["G" * 100 for _ in range(2)]
        rows = left_only + right_only + spanning   # n = 14, W = 100
        seq = "A" * 30 + "C" * 65
        out = detect_chimera(_rec(rows, seq), _cfg(), depth=0)
        # Left arm (~30 bp) < min_chimera_fragment (50) -> dropped; right arm kept.
        self.assertEqual(len(out), 1)
        self.assertNotIn('_aln', out[0])


# ═══════════════════════════════════════════════════════════════════════
# strip_aln / make_fragment helpers
# ═══════════════════════════════════════════════════════════════════════

class TestHelpers(unittest.TestCase):

    def test_strip_aln_removes_handle(self):
        rec = {'id': 'x', 'sequence': 'ACGT', '_aln': (['ACGT'], None, 0, 3)}
        out = strip_aln(rec)
        self.assertNotIn('_aln', out)
        self.assertIn('_aln', rec, "strip_aln must not mutate the caller's dict")

    def test_strip_aln_noop_when_absent(self):
        rec = {'id': 'x', 'sequence': 'ACGT'}
        self.assertIs(strip_aln(rec), rec)

    def test_make_fragment_fields(self):
        rows = ["ACGTACGT", "ACGTACGT"]
        rec = _rec(rows, "ACGTACGT")
        sub_rows = [r[0:4] for r in rows]
        frag = make_fragment(rec, "ACGT", sub_rows, '_chimfrag1')
        self.assertEqual(frag['id'], 'mdl_RX_chimfrag1')
        self.assertEqual(frag['sequence'], 'ACGT')
        self.assertEqual(frag['actual_length'], 4)
        self.assertTrue(frag['is_chimera_fragment'])
        self.assertEqual(frag['_aln'][2], 0)
        self.assertEqual(frag['_aln'][3], 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
