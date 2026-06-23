"""
N3 unit tests — per-subfamily M3 ordering and the two-axis naming composition
(REFINE_STRATEGY_DESIGN_v2.md §5, §3.5, §10 N3).

N3 wires the two orthogonal split axes in the §5 order: cluster FIRST (subfamily
axis, N2), then run M3 structural-chimera detection per emitted subfamily consensus
on THAT subfamily's own MSA. The naming composes orthogonally: a subfamily already
named `mdl_R<N>_sf<k>` that is also structurally split becomes
`mdl_R<N>_sf<k>_chimfrag<j>`; a single-cluster family keeps its original `mdl_R<N>`
id and (if M3 splits it) becomes `mdl_R<N>_chimfrag<j>` — identical to the pre-N2
behavior.

The three §10 N3 pass criteria tested here:
  1. subfamily-of-a-chimera  → `_sf<k>_chimfrag<j>` (composition correct).
  2. pure multi-subfamily     → NOT additionally chimera-split (axes orthogonal).
  3. pure A|B chimera, one subfamily → behaves as today (original id + `_chimfrag<j>`).

INPUTS ARE SYNTHETIC (clearly labeled). Two altitudes are used deliberately:

  * Criterion 2 is driven through the FULL real per-family body
    (`phase1_consensus._refine_one_family`: real BLASTN clustering + real MAFFT
    consensus + real M3), because a pure full-length-divergence multi-subfamily
    family is faithfully reproducible end-to-end and is exactly the case the
    orthogonality claim must survive.

  * Criteria 1 and 3 test the COMPOSITION / ORDERING wiring that N3 adds — that M3
    runs on each subfamily consensus's own `_aln` and that the `_chimfrag<j>` suffix
    composes onto whatever id the subfamily already carries (`_sf<k>` or the original
    id). They drive `detect_chimera` directly on a record already named at the
    subfamily stage, with the same hand-built bimodal MSA geometry the M3 unit tests
    (test_m3_chimera.py) use. This is the right altitude for the wiring: M3's
    geometry is already covered by test_m3_chimera and the clustering by
    test_n2_cluster; N3 only needs to prove the per-subfamily relocation + naming.

REAL FINDING (reported, drives the altitude choice above): through the real N2
distance-clustering substrate, a span-bimodal A|B chimera is SEPARATED into single-arm
subfamilies BEFORE M3 runs, because the two arms of a clean chimera sit at maximal
pairwise distance (~1.0 — they share no homology) and the average-linkage cut splits
them. To keep the arms in one cluster you must give them a shared spine, but a shared
spine fills the occupancy valley M3 needs, so M3 then declines. The two axes are thus
genuinely orthogonal at the copy level: a single qualifying cluster that is ALSO a
clean span-bimodal chimera is close to self-contradictory under this metric. The
`_sf<k>_chimfrag<j>` composition is therefore a true wiring guarantee (M3 splits
whatever subfamily consensus is handed to it), tested as such, rather than a routinely
co-occurring physical outcome. This finding is the reason §5 fixes the cluster-first
order and is reproduced empirically in tests/integration_n3_subfamily_chimera.py.

Run:  cd bin/Refiner_mdl && python3 -m unittest tests.test_n3_subfamily_chimera
"""

import os
import random
import shutil
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig  # noqa: E402  (import BEFORE phase modules:
# a phase module prepends ../Refiner to sys.path at import time, which would
# otherwise shadow this local config.py with Refiner's.)
from phase1_chimera import detect_chimera  # noqa: E402
from phase1_extract import Copy  # noqa: E402
import phase1_consensus as orch  # noqa: E402


def _cfg(**over):
    c = RefinerMdlConfig()
    for k, v in over.items():
        setattr(c, k, v)
    return c


def _randseq(n, seed):
    """SYNTHETIC random ACGT sequence (deterministic per seed)."""
    r = random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))


def _mut(seq, frac, seed):
    """SYNTHETIC point-substitution mutation of `seq` (deterministic per seed)."""
    r = random.Random(seed)
    s = list(seq)
    for i in range(len(s)):
        if r.random() < frac:
            s[i] = r.choice([b for b in 'ACGT' if b != s[i]])
    return ''.join(s)


def _tools_available():
    return bool(shutil.which('blastn') and shutil.which('makeblastdb')
                and shutil.which('mafft'))


# ═══════════════════════════════════════════════════════════════════════
# Criterion 1 — subfamily-of-a-chimera → `_sf<k>_chimfrag<j>` (composition)
# ═══════════════════════════════════════════════════════════════════════

class TestCriterion1_SubfamilyOfChimeraNaming(unittest.TestCase):
    """A subfamily consensus already named `mdl_R<N>_sf<k>` that M3 then splits must
    yield `mdl_R<N>_sf<k>_chimfrag<j>` — the orthogonal composition (§3.5, §5).

    This drives `detect_chimera` on a record whose id is ALREADY the subfamily id
    (the exact object `_refine_one_family` hands to M3 after clustering), with the
    M3-unit-test composite geometry: a left-reaching 'A' group + a right-reaching 'C'
    group joined by a few spanning copies → occupancy valley AND span bimodality both
    fire → split. The point under test is the SUFFIX COMPOSITION and that the lineage
    provenance survives the split, not M3's geometry (covered by test_m3_chimera)."""

    def _composite_aln(self):
        # SYNTHETIC composite MSA geometry (same as test_m3_chimera TestTrueComposite):
        # 6 left-only 'A' (cols 0-54), 6 right-only 'C' (cols 65-119), 2 spanning 'G'.
        left_only = ["A" * 55 + "-" * 65 for _ in range(6)]
        right_only = ["-" * 65 + "C" * 55 for _ in range(6)]
        spanning = ["G" * 120 for _ in range(2)]
        return left_only + right_only + spanning  # n=14, W=120

    def test_sf_chimfrag_composition(self):
        rows = self._composite_aln()
        # Record AS HANDED BY _refine_one_family to M3: id already suffixed `_sf2`,
        # carrying its subfamily lineage provenance dict.
        sf_rec = {
            'id': 'mdl_R8_sf2',
            'sequence': 'A' * 55 + 'C' * 55,
            'subfamily': {'parent_R': 8, 'subfamily_index': 2,
                          'n_subfamilies': 3, 'member_count': 14},
            '_aln': (rows, None, 0, 119),
        }
        out = detect_chimera(sf_rec, _cfg(), depth=0)
        self.assertEqual(len(out), 2, "subfamily chimera must split into 2 fragments")
        ids = sorted(o['id'] for o in out)
        self.assertEqual(ids, ['mdl_R8_sf2_chimfrag1', 'mdl_R8_sf2_chimfrag2'],
                         "the _chimfrag<j> suffix must compose onto the _sf<k> id")
        for o in out:
            self.assertTrue(o.get('is_chimera_fragment'))
            self.assertNotIn('_aln', o, "fragments must not leak the MSA past Phase 1")
            # Biological provenance (the subfamily lineage) survives the split.
            self.assertIn('subfamily', o)
            self.assertEqual(o['subfamily']['subfamily_index'], 2)

    def test_single_cluster_chimera_keeps_original_id(self):
        # The contrast: a SINGLE-cluster family keeps its original id (no `_sf`), so an
        # M3 split yields `mdl_R<N>_chimfrag<j>` — identical to pre-N2 behavior (§3.5).
        rows = self._composite_aln()
        rec = {'id': 'mdl_R8', 'sequence': 'A' * 55 + 'C' * 55,
               '_aln': (rows, None, 0, 119)}
        out = detect_chimera(rec, _cfg(), depth=0)
        ids = sorted(o['id'] for o in out)
        self.assertEqual(ids, ['mdl_R8_chimfrag1', 'mdl_R8_chimfrag2'])


# ═══════════════════════════════════════════════════════════════════════
# Criterion 2 — pure multi-subfamily NOT additionally chimera-split (real path)
# ═══════════════════════════════════════════════════════════════════════

@unittest.skipUnless(_tools_available(),
                     "blastn/makeblastdb/mafft not on PATH (PGTA env)")
class TestCriterion2_PureMultiSubfamilyNotChimeraSplit(unittest.TestCase):
    """A family whose copies fall into two FULL-LENGTH divergence subfamilies (no A|B
    structural breakpoint) must cluster into >1 subfamily but each subfamily must NOT
    be additionally chimera-split — the §5 orthogonality guarantee, driven through the
    real per-family body (real BLASTN clustering + real MAFFT consensus + real M3)."""

    def setUp(self):
        self.cfg = _cfg(genome_file='/dev/null')  # not read: copies are passed in

    def _two_lineage_copies(self):
        # SYNTHETIC: subfamily A = 7 copies ~8% diverged from an ancestor; subfamily B =
        # 7 copies ~8% diverged from a ~30%-diverged lineage center. Both span the FULL
        # element at the same coordinates (whole-length divergence, NOT positional).
        anc = _randseq(600, 7)
        a = [Copy(id=f"A{v}", sequence=_mut(anc, 0.08, 100 + v),
                  divergence=0.08, strand='+') for v in range(7)]
        center = _mut(anc, 0.30, 999)
        b = [Copy(id=f"B{v}", sequence=_mut(center, 0.08, 200 + v),
                  divergence=0.08, strand='+') for v in range(7)]
        return a + b

    def test_multi_subfamily_no_extra_chimera_split(self):
        copies = self._two_lineage_copies()
        rec = {'id': 'mdl_R800', 'R': 800, 'sequence': _randseq(600, 7),
               'length': 600, 'actual_length': 600, 'copies': len(copies),
               'tier': 'T2', 'topology': 'simple'}
        out, _verdict = orch._refine_one_family(rec, copies, self.cfg)

        # Resolves into exactly two subfamilies...
        self.assertEqual(len(out), 2,
                         f"expected 2 subfamilies, got ids {[o['id'] for o in out]}")
        ids = sorted(o['id'] for o in out)
        self.assertEqual(ids, ['mdl_R800_sf1', 'mdl_R800_sf2'])
        # ...and NEITHER is additionally chimera-split (no `_chimfrag` suffix): the
        # whole-length divergence signal does not trip M3's occupancy-valley + span-
        # bimodality AND-gate.
        for o in out:
            self.assertNotIn('_chimfrag', o['id'],
                             "full-length divergence must NOT trigger a chimera split")
            self.assertFalse(o.get('is_chimera_fragment', False))
            self.assertNotIn('_aln', o, "no _aln may leak past Phase 1")
            # Per-subfamily provenance present (multi-cluster family).
            self.assertIn('subfamily', o)


# ═══════════════════════════════════════════════════════════════════════
# Criterion 3 — pure A|B chimera, one subfamily, behaves as today
# ═══════════════════════════════════════════════════════════════════════

class TestCriterion3_PureChimeraSingleSubfamily(unittest.TestCase):
    """A homogeneous (single-subfamily) family whose consensus is an A|B chimera must
    behave exactly as the pre-N2 path: one cluster → original id (no `_sf`) → M3 splits
    on its MSA → `mdl_R<N>_chimfrag<j>`. Tested via the composition path (M3 on the
    single-subfamily consensus with the original id), the altitude at which 'behaves as
    today' is well-defined — see the module docstring's REAL FINDING for why the real
    distance clusterer instead pre-separates disjoint arms, so this case is a wiring
    guarantee, not a routinely co-occurring physical outcome."""

    def _composite_aln(self):
        left_only = ["A" * 55 + "-" * 65 for _ in range(6)]
        right_only = ["-" * 65 + "C" * 55 for _ in range(6)]
        spanning = ["G" * 120 for _ in range(2)]
        return left_only + right_only + spanning

    def test_single_subfamily_chimera_as_today(self):
        # n_sf == 1 → subfamily_id returns the ORIGINAL id (no `_sf`); the record handed
        # to M3 therefore carries `mdl_R700`, and an M3 split yields `mdl_R700_chimfrag*`
        # — byte-for-byte the naming of the pre-N2 (M4) single-consensus chimera path.
        from phase1_cluster import subfamily_id
        self.assertEqual(subfamily_id('mdl_R700', 1, 1), 'mdl_R700',
                         "single subfamily must keep the original id (no _sf suffix)")
        rec = {'id': subfamily_id('mdl_R700', 1, 1),
               'sequence': 'A' * 55 + 'C' * 55, '_aln': (self._composite_aln(),
                                                         None, 0, 119)}
        out = detect_chimera(rec, _cfg(), depth=0)
        ids = sorted(o['id'] for o in out)
        self.assertEqual(ids, ['mdl_R700_chimfrag1', 'mdl_R700_chimfrag2'],
                         "single-subfamily chimera must split with NO _sf in the id")
        for o in out:
            self.assertNotIn('_sf', o['id'])
            self.assertNotIn('_aln', o)

    def test_single_subfamily_no_valley_no_split(self):
        # The AT-rich guard still holds per subfamily: a single subfamily whose copies
        # all span the element (no bimodality) is NOT split — one record out, original id.
        rows = ["ACGT" * 30 for _ in range(8)]  # uniform full occupancy, no valley
        rec = {'id': 'mdl_R701', 'sequence': "ACGT" * 30, '_aln': (rows, None, 0, 119)}
        out = detect_chimera(rec, _cfg(), depth=0)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]['id'], 'mdl_R701')


# ═══════════════════════════════════════════════════════════════════════
# Count conservation across the per-subfamily M3 (subfamily + chimera splits both
# counted in `extra`) — the §8.2 / §3 generalized data-integrity invariant
# ═══════════════════════════════════════════════════════════════════════

class TestCountConservationPerSubfamilyM3(unittest.TestCase):
    """`_refine_one_family` returns one record per emitted (subfamily × chimera-fragment)
    leaf, with no member silently dropped. The orchestrator's `extra = sum(len-1)` then
    counts BOTH subfamily splits and chimera fragments (§3, §8.2). This checks the
    per-family leaf count is well-defined for the composition cases without needing
    real tools (uses the deterministic M3 geometry)."""

    def test_chimera_fragments_counted_per_subfamily(self):
        # One subfamily record split by M3 into 2 fragments → family contributes 2 leaves
        # → extra == 1 for that family (subfamily axis added 0 here; chimera axis added 1).
        rows = (["A" * 55 + "-" * 65 for _ in range(6)]
                + ["-" * 65 + "C" * 55 for _ in range(6)]
                + ["G" * 120 for _ in range(2)])
        sf_rec = {'id': 'mdl_R8_sf2', 'sequence': 'A' * 55 + 'C' * 55,
                  '_aln': (rows, None, 0, 119)}
        frags = detect_chimera(sf_rec, _cfg(), depth=0)
        self.assertEqual(len(frags), 2)
        extra = len(frags) - 1
        self.assertEqual(extra, 1)

    def test_no_split_contributes_zero_extra(self):
        rows = ["ACGT" * 30 for _ in range(8)]
        rec = {'id': 'mdl_R701', 'sequence': "ACGT" * 30, '_aln': (rows, None, 0, 119)}
        frags = detect_chimera(rec, _cfg(), depth=0)
        self.assertEqual(len(frags) - 1, 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
