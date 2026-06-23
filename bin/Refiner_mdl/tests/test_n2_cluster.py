"""
N2 unit tests — phase1_cluster: distance, clustering, qualifying gate, naming,
member conservation.

All inputs here are SYNTHETIC (clearly labeled). They exercise the clustering logic
on constructed homogeneous-vs-two-divergent-group copy sets. Real-data validation
(chr4 homogeneous-family-degrades-to-today + a SYNTHETIC multi-subfamily mechanism
check) lives in the separate N2 integration script
(tests/integration_n2_cluster.py).

Run:  cd bin/Refiner_mdl && python3 -m unittest tests.test_n2_cluster
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig  # noqa: E402  (import BEFORE phase modules:
# a phase module prepends ../Refiner to sys.path at import time, which would
# otherwise shadow this local config.py with Refiner's.)
from phase1_extract import Copy  # noqa: E402
from phase1_cluster import (  # noqa: E402
    Cluster,
    cluster_copies,
    pairwise_distance,
    qualifying_clusters,
    subfamily_id,
    subfamily_lineage,
)


# ── SYNTHETIC sequence builders ───────────────────────────────────────────
# Group A and group B are deliberately near-disjoint in k-mer content (different
# repeated motifs) so the k-mer Jaccard distance separates them; within a group the
# copies are near-identical (a couple of point substitutions) so intra-group distance
# is low. These are SYNTHETIC TE-like copies, NOT real genomic sequence.

def _seqA(variant=0):
    """SYNTHETIC group-A copy: ~480 bp built from one motif, with a FEW point
    substitutions (a near-identical homogeneous subfamily — ~1-2 substitutions per
    copy, so at k=13 the intra-group Jaccard distance stays well below the cut)."""
    base = ("ACGTTGCAACCGTTAGGCATGACTTGACCATGGTACATGCATGCATTAGGCATCAGT"
            "TTGACCATGAGTCATGCATTAGCATGCATCATGGTACATGCATGCATTAGGCATCAG") * 4
    s = list(base)
    # One substitution per variant at a distinct position (keeps copies near-identical
    # while remaining distinguishable). variant 0 == the unmutated reference.
    if variant > 0:
        p = (variant * 53) % len(s)
        s[p] = 'A' if s[p] != 'A' else 'C'
    return ''.join(s)


def _seqB(variant=0):
    """SYNTHETIC group-B copy: ~480 bp built from a DIFFERENT motif, near-identical
    within the group (a FEW substitutions)."""
    base = ("TTAACCGGTTAACCGGAATTCCGGAATTGGCCTTAAGGCCTTAACCGGAATTCCGGT"
            "AATTGGCCAATTCCGGTTAACCGGAATTGGCCTTAAGGCCAATTCCGGTTAACCGGA") * 4
    s = list(base)
    if variant > 0:
        p = (variant * 53) % len(s)
        s[p] = 'T' if s[p] != 'T' else 'G'
    return ''.join(s)


def _copy(seq, cid, div=0.1):
    """SYNTHETIC Copy helper (strand-co-oriented, '+' as extract emits)."""
    return Copy(id=cid, sequence=seq, divergence=div, strand='+')


class TestPairwiseDistance(unittest.TestCase):
    """Homogeneous copies → low pairwise distance; two groups → high inter-group."""

    def setUp(self):
        self.cfg = RefinerMdlConfig()

    def test_homogeneous_low_distance(self):
        # SYNTHETIC: six near-identical group-A copies.
        copies = [_copy(_seqA(v), f"a{v}") for v in range(6)]
        d = pairwise_distance(copies, self.cfg)
        self.assertEqual(d.shape, (6, 6))
        # all off-diagonal distances small (near-identical)
        offdiag = d[~__import__('numpy').eye(6, dtype=bool)]
        self.assertLess(offdiag.max(), 0.3)

    def test_two_groups_high_inter_distance(self):
        # SYNTHETIC: three group-A + three group-B copies. Intra low, inter high.
        copies = ([_copy(_seqA(v), f"a{v}") for v in range(3)]
                  + [_copy(_seqB(v), f"b{v}") for v in range(3)])
        d = pairwise_distance(copies, self.cfg)
        intra_a = d[0, 1]
        intra_b = d[3, 4]
        inter = d[0, 3]
        self.assertLess(intra_a, 0.3)
        self.assertLess(intra_b, 0.3)
        self.assertGreater(inter, 0.7)
        self.assertGreater(inter, intra_a)
        self.assertGreater(inter, intra_b)

    def test_symmetry_and_zero_diagonal(self):
        copies = [_copy(_seqA(v), f"a{v}") for v in range(4)]
        d = pairwise_distance(copies, self.cfg)
        import numpy as np
        self.assertTrue(np.allclose(d, d.T))
        self.assertTrue(np.allclose(np.diag(d), 0.0))

    def test_single_copy_trivial(self):
        d = pairwise_distance([_copy(_seqA(0), "a0")], self.cfg)
        self.assertEqual(d.shape, (1, 1))
        self.assertEqual(d[0, 0], 0.0)


class TestClusterCopies(unittest.TestCase):
    """Homogeneous → 1 cluster; two divergent groups → 2 clusters; no member lost."""

    def setUp(self):
        self.cfg = RefinerMdlConfig()

    def test_homogeneous_one_cluster(self):
        # SYNTHETIC: eight near-identical group-A copies → exactly one cluster.
        copies = [_copy(_seqA(v), f"a{v}") for v in range(8)]
        clusters = cluster_copies(copies, self.cfg)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0].size, 8)
        self.assertGreater(clusters[0].mean_identity, 0.7)

    def test_two_divergent_groups_two_clusters(self):
        # SYNTHETIC: five group-A + five group-B → two clusters.
        copies = ([_copy(_seqA(v), f"a{v}") for v in range(5)]
                  + [_copy(_seqB(v), f"b{v}") for v in range(5)])
        clusters = cluster_copies(copies, self.cfg)
        self.assertEqual(len(clusters), 2)
        # Member conservation: every copy lands in exactly one cluster.
        total = sum(c.size for c in clusters)
        self.assertEqual(total, 10)
        all_idx = sorted(i for c in clusters for i in c.member_indices)
        self.assertEqual(all_idx, list(range(10)))

    def test_deterministic_order_largest_first(self):
        # SYNTHETIC: 6 group-A + 4 group-B → first cluster is the larger (A).
        copies = ([_copy(_seqA(v), f"a{v}") for v in range(6)]
                  + [_copy(_seqB(v), f"b{v}") for v in range(4)])
        clusters = cluster_copies(copies, self.cfg)
        self.assertEqual(len(clusters), 2)
        self.assertGreaterEqual(clusters[0].size, clusters[1].size)
        self.assertEqual(clusters[0].size, 6)

    def test_empty_and_single(self):
        self.assertEqual(cluster_copies([], self.cfg), [])
        one = cluster_copies([_copy(_seqA(0), "a0")], self.cfg)
        self.assertEqual(len(one), 1)
        self.assertEqual(one[0].size, 1)


class TestQualifyingClusters(unittest.TestCase):
    """Gate + never-drop: small-cluster reassignment, sparse fallback, close-no-split."""

    def setUp(self):
        self.cfg = RefinerMdlConfig()

    def test_small_cluster_members_reassigned_not_lost(self):
        # SYNTHETIC: 6 group-A (qualifies, >=5) + 2 group-B (sub-min). The 2 B copies
        # must be reassigned to the nearest qualifying cluster, never dropped.
        copies = ([_copy(_seqA(v), f"a{v}") for v in range(6)]
                  + [_copy(_seqB(v), f"b{v}") for v in range(2)])
        clusters = cluster_copies(copies, self.cfg)
        quals, audit = qualifying_clusters(clusters, copies, self.cfg)
        # member conservation: nothing lost
        self.assertEqual(audit['total_members_out'], 8)
        self.assertEqual(sum(c.size for c in quals), 8)
        self.assertEqual(audit['n_reassigned'], 2)
        self.assertFalse(audit['fallback_single_cluster'])
        all_idx = sorted(i for c in quals for i in c.member_indices)
        self.assertEqual(all_idx, list(range(8)))

    def test_sparse_family_falls_back_to_single_consensus(self):
        # SYNTHETIC: only 3 group-A copies → no cluster reaches min_members=5 →
        # never-regress fallback to ONE all-copy cluster (today's behavior).
        copies = [_copy(_seqA(v), f"a{v}") for v in range(3)]
        clusters = cluster_copies(copies, self.cfg)
        quals, audit = qualifying_clusters(clusters, copies, self.cfg)
        self.assertEqual(len(quals), 1)
        self.assertTrue(audit['fallback_single_cluster'])
        self.assertEqual(quals[0].size, 3)
        self.assertEqual(audit['total_members_out'], 3)

    def test_two_qualifying_groups_both_emitted(self):
        # SYNTHETIC: 5 A + 5 B, both qualify → two subfamily clusters.
        copies = ([_copy(_seqA(v), f"a{v}") for v in range(5)]
                  + [_copy(_seqB(v), f"b{v}") for v in range(5)])
        clusters = cluster_copies(copies, self.cfg)
        quals, audit = qualifying_clusters(clusters, copies, self.cfg)
        self.assertEqual(len(quals), 2)
        self.assertFalse(audit['fallback_single_cluster'])
        self.assertEqual(audit['total_members_out'], 10)

    def test_homogeneous_stays_one_cluster(self):
        # SYNTHETIC: 8 homogeneous copies → one qualifying cluster, no reassignment.
        copies = [_copy(_seqA(v), f"a{v}") for v in range(8)]
        clusters = cluster_copies(copies, self.cfg)
        quals, audit = qualifying_clusters(clusters, copies, self.cfg)
        self.assertEqual(len(quals), 1)
        self.assertEqual(audit['n_reassigned'], 0)
        self.assertFalse(audit['fallback_single_cluster'])
        self.assertEqual(quals[0].size, 8)


class TestNaming(unittest.TestCase):
    """Single cluster keeps original id; multi-cluster suffixes _sf<k>."""

    def test_single_cluster_keeps_original_id(self):
        self.assertEqual(subfamily_id("mdl_R7", 1, 1), "mdl_R7")

    def test_multi_cluster_sf_suffix(self):
        self.assertEqual(subfamily_id("mdl_R7", 1, 3), "mdl_R7_sf1")
        self.assertEqual(subfamily_id("mdl_R7", 2, 3), "mdl_R7_sf2")
        self.assertEqual(subfamily_id("mdl_R7", 3, 3), "mdl_R7_sf3")

    def test_lineage_record_fields(self):
        cfg = RefinerMdlConfig()
        cl = Cluster(members=[1, 2, 3, 4, 5], mean_identity=0.83,
                     member_indices=[0, 1, 2, 3, 4])
        rec = {'id': 'mdl_R7', 'R': 7}
        lin = subfamily_lineage(rec, 2, 3, cl, cfg)
        self.assertEqual(lin['parent_R'], 7)
        self.assertEqual(lin['subfamily_index'], 2)
        self.assertEqual(lin['n_subfamilies'], 3)
        self.assertEqual(lin['member_count'], 5)
        self.assertAlmostEqual(lin['intra_identity'], 0.83, places=4)
        self.assertEqual(lin['divergence_cut'], cfg.subfamily_divergence_cut)


class TestMemberConservation(unittest.TestCase):
    """End-to-end member conservation across cluster → qualify, multiple mixes."""

    def setUp(self):
        self.cfg = RefinerMdlConfig()

    def test_conservation_various_mixes(self):
        mixes = [
            (6, 0),   # homogeneous
            (5, 5),   # two equal qualifying
            (6, 2),   # one qualifying + reassigned tail
            (3, 0),   # sparse fallback
            (7, 3),   # one qualifying + sub-min B reassigned
        ]
        for na, nb in mixes:
            copies = ([_copy(_seqA(v), f"a{v}") for v in range(na)]
                      + [_copy(_seqB(v), f"b{v}") for v in range(nb)])
            clusters = cluster_copies(copies, self.cfg)
            quals, audit = qualifying_clusters(clusters, copies, self.cfg)
            n_total = na + nb
            self.assertEqual(audit['total_members_out'], n_total,
                             msg=f"mix A={na} B={nb} lost members")
            self.assertEqual(sum(c.size for c in quals), n_total)


# ── REALISTIC-DIVERGENCE distance tests ───────────────────────────────────
# These exercise the BLASTN local-identity distance on copies with BIOLOGICALLY
# realistic substitution divergence (NOT the near-identical / near-disjoint motif copies
# of the classes above), to prove the substrate measures graded identity correctly: an
# ~8%-diverged subfamily reads as LOW pairwise distance, a ~30%-diverged sibling reads as
# HIGH inter-subfamily distance, and clustering resolves the two. Still SYNTHETIC (a
# seeded RNG mutates one ancestral sequence), but the divergence levels are real TE-scale,
# so the test fails if the metric ever regresses to a saturating/near-identical-only one.
# Requires blastn + makeblastdb on PATH (PGTA env) — the real distance shells out.

import random as _random  # noqa: E402


def _ancestral(n=600, seed=7):
    r = _random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))


def _mutate(seq, frac, seed):
    """Substitute `frac` of positions to a different base (point divergence only)."""
    r = _random.Random(seed)
    s = list(seq)
    for i in range(len(s)):
        if r.random() < frac:
            s[i] = r.choice([b for b in 'ACGT' if b != s[i]])
    return ''.join(s)


class TestRealisticDivergenceDistance(unittest.TestCase):
    """BLASTN distance on TE-scale divergence: intra-subfamily LOW, inter-subfamily HIGH.

    Skips cleanly if blastn is unavailable so the suite still runs without the tool; in
    the PGTA env it shells out to the real distance, the same code the production path
    uses."""

    def setUp(self):
        import shutil as _sh
        if not (_sh.which('blastn') and _sh.which('makeblastdb')):
            self.skipTest('blastn/makeblastdb not on PATH')
        self.cfg = RefinerMdlConfig()

    def _two_subfamily_copies(self):
        anc = _ancestral()
        # Subfamily A: 6 copies each ~8% diverged from the ancestor (so pairwise ~15%).
        a = [_copy(_mutate(anc, 0.08, 100 + v), f"A{v}", div=0.08) for v in range(6)]
        # Subfamily B: a ~30%-diverged lineage center, then 6 copies ~8% diverged from IT
        # (so intra-B ~15%, but A↔B ~38% — a genuine subfamily split, not near-identical).
        bcenter = _mutate(anc, 0.30, 999)
        b = [_copy(_mutate(bcenter, 0.08, 200 + v), f"B{v}", div=0.30) for v in range(6)]
        return a, b

    def test_intra_low_inter_high(self):
        a, b = self._two_subfamily_copies()
        copies = a + b
        d = pairwise_distance(copies, self.cfg)
        # Intra-subfamily-A pairs (indices 0..5) are tight; A↔B pairs (0..5 vs 6..11) far.
        intra_a = [d[i, j] for i in range(6) for j in range(i + 1, 6)]
        inter_ab = [d[i, j] for i in range(6) for j in range(6, 12)]
        self.assertLess(max(intra_a), 0.3,
                        msg=f"intra-A distance not low: max={max(intra_a):.3f}")
        self.assertGreater(min(inter_ab), 0.3,
                           msg=f"inter A-B distance not high: min={min(inter_ab):.3f}")
        # The metric must SEPARATE: every inter pair strictly farther than every intra.
        self.assertGreater(min(inter_ab), max(intra_a))

    def test_resolves_two_subfamilies(self):
        a, b = self._two_subfamily_copies()
        copies = a + b
        clusters = cluster_copies(copies, self.cfg)
        quals, audit = qualifying_clusters(clusters, copies, self.cfg)
        self.assertEqual(len(quals), 2,
                         msg=f"expected 2 subfamilies, got {len(quals)} "
                             f"(sizes {[c.size for c in quals]})")
        self.assertEqual(audit['total_members_out'], 12)   # member conservation
        for c in quals:
            # Each emitted subfamily is internally tight (real identity, clears the floor).
            self.assertGreaterEqual(c.mean_identity,
                                    self.cfg.subfamily_intra_homogeneity_floor)

    def test_homogeneous_realistic_stays_one(self):
        # One ~8%-diverged subfamily alone (no second lineage) must stay ONE cluster.
        anc = _ancestral(seed=21)
        copies = [_copy(_mutate(anc, 0.08, 300 + v), f"H{v}", div=0.08) for v in range(8)]
        clusters = cluster_copies(copies, self.cfg)
        quals, audit = qualifying_clusters(clusters, copies, self.cfg)
        self.assertEqual(len(quals), 1)
        self.assertEqual(audit['total_members_out'], 8)


if __name__ == '__main__':
    unittest.main(verbosity=2)
