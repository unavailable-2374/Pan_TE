"""
N6 metric-harness unit tests (phase1_metrics, REFINE_STRATEGY_DESIGN_v2.md §6).

Covers the harness's deterministic pieces WITHOUT a blastn run:
  * parse_bed_truth: BED -> {mdl_R<N>: [TrueMember]} with R= normalization + divergence.
  * family_key: _sf / _chimfrag suffix stripping back to the parent family id.
  * _member_recovered: coverage-thresholded overlap recovery.
  * compute_member_recovery: MULTI vs SINGLE shape, with blast monkeypatched so the
    recovery math is exercised deterministically (a constructed family where the
    divergent member is recovered only by the second subfamily consensus -> multi lifts
    over single).
  * subfamily_count_vs_divergence, fp_side_metric, boundary_metric.

Run:  cd bin/Refiner_mdl && python3 -m unittest tests.test_n6_metrics
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import phase1_metrics as M                       # noqa: E402
from phase1_metrics import Hit, TrueMember        # noqa: E402


class TestBedTruth(unittest.TestCase):
    def test_parse_and_normalize(self):
        tmp = tempfile.mkdtemp(prefix='n6m_bed_')
        bed = os.path.join(tmp, 'inst.bed')
        with open(bed, 'w') as fh:
            fh.write("chr4\t100\t200\tR=8\t980\t+\n")     # div 0.02
            fh.write("chr4\t500\t700\tR=8\t324\t-\n")      # div 0.676
            fh.write("chr4\t10\t60\tR=150\t900\t+\n")
        truth = M.parse_bed_truth(bed)
        self.assertIn('mdl_R8', truth)
        self.assertIn('mdl_R150', truth)
        self.assertEqual(len(truth['mdl_R8']), 2)
        # divergence = 1 - score/1000
        divs = sorted(round(m.divergence, 3) for m in truth['mdl_R8'])
        self.assertEqual(divs, [0.02, 0.676])

    def test_missing_file_returns_empty(self):
        self.assertEqual(M.parse_bed_truth('/no/such/file.bed'), {})


class TestFamilyKey(unittest.TestCase):
    def test_plain_id_unchanged(self):
        self.assertEqual(M.family_key('mdl_R8'), 'mdl_R8')

    def test_subfamily_suffix_stripped(self):
        self.assertEqual(M.family_key('mdl_R8_sf2'), 'mdl_R8')

    def test_chimera_suffix_stripped(self):
        self.assertEqual(M.family_key('mdl_R8_chimfrag1'), 'mdl_R8')

    def test_composed_suffixes_stripped(self):
        self.assertEqual(M.family_key('mdl_R8_sf2_chimfrag1'), 'mdl_R8')


class TestMemberRecovered(unittest.TestCase):
    def test_full_overlap_recovered(self):
        m = TrueMember('chr4', 100, 200, 0.1)
        hits = [Hit('chr4', 90, 210, 95.0)]
        self.assertTrue(M._member_recovered(m, hits, 0.5))

    def test_no_overlap_not_recovered(self):
        m = TrueMember('chr4', 100, 200, 0.1)
        hits = [Hit('chr4', 500, 600, 95.0)]
        self.assertFalse(M._member_recovered(m, hits, 0.5))

    def test_wrong_chrom_not_recovered(self):
        m = TrueMember('chr4', 100, 200, 0.1)
        hits = [Hit('chr1', 90, 210, 95.0)]
        self.assertFalse(M._member_recovered(m, hits, 0.5))

    def test_partial_below_threshold_not_recovered(self):
        m = TrueMember('chr4', 100, 200, 0.1)        # len 100
        hits = [Hit('chr4', 100, 130, 95.0)]         # 30 bp overlap = 0.30 < 0.5
        self.assertFalse(M._member_recovered(m, hits, 0.5))


class TestComputeMemberRecovery(unittest.TestCase):
    """Constructed family where the divergent member is recovered ONLY by the second
    subfamily consensus, so MULTI recovers it and SINGLE (longest rep = the first,
    non-divergent consensus) does not — the v2 net-value direction, on synthetic input
    with blast monkeypatched (no real blastn)."""

    def setUp(self):
        self._orig_blast = M.blast_consensi_to_genome

    def tearDown(self):
        M.blast_consensi_to_genome = self._orig_blast

    def test_multi_lifts_divergent_tail(self):
        # Family mdl_R8 split into two subfamilies. sf1 is the longer (single rep).
        records = [
            {'id': 'mdl_R8_sf1', 'sequence': 'A' * 300,
             'subfamily': {'intra_identity': 0.9}},
            {'id': 'mdl_R8_sf2', 'sequence': 'A' * 100,
             'subfamily': {'intra_identity': 0.85}},
        ]
        # Two true members: one young (recovered by sf1), one divergent (only sf2 hits).
        truth = {'mdl_R8': [
            TrueMember('chr4', 100, 200, 0.03),    # young member
            TrueMember('chr4', 5000, 5100, 0.22),  # divergent member
        ]}
        # Monkeypatch blast: sf1 hits only the young member's locus; sf2 hits only the
        # divergent member's locus.
        def fake_blast(consensi, config, pident_min=70.0, wall_s=3600):
            return {
                'mdl_R8_sf1': [Hit('chr4', 90, 210, 96.0)],
                'mdl_R8_sf2': [Hit('chr4', 4990, 5110, 80.0)],
            }
        M.blast_consensi_to_genome = fake_blast

        res = M.compute_member_recovery(records, truth, config=object(),
                                        min_cov=0.5, divergent_band=(0.12, 0.50))
        fam = res['families'][0]
        # multi recovers BOTH (1.0); single (longest = sf1) recovers only the young (0.5).
        self.assertAlmostEqual(fam['recovery_multi'], 1.0)
        self.assertAlmostEqual(fam['recovery_single'], 0.5)
        # divergent member: multi recovers it, single does not.
        self.assertAlmostEqual(fam['divergent_recovery_multi'], 1.0)
        self.assertAlmostEqual(fam['divergent_recovery_single'], 0.0)
        # aggregate divergent-member recovery reflects the lift.
        dm = res['divergent_member_recovery']
        self.assertEqual(dm['n_divergent_members'], 1)
        self.assertEqual(dm['recovered_multi'], 1)
        self.assertEqual(dm['recovered_single'], 0)

    def test_homogeneous_family_multi_equals_single(self):
        # One-cluster family: multi and single are the SAME record -> identical recovery
        # (degrade-to-today invariant for the metric).
        records = [{'id': 'mdl_R99', 'sequence': 'C' * 200}]
        truth = {'mdl_R99': [TrueMember('chr4', 100, 200, 0.05)]}

        def fake_blast(consensi, config, pident_min=70.0, wall_s=3600):
            return {'mdl_R99': [Hit('chr4', 90, 210, 97.0)]}
        M.blast_consensi_to_genome = fake_blast

        res = M.compute_member_recovery(records, truth, config=object(), min_cov=0.5)
        fam = res['families'][0]
        self.assertEqual(fam['recovery_multi'], fam['recovery_single'])
        self.assertEqual(res['n_split_families'], 0)


class TestSubfamilyCountVsDivergence(unittest.TestCase):
    def test_proxy_and_count(self):
        records = [
            {'id': 'mdl_R8_sf1', 'subfamily': {'intra_identity': 0.9}},
            {'id': 'mdl_R8_sf2', 'subfamily': {'intra_identity': 0.8}},
            {'id': 'mdl_R99'},   # single-cluster, no lineage
        ]
        rows = M.subfamily_count_vs_divergence(records)
        by_fam = {r['family']: r for r in rows}
        self.assertEqual(by_fam['mdl_R8']['n_subfamilies'], 2)
        # proxy = 1 - mean(0.9, 0.8) = 0.15
        self.assertAlmostEqual(by_fam['mdl_R8']['intra_divergence_proxy'], 0.15)
        self.assertEqual(by_fam['mdl_R99']['n_subfamilies'], 1)
        self.assertEqual(by_fam['mdl_R99']['intra_divergence_proxy'], 0.0)


class TestFpSideMetric(unittest.TestCase):
    def test_drop_rate_and_reasons(self):
        pruned = [
            {'id': 'mdl_R1', 'reason': 'PSEUDO: incoherent grab-bag', 'metrics': {}},
            {'id': 'mdl_R2', 'reason': 'PSEUDO: incoherent grab-bag', 'metrics': {}},
        ]
        res = M.fp_side_metric(pruned, n_input_families=100)
        self.assertEqual(res['n_dropped'], 2)
        self.assertAlmostEqual(res['drop_rate'], 0.02)
        self.assertEqual(res['by_reason']['PSEUDO'], 2)
        self.assertEqual(len(res['spot_check_sample']), 2)


class TestBoundaryMetric(unittest.TestCase):
    def test_refined_vs_original(self):
        records = [
            {'id': 'mdl_R8', 'sequence': 'A' * 320, 'consensus_source': 'msa'},
            {'id': 'mdl_R9', 'sequence': 'A' * 100, 'consensus_source': 'original'},
        ]
        original_lengths = {'mdl_R8': 300, 'mdl_R9': 100}
        res = M.boundary_metric(records, original_lengths)
        self.assertEqual(res['n_refined_source'], 1)
        self.assertEqual(res['n_original_source'], 1)
        self.assertEqual(res['length_delta']['n'], 2)
        # deltas: +20 (R8), 0 (R9)
        self.assertEqual(sorted([res['length_delta']['min'],
                                 res['length_delta']['max']]), [0, 20])


if __name__ == '__main__':
    unittest.main()
