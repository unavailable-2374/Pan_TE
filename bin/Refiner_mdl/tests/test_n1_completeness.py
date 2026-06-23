"""
N1 unit tests — phase1_completeness cheap signals + assess_completeness.

All inputs here are SYNTHETIC (clearly labeled), exercising the signal logic on
constructed homogeneous-vs-divergent / high-vs-low-copy families.  Real-data
validation against blastn divergent-tail ground truth lives in the separate N1
integration script (tests/integration_n1_completeness.py).

Run:  cd bin/Refiner_mdl && python3 -m unittest tests.test_n1_completeness
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig  # noqa: E402  (import before phase modules:
# phase0_triage prepends ../Refiner to sys.path at import time, which would
# otherwise shadow this local config.py with Refiner's.)
from phase1_extract import Instance  # noqa: E402
from phase1_completeness import (  # noqa: E402
    CompletenessVerdict,
    assess_completeness,
    bed_divergence_spread,
    copy_count_vs_length_class,
    under_instanced,
    load_presignals,
)


def _inst(div, chrom="chr1", start=0, end=100, strand="+"):
    """SYNTHETIC instance helper."""
    return Instance(chrom=chrom, start=start, end=end, strand=strand,
                    divergence=div)


class TestBedDivergenceSpread(unittest.TestCase):
    """Spread metric: homogeneous → low; mixed-divergence → high."""

    def test_empty_returns_none(self):
        self.assertIsNone(bed_divergence_spread([]))

    def test_single_instance_zero_spread(self):
        self.assertEqual(bed_divergence_spread([_inst(0.3)]), 0.0)

    def test_homogeneous_low_spread(self):
        # SYNTHETIC: four tightly-clustered divergences.
        insts = [_inst(0.20), _inst(0.21), _inst(0.19), _inst(0.20)]
        sp = bed_divergence_spread(insts)
        self.assertLess(sp, 0.05)

    def test_mixed_high_spread(self):
        # SYNTHETIC mirror of the real R=150 case: divs {0.02, 0.676}.
        insts = [_inst(0.02), _inst(0.676)]
        sp = bed_divergence_spread(insts)
        self.assertGreater(sp, 0.15)  # population std = 0.328

    def test_population_std_value(self):
        # Population std of {0.0, 0.4} = 0.2 exactly.
        self.assertAlmostEqual(bed_divergence_spread([_inst(0.0), _inst(0.4)]),
                               0.2, places=6)


class TestCopyCountVsLengthClass(unittest.TestCase):
    """Short low-copy element OK; long low-copy element flagged."""

    def setUp(self):
        self.cfg = RefinerMdlConfig()  # default length_class_copy_expectation

    def test_short_low_copy_not_flagged(self):
        # SYNTHETIC: 150 bp element, 2 instances. Short class expects >=2 → OK.
        rec = {'id': 'mdl_R1', 'actual_length': 150, 'n_instances': 2}
        self.assertFalse(copy_count_vs_length_class(rec, self.cfg))

    def test_long_low_copy_flagged(self):
        # SYNTHETIC: 1500 bp element, only 2 instances. Long class expects >=4.
        rec = {'id': 'mdl_R2', 'actual_length': 1500, 'n_instances': 2}
        self.assertTrue(copy_count_vs_length_class(rec, self.cfg))

    def test_long_adequate_copy_not_flagged(self):
        # SYNTHETIC: 1500 bp, 10 instances → above the >=4 expectation.
        rec = {'id': 'mdl_R3', 'actual_length': 1500, 'n_instances': 10}
        self.assertFalse(copy_count_vs_length_class(rec, self.cfg))

    def test_falls_back_to_copies_when_no_instances(self):
        # SYNTHETIC: no n_instances → use header copies.
        rec = {'id': 'mdl_R4', 'actual_length': 1500, 'copies': 2}
        self.assertTrue(copy_count_vs_length_class(rec, self.cfg))

    def test_empty_table_disables_signal(self):
        cfg = RefinerMdlConfig()
        cfg.length_class_copy_expectation = ()
        rec = {'id': 'mdl_R5', 'actual_length': 1500, 'n_instances': 1}
        self.assertFalse(copy_count_vs_length_class(rec, cfg))


class TestUnderInstanced(unittest.TestCase):
    """BED instance count below min_copies_for_msa."""

    def setUp(self):
        self.cfg = RefinerMdlConfig()  # min_copies_for_msa default = 5

    def test_under_instanced_flagged(self):
        rec = {'id': 'mdl_R1', 'n_instances': 3}
        self.assertTrue(under_instanced(rec, self.cfg))

    def test_well_instanced_not_flagged(self):
        rec = {'id': 'mdl_R2', 'n_instances': 10}
        self.assertFalse(under_instanced(rec, self.cfg))

    def test_boundary_equal_not_flagged(self):
        rec = {'id': 'mdl_R3', 'n_instances': 5}
        self.assertFalse(under_instanced(rec, self.cfg))


class TestAssessCompleteness(unittest.TestCase):
    """Combined verdict: complete family not eligible; divergent/low-copy eligible."""

    def setUp(self):
        self.cfg = RefinerMdlConfig()

    def test_complete_homogeneous_not_eligible(self):
        # SYNTHETIC clean family: long, many tight-divergence instances.
        # Mirrors real R=8 (56 inst, len 535, spread ~0.06).
        rec = {'id': 'mdl_R8', 'actual_length': 535, 'copies': 56}
        insts = [_inst(0.50 + 0.01 * (k % 3), start=k * 1000, end=k * 1000 + 535)
                 for k in range(56)]
        v = assess_completeness(rec, insts, self.cfg)
        self.assertIsInstance(v, CompletenessVerdict)
        self.assertFalse(v.incomplete_recall_eligible)
        self.assertFalse(any(v.flags.values()))
        self.assertEqual(v.n_instances, 56)

    def test_divergent_tail_family_eligible(self):
        # SYNTHETIC mirror of real R=150: 2 inst, divs {0.02, 0.676}, len 472.
        rec = {'id': 'mdl_R150', 'actual_length': 472, 'copies': 2}
        insts = [_inst(0.02), _inst(0.676, start=5000, end=5472)]
        v = assess_completeness(rec, insts, self.cfg)
        self.assertTrue(v.incomplete_recall_eligible)
        # Two distinct signals should fire: wide spread AND under-instanced.
        self.assertTrue(v.flags['bed_divergence_spread'])
        self.assertTrue(v.flags['under_instanced'])
        self.assertTrue(len(v.reasons) >= 2)

    def test_long_low_copy_eligible_via_length_class(self):
        # SYNTHETIC: 1500 bp element with 3 tight-divergence instances.
        # Spread is low and 3 >= min? no (min=5) so under_instanced also fires;
        # construct 6 instances to isolate the length-class signal.
        rec = {'id': 'mdl_R20', 'actual_length': 1500, 'copies': 3}
        insts = [_inst(0.20 + 0.005 * k, start=k * 2000, end=k * 2000 + 1500)
                 for k in range(3)]
        v = assess_completeness(rec, insts, self.cfg)
        self.assertTrue(v.incomplete_recall_eligible)
        self.assertTrue(v.flags['copy_count_vs_length_class'])

    def test_well_instanced_long_homogeneous_not_eligible(self):
        # SYNTHETIC: 1500 bp, 10 tight instances → no signal fires.
        rec = {'id': 'mdl_R21', 'actual_length': 1500, 'copies': 10}
        insts = [_inst(0.20 + 0.005 * (k % 4), start=k * 2000, end=k * 2000 + 1500)
                 for k in range(10)]
        v = assess_completeness(rec, insts, self.cfg)
        self.assertFalse(v.incomplete_recall_eligible)

    def test_zero_instances_eligible_via_under_instanced(self):
        # SYNTHETIC: record with no BED instances (filtered out) — eligible.
        rec = {'id': 'mdl_R99', 'actual_length': 300, 'copies': 4}
        v = assess_completeness(rec, [], self.cfg)
        self.assertTrue(v.incomplete_recall_eligible)
        self.assertTrue(v.flags['under_instanced'])
        self.assertIsNone(v.bed_divergence_spread)


class TestLoadPresignals(unittest.TestCase):
    """Round-trip the presignal TSV reader on a SYNTHETIC file."""

    def test_parse(self):
        import tempfile
        content = (
            "record_id\tn_instances\tbed_divergence_mean\tbed_divergence_spread\t"
            "length\tlength_class_expectation\tlength_class_flag\n"
            "mdl_R8\t56\t0.5100\t0.0600\t535\t3\t0\n"
            "mdl_R150\t2\t0.3480\t0.3280\t472\t3\t1\n"
            "mdl_R99\t0\t\t\t300\t3\t1\n"
        )
        with tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False) as fh:
            fh.write(content)
            path = fh.name
        try:
            d = load_presignals(path)
            self.assertEqual(set(d), {'mdl_R8', 'mdl_R150', 'mdl_R99'})
            self.assertEqual(d['mdl_R8']['n_instances'], 56)
            self.assertAlmostEqual(d['mdl_R150']['bed_divergence_spread'], 0.328)
            self.assertTrue(d['mdl_R150']['length_class_flag'])
            self.assertIsNone(d['mdl_R99']['bed_divergence_spread'])
        finally:
            os.unlink(path)


if __name__ == '__main__':
    unittest.main(verbosity=2)
