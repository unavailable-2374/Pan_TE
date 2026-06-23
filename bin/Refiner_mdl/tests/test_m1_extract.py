"""
M1 unit tests — phase1_extract coordinate math, RC, seqid stripping, subsampling.

Pure-function gate for the BED-seeded copy extractor. Off-by-one in either
direction must be caught here before any real TE library is built.

Run:  python3 -m unittest bin.Refiner_mdl.tests.test_m1_extract  (from repo root)
   or:  cd bin/Refiner_mdl && python3 -m unittest tests.test_m1_extract
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from phase1_extract import (  # noqa: E402
    Instance,
    bed_to_samtools_region,
    compute_pad,
    reverse_complement,
    stratified_subsample,
    strip_seqid_prefix,
)


class TestBedToSamtoolsRegion(unittest.TestCase):
    """BED 0-based half-open -> samtools 1-based inclusive, both directions."""

    def test_no_pad_core(self):
        # BED [100, 200) == 1-based bases 101..200 (length 100).
        region, s1, e1 = bed_to_samtools_region("chr1", 100, 200, 0, 1000)
        self.assertEqual(region, "chr1:101-200")
        self.assertEqual((s1, e1), (101, 200))
        self.assertEqual(e1 - s1 + 1, 100)  # core length preserved

    def test_with_pad(self):
        region, s1, e1 = bed_to_samtools_region("chr1", 100, 200, 50, 1000)
        # s1 = max(1, 100+1-50) = 51 ; e1 = min(1000, 200+50) = 250
        self.assertEqual(region, "chr1:51-250")
        self.assertEqual((s1, e1), (51, 250))

    def test_left_clamp_to_one(self):
        # pad pushes left edge below 1 -> clamp to 1.
        region, s1, e1 = bed_to_samtools_region("chr1", 10, 100, 50, 1000)
        # s1 = max(1, 10+1-50) = max(1, -39) = 1 ; e1 = min(1000, 150) = 150
        self.assertEqual((s1, e1), (1, 150))
        self.assertEqual(region, "chr1:1-150")

    def test_right_clamp_to_chrom_len(self):
        region, s1, e1 = bed_to_samtools_region("chr1", 900, 990, 50, 1000)
        # s1 = max(1, 851) = 851 ; e1 = min(1000, 1040) = 1000
        self.assertEqual((s1, e1), (851, 1000))
        self.assertEqual(region, "chr1:851-1000")

    def test_both_clamps_tiny_chrom(self):
        region, s1, e1 = bed_to_samtools_region("c", 0, 5, 100, 5)
        # s1 = max(1, 0+1-100) = 1 ; e1 = min(5, 105) = 5
        self.assertEqual((s1, e1), (1, 5))
        self.assertEqual(region, "c:1-5")

    def test_end_le_start_raises(self):
        # end < start (corrupt) -> ValueError.
        with self.assertRaises(ValueError):
            bed_to_samtools_region("chr1", 500, 400, 0, 1000)

    def test_empty_interval_raises(self):
        # zero-length BED [500, 500): s1=501, e1=500 -> e1<s1 -> ValueError.
        with self.assertRaises(ValueError):
            bed_to_samtools_region("chr1", 500, 500, 0, 1000)

    def test_start_past_chrom_end_raises(self):
        # start beyond chromosome end -> s1 > chrom_len == e1 -> ValueError.
        with self.assertRaises(ValueError):
            bed_to_samtools_region("chr1", 2000, 2100, 0, 1000)


class TestReverseComplement(unittest.TestCase):

    def test_basic_acgt(self):
        self.assertEqual(reverse_complement("AAAA"), "TTTT")
        self.assertEqual(reverse_complement("ACGT"), "ACGT")  # self-RC palindrome
        self.assertEqual(reverse_complement("AACCGGTT"), "AACCGGTT")
        self.assertEqual(reverse_complement("ATGC"), "GCAT")

    def test_round_trip(self):
        for seq in ["ACGTACGTACGT", "GGGGCCCCAAATTT", "ACGTN", "TACGGTCA"]:
            self.assertEqual(reverse_complement(reverse_complement(seq)), seq)

    def test_n_and_iupac(self):
        self.assertEqual(reverse_complement("ACGTN"), "NACGT")
        self.assertEqual(reverse_complement("R"), "Y")
        self.assertEqual(reverse_complement("Y"), "R")
        self.assertEqual(reverse_complement("W"), "W")
        self.assertEqual(reverse_complement("S"), "S")
        self.assertEqual(reverse_complement("K"), "M")
        self.assertEqual(reverse_complement("B"), "V")

    def test_lowercase_uppercased(self):
        self.assertEqual(reverse_complement("acgt"), "ACGT")

    def test_unknown_to_n(self):
        # Non-IUPAC characters collapse to N.
        self.assertEqual(reverse_complement("Z"), "N")
        self.assertEqual(reverse_complement("AZT"), "ANT")  # rev: T Z A -> A N T


class TestStripSeqidPrefix(unittest.TestCase):

    def test_ncbi_ref_wrapper(self):
        self.assertEqual(strip_seqid_prefix("ref|NC_003070.9|"), "NC_003070.9")

    def test_clean_chr_noop(self):
        self.assertEqual(strip_seqid_prefix("chr1"), "chr1")
        self.assertEqual(strip_seqid_prefix("testA"), "testA")

    def test_gi_style(self):
        self.assertEqual(strip_seqid_prefix("gi|123|ref|NC_1|"), "123")

    def test_empty_token(self):
        # Bare pipe -> no non-empty trailing token -> falls back to parts[0] ("").
        self.assertEqual(strip_seqid_prefix("|"), "")

    def test_trailing_pipe_only(self):
        self.assertEqual(strip_seqid_prefix("NC_1|"), "NC_1")


class _Cfg:
    pad_fraction = 0.2
    pad_min = 50
    pad_cap = 500


class TestComputePad(unittest.TestCase):

    def test_clamp_floor(self):
        # 0.2 * 100 = 20 < pad_min(50) -> 50
        self.assertEqual(compute_pad(100, _Cfg), 50)

    def test_middle(self):
        # 0.2 * 1000 = 200, within [50, 500]
        self.assertEqual(compute_pad(1000, _Cfg), 200)

    def test_clamp_cap(self):
        # 0.2 * 5000 = 1000 > pad_cap(500) -> 500
        self.assertEqual(compute_pad(5000, _Cfg), 500)


def _mk_instances(n):
    out = []
    for i in range(n):
        out.append(Instance(chrom="chr1", start=i * 1000,
                            end=i * 1000 + 100 + (i % 7) * 10,
                            strand="+" if i % 2 == 0 else "-",
                            divergence=(i % 11) / 10.0))
    return out


class TestStratifiedSubsample(unittest.TestCase):

    def test_below_cap_returns_original(self):
        insts = _mk_instances(8)
        out = stratified_subsample(insts, cap=20, seed=42)
        self.assertEqual(len(out), 8)
        self.assertEqual(out, insts)  # same objects, same order

    def test_equal_to_cap_returns_original(self):
        insts = _mk_instances(20)
        out = stratified_subsample(insts, cap=20, seed=42)
        self.assertEqual(len(out), 20)
        self.assertEqual(out, insts)

    def test_caps_count(self):
        insts = _mk_instances(500)
        out = stratified_subsample(insts, cap=100, seed=42)
        self.assertEqual(len(out), 100)

    def test_deterministic_same_seed(self):
        insts = _mk_instances(500)
        a = stratified_subsample(insts, cap=100, seed=42)
        b = stratified_subsample(insts, cap=100, seed=42)
        self.assertEqual([(i.chrom, i.start) for i in a],
                         [(i.chrom, i.start) for i in b])

    def test_different_seed_may_differ_but_same_size(self):
        insts = _mk_instances(500)
        a = stratified_subsample(insts, cap=100, seed=42)
        c = stratified_subsample(insts, cap=100, seed=7)
        self.assertEqual(len(a), len(c))  # both respect the cap deterministically

    def test_subset_of_input(self):
        insts = _mk_instances(300)
        out = stratified_subsample(insts, cap=50, seed=42)
        keys = {(i.chrom, i.start) for i in insts}
        for i in out:
            self.assertIn((i.chrom, i.start), keys)


if __name__ == "__main__":
    unittest.main(verbosity=2)
