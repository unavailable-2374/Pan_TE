import io
import os
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

import build_cellular_db  # noqa: E402
import final_gene_filter as fgf  # noqa: E402
import preflight  # noqa: E402


class CoverageTests(unittest.TestCase):
    def test_multi_exon_hsps_merge_per_protein(self):
        # Three 30% exons of one protein: the single-HSP rule sees 0.3, merged sees 0.9.
        hsps = {"gene": [(1, 300), (351, 650), (701, 1000)]}
        self.assertAlmostEqual(fgf.merged_coverage(hsps, 1000), 0.9)

    def test_coverage_is_not_summed_across_proteins(self):
        hsps = {"a": [(1, 300)], "b": [(401, 700)]}
        self.assertAlmostEqual(fgf.merged_coverage(hsps, 1000), 0.3)

    def test_overlaps_are_not_double_counted(self):
        self.assertAlmostEqual(fgf.merged_coverage({"a": [(1, 500), (251, 750)]}, 1000), 0.75)


class ExcisionTests(unittest.TestCase):
    def test_internal_hit_leaves_two_segments_with_margin(self):
        self.assertEqual(fgf.excise_segments(1000, {"a": [(401, 600)]}), [(1, 370), (631, 1000)])

    def test_short_segments_are_discarded(self):
        # Left 69 bp and right 70 bp are both under the 100 bp floor.
        self.assertEqual(fgf.excise_segments(1000, {"a": [(100, 900)]}), [])
        self.assertEqual(fgf.excise_segments(1000, {"a": [(150, 900)]}), [(1, 119)])

    def test_terminal_hit_leaves_one_segment(self):
        self.assertEqual(fgf.excise_segments(1000, {"a": [(1, 200)]}), [(231, 1000)])

    def test_segment_id_keeps_class(self):
        self.assertEqual(fgf.segment_id("TE_1#LTR/Gypsy", 2), "TE_1_s2#LTR/Gypsy")
        self.assertEqual(fgf.segment_id("TE_1", 1), "TE_1_s1")


class RuleTests(unittest.TestCase):
    def rec(self, rid, n=1000):
        return {"id": rid, "sequence": "A" * n}

    def test_thresholds_depend_on_te_evidence(self):
        records = [self.rec("noev#X"), self.rec("ev#X"), self.rec("evhigh#X"), self.rec("clean#X")]
        signals = {"ev#X": ["ltr"], "evhigh#X": ["ltr"]}
        hsps = {"noev#X": {"p": [(1, 100)]}, "ev#X": {"p": [(401, 600)]}, "evhigh#X": {"p": [(1, 500)]}}
        library, dropped, excised = fgf.apply_rules(records, signals, set(), hsps)
        self.assertEqual({r[0] for r in dropped}, {"noev#X", "evhigh#X"})
        self.assertEqual([e[0] for e in library], ["ev_s1#X", "ev_s2#X", "clean#X"])
        self.assertEqual(excised, [("ev#X", 1000, 2, 740)])

    def test_single_hsp_drop_is_kept(self):
        _, dropped, _ = fgf.apply_rules([self.rec("a#X")], {"a#X": ["ltr"]}, {"a#X"}, {})
        self.assertEqual(dropped, [("a#X", 0.0, "single_hsp")])


class NoSkipPathTests(unittest.TestCase):
    def test_missing_db_aborts(self):
        with tempfile.TemporaryDirectory() as d:
            lib = Path(d, "TEs.fa")
            lib.write_text(">a\nACGT\n")
            with self.assertRaises(FileNotFoundError):
                fgf.filter_library(str(lib), str(Path(d, "missing")), d, 1)
            self.assertFalse(Path(d, "TEs.fa.tmp").exists())

    def test_failed_blastx_raises(self):
        with tempfile.TemporaryDirectory() as d:
            exe = Path(d, "blastx")
            exe.write_text("#!/bin/sh\necho boom >&2\nexit 3\n")
            exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
            with self.assertRaisesRegex(RuntimeError, "blastx failed"):
                fgf.run_hsp_blastx([{"id": "a", "sequence": "ACGT"}], "db", 1, d, str(exe))

    def test_exclusion_blastx_failure_raises_instead_of_keeping_all(self):
        from Refiner_mdl.phase2_gene_exclusion import _blastx_cellular
        with tempfile.TemporaryDirectory() as d:
            fa = Path(d, "q.fa")
            fa.write_text(">a\nACGT\n")
            exe = Path(d, "blastx")
            exe.write_text("#!/bin/sh\nexit 3\n")
            exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
            cfg = mock.Mock(gene_excl_jobs=1, threads=1, blastx_exe=str(exe), gene_excl_evalue=1e-5,
                            gene_excl_max_target_seqs=3, gene_excl_min_pident=35.0)
            with self.assertRaises(RuntimeError):
                _blastx_cellular(str(fa), "db", cfg, d)

    def test_preflight_blocks_when_db_missing(self):
        report = preflight.Report()
        with mock.patch.dict(os.environ, {"PAN_TE_CELLULAR_DB": "/nonexistent/db"}):
            preflight._check_cellular_db(report, {"combine"})
        self.assertFalse(report.ok)
        self.assertIn("pan_te-setup-data --cellular-db", report.failures[0])

    def test_preflight_ignores_db_when_combine_is_done(self):
        report = preflight.Report()
        with mock.patch.dict(os.environ, {"PAN_TE_CELLULAR_DB": "/nonexistent/db"}):
            preflight._check_cellular_db(report, set())
        self.assertTrue(report.ok)

    def test_build_mdl_has_no_private_path_or_silent_filter(self):
        text = (BIN / "build_mdl").read_text()
        self.assertNotIn("/scratch/shuoc", text)
        self.assertNotIn("--cellular-protein-db", text)
        self.assertNotIn("--enable-flank-trim", text)


class CellularDbTests(unittest.TestCase):
    def test_te_product_names_are_removed(self):
        self.assertTrue(build_cellular_db.names_te_product(">sp|P1|X Transposase foo OS=Bar OX=1"))
        self.assertFalse(build_cellular_db.names_te_product(">sp|P2|Y F-box protein OS=Copia virus"))

    def test_only_te_families_are_taken_from_dfam(self):
        embl = """ID   DF1; SV 1; linear; DNA; STD; UNC; 5 BP.
CC   RepeatMasker Annotations:
CC        Type: SINE
SQ   Sequence 5 BP;
     acgta                                                         5
//
ID   DF2; SV 1; linear; DNA; STD; UNC; 5 BP.
CC        Type: Satellite
SQ   Sequence 5 BP;
     ttttt                                                         5
//
""".splitlines(True)
        out = io.StringIO()
        self.assertEqual(build_cellular_db.embl_te_consensi(embl, out), 1)
        self.assertEqual(out.getvalue(), ">DF1\nacgta\n")

    def test_te_segments_are_masked_and_the_rest_kept(self):
        self.assertEqual(build_cellular_db.mask_sequence("MKLVAGT", [(2, 3), (7, 99)]), "MXXVAGX")


if __name__ == "__main__":
    unittest.main()
