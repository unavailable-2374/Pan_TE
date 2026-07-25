import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

from gfa_annotator import GfaAnnotator, RepeatHit, _covered_bp, _simple_overlap  # noqa: E402


class GfaAnnotatorTests(unittest.TestCase):
    def make_annotator(self, tmp, gfa_text):
        gfa = Path(tmp) / "graph.gfa"
        library = Path(tmp) / "library.fa"
        output = Path(tmp) / "anno"
        gfa.write_text(gfa_text)
        library.write_text(">TEfam#LTR/Gypsy\nACGTACGTACGT\n")
        return GfaAnnotator(str(gfa), str(output), str(library), threads=2)

    def read_node_table(self, annotator):
        with open(annotator.node_table_path, newline="") as handle:
            return {row["node_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def test_interval_union_does_not_double_count_repeated_occurrences(self):
        self.assertEqual(_covered_bp([(1, 6), (1, 6), (4, 9)]), 8)

    def test_only_exact_match_cigars_are_linearized(self):
        self.assertEqual(_simple_overlap("0M"), 0)
        self.assertEqual(_simple_overlap("4M3=2X"), 9)
        self.assertIsNone(_simple_overlap("4M1I3M"))
        self.assertIsNone(_simple_overlap("*"))

    def test_p_and_w_records_build_sample_presence(self):
        gfa = (
            "H\tVN:Z:1.1\n"
            "S\tn1\tAAAAAAAAAA\n"
            "S\tn2\tCCCCCCCCCC\n"
            "S\tn3\tGGGGGGGG\n"
            "L\tn1\t+\tn2\t+\t0M\n"
            "L\tn2\t+\tn3\t-\t0M\n"
            "P\tsampleA#1#chr1\tn1+,n2+\t0M\n"
            "W\tsampleB\t1\tchr1\t0\t18\t>n2<n3\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            annotator = self.make_annotator(tmp, gfa)
            annotator.parse_gfa()
            annotator.build_contexts()

        self.assertEqual(annotator.samples, {"sampleA", "sampleB"})
        self.assertEqual(annotator.presence["n1"], {"sampleA"})
        self.assertEqual(annotator.presence["n2"], {"sampleA", "sampleB"})
        self.assertEqual(annotator.presence["n3"], {"sampleB"})

    def test_partial_hits_are_projected_to_each_node(self):
        gfa = (
            "S\tn1\tAAAAAAAAAA\n"
            "S\tn2\tCCCCCCCCCC\n"
            "P\tsampleA#1#chr1\tn1+,n2+\t0M\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            annotator = self.make_annotator(tmp, gfa)
            annotator.parse_gfa()
            annotator.build_contexts()
            hit = RepeatHit(
                100, 1.0, 0.0, 0.0, "path_00000001", 5, 15,
                "+", "TEfam", "LTR/Gypsy")
            projected = annotator.project_hits([hit])
            annotator.write_outputs(projected)
            rows = self.read_node_table(annotator)

        self.assertEqual((projected[0]["node_start"], projected[0]["node_end"]), (5, 10))
        self.assertEqual((projected[1]["node_start"], projected[1]["node_end"]), (0, 5))
        self.assertEqual(rows["n1"]["annotation_status"], "PartialTE")
        self.assertEqual(rows["n1"]["te_covered_bp"], "5")
        self.assertEqual(rows["n1"]["te_uncovered_bp"], "5")
        self.assertEqual(rows["n1"]["te_fraction"], "0.500000")
        self.assertEqual(rows["n1"]["te_annotations"], "TEfam#LTR/Gypsy")

    def test_reverse_walk_projects_to_forward_node_coordinates(self):
        gfa = (
            "S\tn1\tAAAAAAAAAA\n"
            "S\tn2\tCCCCCCCCCC\n"
            "W\tsampleB\t1\tchr1\t0\t20\t>n1<n2\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            annotator = self.make_annotator(tmp, gfa)
            annotator.parse_gfa()
            annotator.build_contexts()
            hit = RepeatHit(
                100, 1.0, 0.0, 0.0, "path_00000001", 12, 16,
                "+", "TEfam", "LTR/Gypsy")
            projected = annotator.project_hits([hit])

        self.assertEqual(len(projected), 1)
        self.assertEqual(projected[0]["node_id"], "n2")
        self.assertEqual((projected[0]["node_start"], projected[0]["node_end"]), (4, 8))

    def test_complex_overlap_splits_context_instead_of_approximating(self):
        gfa = (
            "S\tn1\tAAAAAAAAAA\n"
            "S\tn2\tCCCCCCCCCC\n"
            "P\tsampleA#1#chr1\tn1+,n2+\t4M1I3M\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            annotator = self.make_annotator(tmp, gfa)
            annotator.parse_gfa()
            annotator.build_contexts()

        self.assertEqual(len(annotator.query_sequences), 2)
        self.assertEqual(annotator._unsupported_edges, 1)

    def test_later_paths_only_emit_windows_around_novel_nodes(self):
        shared = "A" * 20000
        novel = "C" * 100
        gfa = (
            f"S\tshared\t{shared}\n"
            f"S\tnovel\t{novel}\n"
            "P\tsampleA#1#chr1\tshared+\t*\n"
            "P\tsampleB#1#chr1\tshared+,novel+\t0M\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            annotator = self.make_annotator(tmp, gfa)
            annotator.parse_gfa()
            annotator.build_contexts()

        self.assertEqual(sorted(map(len, annotator.query_sequences.values())), [10100, 20000])
        self.assertEqual(annotator.presence["shared"], {"sampleA", "sampleB"})
        novel_occurrences = [
            occurrence for occurrence in annotator.occurrences
            if occurrence.node_id == "novel"
        ]
        self.assertEqual(len(novel_occurrences), 1)
        self.assertEqual(
            (novel_occurrences[0].node_start, novel_occurrences[0].node_end),
            (0, 100))

    def test_unavailable_sequence_is_distinct_from_none(self):
        gfa = "S\tn1\t*\tLN:i:42\nP\tsampleA#1#chr1\tn1+\t*\n"
        with tempfile.TemporaryDirectory() as tmp:
            annotator = self.make_annotator(tmp, gfa)
            annotator.parse_gfa()
            annotator.build_contexts()
            annotator.write_outputs([])
            row = self.read_node_table(annotator)["n1"]

        self.assertEqual(row["length"], "42")
        self.assertEqual(row["annotation_status"], "SequenceUnavailable")
        self.assertEqual(row["te_annotations"], "None")
        self.assertEqual(row["sample:sampleA"], "1")

    def test_end_to_end_with_mock_repeatmasker_preserves_partial_coverage(self):
        gfa = (
            "S\tn1\tAAAAAAAAAA\n"
            "S\tn2\tCCCCCCCCCC\n"
            "P\tsampleA#1#chr1\tn1+,n2+\t0M\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            annotator = self.make_annotator(tmp, gfa)
            repeat_dir = Path(annotator.output_dir) / "repeatmasker"
            repeat_dir.mkdir(parents=True)
            repeat_out = repeat_dir / "gfa_annotation_context.fa.out"
            repeat_out.write_text(
                "  100  1.0  0.0  0.0  path_00000001  6  15  (5)  +  "
                "TEfam  LTR/Gypsy  1  10  (0)  1\n"
                "   20  0.0  0.0  0.0  path_00000001  1  20  (0)  +  "
                "(AT)n  Simple_repeat  1  20  (0)  2\n")
            completed = subprocess.CompletedProcess([], 0, "", "")
            with mock.patch("gfa_annotator.shutil.which", return_value="/fake/RepeatMasker"), \
                    mock.patch("gfa_annotator.subprocess.run", return_value=completed) as run:
                result = annotator.annotate()
            rows = self.read_node_table(annotator)

        self.assertEqual(result, annotator.node_table_path)
        self.assertEqual(rows["n1"]["annotation_status"], "PartialTE")
        self.assertEqual(rows["n2"]["annotation_status"], "PartialTE")
        command = run.call_args.args[0]
        self.assertIn("-lib", command)
        self.assertIn("-s", command)
        self.assertNotIn("-nolow", command)


if __name__ == "__main__":
    unittest.main()
