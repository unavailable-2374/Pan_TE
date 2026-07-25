import csv
import tempfile
import unittest
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

from gfa_alleles import (  # noqa: E402
    build_sequence,
    deduplicate_locus_alleles,
    extract_pggb_alleles,
    parse_p_steps,
    window_starts,
)
from gfa_processor import GfaProcessor  # noqa: E402


def diverse(length, phase=0):
    motif = "ACGTGCACTGATCGTA"
    sequence = (motif * ((length // len(motif)) + 2))[phase:phase + length]
    return sequence


class PggbAlleleExtractionTests(unittest.TestCase):
    def test_repeated_reference_node_remains_inside_one_allele(self):
        sequences = {
            "left": diverse(50),
            "repeat": diverse(40, 1),
            "right": diverse(50, 2),
            "alt_a": diverse(60, 3),
            "alt_b": diverse(60, 4),
        }
        paths = {
            "Ref#1#Chr1": [
                ("left", "+"), ("repeat", "+"),
                ("repeat", "+"), ("right", "+")],
            "Alt#1#Chr1": [
                ("left", "+"), ("alt_a", "+"),
                ("repeat", "+"), ("alt_b", "+"),
                ("right", "+")],
        }

        records = extract_pggb_alleles(
            paths, ["Ref#1#Chr1"], sequences,
            min_size=20, window_size=300, window_overlap=50)

        self.assertEqual(len(records), 1)
        self.assertEqual(
            records[0].sequence,
            sequences["alt_a"] + sequences["repeat"] + sequences["alt_b"])
        self.assertEqual(records[0].ref_span_bp, 80)
        self.assertEqual(records[0].allele_type, "insertion_like")

    def test_long_allele_is_windowed_without_reset_or_loss(self):
        sequences = {
            "left": diverse(40),
            "long": diverse(650, 2),
            "right": diverse(40, 4),
        }
        paths = {
            "Ref#1#Chr1": [("left", "+"), ("right", "+")],
            "Alt#1#Chr1": [
                ("left", "+"), ("long", "+"), ("right", "+")],
        }
        stats = {}

        records = extract_pggb_alleles(
            paths, ["Ref#1#Chr1"], sequences, min_size=20,
            window_size=300, window_overlap=50, stats=stats)

        self.assertEqual(
            [(record.chunk_start, record.chunk_end) for record in records],
            [(0, 300), (250, 550), (500, 650)])
        covered = set()
        for record in records:
            covered.update(range(record.chunk_start, record.chunk_end))
        self.assertEqual(covered, set(range(650)))
        self.assertEqual(stats["long_alleles_windowed"], 1)

    def test_zero_overlap_does_not_emit_or_lose_a_short_tail(self):
        starts = window_starts(30001, 30000, 0, min_size=100)
        windows = [(start, min(start + 30000, 30001)) for start in starts]
        self.assertEqual(windows, [(0, 30000), (1, 30001)])
        covered = set()
        for start, end in windows:
            covered.update(range(start, end))
        self.assertEqual(covered, set(range(30001)))

    def test_same_locus_and_sequence_collapses_with_sample_support(self):
        insertion = diverse(160)
        sequences = {
            "left": diverse(50), "right": diverse(50, 2),
            "a": insertion, "b": insertion,
        }
        paths = {
            "Ref#1#Chr1": [("left", "+"), ("right", "+")],
            "A#1#Chr1": [("left", "+"), ("a", "+"), ("right", "+")],
            "B#1#Chr1": [("left", "+"), ("b", "+"), ("right", "+")],
        }

        occurrences = extract_pggb_alleles(
            paths, ["Ref#1#Chr1"], sequences,
            min_size=20, window_size=300, window_overlap=50)
        representatives = deduplicate_locus_alleles(occurrences)

        self.assertEqual(len(occurrences), 2)
        self.assertEqual(len(representatives), 1)
        self.assertEqual(representatives[0].occurrence_count, 2)
        self.assertEqual(representatives[0].support_sample_count, 2)
        self.assertEqual(representatives[0].support_samples, "A#1,B#1")
        self.assertEqual(
            {record.representative_id for record in occurrences},
            {representatives[0].candidate_id})
        self.assertEqual(
            {record.support_sample_count for record in occurrences}, {2})

    def test_identical_sequence_at_different_loci_preserves_copy_number(self):
        insertion = diverse(160)
        sequences = {
            "left": diverse(50), "middle": diverse(50, 2),
            "right": diverse(50, 4), "a": insertion, "b": insertion,
        }
        paths = {
            "Ref#1#Chr1": [
                ("left", "+"), ("middle", "+"), ("right", "+")],
            "Alt#1#Chr1": [
                ("left", "+"), ("a", "+"), ("middle", "+"),
                ("b", "+"), ("right", "+")],
        }

        occurrences = extract_pggb_alleles(
            paths, ["Ref#1#Chr1"], sequences,
            min_size=20, window_size=300, window_overlap=50)
        representatives = deduplicate_locus_alleles(occurrences)

        self.assertEqual(len(occurrences), 2)
        self.assertEqual(len(representatives), 2)
        self.assertNotEqual(occurrences[0].locus_key, occurrences[1].locus_key)

    def test_unchanged_repeated_terminal_sequence_is_not_emitted(self):
        sequences = {"repeat": diverse(120), "anchor": diverse(80, 2)}
        same_steps = [
            ("repeat", "+"), ("repeat", "+"), ("anchor", "+")]
        paths = {
            "Ref#1#Chr1": same_steps,
            "Alt#1#Chr1": list(same_steps),
        }
        records = extract_pggb_alleles(
            paths, ["Ref#1#Chr1"], sequences,
            min_size=20, window_size=300, window_overlap=50)
        self.assertEqual(records, [])

    def test_unchanged_path_without_eligible_anchors_is_not_emitted(self):
        sequences = {"repeat": diverse(120)}
        same_steps = [("repeat", "+")] * 5
        paths = {
            "Ref#1#Chr1": same_steps,
            "Alt#1#Chr1": list(same_steps),
        }
        records = extract_pggb_alleles(
            paths, ["Ref#1#Chr1"], sequences,
            min_size=20, window_size=300, window_overlap=50,
            max_anchor_occurrences=2)
        self.assertEqual(records, [])

    def test_maximum_collinear_chain_ignores_early_rearranged_anchor(self):
        sequences = {
            name: diverse(80, index)
            for index, name in enumerate(("A", "B", "C", "D", "E", "Z"))
        }
        paths = {
            "Ref#1#Chr1": [(name, "+") for name in "ABCDE"],
            "Alt#1#Chr1": [(name, "+") for name in "EAZBCD"],
        }
        records = extract_pggb_alleles(
            paths, ["Ref#1#Chr1"], sequences,
            min_size=20, window_size=300, window_overlap=50)
        inserted = [record for record in records
                    if record.left_anchor == "A" and record.right_anchor == "B"]
        self.assertEqual(len(inserted), 1)
        self.assertEqual(inserted[0].sequence, sequences["Z"])
        self.assertEqual(inserted[0].ref_span_bp, 0)

    def test_casefold_contig_match_and_unmatched_contig_failure(self):
        sequences = {
            "left": diverse(80), "right": diverse(80, 2),
            "alt": diverse(120, 4),
        }
        casefold_paths = {
            "Ref#1#Chr1": [("left", "+"), ("right", "+")],
            "Alt#1#chr1": [
                ("left", "+"), ("alt", "+"), ("right", "+")],
        }
        records = extract_pggb_alleles(
            casefold_paths, ["Ref#1#Chr1"], sequences,
            min_size=20, window_size=300, window_overlap=50)
        self.assertEqual(len(records), 1)

        unmatched_paths = dict(casefold_paths)
        unmatched_paths.pop("Alt#1#chr1")
        unmatched_paths["Alt#1#Chr2"] = [
            ("left", "+"), ("alt", "+"), ("right", "+")]
        with self.assertRaisesRegex(ValueError, "no matching reference contig"):
            extract_pggb_alleles(
                unmatched_paths, ["Ref#1#Chr1"], sequences,
                min_size=20, window_size=300, window_overlap=50)

    def test_nonzero_path_overlap_and_missing_segment_fail_fast(self):
        with self.assertRaisesRegex(ValueError, "requires 0M"):
            parse_p_steps(["P", "sample#1#Chr1", "a+,b+", "10M"])
        with self.assertRaisesRegex(ValueError, "missing segment"):
            build_sequence([("missing", "+")], {})

    def test_full_processor_emits_independent_records_and_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            gfa = tmp_path / "graph.gfa"
            sequence = diverse(200)
            gfa.write_text(
                "H\tVN:Z:1.0\n"
                f"S\tleft\t{diverse(100)}\n"
                f"S\tright\t{diverse(100, 3)}\n"
                f"S\ta\t{sequence}\n"
                f"S\tb\t{sequence}\n"
                "P\tRef#1#Chr1\tleft+,right+\t0M\n"
                "P\tA#1#Chr1\tleft+,a+,right+\t0M,0M\n"
                "P\tB#1#Chr1\tleft+,b+,right+\t0M,0M\n"
            )
            processor = GfaProcessor(
                str(gfa), str(tmp_path / "out"), threads=1,
                gfa_ref="Ref#1", spacer_length=0)
            output = Path(processor.process())
            scaffold = tmp_path / "out" / "gfa_scaffolds.fa"
            manifest = tmp_path / "out" / "gfa_allele_manifest.tsv"
            temporary_store = tmp_path / "out" / ".gfa_alleles.sqlite.tmp"

            scaffold_text = scaffold.read_text()
            with manifest.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertFalse(temporary_store.exists())
        self.assertEqual(output.name, "gfa_pan_genome.fa")
        self.assertEqual(scaffold_text.count(">"), 1)
        self.assertNotIn("N" * 20, scaffold_text)
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["support_sample_count"] for row in rows}, {"2"})
        self.assertEqual({row["qc_status"] for row in rows}, {"passed"})


if __name__ == "__main__":
    unittest.main()
