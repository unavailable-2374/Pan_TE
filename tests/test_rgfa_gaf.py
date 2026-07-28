import csv
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

from gfa_processor import GfaProcessor
from rgfa_gaf import (RgfaGafCallStore, RgfaGafError,
                      derive_rgfa_node_presence, extract_rgfa_gaf_alleles,
                      populate_rgfa_gaf_locus_calls)


class RgfaGafTests(unittest.TestCase):
    def make_inputs(self, root: Path, mapq: int = 60):
        graph = root / "graph.gfa"
        fasta = root / "query.fa"
        gaf = root / "query.gaf"
        graph.write_text(
            "S\tr1\tAAA\tSR:i:0\tSN:Z:Ref#1#Chr1\tSO:i:0\n"
            "S\tr2\tCCC\tSR:i:0\tSN:Z:Ref#1#Chr1\tSO:i:3\n"
            "S\tr3\tGGG\tSR:i:0\tSN:Z:Ref#1#Chr1\tSO:i:6\n"
            "S\ta1\tTTT\tSR:i:1\tSN:Z:SampleA#1#Chr1\tSO:i:3\n"
            "L\tr1\t+\tr2\t+\t0M\n"
            "L\tr2\t+\tr3\t+\t0M\n"
            "L\tr1\t+\ta1\t+\t0M\n"
            "L\ta1\t+\tr3\t+\t0M\n")
        fasta.write_text(
            ">SampleA#1#Chr1\nAAATTTGGG\n"
            ">SampleB#1#Chr1\nAAACCCGGG\n")
        path_a = ">Ref#1#Chr1:0-3>SampleA#1#Chr1:3-6>Ref#1#Chr1:6-9"
        # minigraph may encode this same source path as a bare SN name.
        path_b = "Ref#1#Chr1"
        gaf.write_text(
            f"SampleA#1#Chr1\t9\t0\t9\t+\t{path_a}\t9\t0\t9\t9\t9\t{mapq}\ttp:A:P\tcg:Z:9=\n"
            f"SampleB#1#Chr1\t9\t0\t9\t+\t{path_b}\t9\t0\t9\t9\t9\t{mapq}\ttp:A:P\tcg:Z:9=\n")
        return graph, gaf, fasta

    def test_candidate_and_outer_anchor_presence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph, gaf, fasta = self.make_inputs(root)
            calls = RgfaGafCallStore(str(root / "calls.sqlite"))
            candidates = extract_rgfa_gaf_alleles(
                str(graph), str(gaf), min_size=3, window_size=30,
                query_fasta=str(fasta), call_sink=calls.add)
            self.assertEqual(len(candidates), 1)
            candidate = candidates[0]
            self.assertEqual(candidate.sequence, "TTT")
            self.assertEqual(candidate.sample, "SampleA")
            self.assertEqual(candidate.sequence_source, "graph_complete_nodes")
            self.assertEqual(candidate.locus_key, "Chr1:0:2")
            populate_rgfa_gaf_locus_calls(
                str(graph), str(gaf), candidates, calls.add,
                min_size=3, query_fasta=str(fasta))
            matrix = root / "presence.tsv"
            calls.write_presence_matrix(str(matrix), [candidate.locus_key])
            calls.close(remove=True)
            with matrix.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["sample:SampleA"], "1")
            self.assertEqual(rows[0]["sample:SampleB"], "0")

    def test_reference_backbone_is_forced_to_reference_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph, gaf, fasta = self.make_inputs(root)
            # Deliberately give the declared reference sample a mapper route
            # through the alternate node. SR:i:0, not that route, is truth.
            fasta.write_text(fasta.read_text() + ">Ref#1#Chr1\nAAATTTGGG\n")
            gaf.write_text(gaf.read_text() + (
                "Ref#1#Chr1\t9\t0\t9\t+\t"
                ">Ref#1#Chr1:0-3>SampleA#1#Chr1:3-6>Ref#1#Chr1:6-9\t"
                "9\t0\t9\t9\t9\t60\ttp:A:P\tcg:Z:9=\n"))
            calls = RgfaGafCallStore(str(root / "calls.sqlite"))
            candidates = extract_rgfa_gaf_alleles(
                str(graph), str(gaf), min_size=3, window_size=30,
                query_fasta=str(fasta), call_sink=calls.add)
            self.assertEqual([candidate.sample for candidate in candidates], ["SampleA"])
            populate_rgfa_gaf_locus_calls(
                str(graph), str(gaf), candidates, calls.add,
                min_size=3, query_fasta=str(fasta))
            matrix = root / "presence.tsv"
            calls.write_presence_matrix(
                str(matrix), [candidates[0].locus_key], reference_samples={"Ref"})
            calls.close(remove=True)
            with matrix.open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["sample:Ref"], "0")

    def test_rgfa_default_mode_requires_gaf(self):
        with tempfile.TemporaryDirectory() as tmp:
            graph, _gaf, _fasta = self.make_inputs(Path(tmp))
            processor = GfaProcessor(str(graph), str(Path(tmp) / "out"))
            with self.assertRaisesRegex(ValueError, "requires --gfa-gaf"):
                processor.process()

    def test_low_quality_sample_is_retained_as_na_in_locus_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph, gaf, fasta = self.make_inputs(root)
            fasta.write_text(fasta.read_text() + ">SampleLow#1#Chr1\nAAACCCGGG\n")
            gaf.write_text(gaf.read_text() + (
                "SampleLow#1#Chr1\t9\t0\t9\t+\tRef#1#Chr1\t9\t0\t9\t"
                "9\t9\t5\ttp:A:P\tcg:Z:9=\n"))
            calls = RgfaGafCallStore(str(root / "calls.sqlite"))
            candidates = extract_rgfa_gaf_alleles(
                str(graph), str(gaf), min_size=3, window_size=30,
                query_fasta=str(fasta), call_sink=calls.add,
                sample_sink=calls.register_sample)
            populate_rgfa_gaf_locus_calls(
                str(graph), str(gaf), candidates, calls.add,
                min_size=3, query_fasta=str(fasta),
                sample_sink=calls.register_sample)
            matrix = root / "presence.tsv"
            calls.write_presence_matrix(str(matrix), [candidates[0].locus_key])
            calls.close(remove=True)
            with matrix.open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["sample:SampleLow"], "NA")

    def test_rgfa_failure_removes_temporary_sqlite_stores(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph, gaf, fasta = self.make_inputs(root)
            gaf.write_text("malformed\n")
            outdir = root / "out"
            processor = GfaProcessor(
                str(graph), str(outdir), gfa_gaf=str(gaf),
                gfa_gaf_fasta=str(fasta))
            with self.assertRaisesRegex(RgfaGafError, "malformed GAF"):
                processor.process()
            self.assertFalse((outdir / ".gfa_alleles.sqlite.tmp").exists())
            self.assertFalse((outdir / ".rgfa_gaf_calls.sqlite.tmp").exists())

    def test_low_quality_is_not_presence_or_absence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph, gaf, fasta = self.make_inputs(root, mapq=5)
            candidates = extract_rgfa_gaf_alleles(
                str(graph), str(gaf), min_size=3, window_size=30,
                min_mapq=20, query_fasta=str(fasta))
            self.assertEqual(candidates, [])
            samples, states, _stats = derive_rgfa_node_presence(
                str(graph), str(gaf), min_mapq=20, query_fasta=str(fasta))
            self.assertEqual(samples, {"SampleA", "SampleB"})
            self.assertEqual(states, {})

    def test_gaf_length_mismatch_fails_fast(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph, gaf, fasta = self.make_inputs(root)
            fasta.write_text(">SampleA#1#Chr1\nAAATTTGG\n>SampleB#1#Chr1\nAAACCCGGG\n")
            with self.assertRaisesRegex(RgfaGafError, "length differs"):
                extract_rgfa_gaf_alleles(
                    str(graph), str(gaf), min_size=3, window_size=30,
                    query_fasta=str(fasta))

if __name__ == "__main__":
    unittest.main()
