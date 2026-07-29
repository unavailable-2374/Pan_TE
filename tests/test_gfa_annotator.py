import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

from gfa_annotator import GfaAnnotator, _covered_bp, _simple_overlap  # noqa: E402


class GfaAnnotatorTests(unittest.TestCase):
    @staticmethod
    def read_table(path, key):
        with open(path, newline="") as handle:
            return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}

    @staticmethod
    def write_manifest(root, rows):
        path = root / "annotations.tsv"
        fields = [
            "sample", "haplotype", "fasta", "gff3", "te_annotation", "seqid_map",
            "path_prefix", "gaf", "annotation_complete",
        ]
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_interval_helpers(self):
        self.assertEqual(_covered_bp([(1, 6), (1, 6), (4, 9)]), 8)
        self.assertEqual(_simple_overlap("0M"), 0)
        self.assertEqual(_simple_overlap("4M3=2X"), 9)
        self.assertIsNone(_simple_overlap("4M1I3M"))
        self.assertIsNone(_simple_overlap("*"))

    def test_repeatmasker_style_tbl_reports_union_gc_and_machine_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text(
                "S\tn1\tAAAAANNNNN\n"
                "S\tn2\tCCCCCCCCCC\n"
                "S\tn3\tGGGGGGGGGG\n"
                "P\tsample#1#chr1\tn1+,n2+,n3+\t0M,0M\n")
            fasta = root / "sample.fa"
            fasta.write_text(
                ">chr1\nAAAAANNNNNCCCCCCCCCCGGGGGGGGGG\n")
            repeats = root / "sample.out"
            repeats.write_text(
                " 100 1.0 0.0 0.0 chr1 1 5 (25) + Gypsy1 LTR/Gypsy 1 5 (0) 1\n"
                " 100 1.0 0.0 0.0 chr1 4 5 (25) + Gypsy1 LTR/Gypsy 4 5 (0) 1\n"
                " 100 1.0 0.0 0.0 chr1 8 15 (15) + Mu1 DNA/MuDR 1 8 (0) 2\n"
                " 100 1.0 0.0 0.0 chr1 21 25 (5) + (AT)n Simple_repeat 1 5 (0) 3\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1", "fasta": fasta.name,
                "te_annotation": repeats.name, "annotation_complete": "true",
            }])
            annotator = GfaAnnotator(
                str(graph), str(root / "anno"), str(manifest))
            annotator.annotate()

            with open(annotator.repeat_summary_path) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            indexed = {(row["scope"], row["category_key"]): row for row in rows}
            graph_rows = {key[1]: value for key, value in indexed.items()
                          if key[0] == "graph"}
            sample_rows = {row["category_key"]: row for row in rows
                           if row["scope"] == "sample" and row["sample"] == "sample"}

            self.assertEqual(graph_rows["interspersed"]["occupied_bp"], "13")
            self.assertEqual(graph_rows["masked"]["occupied_bp"], "18")
            self.assertEqual(graph_rows["ltr_gypsy"]["occupied_bp"], "5")
            self.assertEqual(graph_rows["ltr_gypsy"]["element_count"], "1")
            self.assertEqual(graph_rows["dna_mule"]["occupied_bp"], "8")
            self.assertEqual(graph_rows["simple_repeat"]["occupied_bp"], "5")
            self.assertEqual(sample_rows["masked"]["total_length"], "30")
            self.assertEqual(sample_rows["masked"]["non_nx_length"], "25")
            self.assertEqual(sample_rows["masked"]["gc_percent"], "80.000000")
            self.assertEqual(sample_rows["masked"]["denominator_kind"], "source_fasta")
            with open(annotator.repeat_table_path) as handle:
                table = handle.read()
            self.assertIn("bases masked:", table)
            self.assertIn("Gypsy/DIRS1", table)
            self.assertIn("summary, not a RepeatMasker run", table)

    def test_graph_tbl_counts_a_shared_node_once_across_sample_classes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text(
                "S\tn1\tAAAAAAAAAA\n"
                "P\tsampleA#1#chr1\tn1+\t*\n"
                "P\tsampleB#1#chr1\tn1+\t*\n")
            repeat_a = root / "a.out"
            repeat_b = root / "b.out"
            repeat_a.write_text(
                " 100 1.0 0.0 0.0 chr1 1 10 (0) + Gypsy1 LTR/Gypsy 1 10 (0) 1\n")
            repeat_b.write_text(
                " 100 1.0 0.0 0.0 chr1 1 10 (0) + Mu1 DNA/MuDR 1 10 (0) 1\n")
            manifest = self.write_manifest(root, [
                {"sample": "sampleA", "haplotype": "1",
                 "te_annotation": repeat_a.name, "annotation_complete": "true"},
                {"sample": "sampleB", "haplotype": "1",
                 "te_annotation": repeat_b.name, "annotation_complete": "true"},
            ])
            annotator = GfaAnnotator(
                str(graph), str(root / "anno"), str(manifest))
            annotator.annotate()

            with open(annotator.repeat_summary_path) as handle:
                graph_rows = {row["category_key"]: row
                              for row in csv.DictReader(handle, delimiter="\t")
                              if row["scope"] == "graph"}
            self.assertEqual(graph_rows["interspersed"]["occupied_bp"], "10")
            self.assertEqual(graph_rows["masked"]["occupied_bp"], "10")
            self.assertEqual(graph_rows["retroelements"]["occupied_bp"], "10")
            self.assertEqual(graph_rows["dna"]["occupied_bp"], "10")
            self.assertEqual(graph_rows["retroelements"]["element_count"], "1")
            self.assertEqual(graph_rows["dna"]["element_count"], "1")

    def test_pggb_paths_project_all_genome_te_and_gene_annotations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text(
                "H\tVN:Z:1.0\n"
                "S\tn1\tAAAAAAAAAA\n"
                "S\tn2\tCCCCCCCCCC\n"
                "S\tn3\tGGGGGGGGGG\n"
                "P\tsampleA#1#chr1\tn1+,n2+\t0M\n"
                "P\tsampleB#1#chr1\tn2-,n3+\t0M\n")
            gff_a = root / "sampleA.gff3"
            gff_a.write_text(
                "##gff-version 3\n"
                "chr1\tcurated\ttransposable_element\t1\t10\t.\t+\t.\t"
                "ID=teA;Name=Athila;Classification=LTR/Gypsy\n"
                "chr1\tgenes\tgene\t11\t20\t.\t+\t.\tID=geneA\n"
                "chr1\tgenes\tmRNA\t11\t20\t.\t+\t.\tID=txA;Parent=geneA\n"
                "chr1\tgenes\texon\t11\t13\t.\t+\t.\tID=ex1;Parent=txA\n"
                "chr1\tgenes\tCDS\t11\t13\t.\t+\t0\tID=cds1;Parent=txA\n"
                "chr1\tgenes\texon\t18\t20\t.\t+\t.\tID=ex2;Parent=txA\n")
            gff_b = root / "sampleB.gff3"
            gff_b.write_text(
                "##gff-version 3\n"
                "chr1\tcurated\ttransposable_element\t1\t10\t.\t+\t.\t"
                "ID=teB;Name=Angela;Classification=LTR/Copia\n")
            manifest = self.write_manifest(root, [
                {"sample": "sampleA", "haplotype": "1", "gff3": gff_a.name,
                 "annotation_complete": "true"},
                {"sample": "sampleB", "haplotype": "1", "gff3": gff_b.name,
                 "annotation_complete": "true"},
            ])
            annotator = GfaAnnotator(str(graph), str(root / "anno"), str(manifest))
            result = annotator.annotate()

            nodes = self.read_table(result, "node_id")
            with open(annotator.interval_table_path) as handle:
                intervals = list(csv.DictReader(handle, delimiter="\t"))
            with open(annotator.provenance_path) as handle:
                provenance = json.load(handle)

            self.assertEqual(nodes["n1"]["identity"], "TE")
            self.assertEqual(nodes["n1"]["sample:sampleA"], "TE")
            self.assertEqual(nodes["n1"]["sample:sampleB"], "Absent")
            self.assertEqual(nodes["n2"]["identity"], "Mixed")
            self.assertEqual(nodes["n2"]["sample:sampleA"], "NonTE")
            self.assertEqual(nodes["n2"]["sample:sampleB"], "TE")
            self.assertIn("intron", nodes["n2"]["gene_structure"])
            self.assertEqual(nodes["n3"]["identity"], "NonTE")
            reverse_te = next(row for row in intervals if row["feature_id"] == "teB")
            self.assertEqual((reverse_te["node_start"], reverse_te["node_end"]), ("0", "10"))
            self.assertEqual(reverse_te["strand"], "-")
            derived = next(row for row in intervals if row["feature_type"] == "intron")
            self.assertEqual((derived["node_start"], derived["node_end"]), ("3", "7"))
            self.assertEqual(provenance["contract"], "pan-te-gfa-all-genome-annotation-v1")
            annotated = Path(annotator.annotated_gfa_path).read_text()
            self.assertIn("S\tn1\tAAAAAAAAAA\tai:Z:TE\tit:Z:1", annotated)
            self.assertIn("fm:Z:Athila", annotated)

    def test_same_occurrence_te_cds_overlap_is_conflict(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text("S\tn1\tAAAAA\nP\tsample#1#Chr1\tn1+\t*\n")
            gff = root / "annotations.gff3"
            gff.write_text(
                "Chr1\trepeats\ttransposable_element\t1\t5\t.\t+\t.\t"
                "ID=te1;Name=TE1;Classification=DNA/TIR\n"
                "Chr1\tgenes\tCDS\t2\t4\t.\t+\t0\tID=cds1\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1", "gff3": gff.name,
                "annotation_complete": "true",
            }])
            annotator = GfaAnnotator(str(graph), str(root / "anno"), str(manifest))
            annotator.annotate()
            node = self.read_table(annotator.node_table_path, "node_id")["n1"]

            self.assertEqual(node["identity"], "Conflict")
            self.assertEqual(node["conflict"], "1")
            self.assertIn("TE_CDS_overlap\t3", Path(annotator.conflict_path).read_text())

    def test_rgfa_gaf_projects_gapless_query_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text("S\tn1\tAAAAA\nS\tn2\tCCCCC\n")
            gff = root / "sample.gff3"
            gff.write_text(
                "chr1\tcurated\ttransposable_element\t6\t10\t.\t+\t.\t"
                "ID=te2;Name=TE2;Classification=DNA/TIR\n")
            gaf = root / "sample.gaf"
            gaf.write_text(
                "sample#1#chr1\t10\t0\t10\t+\t>n1>n2\t10\t0\t10\t"
                "10\t10\t60\ttp:A:P\tcg:Z:10=\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1", "gff3": gff.name,
                "annotation_complete": "true",
            }])
            annotator = GfaAnnotator(
                str(graph), str(root / "anno"), str(manifest), gfa_gaf=str(gaf))
            annotator.annotate()
            nodes = self.read_table(annotator.node_table_path, "node_id")

            self.assertEqual(nodes["n1"]["identity"], "NonTE")
            self.assertEqual(nodes["n2"]["identity"], "TE")
            self.assertEqual(nodes["n2"]["sample:sample"], "TE")
            self.assertEqual(annotator.graph_mode, "rGFA_GAF")

    def test_minigraph_coordinate_walk_projects_across_cigar_insertion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text(
                "S\tn1\tAAAAA\tSN:Z:reference#1#Chr1\tSO:i:0\tSR:i:0\n"
                "S\tn2\tCCCCC\tSN:Z:reference#1#Chr1\tSO:i:5\tSR:i:0\n"
                "L\tn1\t+\tn2\t+\t0M\n")
            gff = root / "sample.gff3"
            gff.write_text(
                "Chr1\tcurated\ttransposable_element\t7\t12\t.\t+\t.\t"
                "ID=te;Name=TE;Classification=DNA/TIR\n")
            gaf = root / "sample.gaf"
            gaf.write_text(
                "sample#1#Chr1\t12\t0\t12\t+\t"
                ">reference#1#Chr1:0-10\t10\t0\t10\t10\t12\t60\t"
                "tp:A:P\tcg:Z:4=2I6=\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1", "gff3": gff.name,
                "annotation_complete": "true",
            }])
            annotator = GfaAnnotator(
                str(graph), str(root / "anno"), str(manifest), gfa_gaf=str(gaf),
                gfa_gaf_min_identity=0.80)
            annotator.annotate()
            nodes = self.read_table(annotator.node_table_path, "node_id")
            with open(annotator.interval_table_path) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            projected = {(row["node_id"], row["node_start"], row["node_end"])
                         for row in rows}
            self.assertEqual(projected, {("n1", "4", "5"), ("n2", "0", "5")})
            self.assertEqual(nodes["n1"]["identity"], "PartialTE")
            self.assertEqual(nodes["n2"]["identity"], "TE")
            self.assertEqual(nodes["n1"]["total_occurrences"], "1")
            self.assertTrue(all(row["genome_start"] != "4" for row in rows))

    def test_minigraph_named_source_path_and_reverse_query_are_resolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text(
                "S\tn1\tAAAAA\tSN:Z:reference#1#Chr1\tSO:i:0\tSR:i:0\n"
                "S\tn2\tCCCCC\tSN:Z:reference#1#Chr1\tSO:i:5\tSR:i:0\n"
                "L\tn1\t+\tn2\t+\t0M\n")
            gff = root / "sample.gff3"
            gff.write_text(
                "Chr1\tcurated\ttransposable_element\t1\t3\t.\t+\t.\t"
                "ID=te;Name=TE\n")
            gaf = root / "sample.gaf"
            gaf.write_text(
                "sample#1#Chr1\t6\t0\t6\t-\treference#1#Chr1\t"
                "10\t2\t8\t6\t6\t60\ttp:A:P\tcg:Z:6=\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1", "gff3": gff.name,
                "annotation_complete": "true",
            }])
            annotator = GfaAnnotator(
                str(graph), str(root / "anno"), str(manifest), gfa_gaf=str(gaf))
            annotator.annotate()
            with open(annotator.interval_table_path) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["node_id"], "n2")
            self.assertEqual((rows[0]["node_start"], rows[0]["node_end"]),
                             ("0", "3"))
            self.assertEqual(rows[0]["orientation"], "-")
            self.assertEqual(rows[0]["strand"], "-")

    def test_gaf_nonzero_match_overlap_is_included_in_coordinate_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text(
                "S\tn1\tAAAAA\nS\tn2\tAACCC\nL\tn1\t+\tn2\t+\t2M\n")
            gff = root / "sample.gff3"
            gff.write_text(
                "chr1\tcurated\ttransposable_element\t4\t6\t.\t+\t.\t"
                "ID=te;Name=TE\n")
            gaf = root / "sample.gaf"
            gaf.write_text(
                "sample#1#chr1\t8\t0\t8\t+\t>n1>n2\t8\t0\t8\t"
                "8\t8\t60\tcg:Z:8=\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1", "gff3": gff.name,
                "annotation_complete": "true",
            }])
            annotator = GfaAnnotator(
                str(graph), str(root / "anno"), str(manifest), gfa_gaf=str(gaf))
            annotator.annotate()
            with open(annotator.interval_table_path) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            projected = {(row["node_id"], row["node_start"], row["node_end"])
                         for row in rows}
            self.assertEqual(projected, {("n1", "3", "5"), ("n2", "0", "3")})

    def test_repeatmasker_pansn_seqids_and_non_te_classes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text(
                "S\tn1\tAAAAA\nS\tn2\tCCCCC\n"
                "P\tsample#1#Chr1\tn1+,n2+\t0M\n")
            repeat_out = root / "sample.fa.out"
            repeat_out.write_text(
                "  100  1.0  0.0  0.0  sample#1#Chr1  1  5  (5)  +  "
                "TEfam  DNA/TIR  1  5  (0)  1\n"
                "   20  0.0  0.0  0.0  sample#1#Chr1  6  10 (0)  +  "
                "(AT)n  Simple_repeat  1  5  (0)  2\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1",
                "te_annotation": repeat_out.name, "annotation_complete": "true",
            }])
            annotator = GfaAnnotator(str(graph), str(root / "anno"), str(manifest))
            annotator.annotate()
            nodes = self.read_table(annotator.node_table_path, "node_id")

            self.assertEqual(nodes["n1"]["identity"], "TE")
            self.assertEqual(nodes["n2"]["identity"], "NonTE")
            self.assertEqual(nodes["n2"]["feature_types"], "Simple_repeat")

    def test_rgfa_origin_and_incomplete_annotation_remain_conservative(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text(
                "S\tr1\tAAAA\tSN:Z:chr1\tSO:i:0\n"
                "S\tr2\tCCCC\tSN:Z:chr1\tSO:i:4\n")
            gff = root / "sample.gff3"
            gff.write_text(
                "chr1\tcurated\ttransposable_element\t2\t3\t.\t+\t.\t"
                "ID=te1;Name=TE1\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1", "gff3": gff.name,
                "annotation_complete": "false",
            }])
            annotator = GfaAnnotator(str(graph), str(root / "anno"), str(manifest))
            annotator.annotate()
            nodes = self.read_table(annotator.node_table_path, "node_id")

            self.assertEqual(nodes["r1"]["identity"], "PartialTE")
            self.assertEqual(nodes["r1"]["te_covered_bp"], "2")
            self.assertEqual(nodes["r2"]["identity"], "Unknown")
            self.assertEqual(nodes["r2"]["is_te"], "NA")
            self.assertEqual(nodes["r2"]["sample:sample"], "Unknown")

    def test_explicit_seqid_map_connects_gene_gff_to_graph_contig(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text("S\tn1\tAAAAA\nP\tsample#1#Chr1\tn1+\t*\n")
            gff = root / "genes.gff3"
            gff.write_text("NC_000001\tgenes\tCDS\t1\t5\t.\t+\t0\tID=cds1\n")
            seqid_map = root / "seqids.tsv"
            seqid_map.write_text("annotation_seqid\tgraph_seqid\nNC_000001\tChr1\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1", "gff3": gff.name,
                "seqid_map": seqid_map.name, "annotation_complete": "true",
            }])
            annotator = GfaAnnotator(str(graph), str(root / "anno"), str(manifest))
            annotator.annotate()
            node = self.read_table(annotator.node_table_path, "node_id")["n1"]

            self.assertEqual(node["identity"], "GeneStructure")
            self.assertEqual(node["gene_structure"], "CDS")

    def test_gene_only_in_incomplete_annotation_is_not_te_negative(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text("S\tn1\tAAAAA\nP\tsample#1#Chr1\tn1+\t*\n")
            gff = root / "genes.gff3"
            gff.write_text("Chr1\tgenes\tCDS\t1\t5\t.\t+\t0\tID=cds1\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1", "gff3": gff.name,
                "annotation_complete": "false",
            }])
            annotator = GfaAnnotator(str(graph), str(root / "anno"), str(manifest))
            annotator.annotate()
            node = self.read_table(annotator.node_table_path, "node_id")["n1"]
            occurrence = self.read_table(
                annotator.occurrence_table_path, "occurrence_id")
            occurrence = next(iter(occurrence.values()))

            self.assertEqual(node["identity"], "GeneStructure")
            self.assertEqual(node["is_te"], "NA")
            self.assertEqual(node["sample:sample"], "Unknown")
            self.assertEqual(occurrence["annotation_callable"], "0")
            self.assertEqual(occurrence["identity"], "Unknown")

    def test_complex_p_overlap_is_audited_not_approximated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text(
                "S\tn1\tAAAAAAAAAA\nS\tn2\tCCCCCCCCCC\n"
                "P\tsample#1#chr1\tn1+,n2+\t4M1I3M\n")
            gff = root / "sample.gff3"
            gff.write_text("chr1\tgenes\tgene\t1\t20\t.\t+\t.\tID=g1\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1", "gff3": gff.name,
                "annotation_complete": "true",
            }])
            annotator = GfaAnnotator(str(graph), str(root / "anno"), str(manifest))
            annotator.annotate()

            self.assertEqual(annotator.occurrences, [])
            unresolved = Path(annotator.unresolved_path).read_text()
            self.assertIn("complex_or_missing_overlap", unresolved)

    def test_absence_requires_a_complete_embedded_path_for_same_contig(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text(
                "S\tn1\tAAAAA\nS\tn2\tCCCCC\n"
                "P\tsampleA#1#Chr1\tn1+\t*\n"
                "P\tsampleB#1#Chr2\tn2+\t*\n")
            gff_a = root / "a.gff3"
            gff_b = root / "b.gff3"
            gff_a.write_text("Chr1\tgenes\tgene\t1\t5\t.\t+\t.\tID=a\n")
            gff_b.write_text("Chr2\tgenes\tgene\t1\t5\t.\t+\t.\tID=b\n")
            manifest = self.write_manifest(root, [
                {"sample": "sampleA", "haplotype": "1", "gff3": gff_a.name,
                 "annotation_complete": "true"},
                {"sample": "sampleB", "haplotype": "1", "gff3": gff_b.name,
                 "annotation_complete": "true"},
            ])
            annotator = GfaAnnotator(str(graph), str(root / "anno"), str(manifest))
            annotator.annotate()
            nodes = self.read_table(annotator.node_table_path, "node_id")

            self.assertEqual(nodes["n1"]["sample:sampleB"], "NA")
            self.assertEqual(nodes["n2"]["sample:sampleA"], "NA")

    def test_pan_te_anno_only_skips_discovery_and_writes_new_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph = root / "graph.gfa"
            graph.write_text("S\tn1\tAAAAA\nP\tsample#1#Chr1\tn1+\t*\n")
            repeat_out = root / "sample.out"
            repeat_out.write_text(
                "100 1.0 0.0 0.0 sample#1#Chr1 1 5 (0) + TE DNA/TIR 1 5 (0) 1\n")
            manifest = self.write_manifest(root, [{
                "sample": "sample", "haplotype": "1",
                "te_annotation": repeat_out.name, "annotation_complete": "true",
            }])
            output = root / "run"
            completed = subprocess.run([
                sys.executable, str(BIN / "Pan_TE"), "--gfa", str(graph),
                "--anno-only", "--anno-manifest", str(manifest),
                "--out-dir", str(output),
            ], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            self.assertEqual(completed.returncode, 0, completed.stderr)
            with open(output / "anno" / "gfa_annotation_provenance.json") as handle:
                provenance = json.load(handle)
            self.assertEqual(provenance["contract"],
                             "pan-te-gfa-all-genome-annotation-v1")
            self.assertFalse((output / "combine").joinpath("TEs.fa").exists())


if __name__ == "__main__":
    unittest.main()
