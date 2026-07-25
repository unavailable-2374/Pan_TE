import gzip
import io
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

from paf_processor import PAFProcessor, PAFRecord  # noqa: E402
from vcf_processor import VCFFile, VCFProcessor  # noqa: E402


class PAFProcessorTests(unittest.TestCase):
    def make_processor(self, output_dir):
        return PAFProcessor("unused.paf", "unused.fa", output_dir)

    def test_reverse_cigar_uses_forward_query_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            processor = self.make_processor(tmp)
            record = PAFRecord(
                "query", 1000, 100, 400, "-",
                "target", 1000, 200, 400,
                200, 300, 60, "50M100I150M")
            candidates = processor._extract_cigar_insertions(
                {("query", "target"): [record]})

        self.assertEqual(
            [(candidate.start, candidate.end) for candidate in candidates],
            [(250, 350)])

    def test_forward_cigar_coordinate_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            processor = self.make_processor(tmp)
            record = PAFRecord(
                "query", 1000, 100, 400, "+",
                "target", 1000, 200, 400,
                200, 300, 60, "50M100I150M")
            candidates = processor._extract_cigar_insertions(
                {("query", "target"): [record]})

        self.assertEqual(
            [(candidate.start, candidate.end) for candidate in candidates],
            [(150, 250)])

    def test_cigar_tag_is_detected_without_tp_tag(self):
        with tempfile.TemporaryDirectory() as tmp:
            paf = Path(tmp) / "input.paf"
            paf.write_text(
                "q\t7000\t0\t6000\t+\tt\t7000\t0\t6000\t"
                "6000\t6000\t60\tcg:Z:6000M\n")
            processor = PAFProcessor(str(paf), "unused.fa", tmp)
            grouped, aligner = processor._parse_and_filter_paf()

        self.assertEqual(aligner, "minimap2")
        self.assertEqual(len(grouped), 1)

    def test_empty_paf_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            paf = Path(tmp) / "empty.paf"
            paf.write_text("# no records\n")
            processor = PAFProcessor(str(paf), "unused.fa", tmp)
            with self.assertRaisesRegex(ValueError, "no alignment records"):
                processor._validate_paf_format()


class VCFProcessorTests(unittest.TestCase):
    def test_more_files_than_threads_never_allocates_zero(self):
        with tempfile.TemporaryDirectory() as vcf_dir, tempfile.TemporaryDirectory() as out:
            processor = VCFProcessor(vcf_dir, out, 2)
            files = [VCFFile(Path(f"{index}.vcf"), 100) for index in range(3)]
            allocated = processor.allocate_threads(files)

        self.assertEqual([item.allocated_threads for item in allocated], [1, 1, 1])

    def test_proportional_allocation_is_positive_and_bounded(self):
        with tempfile.TemporaryDirectory() as vcf_dir, tempfile.TemporaryDirectory() as out:
            processor = VCFProcessor(vcf_dir, out, 8)
            files = [
                VCFFile(Path("large.vcf"), 100),
                VCFFile(Path("medium.vcf"), 50),
                VCFFile(Path("small.vcf"), 25),
            ]
            allocated = processor.allocate_threads(files)

        counts = [item.allocated_threads for item in allocated]
        self.assertTrue(all(count >= 1 for count in counts))
        self.assertEqual(sum(counts), 8)

    def test_discovers_plain_and_gzip_vcf(self):
        with tempfile.TemporaryDirectory() as vcf_dir, tempfile.TemporaryDirectory() as out:
            Path(vcf_dir, "a.vcf").write_text("")
            with gzip.open(Path(vcf_dir, "b.vcf.gz"), "wt") as handle:
                handle.write("")
            processor = VCFProcessor(vcf_dir, out, 1)
            names = {item.path.name for item in processor.get_vcf_files()}

        self.assertEqual(names, {"a.vcf", "b.vcf.gz"})

    def test_decoder_failure_is_propagated(self):
        class FailedProcess:
            def __init__(self):
                self.stderr = io.StringIO("decoder failed\n")

            def poll(self):
                return 7

            def wait(self):
                return 7

        with tempfile.TemporaryDirectory() as vcf_dir, tempfile.TemporaryDirectory() as out:
            Path(vcf_dir, "bad.vcf").write_text("##fileformat=VCFv4.2\n")
            processor = VCFProcessor(vcf_dir, out, 1)
            with mock.patch("vcf_processor.subprocess.Popen", return_value=FailedProcess()):
                with self.assertRaisesRegex(RuntimeError, "exit 7"):
                    processor.process_all_vcfs()


class DecodeVCFIntegrationTests(unittest.TestCase):
    HEADER = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    def run_decoder(self, vcf_path, output_dir, threads="2"):
        return subprocess.run(
            [str(BIN / "decode_gfa.pl"), "--vcf", str(vcf_path),
             "--threads", threads, "--out", str(output_dir)],
            capture_output=True, text=True, check=False)

    def test_explicit_symbolic_breakend_and_gzip(self):
        insertion = "C" * 60
        records = (
            f"chr1\t10\t.\tA\tA{insertion}\t.\tPASS\t.\n"
            f"chr1\t20\t.\tA\t<INS>\t.\tPASS\tSVINSSEQ={insertion}\n"
            f"chr1\t30\t.\tA\tA{insertion}]chr2:40]\t.\tPASS\tSVTYPE=BND\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "sample.vcf.gz"
            output = Path(tmp) / "out"
            with gzip.open(vcf, "wt") as handle:
                handle.write(self.HEADER + records)
            result = self.run_decoder(vcf, output)
            fasta = (output / "sample.fa").read_text()

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(fasta.count(">"), 3)
        self.assertIn(">SV_sample_chr1_10_R1_A1", fasta)
        self.assertIn(">SV_sample_chr1_20_R2_A1", fasta)
        self.assertIn(">SV_sample_chr1_30_R3_A1", fasta)
        self.assertNotIn("chr2:40", fasta)
        self.assertNotIn("[", fasta)
        self.assertNotIn("]", fasta)

    def test_header_only_vcf_produces_empty_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "empty.vcf"
            output = Path(tmp) / "out"
            vcf.write_text(self.HEADER)
            result = self.run_decoder(vcf, output)
            fasta = output / "empty.fa"

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(fasta.exists())
            self.assertEqual(fasta.stat().st_size, 0)

    def test_malformed_data_record_fails_parent_process(self):
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "bad.vcf"
            output = Path(tmp) / "out"
            vcf.write_text(self.HEADER + "chr1\t10\t.\tA\tACGT\n")
            result = self.run_decoder(vcf, output)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Worker process failure", result.stderr)


if __name__ == "__main__":
    unittest.main()
