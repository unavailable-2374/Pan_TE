import tempfile
import unittest
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

from gfa_processor import GfaProcessor  # noqa: E402


class GfaProcessorPggbTests(unittest.TestCase):
    def make_processor(self, tmp, gfa_ref=None):
        gfa = Path(tmp) / "multi_chr.gfa"
        gfa.write_text(
            "H\tVN:Z:1.0\n"
            f"S\tr1\t{'A' * 200}\n"
            f"S\ta1\t{'C' * 120}\n"
            f"S\tr2\t{'G' * 200}\n"
            f"S\ta2\t{'T' * 120}\n"
            "P\tCol-0_T2T#1#Chr1\tr1+\t*\n"
            "P\tAlt#1#Chr1\tr1+,a1+\t0M\n"
            "P\tCol-0_T2T#1#Chr2\tr2+\t*\n"
            "P\tAlt#1#Chr2\tr2+,a2+\t0M\n"
        )
        return GfaProcessor(str(gfa), str(Path(tmp) / "out"),
                            gfa_ref=gfa_ref)

    def test_all_reference_chromosomes_are_selected_and_reconstructed(self):
        with tempfile.TemporaryDirectory() as tmp:
            processor = self.make_processor(tmp)
            candidates, ref_bp = processor._process_pggb()
            reference = Path(tmp) / "out" / "gfa_reference.fa"
            headers = [line.strip() for line in reference.read_text().splitlines()
                       if line.startswith(">") ]

        self.assertEqual(processor._stats["ref_sample"], "Col-0_T2T#1")
        self.assertEqual(processor._stats["ref_paths"], [
            "Col-0_T2T#1#Chr1", "Col-0_T2T#1#Chr2"])
        self.assertEqual(processor._stats["ref_segments"], 2)
        self.assertEqual(ref_bp, 400)
        self.assertEqual(headers, [">Chr1", ">Chr2"])
        self.assertEqual(len(candidates), 2)
        self.assertEqual({len(record.sequence) for record in candidates}, {120})

    def test_explicit_sample_selector_selects_all_its_chromosomes(self):
        with tempfile.TemporaryDirectory() as tmp:
            processor = self.make_processor(tmp, gfa_ref="Alt")
            selected = processor._select_pggb_references({
                "Col-0_T2T#1#Chr1": [],
                "Alt#1#Chr1": [],
                "Col-0_T2T#1#Chr2": [],
                "Alt#1#Chr2": [],
            })

        self.assertEqual(selected, ["Alt#1#Chr1", "Alt#1#Chr2"])


if __name__ == "__main__":
    unittest.main()
