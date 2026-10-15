import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import drop_contained_fragments as dcf  # noqa: E402

CLSTR = """>Cluster 0
0\t9000nt, >TE_1#LTR/Gypsy... *
1\t3000nt, >TE_2#LTR/Copia... at 1:3000:10:3010/+/95.00%
2\t500nt, >TE_3#SINE/tRNA... at 1:500:20:520/+/90.00%
3\t800nt, >TE_4... at 1:800:30:830/-/85.00%
4\t700nt, >TE_5#Unknown... at 1:700:40:740/+/88.00%
>Cluster 1
0\t1200nt, >TE_6#DNA/MULE-MuDR... *
"""


class TopClassTests(unittest.TestCase):
    def test_top_class(self):
        self.assertEqual(dcf.top_class("TE_1#LTR/Gypsy"), "LTR")
        self.assertEqual(dcf.top_class("TE_1#DNA"), "DNA")
        self.assertEqual(dcf.top_class("TE_1"), "Unknown")
        self.assertEqual(dcf.top_class("TE_1#"), "Unknown")


class SelectDropsTests(unittest.TestCase):
    def test_same_class_and_unknown_dropped_cross_class_kept(self):
        clusters = list(dcf.parse_clusters(CLSTR.splitlines(True)))
        self.assertEqual(clusters[0][0], "TE_1#LTR/Gypsy")
        self.assertEqual(clusters[1], ("TE_6#DNA/MULE-MuDR", []))
        drops, kept = dcf.select_drops(clusters)
        self.assertEqual(set(drops), {"TE_2#LTR/Copia", "TE_4", "TE_5#Unknown"})
        self.assertTrue(all(rep == "TE_1#LTR/Gypsy" for rep in drops.values()))
        self.assertEqual(kept, [("TE_3#SINE/tRNA", "TE_1#LTR/Gypsy")])

    def test_filter_fasta_keeps_order_and_ignores_header_tail(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp, "in.fa"), Path(tmp, "out.fa")
            src.write_text(">TE_1#LTR/Gypsy \nACGT\nACGT\n>TE_2#LTR/Copia\nGGGG\n>TE_3#SINE/tRNA desc\nTTTT\n")
            total, kept = dcf.filter_fasta(src, dst, {"TE_2#LTR/Copia": "TE_1#LTR/Gypsy"})
            self.assertEqual((total, kept), (3, 2))
            self.assertEqual(dst.read_text(), ">TE_1#LTR/Gypsy \nACGT\nACGT\n>TE_3#SINE/tRNA desc\nTTTT\n")


@unittest.skipUnless(shutil.which("cd-hit-est"), "cd-hit-est not installed")
class EndToEndTests(unittest.TestCase):
    def test_contained_same_class_fragment_is_dropped(self):
        import random
        rnd = random.Random(7)
        long_seq = "".join(rnd.choice("ACGT") for _ in range(3000))
        other = "".join(rnd.choice("ACGT") for _ in range(1500))
        with tempfile.TemporaryDirectory() as tmp:
            src, dst, rep = Path(tmp, "in.fa"), Path(tmp, "out.fa"), Path(tmp, "r.tsv")
            src.write_text(f">L#LTR/Gypsy\n{long_seq}\n>F#LTR/Gypsy\n{long_seq[500:1500]}\n"
                           f">S#SINE/tRNA\n{long_seq[2000:2600]}\n>O#DNA\n{other}\n")
            stats = dcf.run(str(src), str(dst), threads=1, report=str(rep))
            ids = [l[1:].strip() for l in dst.read_text().splitlines() if l.startswith(">")]
            self.assertEqual(ids, ["L#LTR/Gypsy", "S#SINE/tRNA", "O#DNA"])
            self.assertEqual(stats, {"input": 4, "kept": 3, "dropped": 1, "kept_cross_class": 1})


if __name__ == "__main__":
    unittest.main()
