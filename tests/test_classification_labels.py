import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROCESS = ROOT / "bin" / "process_for_classify.py"
COMBINE = ROOT / "bin" / "Combine_for_Two"
# Every node of ClassifyTE's nodes/tree.txt.
TREE_NODES = [
    "Retrotransposon", "DNA transposon", "LTR", "Copia", "gypsy", "Bel-Pao", "DIRS",
    "LINE", "R2", "RTE", "Jockey", "L1", "I", "SINE", "tRNA", "7SL", "5S",
    "SubclassI", "TIRS", "Tc1-Mariner", "hAT", "Mutator", "Merlin", "Transib", "P",
    "PiggyBac", "PIF-Harbinger", "CACTA",
]
RM_CLASS = r"^(LTR|LINE|SINE|DNA|RC|Unknown)(/[^\s#]+)?$"


def classify(rc_headers, predictions):
    """Run process_for_classify.py then Combine_for_Two; return {id: label} and logs."""
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        (work / "data").mkdir()
        (work / "data" / "re.fa.classified").write_text("".join(f">{h}\nACGT\n" for h in rc_headers))
        csv_text = "Sequence ID,Predicted label\n" + "".join(f"{i},{label}\n" for i, label in predictions)
        (work / "predicted.csv").write_text(csv_text)
        proc = subprocess.run([sys.executable, str(PROCESS), "predicted.csv", "result.txt"],
                              cwd=work, capture_output=True, text=True, check=True)
        comb = subprocess.run(["perl", str(COMBINE)], cwd=work, capture_output=True, text=True, check=True)
        labels = {}
        for line in (work / "TEs.fa").read_text().splitlines():
            if line.startswith(">"):
                seq_id, label = line[1:].split("#", 1)
                labels[seq_id] = label
        return labels, proc.stderr + comb.stderr


class ClassificationLabelTests(unittest.TestCase):
    def test_every_unknown_gets_its_classifyte_label(self):
        # Regression: the prediction CSV was read two lines at a time, so every
        # second sequence lost its ClassifyTE label and stayed Unknown.
        n = 6
        ids = [f"TE_{i:08d}" for i in range(1, n + 1)]
        labels, _ = classify([f"{i}#Unknown" for i in ids],
                             [(f"{i}#Unknown", "hAT") for i in ids])
        self.assertEqual(labels, {i: "DNA/hAT" for i in ids})

    def test_labels_use_repeatmasker_vocabulary(self):
        labels, _ = classify(
            ["TE_00000001#Unknown", "TE_00000002#Unknown", "TE_00000003#Unknown"],
            [("TE_00000001#Unknown", "gypsy"), ("TE_00000002#Unknown", "SubclassI"),
             ("TE_00000003#Unknown", "tRNA")],
        )
        self.assertEqual(labels["TE_00000001"], "LTR/Gypsy")
        # SubclassI is the DNA-transposon (TIR) subclass, not retrotransposon Class I.
        self.assertEqual(labels["TE_00000002"], "DNA")
        self.assertEqual(labels["TE_00000003"], "SINE/tRNA")

    def test_every_tree_node_is_mapped(self):
        ids = [f"TE_{i:08d}" for i in range(1, len(TREE_NODES) + 1)]
        labels, logs = classify([f"{i}#Unknown" for i in ids],
                                [(f"{i}#Unknown", node) for i, node in zip(ids, TREE_NODES)])
        self.assertNotIn("Unmapped", logs)
        for seq_id, node in zip(ids, TREE_NODES):
            self.assertRegex(labels[seq_id], RM_CLASS, f"{node} -> {labels[seq_id]}")

    def test_repeatclassifier_label_is_kept(self):
        # Positive control: RepeatClassifier already speaks RepeatMasker vocabulary.
        labels, _ = classify(["TE_00000001#LTR/Copia", "TE_00000002#Unknown"],
                             [("TE_00000001#Unknown", "hAT"), ("TE_00000002#Unknown", "L1")])
        self.assertEqual(labels["TE_00000001"], "LTR/Copia")
        self.assertEqual(labels["TE_00000002"], "LINE/L1")

    def test_unrecognised_label_is_reported_not_guessed(self):
        labels, logs = classify(["TE_00000001#Unknown"], [("TE_00000001#Unknown", "Crypton")])
        self.assertEqual(labels["TE_00000001"], "Unknown")
        self.assertIn("Unmapped ClassifyTE label 'Crypton'", logs)


if __name__ == "__main__":
    unittest.main()
