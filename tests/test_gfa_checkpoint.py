import importlib.machinery
import importlib.util
import tempfile
import unittest
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))
loader = importlib.machinery.SourceFileLoader(
    "pan_te_checkpoint_module", str(BIN / "Pan_TE"))
spec = importlib.util.spec_from_loader(loader.name, loader)
pan_te = importlib.util.module_from_spec(spec)
loader.exec_module(pan_te)
PanTEPipeline = pan_te.PanTEPipeline


class GfaCheckpointFingerprintTests(unittest.TestCase):
    def make_checkpointed_pipeline(self, root):
        root = Path(root)
        gfa = root / "graph.gfa"
        gfa.write_text(
            "H\tVN:Z:1.0\n"
            "S\ta\tACGTACGT\n"
            "P\tRef#1#Chr1\ta+\t*\n")
        output = root / "output"
        pipeline = PanTEPipeline(
            None, str(output), threads=1, gfa_file=str(gfa),
            gfa_ref="Ref#1", gfa_mode="alleles")
        processed = Path(pipeline.processed_genome)
        processed.parent.mkdir(parents=True, exist_ok=True)
        processed.write_text(">Chr1\nACGTACGT\n")
        (output / "genome.ok").touch()
        return pipeline

    def test_run_rejects_checkpoint_without_fingerprint(self):
        with tempfile.TemporaryDirectory() as tmp:
            pipeline = self.make_checkpointed_pipeline(tmp)
            with self.assertRaisesRegex(RuntimeError, "has no GFA.*fingerprint"):
                pipeline.run()

    def test_run_rejects_mismatched_gfa_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            pipeline = self.make_checkpointed_pipeline(tmp)
            pipeline._write_gfa_fingerprint()
            pipeline.gfa_window_overlap = 0
            with self.assertRaisesRegex(RuntimeError, "different GFA"):
                pipeline.run()


if __name__ == "__main__":
    unittest.main()
