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
    "pan_te_ltr_handoff_module", str(BIN / "Pan_TE"))
spec = importlib.util.spec_from_loader(loader.name, loader)
pan_te = importlib.util.module_from_spec(spec)
loader.exec_module(pan_te)


class LtrBuilderHandoffTests(unittest.TestCase):
    def test_all_rtr_files_are_sent_to_one_streamed_builder_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "genome").mkdir()
            (root / "genome" / "genome.fa").write_text(">chr1\nACGT\n")
            rtr_dir = root / "Look4LTRs" / "Rtr"
            rtr_dir.mkdir(parents=True)
            rtr_b = rtr_dir / "b.rtr"
            rtr_a = rtr_dir / "a.rtr"
            rtr_b.write_text("header\n")
            rtr_a.write_text("header\n")
            builder = root / "builder.py"
            builder.write_text("pass\n")

            pipeline = pan_te.PanTEPipeline.__new__(pan_te.PanTEPipeline)
            pipeline.output_dir = str(root)
            pipeline.ltr_threads = 3
            calls = []
            pipeline.run_cmd = lambda command, **kwargs: calls.append(
                (command, kwargs))

            old_resolver = pan_te.resolve_ltr_library_builder
            pan_te.resolve_ltr_library_builder = lambda: str(builder)
            try:
                pipeline.build_ltr_consensus()
            finally:
                pan_te.resolve_ltr_library_builder = old_resolver

            self.assertEqual(len(calls), 1)
            command, kwargs = calls[0]
            self.assertIn(f"--rtr {rtr_a} {rtr_b}", command)
            self.assertIn("--threads 3", command)
            self.assertTrue(kwargs["stream"])

    def test_stream_mode_returns_successful_completed_process(self):
        pipeline = pan_te.PanTEPipeline.__new__(pan_te.PanTEPipeline)
        result = pipeline.run_cmd(
            f"{sys.executable} -c 'import sys; print(123); "
            "print(456, file=sys.stderr)'",
            stream=True)
        self.assertEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
