import importlib.machinery
import importlib.util
import json
import shlex
import shutil
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))
loader = importlib.machinery.SourceFileLoader("pan_te_entry", str(BIN / "Pan_TE"))
spec = importlib.util.spec_from_loader(loader.name, loader)
module = importlib.util.module_from_spec(spec)
loader.exec_module(module)


class TeLookerHandoffTests(unittest.TestCase):
    def make_pipeline(self):
        pipeline = module.PanTEPipeline.__new__(module.PanTEPipeline)
        pipeline.threads = 1

        def fake_run(command, work_dir, **_kwargs):
            args = shlex.split(command)
            source = Path(args[args.index("-i") + 1])
            output = Path(args[args.index("-o") + 1])
            shutil.copyfile(source, output)

        pipeline.run_cmd = fake_run
        return pipeline

    def test_discovery_gate_records_copy_length_and_merge_audit(self):
        pipeline = self.make_pipeline()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            families = tmp_path / "families.fasta"
            families.write_text(
                ">too_few members=4\n" + "A" * 120 + "\n"
                ">kept members=5\n" + "C" * 120 + "\n"
                ">missing_members\n" + "G" * 120 + "\n"
            )
            merged, audit = pipeline._gate_and_merge_te_looker(str(families), str(tmp_path))

            self.assertEqual(audit["copy_length_gate"]["kept"], 2)
            self.assertEqual(audit["copy_length_gate"]["dropped"], 1)
            self.assertEqual(audit["copy_length_gate"]["headers_missing_members_assumed_minimum"], 1)
            self.assertEqual(audit["input"]["families"], 3)
            self.assertEqual(audit["output"]["families"], 2)
            self.assertTrue(Path(merged).is_file())

    def test_only_explicit_gated_contract_bypasses_host_gate(self):
        pipeline = self.make_pipeline()
        with tempfile.TemporaryDirectory() as tmp:
            provenance = Path(tmp) / "provenance.json"
            provenance.write_text(json.dumps({"mode": "gated"}))
            self.assertFalse(pipeline._is_gated_dtr_output(str(provenance)))

            provenance.write_text(json.dumps({
                "contract": module.TE_LOOKER_GATED_CONTRACT,
                "mode": "gated",
                "final_library": "families.fasta",
            }))
            self.assertTrue(pipeline._is_gated_dtr_output(str(provenance)))

    def test_handoff_provenance_is_valid_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / module.TE_LOOKER_HANDOFF_PROVENANCE
            module.PanTEPipeline._write_te_looker_handoff(
                str(out), {"schema": "pan-te-te-looker-handoff-v1", "mode": "test"})
            self.assertEqual(json.loads(out.read_text())["mode"], "test")

    def test_final_gated_capability_is_detected_from_help_not_version(self):
        pipeline = self.make_pipeline()
        pipeline._dtr_caps = None

        def fake_run(args, **_kwargs):
            if args[1:] == ['--version']:
                return SimpleNamespace(returncode=0, stdout='dtr 0.3.1\n', stderr='')
            self.assertEqual(args[1:], ['run', '--help'])
            return SimpleNamespace(returncode=0, stdout='usage ... --final-gated\n', stderr='')

        with patch.object(module.subprocess, 'run', side_effect=fake_run):
            caps = pipeline.probe_dtr_capabilities('/fake/dtr', {})
        self.assertTrue(caps['final_gated'])

    def test_final_gated_is_not_inferred_from_version_alone(self):
        pipeline = self.make_pipeline()
        pipeline._dtr_caps = None

        def fake_run(args, **_kwargs):
            if args[1:] == ['--version']:
                return SimpleNamespace(returncode=0, stdout='dtr 9.9.9\n', stderr='')
            return SimpleNamespace(returncode=0, stdout='usage without feature\n', stderr='')

        with patch.object(module.subprocess, 'run', side_effect=fake_run):
            caps = pipeline.probe_dtr_capabilities('/fake/dtr', {})
        self.assertFalse(caps['final_gated'])


if __name__ == "__main__":
    unittest.main()
