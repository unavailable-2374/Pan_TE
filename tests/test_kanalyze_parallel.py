import os
import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "submodule" / "ClassifyTE" / "features" / "kanalyze-2.0.0" / "code"
N_SEQS = 120


@unittest.skipUnless(shutil.which("java"), "java is required for KAnalyze")
class KanalyzeParallelTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        base = Path(self.tmp.name) / "kanalyze-2.0.0"
        self.code = base / "code"
        shutil.copytree(CODE, self.code)
        self.inp = base / "input_data"
        self.out = base / "output_data"
        self.inp.mkdir()
        for k in (2, 3, 4):
            (self.out / f"{k}mer").mkdir(parents=True)
        rng = random.Random(7)
        names = []
        for i in range(1, N_SEQS + 1):
            name = f"seq{i}.fasta"
            seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(200, 3000)))
            (self.inp / name).write_text(f">TE_{i:08d}\n{seq}\n")
            names.append(name)
        (self.code / "list.txt").write_text("\n".join(sorted(names)) + "\n")

    def tearDown(self):
        self.tmp.cleanup()

    def run_parallel(self):
        env = dict(os.environ, TMPDIR=self.tmp.name)
        return subprocess.run(["bash", "runKanalyzer_parallel"], cwd=self.code, env=env,
                              capture_output=True, text=True)

    def serial_count(self, k, name):
        out = Path(self.tmp.name) / f"serial.{k}.{name}.txt"
        subprocess.run(["java", "-jar", "kanalyze.jar", "count", "-k", str(k), "-f", "fasta",
                        "--temploc", self.tmp.name, "-o", str(out), str(self.inp / name)],
                       cwd=self.code, check=True, capture_output=True)
        return out.read_text()

    def test_concurrent_counts_are_complete_and_match_serial(self):
        # Regression: jobs shared segment_1.bin in the output directory and about
        # half of the k-mer tables came back empty.
        proc = self.run_parallel()
        self.assertEqual(proc.returncode, 0, proc.stderr[-2000:])
        empty = [p.name for k in (2, 3, 4) for p in (self.out / f"{k}mer").iterdir() if p.stat().st_size == 0]
        self.assertEqual(empty, [])
        for k in (2, 3, 4):
            self.assertEqual(len(list((self.out / f"{k}mer").glob("*.txt"))), N_SEQS)
        for name in ("seq1.fasta", "seq57.fasta", f"seq{N_SEQS}.fasta"):
            for k in (2, 3, 4):
                self.assertEqual((self.out / f"{k}mer" / f"{name}.txt").read_text(),
                                 self.serial_count(k, name), f"{k}-mer {name}")
        self.assertFalse(list(self.out.rglob("segment_*.bin")), "segment files leaked into output dirs")

    def test_failed_count_fails_the_step(self):
        # A sequence KAnalyze cannot read must not be silently left as an empty table.
        (self.inp / "seq5.fasta").unlink()
        proc = self.run_parallel()
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("seq5.fasta", proc.stderr)


if __name__ == "__main__":
    unittest.main()
