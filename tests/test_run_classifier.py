import os
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_CLASSIFIER = ROOT / "bin" / "run_Classifier"


class RunClassifierEnvironmentTests(unittest.TestCase):
    @staticmethod
    def make_executable(path, content):
        path.write_text(textwrap.dedent(content).lstrip())
        path.chmod(0o755)

    def test_uses_target_python_without_nested_conda_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            work = tmp_path / "work"
            classifyte = tmp_path / "ClassifyTE"
            env_prefix = tmp_path / "ClassifyTE_env"
            fake_bin = tmp_path / "bin"
            trace = tmp_path / "python.trace"

            work.mkdir()
            fake_bin.mkdir()
            (env_prefix / "bin").mkdir(parents=True)
            for dirname in ("features", "models", "nodes"):
                (classifyte / dirname).mkdir(parents=True)
            (classifyte / "generate_feature_file_parallel.py").write_text("# fake\n")
            (classifyte / "evaluate.py").write_text("# fake\n")
            (work / "raw_te_library.fa").write_text(">original\nACGT\n")

            self.make_executable(
                env_prefix / "bin" / "python",
                """
                #!/usr/bin/env bash
                printf '%s\n' "$*" >> "$PYTHON_TRACE"
                case "$1" in
                    *generate_feature_file_parallel.py)
                        mkdir -p data
                        printf 'feature\n' > data/sup.csv
                        ;;
                    *evaluate.py)
                        mkdir -p output
                        printf 'prediction\n' > output/predicted_out_sup_features.csv
                        ;;
                esac
                """,
            )
            self.make_executable(
                fake_bin / "renameTE",
                """
                #!/usr/bin/env bash
                cp "$2" "$3"
                """,
            )
            self.make_executable(
                fake_bin / "RepeatClassifier",
                """
                #!/usr/bin/env bash
                touch "$2.classified"
                """,
            )
            self.make_executable(
                fake_bin / "process_for_classify.py",
                """
                #!/usr/bin/env bash
                printf 'result\n' > "$2"
                """,
            )
            self.make_executable(
                fake_bin / "Combine_for_Two",
                """
                #!/usr/bin/env bash
                printf '>TE#I\ngypsy#Unknown#\n' > TEs.fa
                """,
            )
            self.make_executable(
                fake_bin / "conda",
                """
                #!/usr/bin/env bash
                touch "$CONDA_CALLED"
                exit 99
                """,
            )

            env = os.environ.copy()
            env.update({
                "PATH": f"{fake_bin}:{env['PATH']}",
                "CLASSIFYTE_ENV_PREFIX": str(env_prefix),
                "CONDA_EXE": str(fake_bin / "conda"),
                "CONDA_PREFIX": "/broken/active/environment",
                "CONDA_DEFAULT_ENV": "PGTA",
                "CONDA_SHLVL": "1",
                "PYTHON_TRACE": str(trace),
                "CONDA_CALLED": str(tmp_path / "conda.called"),
            })
            result = subprocess.run(
                [str(RUN_CLASSIFIER), str(classifyte), "2"],
                cwd=work,
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse((tmp_path / "conda.called").exists())
            calls = trace.read_text().splitlines()
            self.assertEqual(len(calls), 2)
            self.assertIn("generate_feature_file_parallel.py -j 1", calls[0])
            self.assertIn("evaluate.py -f sup.csv", calls[1])
            self.assertIn(">TE#LINE", (work / "TEs.fa").read_text())

    def test_missing_environment_does_not_fall_back_to_system_python(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            work = tmp_path / "work"
            classifyte = tmp_path / "ClassifyTE"
            fake_bin = tmp_path / "bin"
            work.mkdir()
            fake_bin.mkdir()
            for dirname in ("features", "models", "nodes"):
                (classifyte / dirname).mkdir(parents=True)
            (work / "raw_te_library.fa").write_text(">original\nACGT\n")

            self.make_executable(
                fake_bin / "conda",
                f"""
                #!/usr/bin/env bash
                if [ "$1 $2" = "info --base" ]; then
                    printf '%s\n' '{tmp_path / "conda-base"}'
                elif [ "$1 $2" = "env list" ]; then
                    printf '%s\n' '# conda environments:'
                fi
                """,
            )
            env = os.environ.copy()
            env.pop("CLASSIFYTE_ENV_PREFIX", None)
            env.update({
                "PATH": f"{fake_bin}:{env['PATH']}",
                "CONDA_EXE": str(fake_bin / "conda"),
            })
            result = subprocess.run(
                [str(RUN_CLASSIFIER), str(classifyte), "2"],
                cwd=work,
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Unable to locate Conda environment", result.stderr)


if __name__ == "__main__":
    unittest.main()
