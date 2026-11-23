"""CLI contract checks for the Typer-based parser adapters."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
PARSERS = (
    "annotapipeline.py",
    "blastp_parser.py",
    "fasta_simple.py",
    "fastatogff.py",
    "funcannotation_parser.py",
    "gfftofasta_parser.py",
    "info_parser.py",
    "kallisto_parser.py",
    "percolator_parser.py",
    "summary_parser.py",
)


class CliTest(unittest.TestCase):
    def run_script(self, script, *args, cwd=None):
        return subprocess.run(
            [sys.executable, str(SCRIPTS / script), *(str(a) for a in args)],
            cwd=cwd,
            capture_output=True,
            text=True,
        )

    def test_help_is_available_for_every_parser(self):
        for script in PARSERS:
            with self.subTest(script=script):
                proc = self.run_script(script, "-h")
                self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_blastp_accepts_legacy_flags(self):
        with tempfile.TemporaryDirectory(prefix="annota-cli-") as tmp:
            root = Path(tmp)
            (root / "query.fa").write_text("")
            proc = self.run_script(
                "blastp_parser.py",
                "-s",
                "query.fa",
                "-sp",
                "swiss.fa",
                "-basename",
                "fix",
                "-spdb",
                "custom.fa",
                "-id",
                40,
                "-pos",
                60,
                "-cov",
                30,
                "-kw",
                "hypothetical,unknown",
                "-t",
                2,
                "-max_target_seqs",
                10,
                "-evalue",
                "1e-5",
                "-blastp",
                "/bin/true",
                "-customsep",
                "|",
                "-customcolumn",
                5,
                cwd=root,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            for name in (
                "fix_annotated_products.txt",
                "fix_hypothetical_products.txt",
                "fix_no_hit_products.txt",
            ):
                with self.subTest(output=name):
                    self.assertTrue((root / name).is_file())

    def test_legacy_runner_requires_exactly_one_input(self):
        for args in (
            ("-c", "config.yaml"),
            ("-s", "seq.fa", "-p", "prot.fa", "-c", "config.yaml"),
        ):
            with self.subTest(args=args):
                proc = self.run_script("annotapipeline.py", *args)
                self.assertNotEqual(proc.returncode, 0)
                self.assertIn("exactly one of", proc.stderr)

    def test_kallisto_rejects_two_thresholds(self):
        proc = self.run_script(
            "kallisto_parser.py",
            "-ktfile",
            "abundance.tsv",
            "-basename",
            "fix",
            "-tpmavg",
            "-tpmmd",
        )
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("exactly one of", proc.stderr)


if __name__ == "__main__":
    unittest.main()
