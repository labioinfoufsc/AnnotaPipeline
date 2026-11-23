import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tool_fixture import PROTEINS

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_CMD = shlex.split(os.environ.get("NEXTFLOW_CMD", "nextflow"))


@unittest.skipUnless(
    NEXTFLOW_CMD and shutil.which(NEXTFLOW_CMD[0]), "requires Nextflow"
)
class WorkflowTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="annota-test-")
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.protein = self.root / "plant protein.fasta"
        self.protein.write_text(PROTEINS)
        self.genome = self.root / "plant genome.fasta"
        self.genome.write_text(">scaffold\n" + "ATG" * 32 + "\n")
        for mate in ("z-R1", "a-R2"):
            (self.root / f"{mate}.fq").write_text("@record\nAC\n+\nII\n")

        tool_dir = self.root / "tools"
        tool_dir.mkdir()
        self.tools = {}
        for tool in (
            "blastp",
            "rpsblast",
            "hmmscan",
            "augustus",
            "getannofasta",
            "kallisto",
            "comet",
            "percolator",
        ):
            self.tools[tool] = str(self.wrapper(tool_dir / tool, tool))

        refs = {}
        for key, entry in (
            ("swissprot", "swiss"),
            ("secondary", "secondary"),
            ("pfam", "pfam.hmm"),
            ("cdd", "cdd"),
            ("interproscan", "interproscan.sh"),
            ("augustus", "."),
        ):
            bundle = self.root / key
            bundle.mkdir()
            if key == "interproscan":
                self.wrapper(bundle / entry, "interproscan")
            elif entry != ".":
                (bundle / entry).write_text("test database")
            (bundle / "SHA256SUMS").write_text(f"{key} fixture reference\n")
            refs[key] = {"directory": str(bundle), "entry": entry}

        self.config = {
            "basename": "fixture",
            "organism": "Model plant",
            "software_id": "fixture-tools-v1",
            "protein": self.protein.name,
            "tools": self.tools,
            "references": refs,
        }

    def wrapper(self, output, tool):
        command = shlex.join(
            [sys.executable, str(ROOT / "tests" / "tool_fixture.py"), tool]
        )
        output.write_text(f'#!/bin/bash\nexec {command} "$@"\n')
        output.chmod(0o755)
        return output

    def invoke(self, name, *, resume=False, success=True):
        params = self.root / "params.json"
        params.write_text(json.dumps(self.config))
        trace = self.root / f"{name}.tsv"
        command = [
            *NEXTFLOW_CMD,
            "run",
            str(ROOT),
            "-ansi-log",
            "false",
            "-params-file",
            str(params),
            "-with-trace",
            str(trace),
            "-output-dir",
            str(self.root / name),
        ]
        if resume:
            command.append("-resume")
        env = dict(os.environ, NXF_OFFLINE="true")
        env["PATH"] = (
            str(Path(sys.executable).parent) + os.pathsep + env["PATH"]
        )
        proc = subprocess.run(
            command,
            cwd=self.root,
            env=env,
            capture_output=True,
            text=True,
            timeout=150,
        )
        report = proc.stdout + proc.stderr
        metadata = self.root / name / "fixture" / "run-metadata.json"
        if success:
            assert proc.returncode == 0, report
            assert metadata.is_file(), report
        else:
            assert proc.returncode != 0, report
            assert not metadata.exists(), report
        return self.root / name / "fixture", trace, report

    def trace_status(self, trace, process):
        with trace.open() as handle:
            rows = csv.DictReader(handle, delimiter="\t")
            return [
                row["status"]
                for row in rows
                if row["name"].split(" ")[0].endswith(":" + process)
            ]

    def test_protein_only_and_reference_resume(self):
        out, _, _ = self.invoke("original")
        annotations = (out / "All_Annotated_Products.txt").read_text()
        assert (
            "geneA\tFirst enzyme (InterPro:IPR000001,GO:0000001)" in annotations
        )
        assert "geneB\tCatalytic B\n" in annotations
        assert "geneC\thypothetical protein" in annotations
        assert "geneD\thypothetical protein\n" in annotations
        assert not list(out.glob("*.gff"))
        assert not list(out.glob("*transcripts.fasta"))

        repeat, trace, _ = self.invoke("repeat", resume=True)
        assert self.trace_status(trace, "BLAST_SWISSPROT") == ["CACHED"]
        assert self.trace_status(trace, "SUMMARY") == ["CACHED"]
        assert (
            annotations == (repeat / "All_Annotated_Products.txt").read_text()
        )

        (self.root / "cdd" / "SHA256SUMS").write_text("updated manifest\n")
        _, trace, _ = self.invoke("updated", resume=True)
        assert self.trace_status(trace, "RPSBLAST") == ["COMPLETED"]
        assert self.trace_status(trace, "HMMSCAN") == ["CACHED"]
        assert self.trace_status(trace, "SUMMARY") == ["CACHED"]

    def test_genome_and_both_omics(self):
        del self.config["protein"]
        self.config["genome"] = self.genome.name
        self.config["augustus"] = {"species": "model_plant"}
        self.config["kallisto"] = {"reads": ["z-R1.fq", "a-R2.fq"]}
        spectra = ["sample space.mgf", "no_hits.mgf"]
        for name in spectra:
            (self.root / name).write_text("fixture spectra\n")
        src = self.root / "input comet.params"
        original = "decoy_search = 0\noutput_percolatorfile = 0\n"
        src.write_text(original)
        self.config["proteomics"] = {
            "spectra": spectra,
            "comet_params": src.name,
        }

        out, _, _ = self.invoke("genome")
        proteins = (out / "AnnotaPipeline_fixture_proteins.fasta").read_text()
        assert "Strand: -" in proteins
        transcripts = (
            out / "AnnotaPipeline_fixture_transcripts.fasta"
        ).read_text()
        assert "ATGGCTGCTGCT" in transcripts
        assert (out / "fixture_Annotated_GFF.gff").is_file()
        with (out / "AnnotaPipeline_fixture_Summary.tsv").open() as handle:
            rows = {
                row["ProteinID"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        assert rows["geneA"]["TPM"] == "9.0"
        assert rows["geneA"]["Total_Peptide"] == "2"
        assert rows["geneA"]["Unique_Peptide"] == "2"
        assert rows["geneB"]["Transcript"] == "False"
        assert src.read_text() == original
        assert not list(self.root.glob("*.pin"))
        parsed = (
            out
            / "5_PeptideIdentification_fixture"
            / "s_no_hits"
            / "fixture_s_no_hits_parsed.tsv"
        )
        assert parsed.is_file()

    def test_genome_single_end_kallisto(self):
        del self.config["protein"]
        self.config["genome"] = self.genome.name
        self.config["augustus"] = {"species": "model_plant"}
        self.config["kallisto"] = {
            "reads": ["z-R1.fq"],
            "l": 180,
            "s": 20,
        }
        out, _, _ = self.invoke("single")
        with (out / "AnnotaPipeline_fixture_Summary.tsv").open() as handle:
            rows = {
                row["ProteinID"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        assert rows["geneA"]["TPM"] == "9.0"

    def test_protein_with_gff_and_no_interpro_matches(self):
        result = subprocess.check_output([self.tools["augustus"]], text=True)
        gff = self.root / "supplied genome.gff"
        gff.write_text(result)
        self.config["gff"] = gff.name
        (self.root / "interproscan" / "no-matches").touch()
        out, _, _ = self.invoke("with-gff")
        assert (out / "fixture_Annotated_GFF.gff").is_file()
        annotations = (out / "All_Annotated_Products.txt").read_text()
        assert "geneA\tFirst enzyme\n" in annotations
        assert not list(out.glob("*transcripts.fasta"))

    def test_interpro_failure_is_not_an_empty_result(self):
        (self.root / "interproscan" / "fail").touch()
        _, _, report = self.invoke("failure", success=False)
        assert "fixture InterProScan failure" in report

    def test_reject_transcriptomics_without_coding_sequence(self):
        self.config["kallisto"] = {"reads": ["z-R1.fq"]}
        _, _, report = self.invoke("invalid", success=False)
        assert "Kallisto requires genome mode" in report


if __name__ == "__main__":
    unittest.main()
