"""Internal task adapters. Nextflow supplies staged files and task settings."""

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from base64 import urlsafe_b64decode
from pathlib import Path


def run(*argv, stdout=None, cwd=None, allowed=(0,)):
    command = list(map(str, argv))
    print(f"$ {shlex.join(command)}", file=sys.stderr, flush=True)
    result = subprocess.run(command, stdout=stdout, cwd=cwd)
    if result.returncode not in allowed:
        raise subprocess.CalledProcessError(result.returncode, command)
    return result.returncode


def parser(script, *argv, allowed=(0,)):
    return run(
        sys.executable, Path(__file__).parent / script, *argv, allowed=allowed
    )


def required(path, *, empty=False):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"missing output: {path}")
    if not empty and path.stat().st_size == 0:
        raise ValueError(f"empty output: {path}")
    return path


def fasta_ids(path):
    from Bio.SeqIO import parse

    with open(path) as fasta:
        ids = [record.id for record in parse(fasta, "fasta")]
    if len(ids) != len(set(ids)):
        raise ValueError(f"duplicate FASTA identifiers: {path}")
    return set(ids)


def predict(options, _):
    gff = f"AUGUSTUS_{options['basename']}.gff"
    flags = [f"--{key}={value}" for key, value in options["flags"].items()]
    for key, name in (
        ("hintsfile", "hints.gff"),
        ("extrinsicCfgFile", "extrinsic.cfg"),
    ):
        if Path(name).is_file():
            flags.append(f"--{key}={name}")
    with open(gff, "w") as output:
        run(
            options["exe"],
            f"--AUGUSTUS_CONFIG_PATH=reference/{options['entry']}",
            *flags,
            "genome.fasta",
            stdout=output,
        )
    required(gff)


def extract(options, _):
    prefix = f"AUGUSTUS_{options['basename']}"
    run(options["exe"], f"{prefix}.gff", "--seqfile=genome.fasta")
    for extension in ("aa", "codingseq", "cdsexon"):
        required(f"{prefix}.{extension}")


def clean(options, _):
    from annotation_io import sequence_cleaner

    output = f"Clear_{options['basename']}.fasta"
    sequence_cleaner("proteins.fasta", options["min_length"], output)
    if not fasta_ids(required(output)):
        raise ValueError("no protein sequences passed the length filter")


def blast_search(options, cpus):
    from blast_classification import blast

    out = f"{options['basename']}_BLASTp_AAvs{options['database']}.outfmt6"
    blast(
        options["exe"],
        "query.fasta",
        out,
        f"reference/{options['entry']}",
        options["max_targets"],
        options["evalue"],
        cpus,
    )


def parse_swiss(options, _):
    from blast_classification import process_swiss

    process_swiss(options["basename"], "clean.fasta", "swiss.tsv", options)


def parse_secondary(options, _):
    from blast_classification import process_secondary

    basename = options["basename"]
    shutil.copyfile("swiss_accepted.txt", f"{basename}_annotated_products.txt")
    process_secondary(
        basename,
        "unresolved.fasta",
        "secondary.tsv",
        options["db_format"],
        options,
    )


def partition(_, __):
    from annotation_io import fasta_fetcher

    query = "clean.fasta"
    groups = [
        {
            line.split("\t", 1)[0]
            for line in Path(f"{kind}.txt").read_text().splitlines()
            if line
        }
        for kind in ("annotated", "hypothetical", "no_hit")
    ]
    union = set().union(*groups)
    if union != fasta_ids(query) or sum(map(len, groups)) != len(union):
        raise ValueError(
            "BLAST classes do not partition the cleaned protein IDs; "
            "check ambiguous identifiers and BLAST query ID formatting"
        )
    for label, ids in (
        ("Annotated", groups[0]),
        ("Hypothetical", groups[1] | groups[2]),
    ):
        output = f"{label}_Products.fasta"
        fasta_fetcher(query, ids, output)
        if fasta_ids(output) != ids:
            raise ValueError(f"{label} FASTA has ambiguous protein IDs")


def interproscan(options, cpus):
    ref = Path("reference") / options["entry"]
    out = f"{options['basename']}_interproscan_{options['group']}_output.gff3"
    if not Path("query.fasta").stat().st_size:
        Path(out).write_text("")
        return
    scratch = Path("ips-tmp").absolute()
    scratch.mkdir()
    run(
        "bash",
        ref,
        "-i",
        "query.fasta",
        "-o",
        out,
        "-f",
        "GFF3",
        "-t",
        "p",
        "-goterms",
        "-iprlookup",
        "-cpu",
        cpus,
        "-T",
        scratch,
        *options["args"],
    )
    required(out, empty=True)
    with open(out) as output:
        header = output.readline()
    if header and not header.startswith("##gff-version 3"):
        raise ValueError(f"missing GFF3 header: {out}")


def hmmscan(options, cpus):
    out = f"{options['basename']}_hmmscan_output.txt"
    if not Path("query.fasta").stat().st_size:
        Path(out).write_text("")
        return
    # hmmer starts worker threads in addition to the main thread
    run(
        options["exe"],
        "--cpu",
        max(cpus - 1, 0),
        "--tblout",
        out,
        "--noali",
        "-E",
        options["evalue"],
        "--domE",
        options["domain_evalue"],
        f"reference/{options['entry']}",
        "query.fasta",
    )
    required(out, empty=True)


def rpsblast(options, cpus):
    from blast_classification import BLAST_FIELDS

    out = f"{options['basename']}_rpsblast_output.outfmt6"
    if not Path("query.fasta").stat().st_size:
        Path(out).write_text("")
        return
    run(
        options["exe"],
        "-query",
        "query.fasta",
        "-out",
        out,
        "-db",
        f"reference/{options['entry']}",
        "-outfmt",
        BLAST_FIELDS,
        "-num_threads",
        cpus,
        "-evalue",
        options["evalue"],
        "-max_target_seqs",
        options["max_targets"],
    )
    required(out, empty=True)


def functional_report(options, _):
    parser(
        "funcannotation_parser.py",
        "-ipr_annot",
        "annotated.gff3",
        "-ipr_hyp",
        "hypothetical.gff3",
        "-hmm",
        "hmmer.txt",
        "-rpsblast",
        "rps.tsv",
        "-basename",
        options["basename"],
    )


def merge_annotation(_, __):
    parser(
        "info_parser.py",
        "-ipr1",
        "annotated.gff3",
        "-ipr2",
        "hypothetical.gff3",
        "-a",
        "annotated.txt",
        "-hy",
        "hypothetical.txt",
        "-nh",
        "no_hit.txt",
    )
    version_sort("All_annotation_products.txt", "All_Annotated_Products.txt")


def version_sort(source, target):
    with open(target, "w") as output:
        run("sort", "-V", required(source), stdout=output)


def export_proteins(options, _):
    basename = options["basename"]
    gff = Path("prediction.gff").is_file()
    args = ["-gff", "prediction.gff"] if gff else []
    organism = options["organism"]
    parser(
        "gfftofasta_parser.py" if gff else "fasta_simple.py",
        *args,
        "-annot",
        "annotations.txt",
        "-b",
        basename,
        "-faf",
        "clean.fasta",
        "-org",
        f'"{organism}"',
    )
    output = required(f"AnnotaPipeline_{basename}_proteins.fasta")
    if fasta_ids(output) != fasta_ids("clean.fasta"):
        raise ValueError("not all protein IDs matched the annotations/GFF")


def export_gff(options, _):
    parser(
        "fastatogff.py",
        "-gff",
        "prediction.gff",
        "-all",
        "annotations.txt",
        "-b",
        options["basename"],
    )


def export_transcripts(options, _):
    from annotation_io import annotate_codingseq

    basename = options["basename"]
    annotate_codingseq("annotated.fasta", "coding.fasta", basename)
    output = required(f"AnnotaPipeline_{basename}_transcripts.fasta")
    if fasta_ids(output) != fasta_ids("annotated.fasta"):
        raise ValueError("protein and annotated transcript IDs do not match")


def kallisto_index(options, _):
    run(
        options["exe"],
        "index",
        "-i",
        f"{options['basename']}_kallisto_index.idx",
        "transcripts.fasta",
    )


def kallisto_quant(options, cpus):
    args = ["--single"] if options["single_end"] else []
    for key in ("l", "s"):
        if options.get(key) is not None:
            args.extend([f"-{key}", options[key]])
    if options.get("seed") is not None:
        args.extend(["--seed", options["seed"]])
    out = f"{options['basename']}_kallisto_output"
    reads = sorted(Path("reads").glob("mate*/*"))
    run(
        options["exe"],
        "quant",
        "-i",
        "transcripts.idx",
        "-o",
        out,
        "-b",
        options["bootstrap"],
        "-t",
        cpus,
        *args,
        *reads,
    )
    required(Path(out) / "abundance.tsv")


def filter_tpm(options, _):
    threshold = options["threshold"]
    if threshold == "mean":
        flags = ["-tpmavg"]
    elif threshold == "median":
        flags = ["-tpmmd"]
    else:
        flags = ["-tpmval", threshold]
    parser(
        "kallisto_parser.py",
        "-ktfile",
        "abundance.tsv",
        "-basename",
        options["basename"],
        *flags,
    )


def rewrite_comet_params(source, target, settings):
    remaining = dict(settings)
    lines = []
    assignment = re.compile(r"^(\s*)(\w+)\s*=[^#\n]*(.*)$")
    for line in Path(source).read_text().splitlines():
        match = assignment.match(line)
        if match and match[2] in settings:
            key = match[2]
            if key not in remaining:
                raise ValueError(f"duplicate Comet parameter: {key}")
            comment = f" {match[3]}" if match[3] else ""
            line = f"{match[1]}{key} = {remaining.pop(key)}{comment}"
        lines.append(line)
    for key, value in remaining.items():
        lines.append(f"{key} = {value}")
    Path(target).write_text("\n".join(lines) + "\n")


def comet_params(_, __):
    rewrite_comet_params(
        "original.params",
        "prepared.params",
        {
            "decoy_search": 1,
            "output_pepxmlfile": 0,
            "output_percolatorfile": 1,
            "decoy_prefix": "DECOY_",
        },
    )


def comet(options, cpus):
    rewrite_comet_params(
        "prepared.params", "comet.params", {"num_threads": cpus}
    )
    staged = next(Path("source").iterdir())
    spectrum = Path(staged.name)
    shutil.copyfile(staged, spectrum)
    flags = []
    if options.get("first") is not None:
        flags = [f"-F{options['first']}", f"-L{options['last']}"]
    run(
        options["exe"],
        "-Pcomet.params",
        "-Dproteins.fasta",
        *flags,
        str(spectrum.absolute()),
    )
    pins = list(Path(".").glob("*.pin"))
    if len(pins) != 1:
        raise ValueError(f"expected one Comet PIN for {spectrum}, got {pins}")
    # use a uniform name internally; the workflow publishes by spectrum ID
    required(pins[0]).rename("result.pin")


def percolator(options, cpus):
    args = ["--num-threads", min(cpus, 128)]
    if options.get("seed") is not None:
        args.extend(["--seed", options["seed"]])
    run(
        options["exe"],
        *args,
        "-r",
        "peptide_output.tsv",
        "-m",
        "percolator_output.tsv",
        "-B",
        "decoy_output.tsv",
        "input.pin",
    )
    for out in (
        "percolator_output.tsv",
        "peptide_output.tsv",
        "decoy_output.tsv",
    ):
        required(out)


def parse_psm(options, _):
    basename = options["basename"]
    code = parser(
        "percolator_parser.py",
        "-p",
        "percolator.tsv",
        "-qv",
        options["qvalue"],
        "-b",
        basename,
        allowed=(0, 2),
    )
    if code == 2:
        filtered = required(f"{basename}_qvalue_filtered.tsv")
        with filtered.open() as table:
            header = table.readline()
            if "q-value" not in header.split("\t") or table.readline():
                raise ValueError("unexpected Percolator parser exit 2")
        Path(f"{basename}_parsed.tsv").write_text(
            "ProteinID\tPeptide\tSpectrum\n"
        )


def aggregate(options, _):
    from annotation_io import aggregate_proteomics

    files = sorted(Path("observations").glob("*_parsed.tsv"))
    if not files:
        raise ValueError("missing Percolator observations")
    # sort data, not the header, so summary can consume header-only output
    aggregate_proteomics(files, "unsorted.tsv")
    with open("unsorted.tsv") as table:
        header = table.readline()
        with open("counts.tsv", "w") as counts:
            shutil.copyfileobj(table, counts)
    with open(
        f"{options['basename']}_Total_Proteomics_Quantification.tsv", "w"
    ) as output:
        output.write(header)
        output.flush()
        run("sort", "-V", "counts.tsv", stdout=output)


def summary(options, _):
    args = []
    for flag, name in (("-tr", "tpm.tsv"), ("-proteomics", "proteomics.tsv")):
        if Path(name).is_file():
            args.extend([flag, name])
    parser(
        "summary_parser.py",
        "-b",
        options["basename"],
        "-annot",
        "annotations.txt",
        "-ipr_annot",
        "annotated.gff3",
        "-ipr_hyp",
        "hypothetical.gff3",
        *args,
    )


def metadata(options, _):
    Path("run-metadata.json").write_text(
        json.dumps(options, indent=2, sort_keys=True) + "\n"
    )


OPERATIONS = {
    "predict": predict,
    "extract": extract,
    "clean": clean,
    "blast": blast_search,
    "parse-swiss": parse_swiss,
    "parse-secondary": parse_secondary,
    "partition": partition,
    "interproscan": interproscan,
    "hmmscan": hmmscan,
    "rpsblast": rpsblast,
    "functional-report": functional_report,
    "merge": merge_annotation,
    "export-proteins": export_proteins,
    "export-gff": export_gff,
    "export-transcripts": export_transcripts,
    "kallisto-index": kallisto_index,
    "kallisto-quant": kallisto_quant,
    "filter-tpm": filter_tpm,
    "comet-params": comet_params,
    "comet": comet,
    "percolator": percolator,
    "parse-psm": parse_psm,
    "aggregate": aggregate,
    "summary": summary,
    "metadata": metadata,
}


if __name__ == "__main__":
    operation, payload, threads = sys.argv[1:]
    os.environ["LC_ALL"] = "C"
    OPERATIONS[operation](json.loads(urlsafe_b64decode(payload)), int(threads))
