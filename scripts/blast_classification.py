"""BLAST invocation and annotation classification shared by both runners."""

import re
import subprocess
from itertools import groupby
from pathlib import Path

BLAST_FIELDS = "6 qseqid sseqid sacc bitscore evalue ppos pident qcovs stitle"


def blast(executable, query, output, database, max_targets, evalue, cpus):
    if not Path(query).stat().st_size:
        Path(output).write_text("")
        return
    subprocess.run(
        [
            executable,
            "-query",
            query,
            "-out",
            output,
            "-db",
            database,
            "-evalue",
            str(evalue),
            "-outfmt",
            BLAST_FIELDS,
            "-max_target_seqs",
            str(max_targets),
            "-num_threads",
            str(cpus),
        ],
        check=True,
    )
    if not Path(output).is_file():
        raise FileNotFoundError(output)


def description(title, db_format, separator="|", column=5):
    if db_format == "nrdb":
        match = re.search(r"\s(.*?)\s\[.*", title)
        return (match.group(1) if match else title).strip()
    if db_format in ("swissprot", "trembldb"):
        return title.split(" ", 1)[-1].split("OS=")[0].strip()
    if db_format in ("eupathdb", "custom"):
        return (
            title.split(separator)[column]
            .replace("transcript_product=", "")
            .strip()
        )
    raise ValueError(f"unsupported database format: {db_format}")


def candidates(hits, db_format, settings):
    keywords = [word.lower() for word in settings["keywords"]]

    def rows(source):
        for number, line in enumerate(source, 1):
            fields = line.rstrip("\r\n").split("\t", 8)
            if len(fields) != 9 or not fields[0]:
                raise ValueError(f"{hits}:{number}: invalid BLAST format 6 row")
            yield fields

    with open(hits) as source:
        for query, group in groupby(rows(source), key=lambda fields: fields[0]):
            accepted = []
            for fields in group:
                title = description(
                    fields[8],
                    db_format,
                    settings.get("customsep", "|"),
                    settings.get("customcolumn", 5),
                )
                if (
                    float(fields[7]) > float(settings["coverage"])
                    and not any(word in title.lower() for word in keywords)
                    and float(fields[5]) >= float(settings["positivity"])
                    and float(fields[6]) >= float(settings["identity"])
                ):
                    accepted.append((float(fields[3]), title))
            yield query, accepted


def write_candidate(query, accepted, products, report):
    # max retains the first hit for equal bit scores
    best = max(accepted, key=lambda hit: hit[0])[1]
    titles = ";".join(title for _, title in accepted)
    products.write(f"{query}\t{best}\n")
    report.write(f"{query}\t{len(accepted)} Annotation(s): [{titles}]\n")


def process_swiss(basename, protein_seq, hits, settings):
    accepted_ids = []
    with (
        open(f"{basename}_annotated_products.txt", "w") as products,
        open(f"{basename}_SwissProt_annotations.txt", "w") as report,
    ):
        for query, accepted in candidates(hits, "swissprot", settings):
            if accepted:
                write_candidate(query, accepted, products, report)
                accepted_ids.append(query)

    # retain legacy matching until the identity contract has been reviewed
    fasta = Path(protein_seq).read_text().split(">")
    for query in accepted_ids:
        for record in fasta:
            if query in record:
                fasta.remove(record)
    Path(f"{basename}_BLASTp_AA_SwissProted.fasta").write_text(">".join(fasta))


def process_secondary(basename, unresolved, hits, db_format, settings):
    seen = set()
    with (
        open(f"{basename}_hypothetical_products.txt", "w") as hypothetical,
        open(f"{basename}_annotated_products.txt", "a") as products,
        open(f"{basename}_SpecifiedDB_annotations.txt", "w") as report,
    ):
        for query, accepted in candidates(hits, db_format, settings):
            seen.add(query)
            if accepted:
                write_candidate(query, accepted, products, report)
            else:
                hypothetical.write(f"{query}\n")

    with (
        open(unresolved) as fasta,
        open(f"{basename}_no_hit_products.txt", "w") as no_hits,
    ):
        for line in fasta:
            if line.startswith(">") and line[1:].strip() not in seen:
                no_hits.write(line[1:].strip() + "\n")
