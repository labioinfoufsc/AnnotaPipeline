"""Legacy CLI for writing GFF information into FASTA headers."""

import re
import warnings

import typer

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    gff: str = typer.Option(..., "-gff"),
    annot: str = typer.Option(..., "-annot", "--annotations"),
    basename: str = typer.Option(..., "-b", "--basename"),
    seqs: str = typer.Option(..., "-faf", "--fasta"),
    org: str = typer.Option(..., "-org", "--organism"),
) -> None:
    def fasta_records(path):
        header = None
        lines = []
        with open(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    if header is not None:
                        yield header, lines
                    header = line[1:].rstrip("\r\n")
                    lines = []
                else:
                    lines.append(line)
        if header is not None:
            yield header, lines

    # Info from transcript products will be stored here
    trans = {}
    # Info from annotation will be stored here
    ant = {}

    # ============================= Extract info from gff =========================
    with open(str(gff)) as source:
        for line in source:
            if "\ttranscript\t" not in line:
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 9:
                continue
            # use regex to remove ID= and ;PARENT from annotation ID in gff AUGUSTUS output
            id_store = re.sub(r";.*", "", fields[-1])
            id_store = re.sub(r".*=", "", id_store)
            trans[id_store] = {
                "scaff": str(fields[0]),
                "start": str(fields[3]),
                "end": str(fields[4]),
                "strand": str(fields[6]),
            }

    # ============================= Extract annotations ==========================

    with open(str(annot)) as source:
        for line in source:
            fields = line.rstrip("\r\n").split("\t", 1)
            if len(fields) == 2:
                ant[fields[0]] = fields[1]

    # ============== Write fasta file with new header ==========================

    # Remove quotes from organism
    organism = str(org).replace('"', "").replace("'", "")

    with open(f"AnnotaPipeline_{basename!s}_proteins.fasta", "w") as annotafile:
        ids_warn = []  # store ids with no annotation
        for header, sequence in fasta_records(str(seqs)):
            try:
                annotafile.write(
                    f">{header} | "
                    f"Organism: {organism} | "
                    f"Location: {trans[header]['scaff']!s} | "
                    f"Start: {trans[header]['start']!s} | "
                    f"End: {trans[header]['end']!s} | "
                    f"Strand: {trans[header]['strand']!s} | "
                    f"Description: {ant[header]!s}\n"
                )
                annotafile.write("".join(sequence))
            except Exception:
                # warn ignored, user will be informed all at once
                ids_warn.append(header)

    # =============== Warning message, if something goes wrong =====================

    if len(ids_warn) > 0:
        warnings.warn(
            "WARNING: Not all sequences from fasta file were in annotation file",
            stacklevel=2,
        )
        with open(f"{basename!s}_ids_with_no_annotations.txt", "w") as warn_seq:
            for header in ids_warn:
                warn_seq.write(str(header) + "\n")
        warnings.warn(
            "INFO: IDs stored in ids_with_no_annotations.txt", stacklevel=2
        )


if __name__ == "__main__":
    app()
