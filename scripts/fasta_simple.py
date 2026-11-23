"""Legacy CLI for joining FASTA sequences with annotations."""

import warnings

import typer

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    annot: str = typer.Option(..., "-annot", "--annotations"),
    basename: str = typer.Option(..., "-b", "--basename"),
    seqs: str = typer.Option(..., "-faf", "--fasta"),
    org: str = typer.Option(..., "-org", "--organism"),
) -> None:
    # Info from annotation will be stored here
    ant = {}

    # ============================= Extract annotations ==========================

    with open(annot) as source:
        annotation_lines = source.read().splitlines()
    for line in annotation_lines:
        line = line.split("\t")
        ant[line[0]] = line[1]

    with open(seqs) as source:
        fasta = source.read().split(">")
    del fasta[0]  # remove empty value

    # ============== Write fasta file with new header ==========================

    ids_warn = []  # store ids with no annotation
    with open(f"AnnotaPipeline_{basename!s}_proteins.fasta", "w") as out:
        for sequence in fasta:
            sequence = sequence.split("\n")
            sequence_id = sequence[0]
            try:
                out.write(
                    f">{sequence_id!s} | Organism: {org!s} | "
                    f"Description: {ant[sequence_id]!s}\n"
                )
                out.write(str("\n".join(sequence[1:])))
            except Exception:
                ids_warn.append(sequence_id)

    # =============== Warning message, if something goes wrong =====================

    if len(ids_warn) > 0:
        warnings.warn(
            "WARNING: Not all sequences from fasta file were in annotation file",
            stacklevel=2,
        )
        with open(f"{basename!s}_ids_with_no_annotations.txt", "w") as warn_seq:
            warn_seq.writelines(str(a) + "\n" for a in ids_warn)
        warnings.warn(
            "INFO: IDs stored in ids_with_no_annotations.txt", stacklevel=2
        )


if __name__ == "__main__":
    app()
