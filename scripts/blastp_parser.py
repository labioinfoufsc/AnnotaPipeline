#!/usr/bin/env python3
"""Run the two-pass similarity analysis from the legacy command line."""

import logging

import typer
from blast_classification import blast, process_secondary, process_swiss

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    seq: str = typer.Option(..., "-s", "--seq"),
    swissprot: str = typer.Option(..., "-sp", "--swissprot"),
    basename: str = typer.Option(..., "-basename"),
    blastp: str = typer.Option(..., "-blastp"),
    nr: str = typer.Option(None, "-nr", "--nrdb"),
    trembl: str = typer.Option(None, "-trbl", "--Trembl"),
    specificdb: str = typer.Option(None, "-spdb", "--specificdb"),
    identity: float = typer.Option(40, "-id", "--identity"),
    positivity: float = typer.Option(60, "-pos", "--positivity"),
    coverage: float = typer.Option(30, "-cov", "--coverage"),
    keywords: str = typer.Option(
        "hypothetical,unspecified,fragment,partial,unknown,fragemnt",
        "-kw",
        "--keywords",
    ),
    threads: int = typer.Option(20, "-t", "--threads"),
    max_target_seqs: int = typer.Option(10, "-max_target_seqs"),
    evalue: float = typer.Option(0.00001, "-evalue"),
    customsep: str = typer.Option("|", "-customsep"),
    customcolumn: int = typer.Option(5, "-customcolumn"),
) -> None:
    """Search SwissProt, then a secondary database, with two BLAST passes."""
    databases = (nr, trembl, specificdb)
    if sum(value is not None for value in databases) != 1:
        raise typer.BadParameter(
            "exactly one of -nr, -trbl or -spdb must be given"
        )

    logging.basicConfig(
        level=logging.INFO,
        filename="SimilarityAnalysis.log",
        format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
    )
    if nr is not None:
        db_format, secondary_db = "nrdb", nr
    elif trembl is not None:
        db_format, secondary_db = "trembldb", trembl
    else:
        db_format, secondary_db = "custom", specificdb
    settings = {
        "identity": identity,
        "positivity": positivity,
        "coverage": coverage,
        "keywords": keywords.split(","),
        "customsep": customsep,
        "customcolumn": customcolumn,
    }
    swiss = f"{basename}_BLASTp_AAvsSwissProt.outfmt6"
    secondary = f"{basename}_BLASTp_AAvsSpecifiedDB.outfmt6"
    unresolved = f"{basename}_BLASTp_AA_SwissProted.fasta"
    logging.info("Searching SwissProt")
    blast(
        blastp,
        seq,
        swiss,
        swissprot,
        max_target_seqs,
        evalue,
        threads,
    )
    process_swiss(basename, seq, swiss, settings)
    logging.info("Searching the secondary database (%s)", db_format)
    blast(
        blastp,
        unresolved,
        secondary,
        secondary_db,
        max_target_seqs,
        evalue,
        threads,
    )
    process_secondary(basename, unresolved, secondary, db_format, settings)


if __name__ == "__main__":
    app()
