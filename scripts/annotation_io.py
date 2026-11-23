"""Sequence preparation and proteomics aggregation shared by both runners."""

import logging
from pathlib import Path


def sequence_cleaner(fasta_file, min_length=0, output_file=None):
    from Bio import SeqIO

    with open(output_file or f"Clear_{fasta_file}", "w") as output:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if len(record.seq) >= min_length:
                output.write(f">{record.id}\n{str(record.seq).upper()}\n")


def fasta_fetcher(input_fasta, id_list, fetcher_output):
    from Bio import SeqIO

    wanted = sorted(set(id_list))
    # keep legacy matching until transcript/protein identity is characterised
    records = (
        seq
        for seq in SeqIO.parse(input_fasta, "fasta")
        for wanted_id in wanted
        if wanted_id in seq.id
    )
    if SeqIO.write(records, fetcher_output, "fasta") < len(wanted):
        logging.getLogger(__name__).info("IDs not found in input FASTA file")


def annotate_codingseq(aa_fasta, codingseq_fasta, basename):
    from Bio import SeqIO

    with open(aa_fasta) as proteins:
        headers = [
            line.strip().replace(">", "") for line in proteins if ">" in line
        ]
    annotations = {header.split("|")[0].strip(): header for header in headers}
    records = SeqIO.to_dict(SeqIO.parse(codingseq_fasta, "fasta"))

    with open(f"AnnotaPipeline_{basename}_transcripts.fasta", "w") as output:
        for key, record in records.items():
            for protein_id, header in annotations.items():
                if protein_id in key:
                    record.id = header
                    record.description = ""
                    record.seq = record.seq.upper()
                    SeqIO.write(record, output, "fasta-2line")


def add_features(feature, data, save):
    """Compute the legacy unique and total counts for one feature."""
    unique_column = f"Unique {feature}"
    total_column = f"Total {feature}"
    selected = data.copy()
    singletons = selected.drop_duplicates(feature, keep=False)
    selected[unique_column] = selected.index.isin(singletons.index)
    selected.replace({False: 0, True: 1}, inplace=True)
    totals = (
        selected.groupby("ProteinID")
        .size()
        .sort_values(ascending=False)
        .reset_index(name=total_column)
    )
    unique = (
        selected.loc[selected[unique_column] == 1]
        .drop(columns=["Peptide", "Spectrum"])
        .groupby("ProteinID")
        .size()
        .reset_index(name=unique_column)
    )
    if "ProteinID" not in save.columns:
        save["ProteinID"] = totals["ProteinID"]
    save = save.set_index("ProteinID").join(unique.set_index("ProteinID"))
    return save.join(totals.set_index("ProteinID")).reset_index()


def aggregate_proteomics(parsed_files, output_file):
    import pandas as pd

    data = pd.DataFrame({"ProteinID": [], "Peptide": [], "Spectrum": []})
    for path in parsed_files:
        data = pd.concat([data, pd.read_csv(path, sep="\t", header=0)])

    # retain legacy counting, including indices and spectrum identifiers
    total = add_features("Peptide", data, pd.DataFrame())
    total = add_features("Spectrum", data, total)
    total = total.fillna(0).astype(
        {"Unique Peptide": int, "Unique Spectrum": int}
    )
    total.sort_values(by="ProteinID", ascending=False).to_csv(
        output_file, sep="\t", index=False
    )


def quantitative_proteomics(path, basename):
    aggregate_proteomics(
        Path(path).glob("*_parsed.tsv"),
        f"{basename}_pre_total_Proteomics_Quantification.tsv",
    )
