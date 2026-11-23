"""Legacy CLI for filtering Percolator output by q-value."""

import sys
import warnings

import typer

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    proteom: str = typer.Option(..., "-p", "--percolator-output"),
    basename: str = typer.Option(..., "-b", "--basename"),
    qvalue: float = typer.Option(0.05, "-qv", "-q-value"),
) -> None:
    # Check qvalue
    if not 0.0 <= qvalue <= 1.0:
        print("-qv value must be between 0-1")
        sys.exit(1)

    # --- PARSE QVALUE -----------------------------

    with open(proteom) as source:
        percolator_file = source.readlines()

    with open(f"{basename}_qvalue_filtered.tsv", "w") as parcial_output:
        parcial_output.write(str(percolator_file[0]))
        del percolator_file[0]
        for line in percolator_file:
            fields = line.split("\t")
            if float(fields[2]) <= qvalue:
                parcial_output.write(str("\t".join(fields)))

    # ---------- Parse information ----------------

    # open file that will be used to extract peptides
    with open(f"{basename}_qvalue_filtered.tsv") as source:
        percolator_out = source.readlines()
    # Check if there are valid queries for qvalue used as cutoff
    if len(percolator_out) < 2:
        warnings.warn(
            f"INFO: file {basename}_qvalue_filtered.tsv return no valid queries for qvalue",
            stacklevel=2,
        )
        sys.exit(2)
    # Remove index
    del percolator_out[0]

    dict_write = {}
    for line in percolator_out:
        line_split = line.rstrip("\n").split("\t")
        # Protein Id > can have more than one
        all_ids = line_split[5:]
        peptide = line_split[4]
        filecontent = line_split[0].split("/")[-1]
        _idname, spectrum, _charge, _num = filecontent.rsplit("_", 3)
        if len(all_ids) > 1:
            for protein_id in all_ids:
                # remove decoy matches
                if "DECOY" in protein_id:
                    pass
                elif protein_id not in dict_write:
                    dict_write[f"{protein_id}"] = [f"{peptide}\t{spectrum}\n"]
                else:
                    dict_write[f"{protein_id}"].append(
                        f"{peptide}\t{spectrum}\n"
                    )
        else:
            # If list contains only one element, get this guy
            only_id = all_ids[0]
            # remove decoy matches
            if "DECOY" in only_id:
                pass
            elif only_id not in dict_write:
                dict_write[only_id] = [f"{peptide}\t{spectrum}\n"]
            else:
                dict_write[only_id].append(f"{peptide}\t{spectrum}\n")

    with open(f"{basename}_parsed.tsv", "w") as output:
        output.write("ProteinID\tPeptide\tSpectrum\n")
        for key in sorted(dict_write.keys()):
            info = dict_write.get(key)
            for peptide in info:
                output.write(f"{key}\t{''.join(peptide)}")


if __name__ == "__main__":
    app()
