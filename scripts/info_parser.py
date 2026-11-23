"""Legacy CLI for merging annotation sources into product tables."""

import os
import os.path

import typer

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    ipr1: str = typer.Option(..., "-ipr1"),
    ipr2: str = typer.Option(..., "-ipr2"),
    annot: str = typer.Option(..., "-a", "--annotated"),
    nohit: str = typer.Option(..., "-nh", "--no_hit"),
    hypo: str = typer.Option(..., "-hy", "--hypothetical"),
) -> None:
    def temporary_query(arq):
        # temporary adding a line to the BLAST analysis
        #   so it accounts for the last query found by the BLAST algorithm
        try:
            temp = arq[-1].split("\t")
            del temp[0]
            temp.insert(0, "QueryTemp")
            arq.append("\t".join(temp))
        # Index error ocours when Interproscan returns no results, as consequence, arq is empty, and split empty file is not allowed
        except IndexError:
            arq.append("\t".join("QueryTemp"))

    # ------------------------------------------

    def parser_interproscan(arq_entrada, arq_ipr):
        with open(arq_entrada) as source:
            entrada = source.read().split("##FASTA")
        interp = entrada[0].split("##sequence-region")
        del interp[0]

        lista = ["coils", "gene3d", "mobidblite"]

        with open(arq_ipr, "a") as ipr:
            for seq_reg in interp:
                seq_reg = seq_reg.splitlines()
                for linha in seq_reg:
                    if linha == seq_reg[0]:
                        linha = linha.split(" ")
                        del linha[0]
                    elif linha == seq_reg[1]:
                        # Ignore the first line of each query group.
                        pass
                    else:
                        ontologia = str(None)
                        anotacao_db = str(None)
                        name = str(None)
                        interpro = str(None)
                        linha = linha.split("\t")
                        nome_subject = linha[0]
                        db = linha[1]
                        try:
                            evalue = str(linha[5]).replace(",", ".").lower()
                        except Exception:
                            # Type of evalue format
                            evalue = str(linha[5])
                        if not any(s in db.lower() for s in lista):
                            db_certo = db
                            anotacao = linha[-1].split(";")
                            for a in anotacao:
                                if "Ontology" in a:
                                    ontologia = a.replace('"', "").replace(
                                        "Ontology_term=", ""
                                    )
                                if "signature_" in a:
                                    anotacao_db = a.replace(
                                        "signature_desc=", ""
                                    )
                                if "Name" in a:
                                    name = a.replace('"', "").replace(
                                        "Name=", ""
                                    )
                                if "Dbxref" in a:
                                    interpro = a.replace('"', "").replace(
                                        "Dbxref=", ""
                                    )
                            ipr.write(
                                f"{nome_subject}\t{db_certo}\t{name}\t"
                                f"{evalue}\t{anotacao_db}\t{interpro}\t"
                                f"{ontologia}\n"
                            )

    # Function to write ids with no ipr and no GO (no result from InterproScan)
    def write_no_ipr(ids_dict):
        # If there's any id without hit in interproscan file
        # It's will write here, without IPR or GO
        if len(ids_dict.keys()) > 0:
            for anot in ids_dict:
                output.write(f"{anot}\t{str(ids_dict.get(anot)).strip()}\n")
        if len(hypo) > 0:
            output.writelines(f"{hyp}\thypothetical protein\n" for hyp in hypo)
        if len(nohit) > 0:
            output.writelines(f"{hyp}\thypothetical protein\n" for hyp in nohit)

    # Function to join IPRs and GOs from InterproScan with IDs hypothetical or without hits
    def intepro_process(ids_dict):
        # File pre-precessed with IDs, IPRs and go_list
        with open("Interpro_out_tmp.txt") as source:
            interpro_out = source.read().splitlines()
        temporary_query(interpro_out)

        old_id = interpro_out[0].split("\t")
        old_id = old_id[0]

        # Initialize lists
        go_list = []
        ipr_list = []
        for query in interpro_out:
            title = query.split("\t")
            new_id = title[0]
            ipr = str(title[5]).split(",")
            go = str(title[6]).split(",")
            if old_id != new_id:
                # Sort IPRs and go_list
                try:
                    ipr_list.sort(key=lambda item: str(item).split("IPR")[1])
                    go_list.sort(key=lambda item: str(item).split("GO")[1])
                except Exception:
                    # No GO or IP for this protein (not a true warning)
                    pass
                # Check if interpro_result is in annotated
                if old_id in ids_dict:
                    if (len(ipr_list) > 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\t{str(ids_dict.get(old_id)).strip()} ({','.join(ipr_list)!s},{','.join(go_list)!s})\n"
                        )
                    elif (len(ipr_list) > 0) and (len(go_list) == 0):
                        output.write(
                            f"{old_id}\t{str(ids_dict.get(old_id)).strip()} ({','.join(ipr_list)!s})\n"
                        )
                    elif (len(ipr_list) == 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\t{str(ids_dict.get(old_id)).strip()} ({','.join(go_list)!s})\n"
                        )
                    else:
                        output.write(
                            f"{old_id}\t{str(ids_dict.get(old_id)).strip()}\n"
                        )
                    del ids_dict[old_id]
                # Else, interpro_result must be in hypothetical
                elif old_id in hypo:
                    if (len(ipr_list) > 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(ipr_list)!s},{','.join(go_list)!s})\n"
                        )
                    elif (len(ipr_list) > 0) and (len(go_list) == 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(ipr_list)!s})\n"
                        )
                    elif (len(ipr_list) == 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(go_list)!s})\n"
                        )
                    else:
                        output.write(f"{old_id}\thypothetical protein\n")
                    hypo.remove(old_id)
                else:
                    if (len(ipr_list) > 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(ipr_list)!s},{','.join(go_list)!s})\n"
                        )
                    elif (len(ipr_list) > 0) and (len(go_list) == 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(ipr_list)!s})\n"
                        )
                    elif (len(ipr_list) == 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(go_list)!s})\n"
                        )
                    else:
                        output.write(f"{old_id}\thypothetical protein\n")
                    nohit.remove(old_id)
                go_list.clear()
                ipr_list.clear()
                # GET FIRST ITEM FOR NEW ID - After save old
                for match in ipr:
                    if (match != "None") and (match not in ipr_list):
                        ipr_list.append(match)
                for match in go:
                    if (match != "None") and (match not in go_list):
                        go_list.append(match)
            else:
                for match in ipr:
                    if (match != "None") and (match not in ipr_list):
                        ipr_list.append(match)
                for match in go:
                    if (match != "None") and (match not in go_list):
                        go_list.append(match)
            old_id = new_id
        write_no_ipr(ids_dict)
        os.remove("Interpro_out_tmp.txt")

    # ---------------------- Create organized file with iprs and gos -----------------------
    if (os.path.getsize(ipr1) == 0) and (os.path.getsize(ipr2) == 0):
        # Can't process empty file
        pass
    else:
        # Almost one file has result - but, if one of them hasn't parser will crash
        if os.path.getsize(ipr1) != 0:
            parser_interproscan(ipr1, "Interpro_out_tmp.txt")
        if os.path.getsize(ipr2) != 0:
            parser_interproscan(ipr2, "Interpro_out_tmp.txt")

    # ---------------------- Pre-parse annotated products -------------------------
    with open(annot) as source:
        annotation_lines = source.read().splitlines()
    ids = []
    desc = []
    for anot in annotation_lines:
        anot = anot.split("\t")
        ids.append(anot[0])
        desc.append(anot[1])

    ids_dict = dict(
        zip(ids, desc, strict=False)
    )  # Create dictionary with id as key and annotation as value

    # ---------------------- Create list with hypothetical ids -----------------
    if os.path.getsize(hypo) != 0:
        with open(hypo) as source:
            hypo = source.read().splitlines()
    else:
        hypo = []
    # ---------------------- Create list with no_hits ids -----------------
    if os.path.getsize(nohit) != 0:
        with open(nohit) as source:
            nohit = source.read().splitlines()
    else:
        nohit = []

    with open("All_annotation_products.txt", "w") as output:
        if os.path.isfile("Interpro_out_tmp.txt") and os.path.getsize(
            "Interpro_out_tmp.txt"
        ):
            intepro_process(ids_dict)
        else:
            write_no_ipr(ids_dict)
            if os.path.isfile("Interpro_out_tmp.txt"):
                os.remove("Interpro_out_tmp.txt")


if __name__ == "__main__":
    app()
