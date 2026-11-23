"""Deterministic command doubles for exercising the workflow, not tool accuracy."""

import sys
from pathlib import Path

PROTEINS = ">geneA label\nmaaa\n>geneB\nmccc\n>geneC\nmeee\n>geneD\nmggg\n"
PSM_HEADER = (
    "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n"
)


def ids(fasta):
    return [
        line[1:].split()[0]
        for line in Path(fasta).read_text().splitlines()
        if line.startswith(">")
    ]


def blast_row(query, title, coverage=90):
    return f"{query}\tref\tref\t50\t1e-20\t60\t40\t{coverage}\t{title}\n"


def blastp(arg):
    queries = ids(arg("-query"))
    database = Path(arg("-db")).name
    hits = []
    if database == "swiss":
        if "geneA" in queries:
            hits.extend(
                [
                    blast_row("geneA", "ref First enzyme OS=Plant"),
                    blast_row("geneA", "ref Equal bitscore OS=Plant"),
                ]
            )
        if "geneB" in queries:
            hits.append(blast_row("geneB", "ref Hypothetical OS=Plant"))
    elif database == "secondary":
        if "geneB" in queries:
            hits.append(blast_row("geneB", "ref Catalytic B [Plant]"))
        if "geneC" in queries:
            hits.append(blast_row("geneC", "ref Low coverage [Plant]", 30))
    else:
        raise ValueError(f"unexpected database: {database}")
    Path(arg("-out")).write_text("".join(hits))


def interproscan(arg):
    if Path("reference/fail").exists():
        sys.exit("fixture InterProScan failure")
    if Path("reference/missing").exists():
        return
    gff = ["##gff-version 3\n"]
    if not Path("reference/no-matches").exists():
        for query in ids(arg("-i")):
            if query not in {"geneA", "geneC"}:
                continue
            gff += [
                f"##sequence-region {query} 1 4\n",
                f"{query}\t.\tpolypeptide\t1\t4\t.\t+\t.\tID={query}\n",
                f"{query}\tPfam\tprotein_match\t1\t4\t.\t+\t.\tName=PF001;"
                'signature_desc=Domain;Dbxref="InterPro:IPR000001";'
                'Ontology_term="GO:0000001"\n',
            ]
    Path(arg("-o")).write_text("".join(gff))


def augustus():
    print("##gff-version 3")
    for idx, query in enumerate(ids_from_proteins(), 1):
        start, end = idx * 12 + 1, idx * 12 + 12
        strand = "-" if idx == 2 else "+"
        print(f"# start gene {query}")
        print(
            f"scaffold\tAUGUSTUS\tgene\t{start}\t{end}\t.\t{strand}\t."
            f"\tID=gene{idx}"
        )
        print(
            f"scaffold\tAUGUSTUS\ttranscript\t{start}\t{end}\t.\t{strand}\t."
            f"\tID={query};Parent=gene{idx}"
        )
        print(
            f"scaffold\tAUGUSTUS\tCDS\t{start}\t{end}\t.\t{strand}\t0"
            f"\tParent={query}"
        )
        print(f"# end gene {query}")


def ids_from_proteins():
    return ["geneA", "geneB", "geneC", "geneD"]


def extract(args):
    prefix = Path(args[0]).with_suffix("")
    prefix.with_suffix(".aa").write_text(PROTEINS)
    dna = "".join(f">{query}\natggctgctgct\n" for query in ids_from_proteins())
    prefix.with_suffix(".codingseq").write_text(dna)
    prefix.with_suffix(".cdsexon").write_text(dna)


def kallisto(args, arg):
    if args[0] == "index":
        assert ids(args[-1]) == ids_from_proteins()
        Path(arg("-i")).write_text("index fixture")
        return
    if "--single" in args:
        assert arg("-l") == "180"
        assert arg("-s") == "20"
        assert Path(args[-1]).name == "z-R1.fq"
    else:
        assert [Path(read).name for read in args[-2:]] == ["z-R1.fq", "a-R2.fq"]
    out = Path(arg("-o"))
    out.mkdir()
    (out / "abundance.tsv").write_text(
        "target_id\tlength\teff_length\test_counts\ttpm\n"
        "geneA\t12\t3\t9\t9\n"
        "geneB\t12\t3\t1\t1\n"
        "geneC\t12\t3\t0\t0\n"
        "geneD\t12\t3\t0\t0\n"
    )


def comet(args):
    spectrum = Path(args[-1])
    assert spectrum.parent.resolve() == Path.cwd()
    parameters = Path(args[0][2:]).read_text()
    assert "output_percolatorfile = 1" in parameters
    assert "decoy_prefix = DECOY_" in parameters
    assert "num_threads = " in parameters
    spectrum.with_suffix(".pin").write_text(spectrum.stem)


def percolator(args, arg):
    source = Path(args[-1]).read_text()
    target = PSM_HEADER
    if source != "no_hits":
        target += (
            f"{source}_1_2_1\t2\t0.05\t0.1\tK.AAA.R\tgeneA\n"
            f"{source}_2_2_1\t2\t0.01\t0.1\tK.CCC.R\tgeneA\tDECOY_geneB\n"
            f"{source}_3_2_1\t2\t0.07\t0.1\tK.EEE.R\tgeneB\n"
        )
    for flag, output in (
        ("-m", target),
        ("-r", PSM_HEADER),
        ("-B", PSM_HEADER),
    ):
        Path(arg(flag)).write_text(output)


def main(tool, args):
    def arg(flag):
        return args[args.index(flag) + 1]

    if tool == "blastp":
        blastp(arg)
    elif tool == "rpsblast":
        hits = "".join(
            blast_row(query, "ref Full CDD description")
            for query in ids(arg("-query"))
        )
        Path(arg("-out")).write_text(hits)
    elif tool == "hmmscan":
        Path(arg("--tblout")).write_text("")
    elif tool == "interproscan":
        interproscan(arg)
    elif tool == "augustus":
        augustus()
    elif tool == "getannofasta":
        extract(args)
    elif tool == "kallisto":
        kallisto(args, arg)
    elif tool == "comet":
        comet(args)
    elif tool == "percolator":
        percolator(args, arg)
    else:
        raise ValueError(f"unknown fixture command: {tool}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
