"""Legacy CLI for joining annotations into a GFF file."""

import re

import typer

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    anot: str = typer.Option(..., "-all", "--allah"),
    gff: str = typer.Option(..., "-gff", "--gff"),
    basename: str = typer.Option(..., "-b", "--basename"),
) -> None:
    def load_annotations(path):
        result = {}
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip("\r\n").split("\t", 1)
                if len(fields) == 2:
                    result[fields[0]] = fields[1]
        return result

    def transcript_id(fields):
        value = re.sub(r";.*", "", fields[-1])
        return re.sub(r".*=", "", value)

    def write_block(block, annotations, outfile):
        identifier = None
        for line in block:
            if "\ttranscript\t" in line:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) >= 2:
                    identifier = transcript_id(fields)
        if identifier not in annotations:
            return
        for line in block:
            if "\ttranscript\t" in line:
                fields = line.rstrip("\r\n").split("\t")
                outfile.write(
                    "\t".join(fields[:8])
                    + "\t"
                    + annotations[identifier]
                    + "\n"
                )
            else:
                outfile.write(line)

    dict_anot = load_annotations(anot)
    with (
        open(str(gff)) as source,
        open(f"{basename!s}_Annotated_GFF.gff", "w") as outfile,
    ):
        block = []
        started = False
        for line in source:
            if line.startswith("# start gene "):
                if block:
                    write_block(block, dict_anot, outfile)
                block = [line]
                started = True
            elif not started:
                outfile.write(line)
            else:
                block.append(line)
        if block:
            write_block(block, dict_anot, outfile)


if __name__ == "__main__":
    app()
