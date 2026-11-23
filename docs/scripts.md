# Python scripts

The `scripts/` directory contains the Python layer shared by both AnnotaPipeline
entry points: the Nextflow task adapter, the sequence and proteomics helpers,
and the command-line parsers. Scripts write their outputs to the current
working directory, which Nextflow sets to the task directory.

## Layout

| File | Role |
| --- | --- |
| `annota_task.py` | Task adapter invoked by Nextflow modules |
| `annotation_io.py` | FASTA cleaning, ID fetching, transcript headers, peptide counts |
| `blast_classification.py` | BLAST hit filtering and product partitioning |
| `annotapipeline.py` | Direct CLI driver for the YAML configuration |
| `blastp_parser.py` | Two-pass BLAST search and classification |
| `fasta_simple.py` | FASTA headers from an annotation table without a GFF |
| `fastatogff.py` | Annotation transfer into an AUGUSTUS GFF3 file |
| `gfftofasta_parser.py` | FASTA headers from an annotation table and GFF |
| `funcannotation_parser.py` | InterProScan, HMMER, and RPS-BLAST tables |
| `info_parser.py` | Product table joining BLAST and InterProScan results |
| `kallisto_parser.py` | TPM thresholding of Kallisto abundance |
| `percolator_parser.py` | q-value filtering and peptide extraction |
| `summary_parser.py` | Final per-protein summary table |

## CLI conventions

Every command accepts `-h`, `-help`, and `--help`. Flags from the original
scripts are retained, including single-dash long options such as `-basename`.
Typer validates required options and exits with status 2 on malformed input.
`kallisto_parser.py` requires exactly one of `-tpmavg`, `-tpmmd`, or `-tpmval`.
`percolator_parser.py` defaults to `-qv 0.05`.

## Parsers

| Script | Required flags | Optional flags | Outputs |
| --- | --- | --- | --- |
| `blastp_parser.py` | `-s`, `-sp`, `-basename`, `-blastp`, one of `-nr`, `-trbl`, `-spdb` | `-id`, `-pos`, `-cov`, `-kw`, `-t`, `-max_target_seqs`, `-evalue`, `-customsep`, `-customcolumn` | `<basename>_annotated_products.txt`, `<basename>_hypothetical_products.txt`, `<basename>_no_hit_products.txt`, two annotation reports, unresolved protein FASTA |
| `fasta_simple.py` | `-annot`, `-b`, `-faf`, `-org` | | `AnnotaPipeline_<basename>_proteins.fasta`, `<basename>_ids_with_no_annotations.txt` |
| `gfftofasta_parser.py` | `-gff`, `-annot`, `-b`, `-faf`, `-org` | | `AnnotaPipeline_<basename>_proteins.fasta`, `<basename>_ids_with_no_annotations.txt` |
| `fastatogff.py` | `-gff`, `-all`, `-b` | | `<basename>_Annotated_GFF.gff` |
| `funcannotation_parser.py` | `-ipr_annot`, `-ipr_hyp`, `-hmm`, `-rpsblast`, `-basename` | | `InterProScan_Out_<basename>.tsv`, `Hmmscan_Out_<basename>.tsv`, `RPSblast_Out_<basename>.tsv`, `<basename>_Grouped_Hypothetical_Information.tsv` |
| `info_parser.py` | `-ipr1`, `-ipr2`, `-a`, `-nh`, `-hy` | | `All_annotation_products.txt` |
| `kallisto_parser.py` | `-ktfile`, `-basename`, one TPM selector | `--tpm-average`, `--tpm-median`, `--tpm-threshold` | `<basename>_Transcript_Quantification.tsv` |
| `percolator_parser.py` | `-p`, `-b` | `-qv` | `<basename>_qvalue_filtered.tsv`, `<basename>_parsed.tsv` |
| `summary_parser.py` | `-b`, `-annot`, `-ipr_annot`, `-ipr_hyp` | `-tr`, `-proteomics` | `AnnotaPipeline_<basename>_Summary.tsv` |

The direct driver uses `annotapipeline.py -c annotapipeline.yaml` with `-s`,
`-p`, or `-p` plus `-gff`. See the [runner instructions](../README.md#how-to-run).

## Task adapter

Nextflow invokes the adapter with three positional arguments:

```bash
python scripts/annota_task.py OPERATION PAYLOAD THREADS
```

`PAYLOAD` is a URL-safe base64-encoded JSON object describing the staged
inputs and settings. `THREADS` comes from the task resource request. The
adapter sets `LC_ALL=C`, dispatches to the named operation, and raises on any
non-zero tool exit that the operation does not explicitly allow.

Shared helpers perform sequence cleaning, protein ID fetching from AUGUSTUS
FASTA files, transcript header transfer, and aggregation of unique and total
peptide or spectrum counts.

## Data contracts

- GFF input uses the AUGUSTUS GFF3 dialect, including `##sequence-region` and
  `##FASTA` sections. `fastatogff.py` expects gene blocks beginning with
  `# start gene `.
- Annotation tables are tab-separated with an identifier in column 1 and a
  description in column 2. `fastatogff.py` splits on the first tab only.
- Protein FASTA headers written by the parsers have the form
  `>ID | Organism: <organism> | Description: <annotation>`.
- Kallisto input is the tab-separated `abundance.tsv` with a `tpm` column.
- Percolator input is tab-separated; q-value is in column 3. The parser emits
  a header row followed by `ProteinID`, peptide, and spectrum fields.
- The summary table has ten columns: `ProteinID`, `Annotation`, `IPR`, `GO`,
  `Superfamily`, `Transcript`, `TPM`, `Expression`, `Total_Peptide`, and
  `Unique_Peptide`.

## Development

Python 3.12 is provided by the Pixi environment. Formatting and linting use
Ruff with the settings in `ruff.toml`:

```bash
pixi run --locked format
pixi run --locked lint
pixi run --locked test
```

Command contracts are covered by `tests/test_cli.py`; workflow behaviour is
covered by `tests/test_workflow.py` with fixture tools. See
[`nextflow.md`](nextflow.md) for reference staging and execution profiles.
