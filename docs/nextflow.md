# Nextflow workflow

This branch replaces the shell-driven runner with a Nextflow DSL2 workflow.
The legacy `scripts/annotapipeline.py` entry point remains available while
the migration is validated.

## Requirements

- Nextflow `25.10.4`
- Python `3.12`
- The external tools and reference databases used by the selected mode
- A nonempty `SHA256SUMS` file in every reference directory

The repository Pixi environment provides Nextflow and the Python
dependencies:

```bash
pixi install --locked
```

## Python adapters

Nextflow tasks call the CLI adapters in `scripts/` with the legacy flags.
They use Typer for parsing and `snake_case` identifiers. Format and lint
them with `pixi run --locked format` and `pixi run --locked lint`. Command
flags, outputs, and data contracts are listed in [`scripts.md`](scripts.md).

## Parameters

The workflow accepts a JSON parameters file. Paths are resolved relative to
the launch directory. Exactly one of `genome` or `protein` is required;
`gff` is only valid with `protein`.

```json
{
  "basename": "my_organism",
  "organism": "My organism",
  "software_id": "annotapipeline-refactor",
  "protein": "proteins.fasta",
  "min_length": 0,
  "tools": {
    "blastp": "blastp",
    "rpsblast": "rpsblast",
    "hmmscan": "hmmscan",
    "augustus": "augustus",
    "getannofasta": "getAnnoFasta.pl",
    "kallisto": "kallisto",
    "comet": "comet",
    "percolator": "percolator"
  },
  "references": {
    "swissprot": {
      "directory": "/databases/swissprot",
      "entry": "uniprot_sprot"
    },
    "secondary": {
      "directory": "/databases/secondary",
      "entry": "secondary"
    },
    "pfam": {
      "directory": "/databases/pfam",
      "entry": "Pfam-A.hmm"
    },
    "cdd": {
      "directory": "/databases/cdd",
      "entry": "Cdd"
    },
    "interproscan": {
      "directory": "/databases/interproscan",
      "entry": "interproscan.sh"
    }
  }
}
```

### Reference contract

Each reference is a directory containing an `entry` and a manifest:

- `directory` is staged as the task-local `reference/` directory.
- `entry` is a path relative to `directory`. It is the BLAST database
  prefix for `swissprot`, `secondary`, and `cdd`, the HMM profile for
  `pfam`, the InterProScan launcher for `interproscan`, and the AUGUSTUS
  configuration directory for `augustus`.
- `manifest` defaults to `SHA256SUMS`. Its SHA-256 digest is part of the
  task inputs, so editing the manifest invalidates dependent tasks on
  `-resume`.

### Genome mode

Set `genome` instead of `protein`. `augustus.species` is required. The
default AUGUSTUS flags match the legacy workflow. The similarity,
functional, and optional omics references from the previous example are
still required; the fragment below only shows the additional AUGUSTUS
entry:

```json
{
  "genome": "genome.fasta",
  "augustus": {
    "species": "my_species",
    "hintsfile": "hints.gff",
    "extrinsicCfgFile": "extrinsic.cfg"
  },
  "references": {
    "augustus": {
      "directory": "/opt/augustus",
      "entry": "config"
    }
  }
}
```

`hintsfile` and `extrinsicCfgFile` are optional and must be regular files.

### Protein mode with a supplied GFF

```json
{
  "protein": "proteins.fasta",
  "gff": "prediction.gff"
}
```

A supplied GFF must use the AUGUSTUS-style GFF3 dialect expected by the
legacy `gfftofasta_parser.py` and `fastatogff.py` scripts.

### Similarity settings

`blast` accepts `identity`, `positivity`, `coverage`, `evalue`,
`max_targets`, `db_format`, `customsep`, `customcolumn`, and `keywords`.
Defaults match the legacy configuration.

### Functional settings

`hmmer` accepts `evalue` and `domain_evalue`. `rpsblast` accepts `evalue`
and `max_targets`. `interproscan` accepts `applications` and `min_size`.

### Optional transcriptomics

Kallisto requires genome mode because it consumes predicted coding
sequences. Supply one read for single-end mode or two reads in mate order
for paired-end mode:

```json
{
  "kallisto": {
    "reads": ["reads_R1.fastq.gz", "reads_R2.fastq.gz"],
    "bootstrap": 0,
    "threshold": "mean"
  }
}
```

Single-end mode requires positive `l` and `s` values.

### Optional proteomics

```json
{
  "proteomics": {
    "spectra": ["sample_a.mgf", "sample_b.mgf"],
    "comet_params": "comet.params",
    "qvalue": 0.05,
    "first": 1,
    "last": 2
  }
}
```

Comet parameters are staged and edited in the task directory; the source
file is never modified. Spectrum names are normalised to unique
`s_<stem>` identifiers before the Comet and Percolator stages.

## Running

```bash
pixi run --locked nextflow run . \
  -params-file params.json \
  -output-dir results \
  -resume
```

The `local` profile is the default. Other profiles are `slurm`, `docker`,
`apptainer`, and `conda`.

## Outputs

Published files are written below `results/<basename>/`:

- `All_Annotated_Products.txt`
- `AnnotaPipeline_<basename>_proteins.fasta`
- `AnnotaPipeline_<basename>_Summary.tsv`
- `<basename>_Annotated_GFF.gff` in genome or supplied-GFF mode
- `AnnotaPipeline_<basename>_transcripts.fasta` in genome mode
- `run-metadata.json`
- Intermediate similarity, functional, and quantification tables

Percolator and Comet outputs are published below
`results/<basename>/5_PeptideIdentification_<basename>/<spectrum>/`.

## Resume and invalidation

`-resume` reuses completed tasks when their inputs and scripts are
unchanged. Editing a reference manifest, tool path, parameter, or staged
input invalidates the affected tasks and their downstream tasks.

## Tests

```bash
pixi run --locked test
```

The integration tests replace external tools with fixtures. They do not
validate scientific parity against real databases; that validation still
requires the legacy output archive and the real toolchain.

## Legacy runner

The command-line flags are unchanged; the parser now uses Typer and also
accepts `-h` and `-help`:

```bash
annotapipeline.py -c annotapipeline.yaml -p proteins.fasta
annotapipeline.py -c annotapipeline.yaml -s genome.fasta
annotapipeline.py -c annotapipeline.yaml -p proteins.fasta -gff prediction.gff
```
