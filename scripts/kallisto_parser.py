"""Legacy CLI for filtering Kallisto output by a TPM threshold."""

from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    ktfile: str = typer.Option(..., "-ktfile", "--kallisto-file"),
    basename: str = typer.Option(..., "-basename"),
    tpm_avg: bool = typer.Option(False, "-tpmavg", "--tpm-average"),
    tpm_median: bool = typer.Option(False, "-tpmmd", "--tpm-median"),
    tpm_value: str = typer.Option(None, "-tpmval", "--tpm-threshold"),
) -> None:
    if sum((tpm_avg, tpm_median, tpm_value is not None)) != 1:
        raise typer.BadParameter(
            "exactly one of -tpmavg, -tpmmd or -tpmval is required"
        )
    # SET FILE LOCATION ------------------------------------------------------------

    kallisto_file = Path(ktfile).absolute()

    # OPEN DATAFRAMES --------------------------------------------------------------

    df_kallisto = pd.read_csv(kallisto_file, sep="\t")

    # PARSE DATAFRAMES -------------------------------------------------------------

    if tpm_avg:
        tpm = df_kallisto["tpm"].mean()
        print(f"Average TPM: {tpm}")
    elif tpm_median:
        tpm = df_kallisto["tpm"].median()
        print(f"Median TPM: {tpm}")
    else:
        tpm = float(tpm_value)
        print(f"Chosen TPM value: {tpm}")

    df_kt_threshold = df_kallisto.loc[df_kallisto["tpm"] >= float(tpm)].copy()
    drop_columns = [
        column
        for column in ("length", "eff_length", "est_counts")
        if column in df_kt_threshold.columns
    ]
    df_kt_threshold = df_kt_threshold.drop(columns=drop_columns)

    df_kt_threshold.to_csv(
        f"{basename}_Transcript_Quantification.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    app()
