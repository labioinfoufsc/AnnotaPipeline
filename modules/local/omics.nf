process KALLISTO_INDEX {
    label 'kallisto'

    input:
    path adapter
    path 'transcripts.fasta'
    val base
    val opts

    output:
    path "${base}_kallisto_index.idx"

    script:
    """
    python3 annota_task.py kallisto-index '${opts}' ${task.cpus}
    """
}

process KALLISTO_QUANT {
    label 'kallisto'

    input:
    path adapter
    path 'transcripts.idx'
    path reads, stageAs: 'reads/mate??/*'
    val base
    val opts

    output:
    path "${base}_kallisto_output", emit: raw
    path "${base}_kallisto_output/abundance.tsv", emit: abundance

    script:
    """
    python3 annota_task.py kallisto-quant '${opts}' ${task.cpus}
    """
}

process FILTER_TPM {
    label 'parsers'

    input:
    path adapter
    path parser
    path 'abundance.tsv'
    val base
    val opts

    output:
    path "${base}_Transcript_Quantification.tsv"

    script:
    """
    python3 annota_task.py filter-tpm '${opts}' ${task.cpus}
    """
}

process PREPARE_COMET {
    label 'parsers'

    input:
    path adapter
    path 'original.params'
    val opts

    output:
    path 'prepared.params'

    script:
    """
    python3 annota_task.py comet-params '${opts}' ${task.cpus}
    """
}

process COMET {
    label 'comet'
    tag "${spectrum}"

    input:
    path adapter
    tuple val(spectrum), path(source, stageAs: 'source/*')
    path 'prepared.params'
    path 'proteins.fasta'
    val opts

    output:
    tuple val(spectrum), path('result.pin'), emit: pin
    tuple val(spectrum), path('comet.params'), emit: parameters

    script:
    """
    python3 annota_task.py comet '${opts}' ${task.cpus}
    """
}

process PERCOLATOR {
    label 'percolator'
    tag "${spectrum}"

    input:
    path adapter
    tuple val(spectrum), path('input.pin')
    val opts

    output:
    tuple val(spectrum), path('percolator_output.tsv'), emit: psms
    tuple val(spectrum), path('*_output.tsv'), emit: raw

    script:
    """
    python3 annota_task.py percolator '${opts}' ${task.cpus}
    """
}

process PARSE_PSM {
    label 'parsers'
    tag "${spectrum}"

    input:
    path adapter
    path parser
    tuple val(spectrum), path('percolator.tsv'), val(opts)
    val base

    output:
    tuple val(spectrum), path("${base}_${spectrum}_parsed.tsv"), emit: parsed
    tuple val(spectrum), path("${base}_${spectrum}_qvalue_filtered.tsv"),
        emit: filtered

    script:
    """
    python3 annota_task.py parse-psm '${opts}' ${task.cpus}
    """
}

process AGGREGATE_PROTEOMICS {
    label 'parsers'

    input:
    path adapter
    path shared
    path observations, stageAs: 'observations/*'
    val base
    val opts

    output:
    path "${base}_Total_Proteomics_Quantification.tsv"

    script:
    """
    python3 annota_task.py aggregate '${opts}' ${task.cpus}
    """
}
