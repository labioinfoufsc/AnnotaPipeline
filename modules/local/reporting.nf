process EXPORT_PROTEINS {
    label 'parsers'

    input:
    path adapter
    path parser
    path 'annotations.txt'
    path 'clean.fasta'
    path gff, stageAs: 'prediction.gff'
    val base
    val opts

    output:
    path "AnnotaPipeline_${base}_proteins.fasta"

    script:
    """
    python3 annota_task.py export-proteins '${opts}' ${task.cpus}
    """
}

process EXPORT_GFF {
    label 'parsers'

    input:
    path adapter
    path parser
    path 'annotations.txt'
    path 'prediction.gff'
    val base
    val opts

    output:
    path "${base}_Annotated_GFF.gff"

    script:
    """
    python3 annota_task.py export-gff '${opts}' ${task.cpus}
    """
}

process EXPORT_TRANSCRIPTS {
    label 'parsers'

    input:
    path adapter
    path shared
    path 'annotated.fasta'
    path 'coding.fasta'
    val base
    val opts

    output:
    path "AnnotaPipeline_${base}_transcripts.fasta"

    script:
    """
    python3 annota_task.py export-transcripts '${opts}' ${task.cpus}
    """
}

process SUMMARY {
    label 'parsers'

    input:
    path adapter
    path parser
    path 'annotations.txt'
    path 'annotated.gff3'
    path 'hypothetical.gff3'
    path tpm, stageAs: 'tpm.tsv'
    path proteomics, stageAs: 'proteomics.tsv'
    val base
    val opts

    output:
    path "AnnotaPipeline_${base}_Summary.tsv"

    script:
    """
    python3 annota_task.py summary '${opts}' ${task.cpus}
    """
}

process RUN_METADATA {
    label 'parsers'

    input:
    path adapter
    path protein_fasta
    path annotated_gff
    path transcripts
    path summary
    path evidence
    val opts

    output:
    path 'run-metadata.json'

    script:
    """
    python3 annota_task.py metadata '${opts}' ${task.cpus}
    """
}
