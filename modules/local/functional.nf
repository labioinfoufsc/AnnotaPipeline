process INTERPROSCAN {
    label 'interproscan'

    input:
    path adapter
    tuple val(group), path('query.fasta'), val(opts)
    path 'reference'
    val base

    output:
    tuple val(group), path("${base}_interproscan_${group}_output.gff3")

    script:
    """
    python3 annota_task.py interproscan '${opts}' ${task.cpus}
    """
}

process HMMSCAN {
    label 'hmmer'

    input:
    path adapter
    path 'query.fasta'
    path 'reference'
    val base
    val opts

    output:
    path "${base}_hmmscan_output.txt"

    script:
    """
    python3 annota_task.py hmmscan '${opts}' ${task.cpus}
    """
}

process RPSBLAST {
    label 'blast'

    input:
    path adapter
    path classification
    path 'query.fasta'
    path 'reference'
    val base
    val opts

    output:
    path "${base}_rpsblast_output.outfmt6"

    script:
    """
    python3 annota_task.py rpsblast '${opts}' ${task.cpus}
    """
}

process FUNCTIONAL_REPORT {
    label 'parsers'

    input:
    path adapter
    path parser
    path 'annotated.gff3'
    path 'hypothetical.gff3'
    path 'hmmer.txt'
    path 'rps.tsv'
    val base
    val opts

    output:
    tuple path("InterProScan_Out_${base}.tsv"), path("Hmmscan_Out_${base}.tsv"),
        path("RPSblast_Out_${base}.tsv"),
        path("${base}_Grouped_Hypothetical_Information.tsv")

    script:
    """
    python3 annota_task.py functional-report '${opts}' ${task.cpus}
    """
}

process MERGE_ANNOTATION {
    label 'parsers'

    input:
    path adapter
    path parser
    path 'annotated.gff3'
    path 'hypothetical.gff3'
    path 'annotated.txt'
    path 'hypothetical.txt'
    path 'no_hit.txt'
    val opts

    output:
    path 'All_Annotated_Products.txt'

    script:
    """
    python3 annota_task.py merge '${opts}' ${task.cpus}
    """
}
