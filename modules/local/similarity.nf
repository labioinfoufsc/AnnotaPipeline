process BLASTP {
    label 'blast'

    input:
    path adapter
    path classification
    path 'query.fasta'
    path 'reference'
    val base
    val database
    val opts

    output:
    path "${base}_BLASTp_AAvs${database}.outfmt6"

    script:
    """
    python3 annota_task.py blast '${opts}' ${task.cpus}
    """
}

process PARSE_SWISSPROT {
    label 'parsers'

    input:
    path adapter
    path classification
    path 'clean.fasta'
    path 'swiss.tsv'
    val base
    val opts

    output:
    path "${base}_annotated_products.txt", emit: accepted
    path "${base}_BLASTp_AA_SwissProted.fasta", emit: unresolved
    path "${base}_SwissProt_annotations.txt", emit: report

    script:
    """
    python3 annota_task.py parse-swiss '${opts}' ${task.cpus}
    """
}

process PARSE_SECONDARY {
    label 'parsers'

    input:
    path adapter
    path classification
    path 'unresolved.fasta'
    path 'secondary.tsv'
    path 'swiss_accepted.txt'
    val base
    val opts

    output:
    path "${base}_annotated_products.txt", emit: annotated
    path "${base}_hypothetical_products.txt", emit: hypothetical
    path "${base}_no_hit_products.txt", emit: no_hit
    path "${base}_SpecifiedDB_annotations.txt", emit: report

    script:
    """
    python3 annota_task.py parse-secondary '${opts}' ${task.cpus}
    """
}

process PARTITION {
    label 'parsers'

    input:
    path adapter
    path shared
    path 'clean.fasta'
    path 'annotated.txt'
    path 'hypothetical.txt'
    path 'no_hit.txt'
    val opts

    output:
    path 'Annotated_Products.fasta', emit: annotated
    path 'Hypothetical_Products.fasta', emit: hypothetical

    script:
    """
    python3 annota_task.py partition '${opts}' ${task.cpus}
    """
}
