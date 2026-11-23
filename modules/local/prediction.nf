process PREDICT {
    label 'prediction'

    input:
    path adapter
    path 'genome.fasta'
    path 'reference'
    path hints, stageAs: 'hints.gff'
    path extrinsic, stageAs: 'extrinsic.cfg'
    val base
    val opts

    output:
    path "AUGUSTUS_${base}.gff"

    script:
    """
    python3 annota_task.py predict '${opts}' ${task.cpus}
    """
}

process EXTRACT {
    label 'prediction'

    input:
    path adapter
    path 'genome.fasta'
    path prediction
    val base
    val opts

    output:
    path "AUGUSTUS_${base}.aa", emit: proteins
    path "AUGUSTUS_${base}.codingseq", emit: coding
    path "AUGUSTUS_${base}.cdsexon", emit: exons

    script:
    """
    python3 annota_task.py extract '${opts}' ${task.cpus}
    """
}

process CLEAN {
    label 'parsers'

    input:
    path adapter
    path shared
    path 'proteins.fasta'
    val base
    val opts

    output:
    path "Clear_${base}.fasta"

    script:
    """
    python3 annota_task.py clean '${opts}' ${task.cpus}
    """
}
