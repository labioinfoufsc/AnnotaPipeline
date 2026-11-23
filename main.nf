include { settings; runMetadata } from './lib/parameters'
include { ANNOTATE } from './workflows/annotapipeline'

workflow {
    main:
    def cfg = settings(params)
    def result = ANNOTATE(cfg, projectDir.resolve('scripts'),
                          runMetadata(params, cfg, nextflow.version.toString()))
    publish:
    prediction = result.prediction
    similarities = result.similarities
    functional = result.functional
    annotations = result.annotations
    rna = result.rna
    ms = result.ms
    spectra = result.spectra
}

output {
    prediction { path "${params.basename}" }
    similarities { path "${params.basename}" }
    functional { path "${params.basename}" }
    annotations { path "${params.basename}" }
    rna { path "${params.basename}" }
    ms { path "${params.basename}" }
    spectra {
        path { item ->
            "${params.basename}/5_PeptideIdentification_${params.basename}/${item.id}"
        }
    }
}
