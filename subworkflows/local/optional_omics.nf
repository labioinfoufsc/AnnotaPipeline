include { payload } from '../../lib/parameters'
include {
    KALLISTO_INDEX
    KALLISTO_QUANT
    FILTER_TPM
    PREPARE_COMET
    COMET
    PERCOLATOR
    PARSE_PSM
    AGGREGATE_PROTEOMICS
} from '../../modules/local/omics'

workflow TRANSCRIPTOMICS {
    take:
    cfg
    src
    transcripts

    main:
    def base = cfg.basename
    def common = [basename: base, software: cfg.software]
    def adapter = src.resolve('annota_task.py')
    def options = common + [exe: cfg.tools.kallisto]
    def k = cfg.kallisto
    KALLISTO_INDEX(adapter, transcripts, base, payload(options))
    def quant = options + (k.options - [threshold: k.options.threshold]) +
        [single_end: k.reads.size() == 1]
    KALLISTO_QUANT(adapter, KALLISTO_INDEX.out, k.reads, base, payload(quant))
    FILTER_TPM(adapter, src.resolve('kallisto_parser.py'),
               KALLISTO_QUANT.out.abundance, base,
               payload(common + [threshold: k.options.threshold]))

    emit:
    tpm = FILTER_TPM.out
    raw = KALLISTO_QUANT.out.raw
    index = KALLISTO_INDEX.out
}

workflow PROTEOMICS {
    take:
    cfg
    src
    proteins

    main:
    def base = cfg.basename
    def common = [basename: base, software: cfg.software]
    def adapter = src.resolve('annota_task.py')
    def ms = cfg.proteomics
    PREPARE_COMET(adapter, ms.parameters, payload(common))
    COMET(adapter, channel.fromList(ms.spectra), PREPARE_COMET.out, proteins,
          payload(common + [exe: cfg.tools.comet,
                            first: ms.options.first, last: ms.options.last]))
    PERCOLATOR(adapter, COMET.out.pin,
               payload(common + [exe: cfg.tools.percolator, seed: ms.options.seed]))
    def psms = PERCOLATOR.out.psms.map { spectrum, psm ->
        def opts = payload(common + [basename: "${base}_${spectrum}",
                                     qvalue: ms.options.qvalue])
        tuple(spectrum, psm, opts)
    }
    PARSE_PSM(adapter, src.resolve('percolator_parser.py'), psms, base)
    def observations = PARSE_PSM.out.parsed.map { _spectrum, file -> file }
        .toSortedList { a, b -> a.name <=> b.name }
    AGGREGATE_PROTEOMICS(adapter, src.resolve('annotation_io.py'),
                        observations, base, payload(common))

    emit:
    quant = AGGREGATE_PROTEOMICS.out
    files = COMET.out.pin.mix(COMET.out.parameters, PERCOLATOR.out.raw,
                            PARSE_PSM.out.parsed, PARSE_PSM.out.filtered)
        .map { spectrum, files -> [id: spectrum, files: files] }
}
