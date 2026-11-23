include { payload } from '../lib/parameters'
include { PREDICT; EXTRACT; CLEAN } from '../modules/local/prediction'
include {
    BLASTP as BLAST_SWISSPROT
    BLASTP as BLAST_SECONDARY
    PARSE_SWISSPROT
    PARSE_SECONDARY
    PARTITION
} from '../modules/local/similarity'
include {
    INTERPROSCAN
    HMMSCAN
    RPSBLAST
    FUNCTIONAL_REPORT
    MERGE_ANNOTATION
} from '../modules/local/functional'
include {
    EXPORT_PROTEINS
    EXPORT_GFF
    EXPORT_TRANSCRIPTS
    SUMMARY
    RUN_METADATA
} from '../modules/local/reporting'
include { TRANSCRIPTOMICS; PROTEOMICS } from '../subworkflows/local/optional_omics'

def refOptions(ref) {
    [entry: ref.entry, reference_manifest: ref.digest]
}

workflow ANNOTATE {
    take:
    cfg
    src
    provenance

    main:
    def base = cfg.basename
    def common = [basename: base, software: cfg.software]
    def basic = payload(common)
    def adapter = src.resolve('annota_task.py')
    def shared = src.resolve('annotation_io.py')
    def classifier = src.resolve('blast_classification.py')
    def refs = cfg.refs

    def rawProtein = cfg.protein
    def predictedGff = channel.empty()
    def genePrediction = channel.empty()
    if (cfg.genome) {
        def options = common + refOptions(refs.augustus) +
            [exe: cfg.tools.augustus, flags: cfg.augustus.flags]
        PREDICT(adapter, cfg.genome, refs.augustus.directory,
                cfg.augustus.hints, cfg.augustus.extrinsic, base, payload(options))
        EXTRACT(adapter, cfg.genome, PREDICT.out, base,
                payload(common + [exe: cfg.tools.getannofasta]))
        rawProtein = EXTRACT.out.proteins
        predictedGff = PREDICT.out
        genePrediction = predictedGff.mix(EXTRACT.out.proteins, EXTRACT.out.coding,
                                          EXTRACT.out.exons)
    }
    def annotationGff = cfg.genome ? predictedGff :
        (cfg.gff ? channel.value(file(cfg.gff)) : [])
    def hasGff = cfg.genome || cfg.gff
    CLEAN(adapter, shared, rawProtein, base,
          payload(common + [min_length: cfg.min_length]))

    def b = cfg.blast
    def search = common + [exe: cfg.tools.blastp, evalue: b.evalue,
                           max_targets: b.max_targets]
    def classify = common + [identity: b.identity, positivity: b.positivity,
                             coverage: b.coverage, keywords: b.keywords]
    BLAST_SWISSPROT(adapter, classifier, CLEAN.out, refs.swissprot.directory, base,
        'SwissProt', payload(search + refOptions(refs.swissprot) +
                             [database: 'SwissProt']))
    PARSE_SWISSPROT(adapter, classifier, CLEAN.out, BLAST_SWISSPROT.out,
                    base, payload(classify))
    BLAST_SECONDARY(adapter, classifier, PARSE_SWISSPROT.out.unresolved,
        refs.secondary.directory, base, 'SpecifiedDB',
        payload(search + refOptions(refs.secondary) + [database: 'SpecifiedDB']))
    PARSE_SECONDARY(adapter, classifier, PARSE_SWISSPROT.out.unresolved,
        BLAST_SECONDARY.out, PARSE_SWISSPROT.out.accepted, base,
        payload(classify + [db_format: b.db_format, customsep: b.customsep,
                            customcolumn: b.customcolumn]))

    def classes = PARSE_SECONDARY.out
    PARTITION(adapter, shared, CLEAN.out, classes.annotated,
              classes.hypothetical, classes.no_hit, basic)
    def ips = common + refOptions(refs.interproscan) + [args: cfg.interproscan]
    def ipsInputs = PARTITION.out.annotated
        .map { fasta -> tuple('annotated', fasta, payload(ips + [group: 'annotated'])) }
        .mix(PARTITION.out.hypothetical.map { fasta ->
            tuple('hypothetical', fasta, payload(ips + [group: 'hypothetical']))
        })
    INTERPROSCAN(adapter, ipsInputs, refs.interproscan.directory, base)
    def ipsAnnotated = INTERPROSCAN.out
        .filter { group, _file -> group == 'annotated' }
        .map { _group, file -> file }.first()
    def ipsHypothetical = INTERPROSCAN.out
        .filter { group, _file -> group == 'hypothetical' }
        .map { _group, file -> file }.first()

    HMMSCAN(adapter, PARTITION.out.hypothetical, refs.pfam.directory, base,
            payload(common + refOptions(refs.pfam) + cfg.hmmer + [exe: cfg.tools.hmmscan]))
    RPSBLAST(adapter, classifier, PARTITION.out.hypothetical, refs.cdd.directory,
             base, payload(common + refOptions(refs.cdd) + cfg.rpsblast +
                           [exe: cfg.tools.rpsblast]))
    FUNCTIONAL_REPORT(adapter, src.resolve('funcannotation_parser.py'),
                      ipsAnnotated, ipsHypothetical, HMMSCAN.out, RPSBLAST.out,
                      base, basic)
    MERGE_ANNOTATION(adapter, src.resolve('info_parser.py'), ipsAnnotated,
                     ipsHypothetical, classes.annotated, classes.hypothetical,
                     classes.no_hit, basic)

    def fastaParser = hasGff ? 'gfftofasta_parser.py' : 'fasta_simple.py'
    EXPORT_PROTEINS(adapter, src.resolve(fastaParser), MERGE_ANNOTATION.out,
                    CLEAN.out, annotationGff, base,
                    payload(common + [organism: cfg.organism]))
    def gffOutput = channel.empty()
    def finalGff = []
    if (hasGff) {
        EXPORT_GFF(adapter, src.resolve('fastatogff.py'), MERGE_ANNOTATION.out,
                   annotationGff, base, basic)
        gffOutput = EXPORT_GFF.out
        finalGff = gffOutput
    }
    def transcriptOutput = channel.empty()
    def finalTranscripts = []
    if (cfg.genome) {
        EXPORT_TRANSCRIPTS(adapter, shared, EXPORT_PROTEINS.out,
                           EXTRACT.out.coding, base, basic)
        transcriptOutput = EXPORT_TRANSCRIPTS.out
        finalTranscripts = transcriptOutput
    }

    def tpm = []
    def rnaOutputs = channel.empty()
    if (cfg.kallisto) {
        TRANSCRIPTOMICS(cfg, src, transcriptOutput)
        tpm = TRANSCRIPTOMICS.out.tpm
        rnaOutputs = tpm.mix(TRANSCRIPTOMICS.out.raw, TRANSCRIPTOMICS.out.index)
    }
    def peptideCounts = []
    def msCounts = channel.empty()
    def msFiles = channel.empty()
    if (cfg.proteomics) {
        PROTEOMICS(cfg, src, EXPORT_PROTEINS.out)
        peptideCounts = PROTEOMICS.out.quant
        msCounts = peptideCounts
        msFiles = PROTEOMICS.out.files
    }
    SUMMARY(adapter, src.resolve('summary_parser.py'), MERGE_ANNOTATION.out,
            ipsAnnotated, ipsHypothetical, tpm, peptideCounts, base, basic)
    RUN_METADATA(adapter, EXPORT_PROTEINS.out, finalGff, finalTranscripts,
                 SUMMARY.out, FUNCTIONAL_REPORT.out, payload(provenance))

    emit:
    prediction = genePrediction.mix(CLEAN.out)
    similarities = BLAST_SWISSPROT.out.mix(BLAST_SECONDARY.out,
        PARSE_SWISSPROT.out.unresolved, PARSE_SWISSPROT.out.report,
        classes.annotated, classes.hypothetical, classes.no_hit, classes.report)
    functional = ipsAnnotated.mix(ipsHypothetical, HMMSCAN.out, RPSBLAST.out,
        FUNCTIONAL_REPORT.out, PARTITION.out.annotated, PARTITION.out.hypothetical)
    annotations = MERGE_ANNOTATION.out.mix(EXPORT_PROTEINS.out, gffOutput,
        transcriptOutput, SUMMARY.out, RUN_METADATA.out)
    rna = rnaOutputs
    ms = msCounts
    spectra = msFiles
}
