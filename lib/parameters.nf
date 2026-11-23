def payload(Map options) {
    groovy.json.JsonOutput.toJson(options)
        .getBytes('UTF-8').encodeBase64().toString()
}

def keys(value, Collection allowed, String name) {
    if (!(value instanceof Map)) {
        error "${name} must be a mapping"
    }
    def unexpected = value.keySet() - allowed
    if (unexpected) {
        error "Unknown ${name} setting(s): ${unexpected.join(', ')}"
    }
    value
}

def text(value, String name) {
    if (!(value instanceof String) || !value.trim()) {
        error "${name} must be a nonempty string"
    }
    value
}

def inputFile(value, String name) {
    def path = file(text(value, name), glob: false, checkIfExists: true)
    if (!java.nio.file.Files.isRegularFile(path)) {
        error "${name} is not a regular file: ${path}"
    }
    path
}

def decimal(value, String name, Number minimum, Number maximum = null) {
    if (!(value instanceof Number) || !Double.isFinite(value.doubleValue()) ||
        value < minimum || (maximum != null && value > maximum)) {
        error "${name} must be a number >= ${minimum}" +
            (maximum == null ? '' : " and <= ${maximum}")
    }
    value
}

def integer(value, String name, int minimum) {
    decimal(value, name, minimum, Integer.MAX_VALUE)
    if (value != value.intValue()) {
        error "${name} must be an integer"
    }
    value.intValue()
}

def reference(value, String name, String defaultEntry = null) {
    def cfg = keys(value, ['directory', 'entry', 'manifest'], name)
    def directory = file(text(cfg.directory, "${name}.directory"),
                         glob: false, checkIfExists: true)
    if (!java.nio.file.Files.isDirectory(directory)) {
        error "${name}.directory is not a directory: ${directory}"
    }
    def entry = text(cfg.entry ?: defaultEntry, "${name}.entry")
    def manifest = text(cfg.manifest ?: 'SHA256SUMS', "${name}.manifest")
    [entry, manifest].each { part ->
        if (part.startsWith('/') || part.tokenize('/').contains('..')) {
            error "${name} requires paths relative to its reference directory"
        }
    }
    def manifestPath = directory.resolve(manifest)
    if (!java.nio.file.Files.isRegularFile(manifestPath) || !java.nio.file.Files.size(manifestPath)) {
        error "${name} has no nonempty reference manifest: ${manifestPath}"
    }
    def digest = java.security.MessageDigest.getInstance('SHA-256')
        .digest(java.nio.file.Files.readAllBytes(manifestPath)).encodeHex().toString()
    [directory: directory, entry: entry, digest: digest, manifest: manifest]
}

def settings(p) {
    def basename = text(p.basename, 'basename')
    if (!(basename ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
        error 'basename must start with a letter/digit and use only A-Z, a-z, 0-9, ., _, -'
    }
    if (!!p.genome == !!p.protein || (p.genome && p.gff)) {
        error 'Supply either genome or protein; gff is only for protein mode'
    }

    def cfg = [
        basename: basename, organism: text(p.organism, 'organism'),
        software: text(p.software_id, 'software_id'),
        min_length: integer(p.min_length, 'min_length', 0),
        genome: p.genome ? inputFile(p.genome, 'genome') : null,
        protein: p.protein ? inputFile(p.protein, 'protein') : null,
        gff: p.gff ? inputFile(p.gff, 'gff') : null,
    ]
    def exe = [
        augustus: 'augustus', getannofasta: 'getAnnoFasta.pl',
        blastp: 'blastp', hmmscan: 'hmmscan', rpsblast: 'rpsblast',
        kallisto: 'kallisto', comet: 'comet', percolator: 'percolator',
    ]
    keys(p.tools, exe.keySet(), 'tools').each { key, value ->
        exe[key] = text(value, "tools.${key}")
        if (exe[key].contains('/') && !exe[key].startsWith('/')) {
            error "tools.${key} must be a PATH command or absolute path"
        }
    }
    cfg.tools = exe

    keys(p.references, [
        'swissprot', 'secondary', 'cdd', 'pfam', 'interproscan', 'augustus',
    ], 'references')
    def refs = [:]
    ['swissprot', 'secondary', 'cdd', 'pfam'].each { name ->
        refs[name] = reference(p.references[name], "references.${name}")
    }
    refs.interproscan = reference(p.references.interproscan,
                                  'references.interproscan', 'interproscan.sh')
    ['pfam', 'interproscan'].each { name ->
        def entryPath = refs[name].directory.resolve(refs[name].entry)
        if (!java.nio.file.Files.isRegularFile(entryPath)) {
            error "Missing ${name} entry: ${refs[name].entry}"
        }
    }
    cfg.refs = refs

    def blastDefaults = [
        identity: 40, positivity: 60, coverage: 30, evalue: 1e-5,
        max_targets: 10, db_format: 'nrdb', customsep: '|', customcolumn: 5,
        keywords: ['fragment', 'hypothetical', 'partial', 'uncharacterized',
                   'unknown', 'unspecified'],
    ]
    def blast = blastDefaults + keys(p.blast, blastDefaults.keySet(), 'blast')
    ['identity', 'positivity', 'coverage'].each { name ->
        decimal(blast[name], "blast.${name}", 0, 100)
    }
    decimal(blast.evalue, 'blast.evalue', 0)
    integer(blast.max_targets, 'blast.max_targets', 1)
    if (!(blast.db_format in ['nrdb', 'trembldb', 'eupathdb', 'custom'])) {
        error 'blast.db_format must be nrdb, trembldb, eupathdb, or custom'
    }
    if (!(blast.keywords instanceof List)) {
        error 'blast.keywords must be a list of strings'
    }
    blast.keywords.each { keyword -> text(keyword, 'blast.keywords item') }
    text(blast.customsep, 'blast.customsep')
    integer(blast.customcolumn, 'blast.customcolumn', 0)
    cfg.blast = blast

    def hmm = [evalue: 1e-5, domain_evalue: 1e-5]
    cfg.hmmer = hmm + keys(p.hmmer, hmm.keySet(), 'hmmer')
    decimal(cfg.hmmer.evalue, 'hmmer.evalue', 0)
    decimal(cfg.hmmer.domain_evalue, 'hmmer.domain_evalue', 0)
    def rps = [evalue: 1e-5, max_targets: 10]
    cfg.rpsblast = rps + keys(p.rpsblast, rps.keySet(), 'rpsblast')
    decimal(cfg.rpsblast.evalue, 'rpsblast.evalue', 0)
    integer(cfg.rpsblast.max_targets, 'rpsblast.max_targets', 1)

    def ips = keys(p.interproscan, ['applications', 'min_size'], 'interproscan')
    cfg.interproscan = []
    if (ips.applications != null) {
        if (!(ips.applications instanceof List) || !ips.applications ||
            !ips.applications.every { application ->
                application instanceof String &&
                    application ==~ /[A-Za-z0-9_-]+/
            }) {
            error 'interproscan.applications must be a nonempty list of names'
        }
        cfg.interproscan.addAll(['-appl', ips.applications.join(',')])
    }
    if (ips.min_size != null) {
        cfg.interproscan.addAll(
            ['-ms', integer(ips.min_size, 'interproscan.min_size', 1)]
        )
    }

    if (cfg.genome) {
        refs.augustus = reference(p.references.augustus, 'references.augustus', '.')
        if (!java.nio.file.Files.isDirectory(refs.augustus.directory.resolve(refs.augustus.entry))) {
            error 'references.augustus.entry must be a configuration directory'
        }
        def flags = [
            strand: 'both', genemodel: 'partial', protein: 'on',
            introns: 'on', start: 'on', stop: 'on', cds: 'on',
        ]
        def aug = keys(p.augustus, flags.keySet() +
            ['species', 'UTR', 'hintsfile', 'extrinsicCfgFile'], 'augustus')
        cfg.augustus = [
            flags: flags + aug.findAll { key, _val ->
                !(key in ['hintsfile', 'extrinsicCfgFile'])
            }.collectEntries { key, val -> [key, text(val, "augustus.${key}")] },
            hints: aug.hintsfile ? [inputFile(aug.hintsfile, 'augustus.hintsfile')] : [],
            extrinsic: aug.extrinsicCfgFile ?
                [inputFile(aug.extrinsicCfgFile, 'augustus.extrinsicCfgFile')] : [],
        ]
        text(aug.species, 'augustus.species')
    }

    if (p.kallisto != null) {
        if (!cfg.genome) {
            error 'Kallisto requires genome mode and predicted coding sequences'
        }
        def k = [bootstrap: 0, threshold: 'mean'] + keys(p.kallisto,
            ['reads', 'l', 's', 'bootstrap', 'threshold', 'seed'], 'kallisto')
        if (!(k.reads instanceof List) || !(k.reads.size() in [1, 2])) {
            error 'kallisto.reads must contain one or two paths, in mate order'
        }
        def reads = k.reads.collect { read ->
            inputFile(read, 'kallisto.reads item')
        }
        if (reads.toSet().size() != reads.size()) {
            error 'kallisto.reads must contain distinct paths'
        }
        if (!(k.threshold in ['mean', 'median'])) {
            decimal(k.threshold, 'kallisto.threshold', 0)
        }
        ['l', 's'].each { name ->
            if (reads.size() == 1 || k[name] != null) {
                if (decimal(k[name], "kallisto.${name}", 0) == 0) {
                    error "kallisto.${name} must be positive"
                }
            }
        }
        integer(k.bootstrap, 'kallisto.bootstrap', 0)
        if (k.seed != null) {
            integer(k.seed, 'kallisto.seed', 0)
        }
        cfg.kallisto = [options: k - [reads: k.reads], reads: reads]
    }

    if (p.proteomics != null) {
        def ms = [qvalue: 0.05] + keys(p.proteomics,
            ['spectra', 'comet_params', 'first', 'last', 'qvalue', 'seed'],
            'proteomics')
        if (!(ms.spectra instanceof List) || !ms.spectra) {
            error 'proteomics.spectra must be a nonempty list of files'
        }
        def spectra = ms.spectra.collect { spectrum ->
            def source = inputFile(spectrum, 'proteomics.spectra item')
            def stem = source.name.replaceFirst(/\.[^.]+$/, '')
                .replaceAll(/[^A-Za-z0-9_.-]/, '_')
            def id = "s_${stem}"
            tuple(id, source)
        }.sort { row -> row[0] }
        if (spectra.collect { row -> row[0] }.toSet().size() != spectra.size()) {
            error 'proteomics.spectra must have unique normalised basenames'
        }
        if (ms.first != null || ms.last != null) {
            integer(ms.last, 'proteomics.last',
                    integer(ms.first, 'proteomics.first', 1) + 1)
        }
        decimal(ms.qvalue, 'proteomics.qvalue', 0, 1)
        if (ms.seed != null) {
            integer(ms.seed, 'proteomics.seed', 1)
            decimal(ms.seed, 'proteomics.seed', 1, 20000)
        }
        cfg.proteomics = [
            options: [first: ms.first, last: ms.last,
                      qvalue: ms.qvalue, seed: ms.seed],
            spectra: spectra,
            parameters: inputFile(ms.comet_params, 'proteomics.comet_params'),
        ]
    }
    cfg
}

def runMetadata(p, Map cfg, String engineVersion) {
    def names = [
        'basename', 'organism', 'genome', 'protein', 'gff', 'min_length',
        'software_id', 'tools', 'augustus', 'blast', 'interproscan', 'hmmer',
        'rpsblast', 'kallisto', 'proteomics',
    ]
    def provided = names.collectEntries { name -> [(name): p[name]] }
    provided.tools = cfg.tools
    provided.blast = cfg.blast
    provided.hmmer = cfg.hmmer
    provided.rpsblast = cfg.rpsblast
    provided.nextflow_version = engineVersion
    provided.references = cfg.refs.collectEntries { name, ref ->
        [(name): [directory: ref.directory.toString(), entry: ref.entry,
                  manifest: ref.manifest, manifest_sha256: ref.digest]]
    }
    provided
}
