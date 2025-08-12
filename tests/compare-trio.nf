params.genome = "${projectDir}/data/genome.fa"
params.sources = file("${projectDir}/data/*.sniffles.vcf")
params.samples = file("${projectDir}/data/samples.tsv")
params.scripts = file("${projectDir}/scripts")

process prepare {
    input:
    tuple val(meta), path(src)

    output:
    tuple val(meta), path("${meta.sample}_renamed.vcf")

    script:
    def vcf_num = task.index
    def out_name = "${meta.sample}_renamed.vcf"
    """
    awk 'BEGIN { OFS="\\t"; n=${vcf_num} * 1000000 } /^#/ { print; next } { \$3 = "x" n; n += 1; print }' < ${src} > ${out_name}
    """

    stub:
    def out_name = "${meta.sample}_renamed.vcf"
    """
    touch ${out_name}
    """
}

process svelt {
    publishDir 'results', mode: 'copy'

    input:
    tuple val(fam), path(src, name: '?/*')

    output:
    tuple val(fam), path("${fam}_merged.svelt.vcf")

    script:
    """
    svelt merge --position-window 75 --reference ${params.genome} -o ${fam}_merged.svelt.vcf --write-merge-table ${fam}_merged.svelt.tsv ${src}
    """

    stub:
    """
    touch ${fam}_merged.svelt.vcf ${fam}_merged.svelt.tsv
    """
}

process jasmine {
    publishDir 'results', mode: 'copy'

    input:
    tuple val(fam), path(src, name: '?/*')

    output:
    tuple val(fam), path("${fam}_merged.jasmine.vcf")

    script:
    """
    ls ${src} > files.txt

    jasmine file_list=files.txt genome_file=${params.genome} out_file=tmp.jasmine.vcf \
        --pre_normalize \
        --output_genotypes \
        --clique_merging \
        --dup_to_ins \
        --normalize_type \
        --default_zero_genotype

    bcftools sort -o ${fam}_merged.jasmine.vcf tmp.jasmine.vcf
    """

    stub:
    """
    touch ${fam}_merged.jasmine.vcf
    """
}

process venn {
    publishDir 'results', mode: 'copy'

    input:
    tuple val(fam), path(s, name: '?/*'), path(j, name: '?/*')

    output:
    tuple val(fam), path("${fam}_venn.pdf")

    script:
    """
    merge-venn-plot ${s} ${j} ${fam}_venn.pdf
    """

    stub:
    """
    touch ${fam}_venn.pdf
    """
}

process difference_summary {
    publishDir 'results', mode: 'copy'

    input:
    tuple val(fam), path(s), path(j), path(src)

    output:
    tuple val(fam), path("${fam}_svelt_summary.tsv"), path("${fam}_jasmine_summary.tsv")

    script:
    """
    python3 ${params.scripts}/differences.py ${fam}_svelt_summary.csv ${fam}_jasmine_summary.csv ${s} ${j} ${src}

    tr ',' '\t' < ${fam}_svelt_summary.csv > ${fam}_svelt_summary.tsv

    tr ',' '\t' < ${fam}_jasmine_summary.csv > ${fam}_jasmine_summary.tsv
    """

    stub:
    """
    touch ${fam}_svelt_summary.tsv ${fam}_jasmine_summary.tsv
    """
}

process positions_summary {
    publishDir 'results', mode: 'copy'

    input:
    tuple val(fam), path(s), path(j), path(src)

    output:
    tuple val(fam), path("${fam}_positions.tsv")

    script:
    """
    python3 ${params.scripts}/positions.py ${fam}_positions.csv ${s} ${j} ${src}

    tr ',' '\t' < ${fam}_positions.csv > ${fam}_positions.tsv
    """

    stub:
    """
    touch ${fam}_positions.tsv
    """
}

process sample_vcf_files {
    publishDir 'results', mode: 'copy'

    input:
    tuple val(fam), path(s), path(j), path(src)

    output:
    path 'sampled/*'

    script:
    """
    mkdir -p sampled

    python3 ${params.scripts}/sample-variants.py ${j} sampled ${src}
    """

    stub:
    """
    mkdir - p sampled

    for f in ${src}
    do
        touch samples/\$f
    done
    """
}

workflow {
    samples_root = params.samples.parent
    samples = Channel.fromPath(params.samples) | \
                splitCsv(header: true, sep: '\t') | map { row ->
        def meta = [:]
        meta.family = row.family_id
        meta.sample = row.sample_id
        meta.source = samples_root / row.vcf_file
        [meta, samples_root / row.vcf_file] }

    prepared = samples | prepare 

    grouped = prepared | map { item -> def meta = item[0]; [meta.family, item] } | \
                groupTuple(sort: {lhs, rhs -> lhs[0].sample.compareTo(rhs[0].sample)}) | \
                map { [it[0], it[1].collect { x -> x[1] }] }

    svelted = grouped | svelt
    jasmined = grouped | jasmine
    paired = svelted.combine(jasmined, by: 0)
    venn(paired) | view()
    extra_paired = paired.combine(grouped, by: 0)
    difference_summary(extra_paired) | view()
    positions_summary(extra_paired) | view()
    sample_vcf_files(extra_paired) | view()
}
