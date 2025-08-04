params.genome = "${projectDir}/data/genome.fa"
params.sources = file("${projectDir}/data/*.sniffles.vcf")
params.samples = file("${projectDir}/data/samples.tsv")

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
    svelt merge --reference ${params.genome} -o ${fam}_merged.svelt.vcf --write-merge-table ${fam}_merged.svelt.tsv ${src}
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
    touch ${fam}_fam.pdf
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
    jasmined = grouped | jasmine | view()
    paired = svelted.combine(jasmined, by: 0)
    venn(paired) | view()
    
    //def n = 0
    //def m = 0
    //trios = Channel.from(params.sources) | map { p -> tuple(++n, p) } | prepare | collate(3) | map { p -> tuple(++m, p) } | view()
    //s = svelt(trios) | view()
    //j = jasmine(trios) | view()
    //c = s.combine(j, by: 0) | view()
    //venn(c) | view()
}
