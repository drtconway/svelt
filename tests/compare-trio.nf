params.genome       = "${projectDir}/data/genome.fa"
params.sources      = file("${projectDir}/data/*.sniffles.vcf")

process prepare {
    input:
    tuple val(vcf_num), path(src)

    output:
    path "out.vcf"

    script:
    """
    awk 'BEGIN { OFS="\\t"; n=$vcf_num * 1000000 } /^#/ { print; next } { \$3 = "x" n; n += 1; print }' < $src > out.vcf
    """
}

process svelt {
    publishDir 'results', mode: 'copy'

    input:
    path(src, name: '?/*')

    output:
    tuple path("out.svelt.vcf"), path("out.svelt.tsv")

    script:
    """
    svelt merge --reference $params.genome -o out.svelt.vcf --write-merge-table out.svelt.tsv $src
    """
}

process jasmine {
    publishDir 'results', mode: 'copy'

    input:
    path(src, name: '?/*')

    output:
    path("out.jasmine.vcf")

    script:
    """
    ls $src > files.txt

    jasmine file_list=files.txt genome_file=$params.genome out_file=tmp.jasmine.vcf \
        --pre_normalize \
        --output_genotypes \
        --clique_merging \
        --dup_to_ins \
        --normalize_type \
        --default_zero_genotype

    bcftools sort -o out.jasmine.vcf tmp.jasmine.vcf
    """
}

workflow {
    def n = 0
    Channel.from(params.sources) | map { p -> tuple(++n, p) } | prepare | collate(3) | view() | (svelt & jasmine)
}