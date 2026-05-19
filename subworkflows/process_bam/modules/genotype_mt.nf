// GENOTYPE_MT module

nextflow.enable.dsl = 2

//

process GENOTYPE_MT {

    tag "${params.individual}"
    label 'pysam'
    publishDir "${params.output_folder}/${params.individual}", mode: 'copy'

    input:
    tuple path(annotated_allelic_table), path(confident_allelic_table)

    output:
    path("${params.individual}.genotyped_allelic_table.tsv.gz"), emit: genotyped

    script:
    """
    python ${baseDir}/bin/genotype.py \\
        --annotated_allelic_table  ${annotated_allelic_table} \\
        --confident_allelic_table  ${confident_allelic_table} \\
        --donor                    ${params.individual} \\
        --gap_th                   ${params.gap_th} \\
        --coverage_th              ${params.min_callable_coverage} \\
        --min_total_counts         ${params.min_total_counts} \\
        --max_sb_pval              ${params.max_sb_pval}
    """

    stub:
    """
    touch ${params.individual}.genotyped_allelic_table.tsv.gz
    """

}


//
