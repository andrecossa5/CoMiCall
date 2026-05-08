// FILTER_MUTS module

nextflow.enable.dsl = 2

//

process FILTER_MUTS {

    tag "${params.individual}"
    label 'pysam'
    publishDir "${params.output_folder}/${params.individual}", mode: 'copy'

    input:
    tuple path(allelic_table), path(coverage_table), path(mt_cn)

    output:
    tuple path("${params.individual}.confident_allelic_table.tsv.gz"),
          path("${params.individual}.final_allelic_table.tsv.gz"),
          path("${params.individual}.metrics.txt"), emit: filtered

    script:
    """
    python ${baseDir}/bin/filter.py \\
        --allelic_table            ${allelic_table} \\
        --coverage_table           ${coverage_table} \\
        --ref                      ${params.ref} \\
        --donor                    ${params.individual} \\
        --min_strand_ratio         ${params.min_strand_ratio} \\
        --max_strand_ratio         ${params.max_strand_ratio} \\
        --af_threshold             ${params.af_threshold} \\
        --min_callable_coverage    ${params.min_callable_coverage} \\
        --max_sb_pval              ${params.max_sb_pval} \\
        --max_prevalence_low_AF    ${params.max_prevalence_low_AF} \\
        --min_gs_count_low_AF      ${params.min_gs_count_low_AF} \\
        --min_overlap_count_low_AF ${params.min_overlap_count_low_AF} \\
        --min_n_mutations          ${params.min_n_mutations} \\
        --max_prevalence_germline  ${params.max_prevalence_germline}
    """

    stub:
    """
    touch ${params.individual}.confident_allelic_table.tsv.gz
    touch ${params.individual}.final_allelic_table.tsv.gz
    touch ${params.individual}.metrics.txt
    """

}


//
