// GATHER_TABLES module

nextflow.enable.dsl = 2

//

process GATHER_TABLES {

    label 'pysam'
    publishDir "${params.output_folder}/${params.individual}", mode: 'copy'

    input:
    path(tables)
    path(mt_cn_tables)

    output:
    tuple path("${params.individual}.allelic_table.tsv.gz"),
          path("${params.individual}.coverage_table.tsv.gz"),
          path("${params.individual}.mt_cn.tsv.gz"), emit: output

    script:
    """
    python ${baseDir}/bin/gather_tables.py --individual ${params.individual}
    """

    stub:
    """
    touch ${params.individual}.allelic_table.tsv.gz
    touch ${params.individual}.coverage_table.tsv.gz
    touch ${params.individual}.mt_cn.tsv.gz
    """

}


//
