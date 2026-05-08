// MAKE_TABLES module

nextflow.enable.dsl = 2

//

process MAKE_TABLES {

    tag "${sample}"
    label 'pysam'

    input:
    tuple val(sample), path(mitobam), path(mitobai)

    output:
    tuple val(sample),
          path("${sample}.allelic_table.tsv.gz"),
          path("${sample}.coverage_table.tsv.gz"), emit: tables

    script:
    """
    # Assume alignment on hg38 reference, with chrM as the mitochondrial chromosome
    python ${baseDir}/bin/make_tables.py \
        --bam      ${mitobam} \
        --sample   ${sample} \
        --ref      ${params.ref} \
        --region   chrM \
        --min_mq       ${params.min_mq} \
        --min_bq       ${params.min_bq} \
        --max_softclip ${params.max_softclip}
    """

    stub:
    """
    touch ${sample}.allelic_table.tsv.gz
    touch ${sample}.coverage_table.tsv.gz
    """

}


//
