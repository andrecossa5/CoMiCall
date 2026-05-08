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
    CONTIG=\$(awk '{print \$1}' ${params.ref}.fai | awk '\$1=="chrM" || \$1=="MT"' | head -n1)
    if [ -z "\$CONTIG" ]; then
        echo "No chrM/MT contig found in reference" >&2
        exit 1
    fi
    python ${baseDir}/bin/make_tables.py \
        --bam      ${mitobam} \
        --sample   ${sample} \
        --ref      ${params.ref} \
        --region   \$CONTIG \
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
