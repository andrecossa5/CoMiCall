// INDEX_BAM module

nextflow.enable.dsl = 2

//

process INDEX_BAM {

    tag "${sample}"
    label 'bcftools'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path(bam), path("${bam}.bai"), emit: indexed_bams

    script:
    """
    samtools index -@ ${task.cpus} ${bam} ${bam}.bai
    """

    stub:
    """
    touch ${bam}.bai
    """

}


//
