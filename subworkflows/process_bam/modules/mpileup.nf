// MPILEUP module

nextflow.enable.dsl = 2

//

process MPILEUP {

    label 'bcftools'
    publishDir "${params.output_folder}/${params.individual}", mode: 'copy'

    input:
    tuple val(samples), path(bams), path(bais)

    output:
    tuple path("${params.individual}.MT.vcf.gz"),
          path("${params.individual}.MT.vcf.gz.tbi"), emit: pileup_vcf

    script:
    """
    # Create BAM list from all staged BAMs
    ls *.bam > bam_list.txt

    # MT pileup across all samples
    bcftools mpileup \
        -f ${params.ref} \
        -r ${params.region} \
        -d ${params.max_depth} \
        -L ${params.max_idepth} \
        -q ${params.min_mq} \
        -Q ${params.min_bq} \
        -a AD,DP \
        -O z \
        --threads ${task.cpus} \
        -b bam_list.txt \
        -o ${params.individual}.MT.vcf.gz

    # Index
    bcftools index --tbi \
        --threads ${task.cpus} \
        ${params.individual}.MT.vcf.gz
    """

    stub:
    """
    touch ${params.individual}.MT.vcf.gz
    touch ${params.individual}.MT.vcf.gz.tbi
    """

}


//
