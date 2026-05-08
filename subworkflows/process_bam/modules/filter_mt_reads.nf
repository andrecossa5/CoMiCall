// FILTER_MT_READS module

nextflow.enable.dsl = 2

//

process FILTER_MT_READS {

    tag "${sample}"
    label 'bcftools'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("mitobam.bam"), path("mitobam.bam.bai"), emit: mitobam

    script:
    """
    CONTIG=\$(samtools view -H ${bam} | awk '\$1=="@SQ"{for(i=2;i<=NF;i++) if(\$i ~ /^SN:/){sub("SN:","",\$i); if(\$i=="chrM" || \$i=="MT") print \$i}}' | head -n1)
    if [ -z "\$CONTIG" ]; then
        echo "No chrM/MT contig found in BAM header" >&2
        exit 1
    fi
    samtools view ${bam} -b -@ ${task.cpus} \$CONTIG > mitobam.bam
    samtools index -@ ${task.cpus} mitobam.bam
    """

    stub:
    """
    touch mitobam.bam
    touch mitobam.bam.bai
    """

}


//
