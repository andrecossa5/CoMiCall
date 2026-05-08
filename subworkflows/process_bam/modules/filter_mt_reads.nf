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
    REGION="${params.region}"
    CONTIG="\${REGION%%:*}"
    HEADER_CONTIGS=\$(samtools view -H ${bam} | awk '\$1=="@SQ"{for(i=2;i<=NF;i++) if(\$i ~ /^SN:/){sub("SN:","",\$i); print \$i}}')
    if ! echo "\$HEADER_CONTIGS" | grep -qx "\$CONTIG"; then
        if [ "\$CONTIG" = "chrM" ]; then
            ALT="MT"
        elif [ "\$CONTIG" = "MT" ]; then
            ALT="chrM"
        else
            echo "Contig \$CONTIG not in BAM header and no fallback available" >&2
            exit 1
        fi
        if echo "\$HEADER_CONTIGS" | grep -qx "\$ALT"; then
            REGION="\${ALT}\${REGION#\$CONTIG}"
        else
            echo "Neither \$CONTIG nor \$ALT found in BAM header" >&2
            exit 1
        fi
    fi
    samtools view ${bam} -b -@ ${task.cpus} \$REGION > mitobam.bam
    samtools index -@ ${task.cpus} mitobam.bam
    """

    stub:
    """
    touch mitobam.bam
    touch mitobam.bam.bai
    """

}


//
