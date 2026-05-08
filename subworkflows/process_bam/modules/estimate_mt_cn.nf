// ESTIMATE_MT_CN module
//
// mtDNA copy number via the PCAWG formula:
//   CN = coverage(mtDNA) / coverage(nDNA) * ploidy   (ploidy = 2).
// Mean depths are computed with bedtools genomecov.

nextflow.enable.dsl = 2

//

process ESTIMATE_MT_CN {

    tag "${sample}"
    label 'bwa'

    input:
    tuple val(sample), path(nuclear_bam), path(nuclear_bai), path(mitobam), path(mitobam_bai)

    output:
    path("${sample}.mt_cn.txt"), emit: mt_cn

    script:
    """
    nuc=\$(bedtools genomecov -ibam ${nuclear_bam} \\
        | awk '\$1=="genome" && \$2>0 {sum+=\$2*\$3} \$1=="genome"{tot=\$4} END{if(tot>0) print sum/tot; else print 0}')

    mt=\$(bedtools genomecov -ibam ${mitobam} \\
        | awk -v r=${params.region} '\$1==r && \$2>0 {sum+=\$2*\$3} \$1==r{tot=\$4} END{if(tot>0) print sum/tot; else print 0}')

    cn=\$(awk -v m=\$mt -v n=\$nuc 'BEGIN{ if(n>0) print (m/n)*2; else print 0 }')

    printf "%s\\t%s\\t%s\\t%s\\n" "${sample}" "\$nuc" "\$mt" "\$cn" > ${sample}.mt_cn.txt
    """

    stub:
    """
    touch ${sample}.mt_cn.txt
    """

}


//
