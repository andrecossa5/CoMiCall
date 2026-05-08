// REALIGN_MT_READS module
//
// Re-aligns the mtDNA-retrieved reads to a NUMT-masked full-genome reference
// (mgatk-style, Lareau et al. Nat Protoc 2022). Reads of NUMT origin lose
// their nuclear landing site (masked out) and fail to produce a confident
// alignment anywhere, so they are dropped by the MAPQ/flag filter. Only
// reads that map confidently to chrM in the full-genome context are kept.
//
// ${params.masked_ref} is the masked full genome produced by code/mask_numts.sh
// (chrM unmasked, nuclear NUMT coords N'd, bwa index beside the fasta).
// Kept separate from ${params.ref} so that the bcftools / chrM-only branch of
// the pipeline does not need a full-genome bwa index.

nextflow.enable.dsl = 2

//

process REALIGN_MT_READS {

    tag "${sample}"
    label 'bwa'

    input:
    tuple val(sample),
          path(mitobam_in,     stageAs: 'input.bam'),
          path(mitobam_in_bai, stageAs: 'input.bam.bai')

    output:
    tuple val(sample), path("mitobam.bam"), path("mitobam.bam.bai"), emit: mitobam

    script:
    """
    samtools collate -O -@ ${task.cpus} ${mitobam_in} | samtools fastq -@ ${task.cpus} -1 r1.fq -2 r2.fq -0 /dev/null -s /dev/null -n -

    bwa mem -t ${task.cpus} ${params.masked_ref} r1.fq r2.fq | samtools sort -@ ${task.cpus} -o realigned.bam -

    samtools index -@ ${task.cpus} realigned.bam

    # Assume alignment on hg38 reference, with chrM as the mitochondrial chromosome
    samtools view -b -@ ${task.cpus} -q ${params.min_mq} realigned.bam chrM > mitobam.bam
    samtools index -@ ${task.cpus} mitobam.bam
    """

    stub:
    """
    touch mitobam.bam
    touch mitobam.bam.bai
    """

}


//
