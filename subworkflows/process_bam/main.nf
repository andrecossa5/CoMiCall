// process_bam subworkflow

nextflow.enable.dsl = 2

include { INDEX_BAM }                 from "./modules/index_bam.nf"
include { FILTER_MT_READS }           from "./modules/filter_mt_reads.nf"
include { REALIGN_MT_READS }          from "./modules/realign_mt_reads.nf"
include { MAKE_TABLES }               from "./modules/make_tables.nf"
include { ESTIMATE_MT_CN }            from "./modules/estimate_mt_cn.nf"
include { GATHER_TABLES }             from "./modules/gather_tables.nf"
include { FILTER_MUTS }               from "./modules/filter_muts.nf"
include { GENOTYPE_MT }               from "./modules/genotype_mt.nf"


//


//----------------------------------------------------------------------------//
// process_bam subworkflow
//----------------------------------------------------------------------------//

workflow process_bam {

    take:
        input_bam_ch   // val(sample), path(bam)  — one element per BAM file

    main:

        INDEX_BAM(input_bam_ch)

        // ----------------------------------------------------------------
        // Single-molecule strand-concordance flow:
        //   1. Filter each colony BAM to MT reads only
        //   2. Scan each mitobam read-by-read, tally A/C/G/T fw+rev
        //      counts and mean qualities at every position (MQ>=30, BQ>=30)
        //   3. Gather all per-colony sparse matrices into a single AnnData
        // ----------------------------------------------------------------

        FILTER_MT_READS(INDEX_BAM.out.indexed_bams)
        REALIGN_MT_READS(FILTER_MT_READS.out.mitobam)
        MAKE_TABLES(REALIGN_MT_READS.out.mitobam)
        mt_cn_input_ch = INDEX_BAM.out.indexed_bams
                         .join(REALIGN_MT_READS.out.mitobam)
        ESTIMATE_MT_CN(mt_cn_input_ch)

        all_tables_ch = MAKE_TABLES.out.tables
                        .flatMap { sample, allelic, coverage -> [allelic, coverage] }
                        .collect()
        mt_cn_ch = ESTIMATE_MT_CN.out.mt_cn.collect()

        GATHER_TABLES(all_tables_ch, mt_cn_ch)
        FILTER_MUTS(GATHER_TABLES.out.output)

        // Genotype using the annotated + confident tables
        genotype_input_ch = FILTER_MUTS.out.filtered
            .map { annotated, confident, metrics -> tuple(annotated, confident) }
        GENOTYPE_MT(genotype_input_ch)

    emit:
        output    = FILTER_MUTS.out.filtered
        genotyped = GENOTYPE_MT.out.genotyped

}


//
