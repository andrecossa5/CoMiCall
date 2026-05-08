// process_bam subworkflow

nextflow.enable.dsl = 2

include { INDEX_BAM }                 from "./modules/index_bam.nf"
include { MPILEUP }                   from "./modules/mpileup.nf"
include { FILTER_MT_READS }           from "./modules/filter_mt_reads.nf"
include { REALIGN_MT_READS }          from "./modules/realign_mt_reads.nf"
include { MAKE_TABLES }               from "./modules/make_tables.nf"
include { ESTIMATE_MT_CN }            from "./modules/estimate_mt_cn.nf"
include { GATHER_TABLES }             from "./modules/gather_tables.nf"


//


//----------------------------------------------------------------------------//
// process_bam subworkflow
//----------------------------------------------------------------------------//

workflow process_bam {

    take:
        input_bam_ch   // val(sample), path(bam)  — one element per BAM file

    main:

        // Index each BAM independently (needed by both branches)
        INDEX_BAM(input_bam_ch)

        // Initialise output channel — populated by whichever branch runs
        output_ch = Channel.empty()

        if ( !params.single_molecule ) {

            // ----------------------------------------------------------------
            // Default branch: joint bcftools MT pileup → VCF
            // ----------------------------------------------------------------

            // Collect all indexed BAMs into a single channel element:
            // tuple( [sample, ...], [bam, ...], [bai, ...] )
            grouped_ch = INDEX_BAM.out.indexed_bams
                .collect()
                .map { items ->
                    def samples = items.collate(3).collect { it[0] }
                    def bams    = items.collate(3).collect { it[1] }
                    def bais    = items.collate(3).collect { it[2] }
                    tuple(samples, bams, bais)
                }

            MPILEUP(grouped_ch)
            output_ch = MPILEUP.out.pileup_vcf

        } else {

            // ----------------------------------------------------------------
            // Strand-concordance branch:
            //   1. Filter each colony BAM to MT reads only
            //   2. Scan each mitobam read-by-read, tally A/C/G/T fw+rev
            //      counts and mean qualities at every position (MQ>=30, BQ>=30)
            //   3. Gather all per-colony sparse matrices into a single AnnData
            // ----------------------------------------------------------------

            FILTER_MT_READS(INDEX_BAM.out.indexed_bams)
            REALIGN_MT_READS(FILTER_MT_READS.out.mitobam)
            MAKE_TABLES(REALIGN_MT_READS.out.mitobam)

            // Estimate mtDNA copy number per colony: needs the nuclear-level
            // bam (from INDEX_BAM) joined with the realigned mitobam.
            mt_cn_input_ch = INDEX_BAM.out.indexed_bams
                .join(REALIGN_MT_READS.out.mitobam)
            ESTIMATE_MT_CN(mt_cn_input_ch)

            // Flatten the per-colony tuple to individual file paths, then
            // collect everything before passing to the aggregation step
            all_tables_ch = MAKE_TABLES.out.tables
                .flatMap { sample, allelic, coverage -> [allelic, coverage] }
                .collect()
            mt_cn_ch = ESTIMATE_MT_CN.out.mt_cn.collect()

            GATHER_TABLES(all_tables_ch, mt_cn_ch)
            output_ch = GATHER_TABLES.out.output

        }

    emit:
        output = output_ch

}


//
