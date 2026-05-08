// nf-mito-pileup: MT pileup from single-cell colony BAM files
nextflow.enable.dsl = 2

include { process_bam } from "./subworkflows/process_bam/main"


//


// Help
if (params.help) {

    println """\
    ==========================================================================
    Single-molecule MT-SNVs from single-cell colony WGS data.
    ==========================================================================

    USAGE:
        nextflow run main.nf -params-file <params.json> -profile <profile>
        nextflow run main.nf -params-file <params.json> -profile <profile> -c config/user.config

    INPUT/OUTPUT OPTIONS:
        --input_folder     Directory containing per-colony deduplicated .bam files (REQUIRED)
        --individual       Donor individual identifier (used in output filenames) (REQUIRED)
        --output_folder    Output directory path (REQUIRED)

    REFERENCE OPTIONS:
        --ref              Faidx-indexed reference FASTA containing chrM
                           (used by make_tables / mpileup)  (REQUIRED)
        --masked_ref       NUMT-masked, bwa-indexed full-genome FASTA — see
                           code/mask_numts.sh. Required only when
                           --single_molecule true (REALIGN_MT_READS target).

    QUALITY / DEPTH OPTIONS (shared by both modes):
        --region           Genomic region to pileup (default: ${params.region})
        --min_mq           Minimum mapping quality (default: ${params.min_mq})
        --min_bq           Minimum base quality (default: ${params.min_bq})
        --max_depth        Max per-file depth, default mode only (default: ${params.max_depth})
        --max_idepth       Max per-file depth for INDELs, default mode only (default: ${params.max_idepth})

    SINGLE-MOLECULE MODE:
        --single_molecule  Enable single-molecule base-calls (default: ${params.single_molecule})
        --max_softclip     Max fraction of softclipped bases per read in
                           make_tables Phase-1 filter (default: ${params.max_softclip})

    EXECUTION PROFILES (choose one):
        -profile conda         Use Conda environments
        -profile singularity   Use Singularity/Apptainer containers
        -profile local         Run locally
        -profile lsf           Submit jobs to LSF

    EXAMPLES:

        # Run locally (stub test)
        nextflow run main.nf \\
            -params-file params/params_test.json \\
            -profile local -stub

        # Run on HPC with LSF + conda
        nextflow run main.nf \\
            -params-file params/params.json \\
            -profile conda,lsf \\
            -c config/user.config

    ==========================================================================
    """.stripIndent()
    exit(0)

}

// Version
if (params.version) {
    println "Version: ${workflow.manifest.version}"
    exit(0)
}

// Parameter validation
if (!params.help && !params.version) {
    if (!params.input_folder)  error "Error: --input_folder is required. Use --help for usage information."
    if (!params.individual)    error "Error: --individual is required. Use --help for usage information."
    if (!params.output_folder) error "Error: --output_folder is required. Use --help for usage information."
    if (!params.ref)           error "Error: --ref is required. Use --help for usage information."
    if (params.single_molecule && !params.masked_ref) {
        error "Error: --masked_ref is required when --single_molecule is true (NUMT-masked, bwa-indexed full genome; see code/mask_numts.sh)."
    }
}


//


//----------------------------------------------------------------------------//
// nf-mito-pileup main workflow
//----------------------------------------------------------------------------//

workflow {

    println "\n"
    println "Single-molecule MT-SNVs from single-cell colony WGS data."
    println "Usage: nextflow run main.nf -params-file <params.json> -profile <profile>"
    println "\n"

    // Discover BAMs from input folder; sample name = filename sans .bam
    input_bam_ch = Channel.fromPath("${params.input_folder}/*.bam")
        .map { bam -> tuple(bam.baseName, bam) }

    // Index BAMs, group all, run joint MT pileup
    process_bam(input_bam_ch)

}


//
