#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/cyptyper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/cyptyper
    Website: https://nf-co.re/cyptyper
    Slack  : https://nfcore.slack.com/channels/cyptyper
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CYPTYPER                } from './workflows/cyptyper'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_cyptyper_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_cyptyper_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_cyptyper_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.bed               = getGenomeAttribute('bed')
params.bwamem2           = getGenomeAttribute('bwamem2')
params.dbsnp             = getGenomeAttribute('dbsnp')
params.dbsnp_vqsr        = getGenomeAttribute('dbsnp_vqsr')
params.fasta             = getGenomeAttribute('fasta')
params.fasta_fai         = getGenomeAttribute('fasta_fai')
params.intervals         = getGenomeAttribute('intervals')
params.known_indels      = getGenomeAttribute('known_indels')
params.known_indels_tbi  = getGenomeAttribute('known_indels_tbi')
params.known_indels_vqsr = getGenomeAttribute('known_indels_vqsr')
params.known_snps        = getGenomeAttribute('known_snps')
params.known_snps_tbi    = getGenomeAttribute('known_snps_tbi')
params.known_snps_vqsr   = getGenomeAttribute('known_snps_vqsr')
params.pypgx_version     = getGenomeAttribute('pypgx_version')
params.snpeff_db         = getGenomeAttribute('snpeff_db')
params.snpeff_genome     = getGenomeAttribute('snpeff_genome')
params.model_url         = getGenomeAttribute('model_url')
params.vep_cache_version = getGenomeAttribute('vep_cache_version')
params.vep_genome        = getGenomeAttribute('vep_genome')
params.vep_species       = getGenomeAttribute('vep_species')


// Initialize fasta file with meta map:
fasta             = params.fasta ? Channel.fromPath(params.fasta).map { it -> [[id: it.baseName], it] }.collect() : Channel.empty()
fasta_fai         = params.fasta_fai ? Channel.fromPath(params.fasta_fai).map { it -> [[id: it.baseName], it] }.collect() : Channel.empty()

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
dbsnp             = params.dbsnp ? Channel.fromPath(params.dbsnp).collect() : Channel.empty()
genome_size       = params.genome_size ? Channel.fromPath(params.genome_size).collect() : Channel.empty()
bed               = params.bed ? Channel.fromPath(params.bed).collect() : Channel.empty()
known_indels      = params.known_indels ? Channel.fromPath(params.known_indels).collect() : Channel.value([])
known_snps        = params.known_snps ? Channel.fromPath(params.known_snps).collect() : Channel.value([])
vcf_header        = params.vcf_header ? Channel.fromPath(params.vcf_header).collect() :  
     Channel.value(
        """
        ##fileformat=VCFv4.2
        ##FILTER=<ID=PASS,Description=\\"All filters passed\\">
        """
    )

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
dbsnp_vqsr        = params.dbsnp_vqsr ? Channel.value(params.dbsnp_vqsr) : Channel.empty()
known_indels_vqsr = params.known_indels_vqsr ? Channel.value(params.known_indels_vqsr) : Channel.empty()
known_snps_vqsr   = params.known_snps_vqsr ? Channel.value(params.known_snps_vqsr) : Channel.empty()
model_url         = params.model_url ?: Channel.empty()
pypgx_version     = params.pypgx_version ? Channel.value(params.pypgx_version) :  Channel.value("0.25.0") 
snpeff_db         = params.snpeff_db ?: Channel.empty()
vep_cache_version = params.vep_cache_version ?: Channel.empty()
vep_genome        = params.vep_genome ?: Channel.empty()
vep_species       = params.vep_species ?: Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_CYPTYPER {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    CYPTYPER (
        samplesheet,
        fasta,
        fasta_fai,
        dbsnp,
        genome_size,
        bed,
        model_url,
        vcf_header,
        pypgx_version
    )

    emit:
    multiqc_report = CYPTYPER.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CYPTYPER (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_CYPTYPER.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
