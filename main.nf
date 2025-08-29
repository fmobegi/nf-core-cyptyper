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
params.fasta             = getGenomeAttribute('fasta')
params.fasta_fai         = getGenomeAttribute('fasta_fai')
params.intervals         = getGenomeAttribute('intervals')
params.pypgx_version     = getGenomeAttribute('pypgx_version')
params.model_url         = getGenomeAttribute('model_url')


// Initialize fasta file with meta map:
fasta             = params.fasta ? Channel.value([[id: file(params.fasta).baseName], file(params.fasta)]) : Channel.empty()
fasta_fai         = params.fasta_fai ? Channel.value([[id: file(params.fasta_fai).baseName], file(params.fasta_fai)]) : Channel.empty()
dbsnp             = params.dbsnp ? Channel.value([[id: file(params.dbsnp).baseName], file(params.dbsnp)]) : Channel.empty()
genome_size       = params.genome_size ? Channel.value(file(params.genome_size))    : Channel.empty()
bed               = params.bed ? Channel.value(file(params.bed))                    : Channel.empty()
vcf_header        = params.vcf_header ? Channel.value(file(params.vcf_header))      :
     Channel.value(
        """
        ##fileformat=VCFv4.2
        ##FILTER=<ID=PASS,Description=\\"All filters passed\\">
        """
    )
model_url         = params.model_url ? Channel.value(params.model_url)                  : Channel.empty()
pypgx_version     = params.pypgx_version ? Channel.value(params.pypgx_version)          :  Channel.value("0.25.0")

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
        bed,
        dbsnp,
        fasta,
        fasta_fai,
        genome_size,
        model_url,
        pypgx_version,
        vcf_header
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
