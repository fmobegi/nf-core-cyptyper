/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                                } from '../modules/nf-core/fastqc/main'
include { FASTQ_ALIGN_BWA                       } from '../subworkflows/nf-core/fastq_align_bwa/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { PREPARE_REFERENCE_INDEXES             } from '../subworkflows/local/prepare_reference_indexes/main'
include { BAM_STATS_BEDTOOLS_DEPTH              } from '../subworkflows/local/bam_stats_samtools_bedtool/main'
include { VARIANT_CALLING                       } from '../subworkflows/local/variant_calling/main'
include { VARIANT_ANNOTATION                    } from '../subworkflows/local/variant_annotation/main'
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_cyptyper_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CYPTYPER {

    take:
    ch_samplesheet  // channel: samplesheet read in from --input

    ch_cyp2d6_bed   // channel: [ val(meta), path(bed) ]
    ch_dbsnp        // channel: [ val(meta), path(dbsnp) ]
    ch_fasta        // channel: [ val(meta), path(fasta) ]
    fasta_fai       // channel: [ val(meta), path(fasta.fai) ]
    genome_size     // channel: [ val(meta), path(genome_size) ]
    model_url       // channel: [ val(meta), path(model_url) ]
    pypgx_version   // channel: [ val(meta), path(pypgx_version) ]
    vcf_header      // channel: [ val(meta), path(vcf_header) ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // /*
    // MODULE: FASTQC
    // */
    // FASTQC (
    //     ch_samplesheet
    // )

    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    /*
    SUBWORKFLOW: PREPARE_REFERENCE_INDEXES - LOCAL
    */
    PREPARE_REFERENCE_INDEXES (
        ch_fasta,                      // channel: [mandatory] fasta
        ch_dbsnp,                      // channel: [mandatory] dbsnp
        model_url,                     // channel: [mandatory] model_url
        pypgx_version                  // channel: [mandatory] pypgx_version
    )
    ch_versions = ch_versions.mix( PREPARE_REFERENCE_INDEXES.out.versions.first() )

    /*
    SUBWORKFLOW: FASTQ_ALIGN_BWA - NF-CORE
    */
    FASTQ_ALIGN_BWA (
        ch_samplesheet,                                  // channel input reads: [ val(meta2), path(index) ]
        PREPARE_REFERENCE_INDEXES.out.bwa_index,         // channel BWA index: [ val(meta2), path(index) ]
        true,                                            // boolean value: true/false for sorting BAM files
        PREPARE_REFERENCE_INDEXES.out.fasta_ref,                                           // channel reference fasta: [ val(meta3), path(fasta) ]
    )
    ch_versions = ch_versions.mix( FASTQ_ALIGN_BWA.out.versions.first() )

    combined_bam_channel  = FASTQ_ALIGN_BWA.out.bam.mix( FASTQ_ALIGN_BWA.out.bai)
    ch_bam_bai = combined_bam_channel
        .groupTuple(by: 0)
        .flatMap { meta, path_list ->
            // Split the path_list into BAM and BAI files
            def bam_files = path_list.findAll { it.name.endsWith('.bam') }
            def bai_files = path_list.findAll { it.name.endsWith('.bai') }
            // Create combinations of BAM and BAI files and return as tuples
            def bam_bai_combinations = bam_files.collectMany { bam ->
                bai_files.collect { bai ->
                    return [meta, bam, bai]
                }
            }

            return bam_bai_combinations
        }

    /*
    SUBWORKFLOW: BAM_STATS_BEDTOOLS_DEPTH
    */
    BAM_STATS_BEDTOOLS_DEPTH (
        ch_bam_bai,                                     // channel: [ val(meta), path(bam), path(bai) ]
        genome_size,
        ch_cyp2d6_bed
    )
    ch_versions = ch_versions.mix(BAM_STATS_BEDTOOLS_DEPTH.out.versions)

    /*
    SUBWORKFLOW: VARIANT_CALLING
    */
    VARIANT_CALLING (
        ch_bam_bai,
        PREPARE_REFERENCE_INDEXES.out.clair3_model_dir,
        ch_cyp2d6_bed,
        PREPARE_REFERENCE_INDEXES.out.dbsnp_tbi,
        PREPARE_REFERENCE_INDEXES.out.dbsnp_vcf,
        PREPARE_REFERENCE_INDEXES.out.dict,
        PREPARE_REFERENCE_INDEXES.out.fasta_fai,
        PREPARE_REFERENCE_INDEXES.out.fasta_ref,
        PREPARE_REFERENCE_INDEXES.out.str_table,
        vcf_header
    )
    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)

    /*
    SUBWORKFLOW: VARIANT_ANNOTATION
    */
    VARIANT_ANNOTATION (
        VARIANT_CALLING.out.combined_vcf_tbi,
        VARIANT_CALLING.out.pypgx_vcf_tbi,
        PREPARE_REFERENCE_INDEXES.out.pypgx_bundle
    )
    ch_versions = ch_versions.mix(VARIANT_ANNOTATION.out.versions)

    /*
    Collate and save software versions
    */
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'cyptyper_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    /*
    MODULE: MultiQC
    */
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
