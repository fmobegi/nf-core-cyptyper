/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                                } from '../modules/nf-core/fastqc/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                      } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_cyptyper_pipeline'
include { PREPARE_GENOME                        } from '../subworkflows/local/prepare_genome_and_map/prep'
include { MAPPING_READS                         } from '../subworkflows/local/prepare_genome_and_map/map'
include { BAMSTATS as BAMSTATS_BWA              } from '../subworkflows/local/bamstats/bamstats'
include { BAMSTATS as BAMSTATS_MINIMAP2         } from '../subworkflows/local/bamstats/bamstats'
include { VARIANT_CALLING_ALL                   } from '../subworkflows/local/variant_calling_main/variant_calling'
include { VARIANT_ANNOTATION                    } from '../subworkflows/local/variant_annotation/variant_annotation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CYPTYPER {

     take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta
    fasta_fai
    dbsnp
    genome_size
    bed
    model_url
    vcf_header
    pypgx_version

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run BWA INDEX; SAMTOOLS SORT; SAMTOOLS INDEX
    //
    PREPARE_GENOME (
        fasta,                         // channel: [mandatory] fasta
        dbsnp,                         // channel: [mandatory] dbsnp
        model_url,                     // channel: [mandatory] model_url
        pypgx_version                  // channel: [mandatory] pypgx_version
    )

    ch_versions = ch_versions.mix( PREPARE_GENOME.out.versions.first() )

    //
    // MODULE: Run BWA MEM; MINIMAP2
    //
    fasta_fai         = PREPARE_GENOME.out.fasta_fai
    bwa               = PREPARE_GENOME.out.bwa
    clair3_model_dir  = PREPARE_GENOME.out.clair3_model_dir

    MAPPING_READS (
        ch_samplesheet,                 // channel for input reads
        bwa,                            // bwa_index chanel from PREPARE_GENOME
        fasta,                          // channel for input FASTA
        sort_bam="sort",                // BAM sorting for BWA
        bam_format="bam",               // BAM index extension for minimap2
        bam_index_extension="bai",      // bam index format for minimap2/samtools
        cigar_paf_format=false,         // CIGAR format for minimap2 PAF
        cigar_bam=false                 // CIGAR BAM for minimap2
    )

    ch_versions = ch_versions.mix( MAPPING_READS.out.versions.first() )

    //
    // MODULE: Run BEDTOOLS_GENOMECOV; SAMTOOLS DEPTH/FLAGSTAT/STATS
    //

    // MAPPING_READS.out.bwa_bam.view()
    BAMSTATS_BWA (
        ch_samplesheet,
        MAPPING_READS.out.bwa_bam,
        MAPPING_READS.out.bwa_bai,
        genome_size,
        bed,
        fasta
    )
    ch_versions = ch_versions.mix(BAMSTATS_BWA.out.versions.first())

    BAMSTATS_MINIMAP2 (
        ch_samplesheet,
        MAPPING_READS.out.minimap2_bam,
        MAPPING_READS.out.minimap2_bai,
        genome_size,
        bed,
        fasta
    )
    ch_versions = ch_versions.mix(BAMSTATS_MINIMAP2.out.versions.first())

    //
    // MODULE: Run VARIANT_CALLING; MEDAKA; CUTESV; SNIFFLES; MEDAKA; CLAIR3; GATK4
    //

    genome_dict = PREPARE_GENOME.out.dict
    dbsnp_tbi   = PREPARE_GENOME.out.dbsnp_tbi
    str_table   = PREPARE_GENOME.out.str_table

    // bwa_bam        = BWA_MEM.out.bam                  // bwa-mem bam file
    // bwa_bai        = BWA_MEM.out.bai                  // bwa-mem bam index file
    // minimap2_bam   = MINIMAP2_ALIGN.out.bam           // minimap2 bam file
    // minimap2_bai   = MINIMAP2_ALIGN.out.index         // minimap2 bam index file

    VARIANT_CALLING_ALL (
        ch_samplesheet,
        MAPPING_READS.out.bwa_bam,
        MAPPING_READS.out.bwa_bai,
        MAPPING_READS.out.minimap2_bam,
        MAPPING_READS.out.minimap2_bai,
        fasta,
        fasta_fai,
        genome_dict,
        dbsnp,
        dbsnp_tbi,
        bed,
        str_table,
        clair3_model_dir,
        vcf_header
    )
    ch_versions = ch_versions.mix(VARIANT_CALLING_ALL.out.versions.first())

    //
    // MODULE: Run VARIANT_ANNOTATE
    //    
    
    VARIANT_ANNOTATION (
        VARIANT_CALLING_ALL.out.vcf,
        VARIANT_CALLING_ALL.out.tbi,
        VARIANT_CALLING_ALL.out.ch_from_pypgx,
        VARIANT_CALLING_ALL.out.gene_ch,
        VARIANT_CALLING_ALL.out.assemby_ch,
        PREPARE_GENOME.out.pypgx_bundle
    )
    ch_versions = ch_versions.mix(VARIANT_ANNOTATION.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
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

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
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
        ch_multiqc_logo.toList()
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
