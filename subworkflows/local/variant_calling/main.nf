/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: VARIANT_CALLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This subworkflow performs comprehensive variant calling using multiple tools
    including GATK4, Clair3, Longshot, and PyPGx. It processes BWA alignments
    and generates variant calls for downstream analysis and annotation.

    Tools included:
    - BCFTOOLS (ANNOTATE, CONCAT, FILTER, NORM, SORT) - VCF processing and filtering
    - CLAIR3 - ONT variant calling
    - CLEAN_VCF_HEADER - VCF header cleanup
    - GATK4 (ADDORREPLACEREADGROUPS, CALIBRATEDRAGSTRMODEL, HAPLOTYPECALLER, MARKDUPLICATES) - Variant calling pipeline
    - LONGSHOT - Long-read variant calling
    - PYPGX_CREATEINPUTVCF - Pharmacogenomic variant calling
    - TABIX_BGZIPTABIX - VCF indexing
    - VCFTOOLS - VCF manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_ANNOTATE               } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_CONCAT                 } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_FILTER                 } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_NORM                   } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_SORT                   } from '../../../modules/nf-core/bcftools/sort/main'
include { CLAIR3                          } from '../../../modules/nf-core/clair3/main'
include { CLEAN_VCF_HEADER                } from '../../../modules/local/cleanvcfheader/main'
include { GATK4_ADDORREPLACEREADGROUPS    } from '../../../modules/nf-core/gatk4/addorreplacereadgroups/main'
include { GATK4_CALIBRATEDRAGSTRMODEL     } from '../../../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { GATK4_HAPLOTYPECALLER           } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MARKDUPLICATES            } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { LONGSHOT                        } from '../../../modules/local/longshot/main'
include { PYPGX_CREATEINPUTVCF            } from '../../../modules/nf-core/pypgx/createinputvcf/main'
include { TABIX_BGZIPTABIX                } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow VARIANT_CALLING {

    take:
    ch_bam_bai
    ch_clair3_model_dir
    ch_cyp2d6_bed
    ch_dbsnp_tbi
    ch_dbsnp_vcf
    ch_genome_dict
    ch_genome_fai
    ch_genome_fasta
    ch_genome_str_table
    vcf_header

    main:
    ch_versions = Channel.empty()

    /*
    MODULE: PYPGX_CREATEINPUTVCF
    */
    PYPGX_CREATEINPUTVCF (
        ch_bam_bai,
        ch_genome_fasta
    )
    ch_versions = ch_versions.mix(PYPGX_CREATEINPUTVCF.out.versions.first())

    /*
    MODULE: CLEAN_VCF_HEADER
    */
    CLEAN_VCF_HEADER (
        PYPGX_CREATEINPUTVCF.out.vcf
    )
    ch_versions = ch_versions.mix(CLEAN_VCF_HEADER.out.versions.first())

    /*
    MODULE: CLAIR3
    */
    ch_clair3_input = ch_bam_bai
        .combine(ch_clair3_model_dir)
        .map { meta, bam, bai, model_dir -> [meta, bam, bai, "", "${model_dir}/r1041_e82_400bps_sup_v420", "ont"] }

    CLAIR3 (
        ch_clair3_input,
        ch_genome_fasta,
        ch_genome_fai
    )
    ch_versions = ch_versions.mix(CLAIR3.out.versions.first())

    /*
    MODULE: LONGSHOT
    */
    LONGSHOT (
        ch_bam_bai.map { meta, bam, bai -> [meta, bam] },
        ch_bam_bai.map { meta, bam, bai -> [meta, bai] },
        ch_genome_fasta,
        ch_genome_fai
    )
    ch_versions = ch_versions.mix(LONGSHOT.out.versions.first())

    /*
    MODULE: TABIX_BGZIPTABIX | Index Longshot VCF files with bgzip and tabix
    */
    TABIX_BGZIPTABIX (
        LONGSHOT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

    /*
    MODULE: GATK4_ADDORREPLACEREADGROUPS
    */
    GATK4_ADDORREPLACEREADGROUPS (
        ch_bam_bai.map { meta, bam, bai -> [meta, bam] },
        ch_genome_fasta,
        ch_genome_fai
    )
    ch_versions = ch_versions.mix(GATK4_ADDORREPLACEREADGROUPS.out.versions.first())

    /*
    MODULE: GATK4_MARKDUPLICATES
    */
    GATK4_MARKDUPLICATES (
        GATK4_ADDORREPLACEREADGROUPS.out.bam,
        ch_genome_fasta.map { meta, fasta -> fasta },
        ch_genome_fai.map { meta, fai -> fai }
    )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())

    /*
    MODULE: GATK4_CALIBRATEDRAGSTRMODEL
    */
    GATK4_CALIBRATEDRAGSTRMODEL (
        GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai),
        ch_genome_fasta.map { meta, fasta -> fasta },
        ch_genome_fai.map { meta, fai -> fai },
        ch_genome_dict.map { meta, dict -> dict },
        ch_genome_str_table
    )
    ch_versions = ch_versions.mix(GATK4_CALIBRATEDRAGSTRMODEL.out.versions.first())

    /*
    MODULE: GATK4_HAPLOTYPECALLER || NF-CORE
    */
    GATK4_HAPLOTYPECALLER (
        GATK4_MARKDUPLICATES.out.bam
            .join(GATK4_MARKDUPLICATES.out.bai, by:0, failOnDuplicate:true, failOnMismatch:true)
            .join(GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model, by:0, failOnDuplicate:true, failOnMismatch:true)
            .combine(ch_cyp2d6_bed)
            .map { meta, bam, bai, dragstr, bed_file -> [ meta, bam, bai, bed_file, dragstr ] },
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dict,
        ch_dbsnp_vcf,
        ch_dbsnp_tbi
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    /*
    Combine all VCF and TBI channels for concatenation
    */
    CLAIR3.out.phased_vcf
        .join(CLAIR3.out.phased_tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        .set { ch_from_clair3 }

    GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        .set { ch_from_gatk4 }

    // Longshot VCF and TBI files
    TABIX_BGZIPTABIX.out.gz_tbi.map { meta, vcf, tbi -> [meta, vcf, tbi] }
        .set { ch_from_longshot }

    CLEAN_VCF_HEADER.out.vcf
        .join(CLEAN_VCF_HEADER.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        .set { ch_from_pypgx }

    ch_from_clair3
        .mix(
            ch_from_longshot,
            ch_from_gatk4,
            ch_from_pypgx
        )
        .groupTuple()
        .set { ch_combined_vcfs }
        // .view()

    /*
    MODULE: BCFTOOLS_CONCAT
        bcftools concat - adds rows to a combined vcf by collating variant calls of same sample. e.g.
        variant calls of same individual using different tools.
        bcftools merge - creates a superset of variant calls across multiple individuals/samples.
    */
    BCFTOOLS_CONCAT ( ch_combined_vcfs )
    ch_versions = ch_versions.mix( BCFTOOLS_CONCAT.out.versions.first() )

    /*
    MODULE: BCFTOOLS_NORM
    */

    BCFTOOLS_NORM (
        BCFTOOLS_CONCAT.out.vcf
            .join(BCFTOOLS_CONCAT.out.tbi),
        ch_genome_fasta
    )
    ch_versions = ch_versions.mix( BCFTOOLS_NORM.out.versions.first() )

    /*
    MODULE: BCFTOOLS_SORT
    */
    BCFTOOLS_SORT (BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix( BCFTOOLS_SORT.out.versions.first() )

    /*
    PROCESS CREATE_VCF_HEADER
    */
    CREATE_VCF_HEADER ( vcf_header )

    /*
    MODULE: BCFTOOLS_ANNOTATE
        Add dbSNP data to the VCF file.
    */
    BCFTOOLS_ANNOTATE (
        BCFTOOLS_SORT.out.vcf
            .join( BCFTOOLS_SORT.out.tbi, failOnDuplicate:true, failOnMismatch:true )
            .combine(ch_dbsnp_vcf.map { meta, vcf -> vcf })
            .combine(ch_dbsnp_tbi.map { meta, tbi -> tbi })
            .map { meta, vcf, tbi, dbsnp_vcf, dbsnp_tbi ->
                [meta, vcf, tbi, dbsnp_vcf, dbsnp_tbi]
            },
        CREATE_VCF_HEADER.out.header,
        []
    )
    ch_versions = ch_versions.mix( BCFTOOLS_ANNOTATE.out.versions.first() )

    /*
    MODULE: BCFTOOLS_FILTER
    */
    BCFTOOLS_FILTER (BCFTOOLS_ANNOTATE.out.vcf.join(BCFTOOLS_ANNOTATE.out.tbi))
    ch_versions = ch_versions.mix( BCFTOOLS_FILTER.out.versions.first() )


    emit:
    clair3_vcf_tbi     = ch_from_clair3                                                    // channel: [val(meta), path(vcf), path(tbi)] - Clair3
    combined_vcf_tbi   = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi)              // channel: [val(meta), path(vcf), path(tbi)] - Merged
    gatk4_vcf_tbi      = GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi)  // channel: [val(meta), path(vcf), path(tbi)] - GATK4
    longshot_vcf_tbi   = TABIX_BGZIPTABIX.out.gz_tbi                                        // channel: [val(meta), path(vcf), path(tbi)] - Longshot
    pypgx_vcf_tbi      = CLEAN_VCF_HEADER.out.vcf.join(CLEAN_VCF_HEADER.out.tbi)        // channel: [val(meta), path(vcf), path(tbi)] - PyPGx
    versions           = ch_versions                                                        // channel: [path(versions.yml)]
}

/*
CUSTOM PROCESSES REQUIRED BY SUBWORKFLOW
    Others were covereted into local modules
*/
process CREATE_VCF_HEADER {
    tag "Create VCF header"
    input:
    val vcf_header

    output:
    path 'vcf_header.txt', emit: header

    script:
    """
    echo "${vcf_header}" | sed 's/^[ \t]*//' > vcf_header.txt
    """
}
