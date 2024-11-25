//
// VARIANT CALLING GATK4
//

include { GATK4_ADDORREPLACEREADGROUPS    } from '../../../modules/nf-core/gatk4/addorreplacereadgroups/main'
include { GATK4_CALIBRATEDRAGSTRMODEL     } from '../../../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { GATK4_COMBINEGVCFS              } from '../../../modules/nf-core/gatk4/combinegvcfs/main'
include { GATK4_CREATESEQUENCEDICTIONARY  } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_HAPLOTYPECALLER           } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MARKDUPLICATES            } from '../../../modules/nf-core/gatk4/markduplicates/main'

workflow GATK4 {

    take:
        bwa_bam
        bwa_bai
        fasta
        fasta_fai
        genome_dict
        dbsnp
        dbsnp_tbi
        bed
        str_table

    main:
    versions = Channel.empty()

    // AddOrReplaceReadGroups
    GATK4_ADDORREPLACEREADGROUPS (
        bwa_bam,
        fasta,
        fasta_fai
    )
    versions = versions.mix( GATK4_ADDORREPLACEREADGROUPS.out.versions.first() )

    // Mark Duplicates
    GATK4_MARKDUPLICATES (
        GATK4_ADDORREPLACEREADGROUPS.out.bam,
        fasta,
        fasta_fai
    )
    versions = versions.mix( GATK4_MARKDUPLICATES.out.versions.first() )

    /* 
    CreateSequenceDictionary
    uses chanel genome_dict gereated by prepare genome GATK4_CREATESEQUENCEDICTIONARY 
    */

    // Mark CalibrateDragstrModel
    GATK4_CALIBRATEDRAGSTRMODEL (
        GATK4_MARKDUPLICATES.out.bam,
        GATK4_MARKDUPLICATES.out.bai,
        fasta,
        fasta_fai,
        genome_dict,
        str_table
    )
    versions = versions.mix( GATK4_CALIBRATEDRAGSTRMODEL.out.versions.first() )

    // gatk haplotypeCaller
    GATK4_HAPLOTYPECALLER (
        GATK4_MARKDUPLICATES.out.bam,
        GATK4_MARKDUPLICATES.out.bai,
        bed,
        GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model,
        fasta,
        fasta_fai,
        genome_dict,
        dbsnp,
        dbsnp_tbi
    )
    versions = versions.mix( GATK4_HAPLOTYPECALLER.out.versions.first() )

    emit:
    vcf      = GATK4_HAPLOTYPECALLER.out.vcf           // channel: [ val(meta), [ bam ] ]
    tbi      = GATK4_HAPLOTYPECALLER.out.tbi           // channel: [ val(meta), [ bai ] ]
    bam      = GATK4_HAPLOTYPECALLER.out.bam           // channel: [ val(meta), [ bai ] ]
    versions
}
