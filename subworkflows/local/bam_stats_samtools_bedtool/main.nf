/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: BAM_STATS_BEDTOOLS_DEPTH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This subworkflow performs additional BAM file statistics that are not covered
    by the FASTQ_ALIGN_BWA:BAM_SORT_STATS_SAMTOOLS subworkflow, specifically coverage
    analysis and depth calculations for quality assessment of aligned sequencing data.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDTOOLS_GENOMECOV    } from '../../../modules/nf-core/bedtools/genomecov/main'
include { SAMTOOLS_AVERAGEDEPTH } from '../../../modules/local/samtools/averagedepth/main'
include { SAMTOOLS_DEPTH        } from '../../../modules/nf-core/samtools/depth/main'

workflow BAM_STATS_BEDTOOLS_DEPTH {

    take:
    ch_bam_bai               // channel: [ val(meta), path(bam), path(bai) ]
    genome_size              // value: genome size
    ch_cyp2d6_bed            // path: bed file

    main:
    ch_versions = Channel.empty()

    // Extract BAM files from the joined channel
    ch_bam = ch_bam_bai.map { meta, bam, bai -> [meta, bam] }

    /*
    MODULE: BEDTOOLS_GENOMECOV
    */
    BEDTOOLS_GENOMECOV (
        ch_bam.map { meta, bam ->
            [meta, bam, 1]
        },
        genome_size,
        ".txt",
        true
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())
    
    /*
    MODULE: SAMTOOLS_DEPTH
    */
    // Prepare the bed file with metadata for SAMTOOLS_DEPTH

    def ch_bed_with_meta = ch_cyp2d6_bed.map { file -> [[id: 'intervals'], file] }
    
    SAMTOOLS_DEPTH (
        ch_bam,
        ch_bed_with_meta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    /*
    MODULE: SAMTOOLS_AVERAGEDEPTH
    */
    SAMTOOLS_AVERAGEDEPTH(SAMTOOLS_DEPTH.out.tsv)
    ch_versions = ch_versions.mix(SAMTOOLS_AVERAGEDEPTH.out.versions.first())

    emit:
    average    = SAMTOOLS_AVERAGEDEPTH.out.avg_depth     // channel: [ val(meta), path(avg_depth) ]
    depth      = SAMTOOLS_DEPTH.out.tsv                  // channel: [ val(meta), path(depth) ]
    genomecov  = BEDTOOLS_GENOMECOV.out.genomecov        // channel: [ val(meta), path(genomecov) ]
    versions   = ch_versions                             // channel: [ path(versions.yml) ]
}