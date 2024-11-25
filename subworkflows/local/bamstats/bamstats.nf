
//
// PREPARE BAMSTATS
//

include { BEDTOOLS_GENOMECOV } from '../../../modules/nf-core/bedtools/genomecov/main'
include { SAMTOOLS_DEPTH     } from '../../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_FLAGSTAT  } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_STATS     } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

workflow BAMSTATS {

    take:
    ch_samplesheet
    ch_bam
    ch_bai
    genome_size
    bed
    fasta
    
    main:

    versions = Channel.empty()

    // Module BEDTOOLS_GENOMECOV 
    BEDTOOLS_GENOMECOV ( 
        ch_bam,
        // ch_bai,
        "1",
        genome_size, 
        extension=".txt",
        sort=true
    )
    versions = versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    // Module SAMTOOLS_DEPTH
    SAMTOOLS_DEPTH (
        ch_bam,
        // ch_bai,
        bed
    )
    versions = versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    // Module SAMTOOLS_FLAGSTAT
    SAMTOOLS_FLAGSTAT (
        ch_bam,
        ch_bai
    )
    versions = versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    // Module SAMTOOLS_STATS
    SAMTOOLS_STATS (
        ch_bam,
        ch_bai,
        fasta
    )
    versions = versions.mix(SAMTOOLS_STATS.out.versions.first())

    emit:
    average    = SAMTOOLS_DEPTH.out.average
    depth      = SAMTOOLS_DEPTH.out.tsv
    flagstat   = SAMTOOLS_FLAGSTAT.out.flagstat
    stats      = SAMTOOLS_STATS.out.stats
    genomecov  = BEDTOOLS_GENOMECOV.out.genomecov

    versions          // channel: [ versions.yml ]
}
