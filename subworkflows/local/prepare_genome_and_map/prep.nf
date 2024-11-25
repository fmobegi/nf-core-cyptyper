//
// PREPARE GENOME
//

include { BWA_INDEX as BWAMEM1_INDEX     } from '../../../modules/nf-core/bwa/index/main'
include { GATK4_COMPOSESTRTABLEFILE      } from '../../../modules/nf-core/gatk4/composestrtablefile/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                 } from '../../../modules/nf-core/samtools/index/main'
include { TABIX_TABIX as TABIX_DBSNP     } from '../../../modules/nf-core/tabix/tabix/main'

workflow PREPARE_GENOME {
    take:
    fasta                // channel: [mandatory] fasta
    dbsnp                // channel: [optional]  dbsnp
    model_url            // chanel: not optional model_url
    pypgx_version

    main:
    versions = Channel.empty()

    BWAMEM1_INDEX ( fasta )     // Index for bwa-mem
    versions = versions.mix(BWAMEM1_INDEX.out.versions)

    GATK4_CREATESEQUENCEDICTIONARY ( fasta )
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    SAMTOOLS_FAIDX ( fasta, [ [ id:'no_fai' ], [] ] )
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    GATK4_COMPOSESTRTABLEFILE (
        fasta,
        SAMTOOLS_FAIDX.out.fai,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
    )
    versions = versions.mix(GATK4_COMPOSESTRTABLEFILE.out.versions)

    DOWNLOAD_CLAIR3_MODEL ( model_url )

    DOWNLOAD_PYPGX_BUNDLE ( pypgx_version )

    // the following are flattened and mapped in case the user supplies more than one value for the param
    // written for KNOWN_INDELS, but preemptively applied to the rest
    // [ file1, file2 ] becomes [ [ meta1, file1 ], [ meta2, file2 ] ]
    // outputs are collected to maintain a single channel for relevant TBI files
    TABIX_DBSNP (
        dbsnp.flatten().map{ it -> [ [ id:it.baseName ], it ] }
    )
    versions = versions.mix(TABIX_DBSNP.out.versions)

    emit:
    bwa                      = BWAMEM1_INDEX.out.index.collect()                                        // path: bwa/*
    dbsnp_tbi                = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()                  // path: dbsnb.vcf.gz.tbi
    dict                     = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()                        // path: genome.fasta.dict
    str_table                = GATK4_COMPOSESTRTABLEFILE.out.str_table.collect()                        // path: str_table.zip
    fasta_fai                = SAMTOOLS_FAIDX.out.fai.collect()                                         // path: genome.fasta.fai
    clair3_model_dir         = DOWNLOAD_CLAIR3_MODEL.out.model_dir                                             // path: clair3 variant calling model
    pypgx_bundle             = DOWNLOAD_PYPGX_BUNDLE.out.pypgx_bundle
    versions                                                                                            // channel: [ versions.yml ]
}

process DOWNLOAD_CLAIR3_MODEL {
    tag { url.tokenize('/').last() }
    input:
    val url

    output:
    path("clair3_model_dir"), emit: model_dir

    script:
    """
    wget -O model.tar.gz ${url}
    mkdir -p clair3_model_dir
    tar -xzvf model.tar.gz -C clair3_model_dir
    """
}

process DOWNLOAD_PYPGX_BUNDLE {
    tag { "PYPGX: ${pypgx_version}" }
    input:
    val(pypgx_version)

    output:
    path("pypgx-bundle"), emit: pypgx_bundle

    script:
    def version = "${pypgx_version}".tokenize('/').last()
    """
    git clone --branch $version --depth 1 https://github.com/sbslee/pypgx-bundle
    """
}
