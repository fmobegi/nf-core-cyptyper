/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PREPARE_REFERENCE_INDEXES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This subworkflow prepares reference genome files and associated indices required
    for variant calling and genomic analysis. It creates BWA indices, GATK sequence
    dictionaries, samtools faidx indices, STR tables, and downloads required models.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_INDEX                      } from '../../../modules/nf-core/bwa/index/main'
include { GATK4_COMPOSESTRTABLEFILE      } from '../../../modules/nf-core/gatk4/composestrtablefile/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx/main'
include { TABIX_TABIX                    } from '../../../modules/nf-core/tabix/tabix/main'

workflow PREPARE_REFERENCE_INDEXES {

    take:
    ch_fasta             // channel: [mandatory] fasta
    ch_dbsnp             // channel: [optional]  dbsnp
    model_url            // channel: not optional model_url
    pypgx_version        // val: pypgx version string

    main:
    ch_versions = Channel.empty()

    /*
    MODULE: BWA_INDEX
    */
    BWA_INDEX ( ch_fasta )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    /*
    MODULE: GATK4_CREATESEQUENCEDICTIONARY
    */
    GATK4_CREATESEQUENCEDICTIONARY ( ch_fasta )
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    /*
    MODULE: SAMTOOLS_FAIDX
    */
    SAMTOOLS_FAIDX (
        ch_fasta,
        [ [ id:'no_fai' ], [] ],
        true
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    /*
    MODULE: GATK4_COMPOSESTRTABLEFILE
    */
    GATK4_COMPOSESTRTABLEFILE (
        ch_fasta.map { meta, file -> file },
        SAMTOOLS_FAIDX.out.fai.map { meta, file -> file },
        GATK4_CREATESEQUENCEDICTIONARY.out.dict.map { meta, file -> file }
    )
    ch_versions = ch_versions.mix(GATK4_COMPOSESTRTABLEFILE.out.versions)

    /*
    MODULE: DOWNLOAD_CLAIR3_MODEL
    */
    DOWNLOAD_CLAIR3_MODEL ( model_url )

    /*
    MODULE: DOWNLOAD_PYPGX_BUNDLE
    */
    DOWNLOAD_PYPGX_BUNDLE ( pypgx_version )

    /*
    MODULE: TABIX_TABIX
    */
    TABIX_TABIX ( ch_dbsnp )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    bwa_index            = BWA_INDEX.out.index.collect()                                    // path: bwa/*
    clair3_model_dir     = DOWNLOAD_CLAIR3_MODEL.out.model_dir
    dbsnp_tbi            = TABIX_TABIX.out.tbi.collect()
    dbsnp_vcf            = ch_dbsnp
    dict                 = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()
    fasta_ref            = ch_fasta
    fasta_fai            = SAMTOOLS_FAIDX.out.fai.collect()
    pypgx_bundle         = DOWNLOAD_PYPGX_BUNDLE.out.pypgx_bundle
    str_table            = GATK4_COMPOSESTRTABLEFILE.out.str_table.collect()
    versions             = ch_versions                                                      
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: DOWNLOAD_CLAIR3_MODEL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downloads and extracts Clair3 variant calling model from specified URL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: DOWNLOAD_PYPGX_BUNDLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downloads PyPGx bundle from GitHub repository for pharmacogenomic analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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