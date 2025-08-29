/*
SUBWORKFLOW VARIANT_ANNOTATION
    This subworkflow is responsible for annotating variants using various tools and methods. It includes
    the following steps:
        1. Running PANNO to annotate variants in the provided VCF files.
        2. Running PYPGX_RUNLONGREADPIPELINE to process long-read sequencing data.
        3. Running PYPGX_RUNNGSPIPELINE to process next-generation sequencing data.
        4. Combining results from PANNO and PYPGX into a comprehensive CYP2D6 dataset.
        It emits the final CYP2D6 annotations and versions of the tools used.
*/

include { COMBINE_CYP2D6_DATA             } from '../../../modules/local/combinecyp2d6data/main'
include { PANNO                           } from '../../../modules/local/panno/main'
include { PYPGX_RUNLONGREADPIPELINE       } from '../../../modules/local/pypgx/runlongreadpipeline/main'
include { PYPGX_RUNNGSPIPELINE            } from '../../../modules/local/pypgx/runngspipeline/main'

workflow VARIANT_ANNOTATION {

    take:
    ch_combined_vcf_tbi
    ch_pypgx_vcf_tbi
    ch_pypgx_bundle


    main:

    ch_versions = Channel.empty()
    gene_ch = Channel.value('CYP2D6')
    assemby_ch = Channel.value('GRCh38')

    /*
    MODULE: PANNO
    */

    PANNO ( ch_combined_vcf_tbi )

    ch_versions = ch_versions.mix(PANNO.out.versions.first())

    /*
    MODULE: PYPGX_RUNLONGREADPIPELINE
    */
    PYPGX_RUNLONGREADPIPELINE (
        ch_pypgx_vcf_tbi,
        gene_ch,
        assemby_ch
    )
    ch_versions = ch_versions.mix(PYPGX_RUNLONGREADPIPELINE.out.versions.first())

    /*
    MODULE: PYPGX_RUNNGSPIPELINE
    */
    PYPGX_RUNNGSPIPELINE (
        ch_pypgx_vcf_tbi,
        gene_ch,
        assemby_ch,
        ch_pypgx_bundle
    )
    ch_versions = ch_versions.mix(PYPGX_RUNNGSPIPELINE.out.versions.first())

    // Get python scripts
    extract_panno = file("${projectDir}/bin/extract_cyp2d6_annotations.py").resolve()
    combine_pypgx = file("${projectDir}/bin/combine_pypgx.py").resolve()

    /*
    PROCESS: COMBINE_PANNO_RESULTS
    */
    COMBINE_PANNO_RESULTS (
        PANNO.out.html.map { meta, file -> file }.collect(),
        extract_panno
    )

    ch_combined_inputs = PYPGX_RUNLONGREADPIPELINE.out.star_alleles
        .mix(PYPGX_RUNNGSPIPELINE.out.star_alleles)
        .map { meta, tsv -> tsv }
        .collect()

    /*
    PROCESS: COMBINE_PYPGX_RESULTS
    */
    COMBINE_PYPGX_RESULTS (
        ch_combined_inputs,
        combine_pypgx
    )

    COMBINE_CYP2D6_DATA (
        COMBINE_PYPGX_RESULTS.out.pypgx_tsv,
        COMBINE_PANNO_RESULTS.out.panno_xls,
        output_prefix = 'FINAL_CYP2D6_RESULTS'
    )

    ch_versions = ch_versions.mix(COMBINE_CYP2D6_DATA.out.versions)


    emit:
    versions                  = ch_versions
}


/*
PROCESS: COMBINE_PANNO_RESULTS
    processes required to summarize Panno CYP2D6 results.   
*/

process COMBINE_PANNO_RESULTS {
    conda "${projectDir}/envs/py_combine.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95ad3cea7cbaeab6e0e8bc1b3556dc8772adaafe73b802e2ded8883a21658c53/data' :
        'community.wave.seqera.io/library/openpyxl_pandas:a7e2c36678fbb02a' }"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path html_files
    path extract_panno

    output:
    path "CYP2D6_PANNO_RESULTS.xlsx", emit: panno_xls
    path "CYP2D6_PANNO_RESULTS.log",  emit: log

    script:
    """
    python3 ${extract_panno} -i ./ -o CYP2D6_PANNO_RESULTS.xlsx | tee CYP2D6_PANNO_RESULTS.log
    """
}

/*
PROCESS: COMBINE_PYPGX_RESULTS
    processes required to summarize PYPGX CYP2D6 results.
*/

process COMBINE_PYPGX_RESULTS {
    conda "${projectDir}/envs/py_combine.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95ad3cea7cbaeab6e0e8bc1b3556dc8772adaafe73b802e2ded8883a21658c53/data' :
        'community.wave.seqera.io/library/openpyxl_pandas:a7e2c36678fbb02a' }"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tsv_files
    path combine_pypgx

    output:
    path "CYP2D6_PYPGX_RESULTS.xlsx", emit: pypgx_xls
    path "CYP2D6_PYPGX_RESULTS.tsv" , emit: pypgx_tsv
    path "CYP2D6_PYPGX_RESULTS.log" , emit: log

    script:
    """
    python3 ${combine_pypgx} -i ./ -o CYP2D6_PYPGX_RESULTS | tee CYP2D6_PYPGX_RESULTS.log
    """
}
