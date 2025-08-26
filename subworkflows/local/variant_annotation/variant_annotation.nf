//
// VARIANT ANNOTATION MAIN
//

include { COMBINE_CYP2D6_DATA             } from '../../../modules/local/combinecyp2d6data/main'
include { PANNO                           } from '../../../modules/local/panno/main'
include { PYPGX_RUNLONGREADPIPELINE            } from '../../../modules/local/pypgx/runlongreadpipeline/main'
include { PYPGX_RUNNGSPIPELINE            } from '../../../modules/local/pypgx/runngspipeline/main'

workflow VARIANT_ANNOTATION {

    take:
    vcf
    tbi
    ch_from_pypgx
    gene_ch
    assemby_ch
    pypgx_bundle


    main:
    
    ch_versions = Channel.empty()

    vcf.join( tbi, failOnDuplicate:true, failOnMismatch:true )
        .set { ch_panno_inputs }
    
    //
    // MODULE: Run PANNO
    // 
    
    PANNO ( ch_panno_inputs, "OCE" )

    ch_versions = ch_versions.mix(PANNO.out.versions.first())

    //
    // MODULE: Run PYPGX_RUNLONGREADPIPELINE
    // 
        PYPGX_RUNLONGREADPIPELINE (
        ch_from_pypgx,
        gene_ch,
        assemby_ch
    )
    ch_versions = ch_versions.mix(PYPGX_RUNLONGREADPIPELINE.out.versions.first())
    
    //
    // MODULE: Run PYPGX_RUNNGSPIPELINE
    // 
    PYPGX_RUNNGSPIPELINE (
        ch_from_pypgx,
        gene_ch,
        assemby_ch,
        pypgx_bundle
    )
    ch_versions = ch_versions.mix(PYPGX_RUNNGSPIPELINE.out.versions.first())

    // Get python scripts
    extract_panno = file("${projectDir}/bin/extract_cyp2d6_annotations.py").resolve()
    combine_pypgx = file("${projectDir}/bin/combine_pypgx.py").resolve()

    //
    // MODULE: Run COMBINE_PANNO_RESULTS
    // 
    COMBINE_PANNO_RESULTS ( 
        PANNO.out.html.map { meta, file -> file }.collect(), 
        extract_panno 
    )

    ch_combined_inputs = PYPGX_RUNLONGREADPIPELINE.out.star_alleles
        .mix(PYPGX_RUNNGSPIPELINE.out.star_alleles) 
        .map { meta, tsv -> tsv } 
        .collect()

    //
    // MODULE: Run COMBINE_PYPGX_RESULTS
    // 
    COMBINE_PYPGX_RESULTS (
        ch_combined_inputs, 
        combine_pypgx
    )

    COMBINE_CYP2D6_DATA (
        COMBINE_PYPGX_RESULTS.out.pypgx_tsv,
        COMBINE_PANNO_RESULTS.out.panno_xls
    )

    ch_versions = ch_versions.mix(COMBINE_CYP2D6_DATA.out.versions.first())


    emit:
    versions                  = ch_versions
}


//
// End of module. See additional processes below
//


// processes required to summarize Panno CYP2D6 results.
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
    path "CYP2D6_annotations.xlsx", emit: panno_xls

    script:
    """
    python3 ${extract_panno} -i ./
    """
}

// processes required to summarize PYPGX CYP2D6 results.
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
    path "combined_pypgx_results.tsv", emit: pypgx_tsv

    script:
    """
    python3 ${combine_pypgx} -i ./
    """
}


