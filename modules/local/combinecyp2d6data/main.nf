process COMBINE_CYP2D6_DATA {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-readxl_r-tidyverse:71c385e0b1be291e' :
        'community.wave.seqera.io/library/r-base_r-readxl_r-tidyverse:bbb57ca864228160' }"

    input:
    path pypgx_results
    path panno_annotations
    val output_prefix

    output:
    path("${output_prefix}.tsv")
    path("${output_prefix}.xlsx")
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript --vanilla ${projectDir}/bin/combine_pypgx_panno.R \\
        --pypgx ${pypgx_results} \\
        --panno ${panno_annotations} \\
        --output ${output_prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep 'R version' | sed 's/R version //g')
        tidyverse: \$(Rscript -e 'cat(as.character(packageVersion("tidyverse")))')
        readxl: \$(Rscript -e 'cat(as.character(packageVersion("readxl")))')
        openxlsx: \$(Rscript -e 'cat(as.character(packageVersion("openxlsx")))')
    END_VERSIONS
    """

    stub:
    """
    touch ${output_prefix}.tsv
    touch ${output_prefix}.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep 'R version' | sed 's/R version //g')
        tidyverse: \$(Rscript -e 'cat(as.character(packageVersion("tidyverse")))')
        readxl: \$(Rscript -e 'cat(as.character(packageVersion("readxl")))')
        openxlsx: \$(Rscript -e 'cat(as.character(packageVersion("openxlsx")))')
    END_VERSIONS
    """
}
