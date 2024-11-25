process PYPGX_RUNLONGREADPIPELINE {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val(pgx_genes)
    val(assembly_version)



    output:
    tuple val(meta), path("*.tsv"), emit: star_alleles
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def assembly = "${assembly_version}" ?: "GRCh38"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pgx_genes = "${pgx_genes}" ?: "CYP2D6"

    """
    pypgx run-long-read-pipeline \\
        ${args} \\
        --assembly ${assembly} \\
        ${pgx_genes} \\
        ${prefix}_alleles \\
        $vcf

    pypgx print-data ${prefix}_alleles/results.zip > ${prefix}_TGS.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_TGS.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
