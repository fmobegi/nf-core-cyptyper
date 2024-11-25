process CUTESV {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutesv:1.0.12--pyhdfd78af_0' :
        'biocontainers/cutesv:2.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz"),   emit: vcf,            optional:true
    tuple val(meta), path("*.tbi"),      emit: tbi,            optional:true
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cuteSV \\
        $args \\
        ${bam} \\
        ${fasta} \\
        ${prefix}_cutesv.vcf \\
        . \\
        --threads $task.cpus
        
    bgzip --threads ${task.cpus} ${prefix}_cutesv.vcf
    tabix ${prefix}_cutesv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuteSV: \$( cuteSV --version 2>&1 | sed 's/cuteSV //g' )
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
