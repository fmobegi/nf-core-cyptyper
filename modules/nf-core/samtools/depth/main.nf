process SAMTOOLS_DEPTH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)
    path(intervals)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.average"), emit: average
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def positions = intervals ? "-b ${intervals}" : ""
    """
    samtools \\
        depth \\
        --threads ${task.cpus-1} \\
        $args \\
        $positions \\
        -o ${prefix}.tsv \\
        $bam

    awk 'NR>1 {sum[\$1] += \$3; count[\$1] += \$3=="" ? 0 : 1} ; \\
        END {for (i in sum) print i, (count[i] > 0 ? sum[i]/count[i] : "-")}' \\
        ${prefix}.tsv | \\
        cut -f 1,2 | sort -nrk2 >  ${prefix}.samtools.depth.average;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
