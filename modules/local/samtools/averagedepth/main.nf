process SAMTOOLS_AVERAGEDEPTH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/gawk:5.3.1--78c56a1b4e4b6534':
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    tuple val(meta), path(depth_file)

    output:
    tuple val(meta), path("*.txt"), emit: avg_depth
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '{sum += \$3; count++} END {if (count > 0) print sum/count; else print 0}' ${depth_file} > ${prefix}.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/^.*gawk //; s/,.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/^.*gawk //; s/,.*\$//')
    END_VERSIONS
    """
}
