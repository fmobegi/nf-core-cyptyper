// MODULE NANOCALLER

process NANOCALLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocaller:3.6.0--h8c52a43_1':
        'biocontainers/nanocaller:3.6.0--h8c52a43_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta1), path(fasta_fai)
    path(bed)


    output:
    tuple val(meta), path("*_nanocaller.vcf.gz"),                 emit: vcf, optional: true
    tuple val(meta), path("*_nanocaller.vcf.gz.tbi"),             emit: tbi, optional: true
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: " "
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_dir = "${prefix}"+"_temp"

    """
    NanoCaller \\
        --bam ${bam} \\
        --ref ${fasta} \\
        --cpu $task.cpus \\
        --bed $bed \\
        --phase \\
        --enable_whatshap \\
        --sample ${prefix} \\
        --output ${output_dir}

    cp ./${output_dir}/variant_calls.vcf.gz ./${prefix}_nanocaller.vcf.gz
    cp ./${output_dir}/variant_calls.vcf.gz.tbi ./${prefix}_nanocaller.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocaller: "v3.4.1"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocaller: "v3.6.0"
    END_VERSIONS
    """
}
