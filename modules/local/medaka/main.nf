// MODULE MEDAKA VARIANT
process MEDAKA_VARIANT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.4.4--py38h130def0_0':
        'biocontainers/medaka:1.4.4--py38h130def0_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta1), path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf.gz"),      emit: vcf, optional: true
    tuple val(meta), path("*.vcf.gz.tbi"),  emit: tbi, optional: true
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_dir = "${prefix}"
    def output_vcf = "${output_dir}"+"/round_1.vcf"

    """
    medaka_variant \\
        -f ${fasta} \\
        -n ${prefix} \\
        -i ${bam} \\
        -o ${output_dir} \\
        -t $task.cpus \\
        -l \\
        -p \\
        -r "chr22:42077656-42253758" \\
        -m r104_e81_sup_g5015 \\
        -s r104_e81_sup_g5015 
    
    cp ${prefix}/round_1_phased.vcf ${prefix}_medaka_phased.vcf
    bgzip ${prefix}_medaka_phased.vcf
    tabix ${prefix}_medaka_phased.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_medaka.vcf.gz
    touch ${prefix}_medaka.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
