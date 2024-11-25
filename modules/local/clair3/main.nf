process CLAIR3 {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.0.10--py39h46983ab_0':
        'biocontainers/clair3:1.0.10--py39h46983ab_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)
    path(bed)
    val(model_path)
    val(model)
    val(sequencer)


    output:
    tuple val(meta), path("*.vcf.gz"),          emit: vcf,            optional:true
    tuple val(meta), path("*.tbi"),             emit: tbi,            optional:true
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    run_clair3.sh \\
        --bam_fn=$input \\
        --ref_fn=$fasta \\
        --platform=$sequencer \\
        --bed_fn=$bed \\
        --threads=$task.cpus \\
        --model_path='$model_path/$model/' \\
        --output=./ \\
        --sample_name=$prefix \\
        --enable_phasing \\
        --include_all_ctgs \\
        --remove_intermediate_dir
    
    rm -f pileup.* merge_* full_align*
    mv phased_merge_output.vcf.gz "${prefix}_clair3.vcf.gz"
    mv phased_merge_output.vcf.gz.tbi "${prefix}_clair3.vcf.gz.tbi"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version)
    END_VERSIONS
    """
}
