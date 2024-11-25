process PANNO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_panno:50ff0f96618804eb':
        'biocontainers/panno' }"

    errorStrategy 'ignore'

    input:
    tuple val(meta), path(vcf), path(tbi)
    val(population)

    output:
    tuple val(meta), path("*.html"), emit: html, optional:true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def population_map = [
        'OCE': 'Oceania',
        'AAC': 'AfricanAmerican',
        'AME': 'American',
        'EAS': 'EasternAsia',
        'EUR': 'European',
        'LAT': 'Latin',
        'NEA': 'NearEastern',
        'SAS': 'SouthASIA',
        'SSA': 'SubSaharanAfrica'
    ]
    def population_name = population_map.get(population, 'Unknown')

    """
    panno -s ${prefix}_population_${population_name} \\
        -i <(zcat $vcf) \\
        -p $population \\
        -o .
        
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panno: \$(panno --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def population_map = [
        'OCE': 'Oceania',
        'AAC': 'AfricanAmerican',
        'AME': 'American',
        'EAS': 'EasternAsia',
        'EUR': 'European',
        'LAT': 'Latin',
        'NEA': 'NearEastern',
        'SAS': 'SouthASIA',
        'SSA': 'SubSaharanAfrica'
    ]
    def population_name = population_map.get(population, 'Unknown')

    """
    touch ${prefix}_population_${population_name}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panno: \$(panno --version)
    END_VERSIONS
    """
}