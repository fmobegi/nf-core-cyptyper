process PANNO {
    tag "$meta.id"
    label 'process_medium'

    conda "/data/conda-envs/PAnno"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_panno:50ff0f96618804eb' :
        'biocontainers/panno' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.html"), emit: html, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def population_arg_match = args.find(/--population\s+['"]?(\w+)['"]?/)
    def population_code = population_arg_match ?
        population_arg_match.replaceAll(/--population\s+['"]?(\w+)['"]?/, '$1') : 'Unknown'

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
    def population_name = population_map.get(population_code, 'Unknown')

    """
    panno \\
        $args \\
        -i <(zcat $vcf) \\
        -o ./$prefix

    cp ./$prefix/${prefix}.PAnno.html ./${prefix}_population_${population_name}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panno: \$(panno --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def population_arg_match = args.find(/--population\s+(\S+)/)
    def population_code = population_arg_match ? population_arg_match.replaceAll(/--population\s+/, '') : 'Unknown'

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
    def population_name = population_map.get(population_code, 'Unknown')

    """
    touch ${prefix}_population_${population_name}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panno: stub
    END_VERSIONS
    """
}
