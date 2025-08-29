process CLEAN_VCF_HEADER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/tabix_gawk_gzip:3331e140ac0a612d' :
        'community.wave.seqera.io/library/tabix_gawk_gzip:aa8115a85b67867b' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.clean.vcf.gz"), emit: vcf
    tuple val(meta), path("*.clean.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -e

    # Decompress the input VCF if it's gzipped
    if [[ $vcf == *.gz ]]; then
        gunzip -c $vcf > input.vcf
    else
        cp $vcf input.vcf
    fi

    # Find the #CHROM line and clean sample names to end at _barcode + 2 digits
    awk -F'\\t' -v OFS='\\t' -v meta_id="${meta.id}" '
    {
        if (\$0 ~ /^#CHROM/) {
            for (i=1; i<=NF; i++) {
                if (i > 9) {
                    if (match(\$i, /_barcode[0-9][0-9]/)) {
                        \$i = substr(\$i, 1, RSTART + RLENGTH - 1)
                    }
                    if (\$i != meta_id) {
                        \$i = meta_id
                    }
                }
            }
        }
        print \$0
    }' input.vcf > ${prefix}.clean.vcf

    # Compress with bgzip
    bgzip ${prefix}.clean.vcf

    # Create tabix index
    tabix -p vcf ${prefix}.clean.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$(echo \$(tabix --version 2>&1) | sed 's/^.*tabix //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clean.vcf.gz
    touch ${prefix}.clean.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$(echo \$(tabix --version 2>&1) | sed 's/^.*tabix //; s/ .*\$//')
    END_VERSIONS
    """
}
