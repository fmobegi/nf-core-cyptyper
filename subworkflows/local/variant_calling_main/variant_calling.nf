//
// VARIANT CALLING ALL
//

include {     BCFTOOLS_ANNOTATE     } from '../../../modules/nf-core/bcftools/annotate/main'
include {     BCFTOOLS_CONCAT       } from '../../../modules/nf-core/bcftools/concat/main'
include {     BCFTOOLS_FILTER       } from '../../../modules/nf-core/bcftools/filter/main'
include {     BCFTOOLS_NORM         } from '../../../modules/nf-core/bcftools/norm/main'
include {     BCFTOOLS_SORT         } from '../../../modules/nf-core/bcftools/sort/main'
include {     CLAIR3                } from '../../../modules/local/clair3/main'
include {     CUTESV                } from '../../../modules/nf-core/cutesv/main'
include {     GATK4                 } from '../variant_calling_gatk/gatk4'
include {     MEDAKA_VARIANT        } from '../../../modules/local/medaka/main'
include {     NANOCALLER            } from '../../../modules/local/nanocaller/main'
include {     PYPGX_CREATEINPUTVCF  } from '../../../modules/nf-core/pypgx/createinputvcf/main'
include {     SNIFFLES              } from '../../../modules/nf-core/sniffles/main'
include {     VCFTOOLS              } from '../../../modules/nf-core/vcftools/main'

workflow VARIANT_CALLING_ALL {

    take:
        ch_samplesheet
        bwa_bam
        bwa_bai
        minimap2_bam
        minimap2_bai
        fasta
        fasta_fai
        genome_dict
        dbsnp
        dbsnp_tbi
        bed
        str_table
        clair3_model_dir
        vcf_header

    main:
    ch_versions = Channel.empty()
    combined_bam_channel  = minimap2_bam.mix( minimap2_bai )
    // combined_bam_channel.flatten().view() // View chanel

    // Join BWA_BAM and BWA_BAI for PYPGX
    bwa_bam.join( bwa_bai , failOnDuplicate:true, failOnMismatch:true)
        .map { meta, bam, bai -> [meta, bam, bai] }
        .set { ch_bwa_bam }

    /*
    Samples at this point need to be processed in an organised manner. Each tuple should contain Sample_ID, BAM, BAM.BAI 
    Files are therefore combined into a single chanel before they are reordered for variant-calling
    */
    combined_bam_channel = combined_bam_channel
        .groupTuple(by: 0) 
        .flatMap { meta, path_list ->
            // Split the path_list into BAM and BAI files
            def bam_files = path_list.findAll { it.name.endsWith('.bam') }
            def bai_files = path_list.findAll { it.name.endsWith('.bai') }
            // Create combinations of BAM and BAI files and return as tuples
            def bam_bai_combinations = bam_files.collectMany { bam ->
                bai_files.collect { bai ->
                    return [meta, bam, bai] 
                }
            }

            return bam_bai_combinations
        }
        // .view()
    
    /*
    MODULE PYPGX || Does region filtering based on the gene
    */
    gene_ch = Channel.value('CYP2D6')
    assemby_ch = Channel.value('GRCh38')

    PYPGX_CREATEINPUTVCF (
        ch_bwa_bam,
        fasta,
        gene_ch,
        assemby_ch
    )
    ch_versions = ch_versions.mix( PYPGX_CREATEINPUTVCF.out.versions.first() )

    CLEAN_VCF_HEADER (
        PYPGX_CREATEINPUTVCF.out.vcf
    )
    ch_versions = ch_versions.mix( CLEAN_VCF_HEADER.out.versions.first() )

    /*
    MODULE MEDAKA VARIANT || Does region filtering chr22:42077656-42253758
    */
    // MEDAKA_VARIANT (
    //     combined_bam_channel,
    //     fasta,
    //     fasta_fai
    // )
    // ch_versions = ch_versions.mix( MEDAKA_VARIANT.out.versions.first() )

    /*
    MODULE CLAIR3 VARIANT || Does region filtering chr22:42077656-42253758 
    */
    CLAIR3 (
        combined_bam_channel,
        fasta,
        fasta_fai,
        bed,
        model_path=clair3_model_dir,
        model="r1041_e82_400bps_sup_v420/", // basename of the params.model_url without the extension  '.tar.gz'
        sequencer="ont"
    )
    ch_versions = ch_versions.mix( CLAIR3.out.versions.first() )

    /*
    MODULE NANOCALLER VARIANT || Does region filtering chr22:42077656-42253758
    */
    NANOCALLER (
        combined_bam_channel,
        fasta,
        fasta_fai,
        bed
    )
    ch_versions = ch_versions.mix( NANOCALLER.out.versions.first() )
    
    /*
    MODULE GATK4 VARIANT || Does region filtering chr22:42077656-42253758 
    */
    GATK4 (
        bwa_bam,
        bwa_bai,
        fasta,
        fasta_fai,
        genome_dict,
        dbsnp,
        dbsnp_tbi,
        bed,
        str_table
    )
    ch_versions = ch_versions.mix( GATK4.out.versions.first() )

    // Combine all VCF and TBI channels for concatenation
    CLAIR3.out.vcf
        .join(CLAIR3.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        .set { ch_from_clair3 }

    GATK4.out.vcf
        .join(GATK4.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        .set { ch_from_gatk4 }
    
    // MEDAKA_VARIANT.out.vcf
    //     .join(MEDAKA_VARIANT.out.tbi, failOnDuplicate:true, failOnMismatch:true)
    //     .map { meta, vcf, tbi -> [meta, vcf, tbi] }
    //     .set { ch_from_medaka }

    NANOCALLER.out.vcf
        .join(NANOCALLER.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        .set { ch_from_nanocaller }    
        
    CLEAN_VCF_HEADER.out.vcf
        .join(CLEAN_VCF_HEADER.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        .set { ch_from_pypgx }
        
    ch_from_clair3.mix(
        ch_from_nanocaller, 
        // ch_from_medaka, 
        ch_from_gatk4,
        ch_from_pypgx
    )
    .groupTuple()
    .set { ch_combined_vcfs }
    // .view()

    /*
    MODULE BCFTOOLS CONCAT
    */
    // bcftools concat - adds rows to a combined vcf by collating variant calls of same sample. e.g.
    // variant calls of same individual using different tools.
    // bcftools merge - creates a superset of variant calls across multiple individuals/samples.
    BCFTOOLS_CONCAT ( ch_combined_vcfs, bed )
    ch_versions = ch_versions.mix( BCFTOOLS_CONCAT.out.versions.first() )

    /*
    MODULE BCFTOOLS NORM
    */
    BCFTOOLS_NORM (
        BCFTOOLS_CONCAT.out.vcf,
        BCFTOOLS_CONCAT.out.tbi,
        fasta
    )
    ch_versions = ch_versions.mix( BCFTOOLS_NORM.out.versions.first() )

    /*
    MODULE BCFTOOLS SORT
    */
    BCFTOOLS_SORT (BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix( BCFTOOLS_SORT.out.versions.first() )

    /*
    MODULE BCFTOOLS ANNOTATE || ADD DBSNP DATA
    */
    BCFTOOLS_SORT.out.vcf
        .join( BCFTOOLS_SORT.out.tbi, failOnDuplicate:true, failOnMismatch:true )
        .map { meta, vcf, tbi -> [meta, vcf, tbi] }
        .set { ch_bcftools_sort }
        
    // Create VCF headers chanel
    vcf_header_ch = CREATE_VCF_HEADER ( vcf_header )

    BCFTOOLS_ANNOTATE (
        ch_bcftools_sort,
        dbsnp,
        dbsnp_tbi,
        vcf_header_ch
    )
    ch_versions = ch_versions.mix( BCFTOOLS_ANNOTATE.out.versions.first() )

    /*
    MODULE BCFTOOLS FILTER || ADD DBSNP DATA
    */
    BCFTOOLS_FILTER (BCFTOOLS_ANNOTATE.out.vcf)
    ch_versions = ch_versions.mix( BCFTOOLS_FILTER.out.versions.first() )

    emit:
        vcf                = BCFTOOLS_FILTER.out.vcf
        tbi                = BCFTOOLS_FILTER.out.tbi
        versions           = ch_versions
        ch_from_pypgx
        gene_ch
        assemby_ch
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~ END OF SUBWORKFLOW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

/*
Processes required in this section
*/
process CREATE_VCF_HEADER {
    tag "Creating Header"
    input:
    val vcf_header

    output:
    path 'vcf_header.txt', emit: vcf_header_ch

    script:
    """
    echo "${vcf_header}" | sed 's/^[ \t]*//' > vcf_header.txt
    """
}

process QC_AND_CLEAN_FINAL_VCF {
    conda 'xxx'
    tag "$name"

    input:
        tuple val(meta), path(input), path(index)

    output:
        tuple val(name), path("*.gz"),  emit: merged
        tuple val(name), path("*.tbi"), emit: mergedindex

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    
    """
    python filter_vcf.py -i ${input}
    """
}

process CLEAN_VCF_HEADER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::htslib=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib:1.17--h81da01d_2' :
        'quay.io/biocontainers/htslib:1.17--h81da01d_2' }"

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

    # Find the #CHROM line, remove .bam extension, and process the rest of the file
    awk -F'\\t' -v OFS='\\t' '
    {
        if (\$0 ~ /^#CHROM/) {
            for (i=1; i<=NF; i++) {
                if (i > 9) {  # Only process sample columns
                    sub(/\\.bam\$/, "", \$i)
                }
            }
        }
        print \$0
    }
    ' input.vcf > ${prefix}.clean.vcf

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
