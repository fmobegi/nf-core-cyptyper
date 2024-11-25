//
// MAP TO GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { BWA_MEM        } from '../../../modules/nf-core/bwa/mem/main'
include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'


workflow MAPPING_READS {
   take:
   ch_samplesheet          // channel for input reads for both workflows
   bwa                     // channel for BWA aligner
   fasta                   // channel for input FASTA for both workflows
   sort_bam                // channel for sorting BAM files after BWA alignment
   bam_format              // format for minimap2 BAM
   bam_index_extension     // BAM index extension for minimap2/samtools
   cigar_paf_format        // CIGAR format for minimap2 PAF
   cigar_bam               // CIGAR BAM for minimap2


   main:
   versions = Channel.empty()

   // Execute BWA mem mapping module
   BWA_MEM(
      ch_samplesheet,
      bwa,
      fasta,
      sort_bam="sort"
   )
   
   // Collect versions from BWA module
   versions = versions.mix(BWA_MEM.out.versions.first())

   // Execute Minimap2 mapping module
   MINIMAP2_ALIGN(
      ch_samplesheet,
      fasta,
      bam_format="bam",
      bam_index_extension="bai",
      cigar_paf_format=false,
      cigar_bam=false
   )

   // Collect versions from Minimap2 module
   versions = versions.mix(MINIMAP2_ALIGN.out.versions.first())

   emit:
   bwa_bam        = BWA_MEM.out.bam                  // bwa-mem bam file
   bwa_bai        = BWA_MEM.out.bai                  // bwa-mem bam index file
   minimap2_bam   = MINIMAP2_ALIGN.out.bam           // minimap2 bam file
   minimap2_bai   = MINIMAP2_ALIGN.out.index         // minimap2 bam index file
   versions
}
