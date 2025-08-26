process COMBINE_CYP2D6_DATA {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-readxl_r-tidyverse:71c385e0b1be291e':
        'community.wave.seqera.io/library/r-base_r-readxl_r-tidyverse:bbb57ca864228160' }"

    input:
    path pypgx_results
    path panno_annotations

    output:
    path("combined_cyp2d6_genotyping.tsv")
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    library(tidyverse)
    library(readxl)

    # Read the data
    data_pypgx <- read_tsv("${pypgx_results}")
    data_panno <- read_excel("${panno_annotations}")

    data_panno_clean <- data_panno %>%
        separate(Oceania, c("PAnno_Genotype_OCE", "Phenotype_panno"), sep = "\\\\|") %>%
        mutate(
            PAnno_Genotype_OCE = trimws(PAnno_Genotype_OCE, "both"),
            Phenotype_panno = paste(trimws(Phenotype_panno, "both"), "Metabolizer")
        )

    # Function to combine haplotypes
    combine_haplotypes <- function(hap1, hap2) {
        hap1_sorted <- sort(unique(strsplit(hap1, ";")[[1]]))
        hap2_sorted <- sort(unique(strsplit(hap2, ";")[[1]]))

        paste(
            paste(hap1_sorted, collapse = ";"),
            paste(hap2_sorted, collapse = ";"),
            sep = " | "
        )
    }

    # Process the data_pypgx
    simplified_data <- data_pypgx %>%
        group_by(SampleID) %>%
        summarize(
            NGS_geno = first(Genotype[Pipeline == "NGS"]),
            TGS_geno = first(Genotype[Pipeline == "TGS"]),
            Phenotype = first(Phenotype),
            NGS_haplotypes = combine_haplotypes(first(Haplotype1[Pipeline == "NGS"]), first(Haplotype2[Pipeline == "NGS"])),
            TGS_haplotypes = combine_haplotypes(first(Haplotype1[Pipeline == "TGS"]), first(Haplotype2[Pipeline == "TGS"])),
            AlternativePhase = paste(unique(unlist(
                strsplit(AlternativePhase, ";")
            )), collapse = ";"),
            VariantData = paste(unique(VariantData), collapse = " | ")
        ) %>%
        mutate(
            AlternativePhase = ifelse(AlternativePhase == "", NA, AlternativePhase),
            across(everything(), ~ na_if(., ""))
        ) %>%
        arrange(SampleID)

    combined <- simplified_data %>%
        full_join(data_panno_clean) %>%
        dplyr::select(SampleID, PAnno_Genotype_OCE, NGS_geno, TGS_geno, Phenotype, Phenotype_panno, everything())

    write_tsv(combined, "combined_cyp2d6_genotyping.tsv")

    # Output versions
    writeLines(
        c(
            "\\"${task.process}\\":",
            "    r-base: \$(R --version | grep 'R version' | sed 's/R version //')",
            "    tidyverse: \$(Rscript -e 'cat(as.character(packageVersion(\\"tidyverse\\")))')",
            "    readxl: \$(Rscript -e 'cat(as.character(packageVersion(\\"readxl\\")))')"
        ),
        "versions.yml"
    )
    """

    stub:
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch combined_cyp2d6_genotyping.tsv
    echo "\\"${task.process}\\":" > versions.yml
    echo "    r-base: 4.2.2" >> versions.yml
    echo "    tidyverse: 1.3.2" >> versions.yml
    echo "    readxl: 1.4.1" >> versions.yml
    """
}
