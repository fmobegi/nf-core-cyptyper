process COMBINE_CYP2D6_DATA {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-readxl_r-tidyverse:71c385e0b1be291e' :
        'community.wave.seqera.io/library/r-base_r-readxl_r-tidyverse:bbb57ca864228160' }"

    input:
    path(pypgx_results)
    path(panno_annotations)
    val(output_prefix)

    output:
    path("${output_prefix}.tsv")
    path("${output_prefix}.xlsx")
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(readxl)
    library(openxlsx)
    library(dplyr)
    library(tidyr)
    library(tidyselect)

    # Read input files
    data_pypgx <- readr::read_tsv("${pypgx_results}")
    data_panno <- readxl::read_excel("${panno_annotations}")

    # Define known population columns
    known_populations <- c("Oceania", "Oceanian", "AfricanAmerican", "American", "EasternAsia", "European",
                           "Latin", "NearEastern", "SouthASIA", "SubSaharanAfrica")

    # Detect which population columns are present
    present_pops <- base::intersect(known_populations, base::colnames(data_panno))
    if (length(present_pops) == 0) {
        stop("No known population columns found in PAnno annotations.")
    }

    # Reshape population columns into long format
    data_panno_long <- data_panno %>%
        dplyr::select(SampleID, tidyselect::all_of(present_pops)) %>%
        tidyr::pivot_longer(cols = tidyselect::all_of(present_pops), names_to = "Population", values_to = "PAnno_Combined") %>%
        tidyr::separate(PAnno_Combined, c("PAnno_Genotype", "Phenotype_panno"), sep = "\\\\|") %>%
        dplyr::mutate(
            PAnno_Genotype = base::trimws(PAnno_Genotype, "both"),
            Phenotype_panno = paste(base::trimws(Phenotype_panno, "both"), "Metabolizer")
        ) %>%
        dplyr::filter(!base::startsWith(SampleID, "Table citation"))

    # Function to combine haplotypes
    combine_haplotypes <- function(hap1, hap2) {
        hap1_sorted <- base::sort(base::unique(base::strsplit(hap1, ";")[[1]]))
        hap2_sorted <- base::sort(base::unique(base::strsplit(hap2, ";")[[1]]))
        paste(paste(hap1_sorted, collapse = ";"), paste(hap2_sorted, collapse = ";"), sep = " | ")
    }

    # Process PyPGx data
    simplified_data <- data_pypgx %>%
        dplyr::group_by(SampleID) %>%
        dplyr::summarize(
            NGS_geno = dplyr::first(Genotype[Pipeline == "NGS"]),
            TGS_geno = dplyr::first(Genotype[Pipeline == "TGS"]),
            Phenotype = dplyr::first(Phenotype),
            NGS_haplotypes = combine_haplotypes(dplyr::first(Haplotype1[Pipeline == "NGS"]), dplyr::first(Haplotype2[Pipeline == "NGS"])),
            TGS_haplotypes = combine_haplotypes(dplyr::first(Haplotype1[Pipeline == "TGS"]), dplyr::first(Haplotype2[Pipeline == "TGS"])),
            AlternativePhase = paste(base::unique(base::unlist(base::strsplit(AlternativePhase, ";"))), collapse = ";"),
            VariantData = paste(base::unique(VariantData), collapse = " | ")
        ) %>%
        dplyr::mutate(
            AlternativePhase = dplyr::if_else(AlternativePhase == "", NA_character_, AlternativePhase),
            dplyr::across(dplyr::everything(), ~ dplyr::na_if(., ""))
        ) %>%
        dplyr::arrange(SampleID)

    # Merge with long-format PAnno annotations
    combined <- simplified_data %>%
        dplyr::full_join(data_panno_long, by = "SampleID") %>%
        dplyr::select(
            SampleID,
            Population,
            PAnno_Genotype,
            NGS_geno,
            TGS_geno,
            Phenotype,
            Phenotype_panno,
            dplyr::everything()
        ) %>%
        dplyr::arrange(SampleID)

    # Write outputs
    readr::write_tsv(combined, "${output_prefix}.tsv")
    openxlsx::write.xlsx(combined, "${output_prefix}.xlsx", overwrite = TRUE)

    # Write version info
    writeLines(
        c(
            "\\"${task.process}\\":",
            "    r-base: \$(R --version | grep 'R version' | sed 's/R version //')",
            "    tidyverse: \$(Rscript -e 'cat(as.character(packageVersion(\\"tidyverse\\")))')",
            "    readxl: \$(Rscript -e 'cat(as.character(packageVersion(\\"readxl\\")))')",
            "    openxlsx: \$(Rscript -e 'cat(as.character(packageVersion(\\"openxlsx\\")))')"
        ),
        "versions.yml"
    )
    """

    stub:
    """
    #!/usr/bin/env bash

    touch ${output_prefix}.tsv
    touch ${output_prefix}.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep 'R version' | sed 's/R version //g')
        tidyverse: \$(Rscript -e 'cat(as.character(packageVersion("tidyverse")))')
        readxl: \$(Rscript -e 'cat(as.character(packageVersion("readxl")))')
        openxlsx: \$(Rscript -e 'cat(as.character(packageVersion("openxlsx")))')
    END_VERSIONS
    """
}
