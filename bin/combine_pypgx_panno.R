#!/usr/bin/env Rscript

suppressPackageStartupMessages({
   library(optparse)
   library(readr)
   library(readxl)
   library(openxlsx)
   library(dplyr)
   library(tidyr)
   library(tidyselect)
})

# Define command-line options
option_list <- list(
   make_option(c("--pypgx"), type = "character", help = "Path to PyPGx results TSV"),
   make_option(c("--panno"), type = "character", help = "Path to PAnno annotations XLSX"),
   make_option(c("--output"), type = "character", help = "Output prefix for result files")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Read input files
data_pypgx <- read_tsv(opt$pypgx, show_col_types = FALSE)
data_panno <- read_excel(opt$panno)

# Define known population columns
known_populations <- c(
   "Oceania", "Oceanian", "AfricanAmerican", "American", "EasternAsia", "European",
   "Latin", "NearEastern", "SouthASIA", "SubSaharanAfrica"
)

# Detect which population columns are present
present_pops <- intersect(known_populations, colnames(data_panno))
if (length(present_pops) == 0) {
   stop("No known population columns found in PAnno annotations.")
}

# Reshape population columns into long format
data_panno_long <- data_panno %>%
   select(SampleID, all_of(present_pops)) %>%
   pivot_longer(cols = all_of(present_pops), names_to = "Population", values_to = "PAnno_Combined") %>%
   separate(PAnno_Combined, c("PAnno_Genotype", "Phenotype_panno"), sep = "\\|") %>%
   mutate(
      PAnno_Genotype = trimws(PAnno_Genotype),
      Phenotype_panno = paste(trimws(Phenotype_panno), "Metabolizer")
   ) %>%
   filter(!startsWith(SampleID, "Table citation"))

# Function to combine haplotypes
combine_haplotypes <- function(hap1, hap2) {
   hap1_sorted <- sort(unique(strsplit(hap1, ";")[[1]]))
   hap2_sorted <- sort(unique(strsplit(hap2, ";")[[1]]))
   paste(paste(hap1_sorted, collapse = ";"), paste(hap2_sorted, collapse = ";"), sep = " | ")
}

# Process PyPGx data
simplified_data <- data_pypgx %>%
   group_by(SampleID) %>%
   summarize(
      NGS_geno = first(Genotype[Pipeline == "NGS"]),
      TGS_geno = first(Genotype[Pipeline == "TGS"]),
      Phenotype = first(Phenotype),
      NGS_haplotypes = combine_haplotypes(first(Haplotype1[Pipeline == "NGS"]), first(Haplotype2[Pipeline == "NGS"])),
      TGS_haplotypes = combine_haplotypes(first(Haplotype1[Pipeline == "TGS"]), first(Haplotype2[Pipeline == "TGS"])),
      AlternativePhase = paste(unique(unlist(strsplit(AlternativePhase, ";"))), collapse = ";"),
      VariantData = paste(unique(VariantData), collapse = " | ")
   ) %>%
   mutate(
      AlternativePhase = if_else(AlternativePhase == "", NA_character_, AlternativePhase),
      across(everything(), ~ na_if(., ""))
   ) %>%
   arrange(SampleID)

# Merge with long-format PAnno annotations
combined <- simplified_data %>%
   full_join(data_panno_long, by = "SampleID") %>%
   select(
      SampleID,
      Population,
      PAnno_Genotype,
      NGS_geno,
      TGS_geno,
      Phenotype,
      Phenotype_panno,
      everything()
   ) %>%
   arrange(SampleID)

# Write outputs
write_tsv(combined, paste0(opt$output, ".tsv"))
write.xlsx(combined, paste0(opt$output, ".xlsx"), overwrite = TRUE)
