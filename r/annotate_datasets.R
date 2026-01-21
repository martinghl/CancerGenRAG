# Load necessary libraries
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", force = TRUE)
library(DBI)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(annotatr)
library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
#' @param input_directory The path to the full dataset CSV file.
#' @param output_directoryThe prefix for the output file name.
#' @return The annotated dataset.
#'
# Define the function to process all files
process_data_folder <- function(input_directory, output_directory) {
  files <- list.files(input_directory, pattern = "\\.csv$", full.names = TRUE)
  
  for (file_path in files) {
    file_name <- basename(file_path)
    output_file <- file.path(output_directory, file_name)
    
    # Read full dataset
    col_types <- c(
      "character",   # rsID
      "character",     # chr_name should actually be numeric if it's a chromosome number
      "numeric",     # chr_position
      "character",   # effect_allele
      "character",   # other_allele
      "numeric",     # effect_weight
      "numeric",     # allelefrequency_effect
      "character",   # locus_name
      "numeric",     # OR
      "character",   # hm_source
      "character",   # hm_rsID
      "character",   # hm_chr
      "numeric",     # hm_pos
      "character",   # hm_inferOtherAllele
      "character",   # ID
      "character",   # variant_description
      "character",     # inclusion_criteria
      "logical",     # hm_match_chr
      "logical",     # hm_match_pos
      "logical"      # is_diplotype
    )
    
    # full_dataset <- read.csv(file_path,colClasses = col_types)
    full_dataset <- read_csv(file_path)
    str(full_dataset)
    # Reshape AD_data to wide format (missing weights will be NA)
    full_dataset_wide <- tidyr::spread(full_dataset, key = ID, value = effect_weight, fill = NA)
    
    # Extract chromosome coordinates
    chrom_coord <- full_dataset_wide %>%
      dplyr::select(hm_chr, hm_pos, hm_rsID) %>%
      dplyr::rename(chrom = hm_chr, chromStart = hm_pos) %>%
      dplyr::mutate(chromEnd = chromStart + 1) %>%
      dplyr::select(chrom, chromStart, chromEnd, hm_rsID) %>%
      dplyr::mutate(chrom = paste("chr", chrom, sep = ""))
    
    # Remove rows where "chrom_coord" column is "NA"
    chrom_coord <- subset(chrom_coord, chromEnd != "NA")
    
    # Save DataFrame as a tab-delimited file without header names
    temp_bed_path <- tempfile(fileext = ".bed")
    write.table(chrom_coord, file = temp_bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")
    
    # Read bed file created above for annotations 
    dm_regions <- read_regions(
      con = temp_bed_path,
      genome = "hg38",
      format = "bed",
      extraCols = c(rsID = "character")
    )
    
    annotations <- build_annotations(genome = "hg38", annotations = c(
      "hg38_genes_1to5kb", "hg38_genes_promoters", "hg38_genes_cds", "hg38_genes_5UTRs", "hg38_genes_exons", 
      "hg38_genes_firstexons", "hg38_genes_introns", "hg38_genes_intronexonboundaries", "hg38_genes_exonintronboundaries", 
      "hg38_genes_3UTRs", "hg38_genes_intergenic", "hg38_enhancers_fantom", "hg38_basicgenes", "hg38_cpgs" 
    ))
    
    # Intersect regions with annotations
    dm_annotated <- annotate_regions(
      regions = dm_regions,
      annotations = annotations,
      minoverlap = 1L,
      ignore.strand = TRUE,
      quiet = TRUE
    )
    
    # Convert GRanges object to data frame
    df_dm_annotated <- data.frame(dm_annotated)
    
    # Save the annotated dataset as a CSV file
    write.csv(df_dm_annotated, file = output_file, row.names = FALSE)
  }
}

# Example usage
input_directory <- '/Users/martinli/Desktop/CancerGeneBot/Notebook/CancerPGS/Out/'
output_directory <- '/Users/martinli/Desktop/CancerGeneBot/Notebook/CancerPGS/Annotated/'
process_data_folder(input_directory, output_directory)
