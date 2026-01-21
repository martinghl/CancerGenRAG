#' Process data
#'
#' This function annotates a merged dataset with genomic information.
#'
#' @param full_dataset_path The path to the full dataset CSV file.
#' @param output_file_prefix The prefix for the output file name.
#' @return The annotated dataset.
#'
#' @import tidyr
#' @importFrom dplyr select mutate group_by arrange distinct left_join
#' @importFrom GenomicRanges read_regions annotate_regions
#' @importFrom annotatr build_annotations
#' @export

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", force=TRUE)
BiocManager::install("annotatr", force=TRUE)
BiocManager::install("org.Hs.eg.db", force=TRUE)

unlink("C:/Users/Gang/AppData/Local/R/win-library/4.5/00LOCK", recursive = TRUE)

install.packages("BiocParallel")
BiocManager::install("Rsamtools")

# ****        process_data function to annotate merged dataset with genomic information
process_data <- function(full_dataset_path, output_file_prefix) {
  options(max.print = 100)
  
  library(DBI)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(annotatr)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  
  # -------- Read --------
  full_dataset <- read.csv(full_dataset_path, stringsAsFactors = FALSE)
  
  # -------- Build GRanges directly (no BED) --------
  # normalize chromosome labels to UCSC style (1..22,X,Y,M)
  normalize_chr_vec <- function(x) {
    x <- as.character(x)
    x <- sub("^chr", "", x, ignore.case = TRUE)
    x <- ifelse(x == "23", "X", x)
    x <- ifelse(x == "24", "Y", x)
    x <- ifelse(x %in% c("MT", "m", "Mt", "mt"), "M", x)
    x
  }
  
  coords <- full_dataset %>%
    dplyr::select(hm_chr, hm_pos, hm_rsID) %>%
    dplyr::mutate(
      chr_std = normalize_chr_vec(hm_chr),
      pos     = as.integer(hm_pos)
    ) %>%
    dplyr::filter(!is.na(chr_std), !is.na(pos)) %>%
    dplyr::distinct(chr_std, pos, hm_rsID, .keep_all = TRUE)
  
  dm_regions <- GenomicRanges::makeGRangesFromDataFrame(
    coords %>% dplyr::transmute(
      chrom = paste0("chr", chr_std),
      start = pos,
      end   = pos,
      rsID  = hm_rsID
    ),
    seqnames.field = "chrom",
    start.field    = "start",
    end.field      = "end",
    keep.extra.columns = TRUE
  )
  GenomeInfoDb::seqlevelsStyle(dm_regions) <- "UCSC"
  GenomeInfoDb::genome(dm_regions) <- "hg38"
  
  # -------- Build annotations --------
  annots <- c(
    "hg38_genes_1to5kb", "hg38_genes_promoters", "hg38_genes_cds",
    "hg38_genes_5UTRs", "hg38_genes_exons", "hg38_genes_firstexons",
    "hg38_genes_introns", "hg38_genes_intronexonboundaries",
    "hg38_genes_exonintronboundaries", "hg38_genes_3UTRs",
    "hg38_genes_intergenic", "hg38_enhancers_fantom",
    "hg38_basicgenes", "hg38_cpgs"
  )
  annotations <- suppressMessages(build_annotations(genome = "hg38", annotations = annots))
  
  # -------- Annotate --------
  dm_annotated <- annotate_regions(
    regions = dm_regions,
    annotations = annotations,
    minoverlap = 1L,
    ignore.strand = TRUE,
    quiet = TRUE
  )
  
  df_dm_annotated <- as.data.frame(dm_annotated)
  print(names(df_dm_annotated))
  print(full_dataset_path)
  
  # prefer annot.symbol if present; otherwise fall back to symbol
  if (!"annot.symbol" %in% names(df_dm_annotated) && "symbol" %in% names(df_dm_annotated)) {
    df_dm_annotated$annot.symbol <- df_dm_annotated$symbol
  }
  
  df_dm_annotated_select <- df_dm_annotated %>%
    dplyr::select(seqnames, start, annot.symbol, annot.type, annot.width) %>%
    dplyr::mutate(
      chr_std = sub("^chr", "", as.character(seqnames)),
      SNP_coord = paste(chr_std, start, sep = "_")
    )
  
  # 以最窄注释为代表（去重）
  df_dm_annotated_filter <- df_dm_annotated_select %>%
    dplyr::group_by(SNP_coord) %>%
    dplyr::filter(!is.na(annot.symbol)) %>%
    dplyr::arrange(annot.width, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # -------- Prepare join key on the original dataset (use the same normalization) --------
  full_dataset <- full_dataset %>%
    dplyr::mutate(
      chr_std = normalize_chr_vec(hm_chr),
      pos     = as.integer(hm_pos),
      SNP_coord = paste(chr_std, pos, sep = "_")
    )
  
  # -------- Merge & post-process --------
  annotated_df <- dplyr::left_join(full_dataset, df_dm_annotated_filter, by = "SNP_coord") %>%
    dplyr::select(
      "SNP_coord", "hm_rsID", "hm_chr", "hm_pos",
      "effect_allele", "effect_weight", "ID",
      "annot.symbol", "annot.type"
    ) %>%
    dplyr::filter(!is.na(SNP_coord) & SNP_coord != "NA_NA")
  
  # 每个 study (ID) 内按 effect_weight 选去重 SNP
  annotated_df <- annotated_df %>%
    dplyr::group_by(ID) %>%
    dplyr::arrange(dplyr::desc(effect_weight), .by_group = TRUE) %>%
    dplyr::distinct(SNP_coord, .keep_all = TRUE) %>%
    dplyr::ungroup()
  
  # 添加 rank（绝对值）
  annotated_df <- annotated_df %>%
    dplyr::arrange(ID, effect_weight) %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(ranks = rank(-abs(as.double(effect_weight)), ties.method = "average")) %>%
    dplyr::ungroup()
  
  # -------- Save --------
  output_file <- paste0(output_file_prefix, "_annotated_dataset.csv")
  write.csv(annotated_df, file = output_file, row.names = FALSE)
  
  return(annotated_df)
}

# ---------- Batch runner ----------
input_directory <- '/Users/martinli/Desktop/CancerGeneBot/Notebook/CancerPGS/Out/'
output_directory <- '/Users/martinli/Desktop/CancerGeneBot/Notebook/CancerPGS/Annotated/'
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

file_list <- list.files(path = input_directory, pattern = "\\.csv$", full.names = TRUE)

for (file_path in file_list) {
  file_name <- basename(file_path)
  output_file_prefix <- file.path(output_directory, tools::file_path_sans_ext(file_name))
  message("Processing: ", file_name)
  invisible(process_data(file_path, output_file_prefix))
}
