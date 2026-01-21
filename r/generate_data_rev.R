## ================================
## PGS hmPOS -> clean CSV (testable)
## ================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(quincunx)
})

## ---- Config ----
trait_term       <- "lung cancer"
directory_path   <- "/Users/martinli/Desktop/CancerGeneBot/Notebook/CancerPGS/PGS_Decompress"
output_directory <- "/Users/martinli/Desktop/CancerGeneBot/Notebook/CancerPGS/Out"
suffix           <- "_hmPOS_GRCh38.txt"

## Required cols you care about (for debug only)
must_have_cols <- c("rsID","chr_name","chr_position","effect_allele","effect_weight")

## Exact standard schema (your list & order)
std_cols <- c(
  "rsID","chr_name","chr_position","effect_allele","other_allele","effect_weight",
  "allelefrequency_effect","locus_name","OR","hm_source","hm_rsID","hm_chr","hm_pos",
  "hm_inferOtherAllele","ID","variant_description","inclusion_criteria",
  "hm_match_chr","hm_match_pos","is_diplotype"
)

enforce_schema <- TRUE  # keep TRUE to force the output header exactly
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

## ---- Helpers ----

# Strip .metadata.* and keep/rename .data.* columns, or extract $data if it's a list
clean_pgs_table <- function(obj) {
  # Case A: list(metadata=..., data=...)
  if (is.list(obj) && all(c("metadata","data") %in% names(obj))) {
    obj <- obj$data
  }
  # Case B: data.frame with path-prefixed names like ...metadata... / ...data...
  if (is.data.frame(obj) && any(grepl("\\.(metadata|data)\\.", names(obj)))) {
    keep <- grepl("\\.data\\.", names(obj))
    obj  <- obj[, keep, drop = FALSE]
    names(obj) <- sub(".*\\.data\\.", "", names(obj), perl = TRUE)
  }
  obj
}

# Fallback: scan file to locate the data header line, re-read as tab-delimited table
fallback_read_hmpos <- function(full_path) {
  lines <- readLines(full_path, warn = FALSE)
  # Find first header line with .data. and tabs
  hdr_idx <- which(grepl("\\.data\\.", lines) & grepl("\t", lines))[1]
  # If not found, try common bare headers (some files may already be clean)
  if (is.na(hdr_idx)) {
    hdr_idx <- which(grepl("(^|\\t)(rsID|chr_name|chr_position)(\\t|$)", lines))[1]
  }
  if (is.na(hdr_idx)) {
    stop("fallback_read_hmpos: could not locate a data header line in: ", full_path)
  }
  df <- utils::read.delim(
    full_path, header = TRUE, sep = "\t",
    skip = hdr_idx - 1, check.names = FALSE, quote = "", comment.char = ""
  )
  # Strip prefixes if present
  if (any(grepl("\\.(metadata|data)\\.", names(df)))) {
    keep <- grepl("\\.data\\.", names(df))
    df   <- df[, keep, drop = FALSE]
    names(df) <- sub(".*\\.data\\.", "", names(df), perl = TRUE)
  }
  df
}

# Enforce exact schema & ensure ID exists
lock_schema <- function(df, std_cols, id_value = NULL) {
  if (!is.data.frame(df)) stop("lock_schema: df is not a data.frame")
  if (!"ID" %in% names(df) && !is.null(id_value)) df$ID <- id_value
  for (cc in setdiff(std_cols, names(df))) df[[cc]] <- NA
  df <- df[, std_cols, drop = FALSE]
  df
}

## ---- Main ----

generate_datasets_per_file <- function(trait_term, directory_path, output_directory) {
  message("Trait term received: ", trait_term)
  
  # Retrieve trait entries and flatten to PGS IDs
  PGS_traits <- quincunx::get_traits(trait_term = trait_term, exact_term = FALSE)
  trait_vector <- unique(as.character(unlist(PGS_traits@pgs_ids, use.names = FALSE)))
  trait_vector <- trait_vector[!is.na(trait_vector) & grepl("^PGS\\d{6}$", trait_vector)]
  
  # Intersect with files actually present on disk
  available_files <- list.files(directory_path, pattern = "_hmPOS_GRCh38\\.txt$", full.names = FALSE)
  available_ids   <- sub("_hmPOS_GRCh38\\.txt$", "", available_files)
  trait_vector    <- intersect(trait_vector, available_ids)
  
  if (!length(trait_vector)) {
    warning("No matching PGS IDs on disk for '", trait_term, "'.")
    return(invisible(NULL))
  }
  
  lapply(trait_vector, function(item) {
    full_path <- file.path(directory_path, paste0(item, suffix))
    message("\n---\n[", item, "] path: ", full_path)
    
    if (!file.exists(full_path)) {
      warning("File not found: ", full_path)
      return(invisible(NULL))
    }
    
    # Try quincunx reader first
    data <- tryCatch(
      read_scoring_file(full_path),
      error = function(e) {
        message("[", item, "] read_scoring_file() error: ", conditionMessage(e))
        return(structure(list(error = TRUE), class = "read_error"))
      }
    )
    if (inherits(data, "read_error")) data <- NULL
    
    # Clean structure & headers if we got something
    if (!is.null(data)) data <- clean_pgs_table(data)
    
    # Fallback if not a proper data.frame with >= 2 columns
    need_fallback <- is.null(data) || !is.data.frame(data) || ncol(data) < 2
    if (need_fallback) {
      message("[", item, "] fallback parse: scanning for '.data.' header")
      data <- fallback_read_hmpos(full_path)
    }
    
    if (!is.data.frame(data)) {
      stop("[", item, "] parsed object is not a data.frame; got class: ", paste(class(data), collapse = ", "))
    }
    
    message("[", item, "] colnames (first 12): ", paste(utils::head(names(data), 12), collapse = ", "))
    
    # ID from filename
    cleaned_name <- sub("_hmPOS_GRCh38.txt$", "", basename(full_path))
    data$ID <- cleaned_name
    
    # Debug: report missing key columns and write a colname dump
    missing <- setdiff(must_have_cols, names(data))
    if (length(missing)) {
      message("[", item, "] MISSING columns: ", paste(missing, collapse = ", "))
      diag_file <- file.path(output_directory, paste0(cleaned_name, "_DEBUG_COLNAMES.txt"))
      writeLines(paste(names(data), collapse = "\n"), diag_file)
      message("[", item, "] wrote colnames to: ", diag_file)
      # Continue; schema lock below will fill NA so downstream still works
    }
    
    # Enforce exact schema if requested
    if (enforce_schema) {
      data <- lock_schema(data, std_cols, id_value = cleaned_name)
    }
    
    # Save
    out_csv <- file.path(output_directory, paste0(cleaned_name, ".csv"))
    write.csv(data, file = out_csv, row.names = FALSE)
    cat("File saved as ", out_csv, "\n", sep = "")
    invisible(out_csv)
  })
  
  invisible(NULL)
}

## ---- Run ----
generate_datasets_per_file(trait_term, directory_path, output_directory)

## ---- Single-file quick test (optional) ----
# test_id <- "PGS000078"
# fp <- file.path(directory_path, paste0(test_id, suffix))
# if (file.exists(fp)) {
#   dd <- try(read_scoring_file(fp), silent = TRUE)
#   if (inherits(dd, "try-error")) dd <- NULL
#   if (!is.null(dd)) dd <- clean_pgs_table(dd)
#   if (is.null(dd) || !is.data.frame(dd) || ncol(dd) < 2) dd <- fallback_read_hmpos(fp)
#   dd$ID <- test_id
#   if (enforce_schema) dd <- lock_schema(dd, std_cols, id_value = test_id)
#   write.csv(dd, file.path(output_directory, paste0(test_id, ".csv")), row.names = FALSE)
#   message("Single-file test saved: ", file.path(output_directory, paste0(test_id, ".csv")))
# } else {
#   message("Single-file test not found: ", fp)
# }
