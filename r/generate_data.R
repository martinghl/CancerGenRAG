#' Generate datasets for each PGS file
#'
#' This function processes each PGS score file separately and saves the output with a cleaned name.
#'
#' @param trait_term The trait term to retrieve PGS scores for.
#' @param directory_path The directory path where the PGS score files are located.
#' @return None. Files are saved directly to the specified directory.
#' @export

generate_datasets_per_file <- function(trait_term, directory_path, output_directory) {
  print(paste("Trait term received:", trait_term))
  library(dplyr)
  library(tidyr)
  library(quincunx)
  
  # Retrieve PGS traits related to the provided trait term
  PGS_traits <- quincunx::get_traits(trait_term = trait_term, exact_term = FALSE)
  # trait_vector <- PGS_traits@pgs_ids
  trait_vector <- unique(as.character(unlist(PGS_traits@pgs_ids, use.names = FALSE)))
  
  # keep ONLY real PGS IDs like "PGS000388"
  trait_vector <- trait_vector[!is.na(trait_vector) & grepl("^PGS\\d{6}$", trait_vector)]
  
  # (optional but robust) intersect with what's actually present on disk
  available_files <- list.files(directory_path, pattern = "_hmPOS_GRCh38\\.txt$", full.names = FALSE)
  available_ids   <- sub("_hmPOS_GRCh38\\.txt$", "", available_files)
  trait_vector    <- intersect(trait_vector, available_ids)
  
  # BEFORE lapply: make sure the output dir exists
  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)
  
  # Define the file suffix
  suffix <- "_hmPOS_GRCh38.txt"
  
  # Process each PGS file separately
  lapply(trait_vector, function(item) {
    full_path <- file.path(directory_path, paste0(item, suffix))
    
    if (file.exists(full_path)) {
      # Read the scoring file
      data <- read_scoring_file(full_path)
      
      # Clean the file name
      cleaned_name <- sub("_hmPOS_GRCh38.txt$", "", basename(full_path))
      
      # Add the ID column
      data$ID <- cleaned_name
      
      # Save the dataset as a CSV file with the cleaned name
      output_file <- file.path(output_directory, paste0(cleaned_name, ".csv"))
      write.csv(data, file = output_file, row.names = FALSE)
      
      cat("File saved as", output_file, "\n")
    } else {
      warning(paste("File not found:", full_path))
    }
  })
}


# Set parameters
directory_path <- 'F:/TulaneOneDrive/OneDrive - Tulane University/Desktop/CancerGeneBot/Notebook/CancerPGS/PGS_Decompress'
output_directory <- 'F:/TulaneOneDrive/OneDrive - Tulane University/Desktop/CancerGeneBot/Notebook/CancerPGS/Out'

# Run the function
generate_datasets_per_file("lung cancer", directory_path, output_directory)
