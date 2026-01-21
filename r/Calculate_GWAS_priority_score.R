# Load necessary libraries
library(gwasrapidd)
library(readr)
library(dplyr)

# Set the directory where you want to save the file
output_directory <- "/Users/martinli/Desktop/CancerGeneBot/data"  # Update this path to your specific directory

# Retrieve GWAS catalog data for Lung Cancer
LC_VariantAssociations <- get_associations(efo_id = "MONDO_0004992")

# Extract associations and risk alleles from the retrieved data
LC_associations <- data.frame(
  association_id = LC_VariantAssociations@associations$association_id,
  pvalue = LC_VariantAssociations@associations$pvalue,
  or_per_copy_number = LC_VariantAssociations@associations$or_per_copy_number
)

LC_riskAlleles <- data.frame(
  association_id = LC_VariantAssociations@risk_alleles$association_id,
  variant_id = LC_VariantAssociations@risk_alleles$variant_id
)

# Merge association data with risk allele data to link variant IDs
LC_GWAS_Summary <- left_join(LC_associations, LC_riskAlleles, by = "association_id") %>%
  group_by(variant_id) %>%
  summarise(
    GWAS_significant_count = sum(pvalue <= 5e-08, na.rm = TRUE),
    mean_log_or = mean(abs(log(or_per_copy_number)), na.rm = TRUE),
    priority_score = mean_log_or * GWAS_significant_count
  ) %>%
  filter(GWAS_significant_count > 0)  # Keep only variants with significant associations

# Save the results
write.csv(LC_GWAS_Summary, file = paste0(output_directory, "/PC_GWAS_Priority_Scores.csv"), row.names = FALSE)
