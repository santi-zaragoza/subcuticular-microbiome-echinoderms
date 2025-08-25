# ========================================================
# Script: 07_phylogenetic_signal_abundance.R
# Phylogenetic Signal Analysis (Pagel’s λ): Average Functional Abundances
# --------------------------------------------------------
# This script is part of an ongoing project.
# The datasets required to run it are NOT publicly available.
# Example input files (not provided):
#   - path_abundance_data.tsv
#   - metadata.tsv
#   - tree.txt
# ========================================================

# === Required libraries ===
library(dplyr)
library(readr)
library(ape)
library(geiger)
library(writexl)

# === Input files ===
abundance_file <- "data/path_abundance_data.tsv"
metadata_file  <- "data/metadata.tsv"
tree_file      <- "data/tree.txt"

# === Parameters ===
sample_type <- "subcuticle"  # Replace with any sample type of interest
functions_to_analyze <- c(
  "1CMET2-PWY", "ARGSYN-PWY", ..., "URSIN-PWY", "VALSYN-PWY"  # Full list here
)

# === Load data ===
abundance_data <- read_tsv(abundance_file, trim_ws = TRUE)
metadata       <- read_tsv(metadata_file, trim_ws = TRUE)
tree           <- read.tree(tree_file)

# === Filter metadata by sample type ===
metadata_filtered <- metadata %>% filter(type == sample_type)

# Check that all filtered samples exist in abundance data
missing_codes <- setdiff(metadata_filtered$code, colnames(abundance_data))
if (length(missing_codes) > 0) {
  stop("The following samples are missing in abundance data: ", paste(missing_codes, collapse = ", "))
}

# Filter abundance data for selected samples + pathway column
abundance_filtered <- abundance_data %>% select(c("pathway", metadata_filtered$code))

# === List to store results ===
results_lambda <- list()

# === Iterate over functions ===
for (func in functions_to_analyze) {
  cat("Processing:", func, "\n")
  
  row_func <- abundance_filtered %>% filter(pathway == func)
  if (nrow(row_func) == 0) {
    warning(paste("Function not found:", func))
    next
  }
  
  # Create vector of abundance per sample
  abundance_vector <- as.numeric(row_func[, -1])
  names(abundance_vector) <- colnames(row_func)[-1]
  
  # Join with metadata 
  df_abundance <- data.frame(
    sample = names(abundance_vector),
    abundance = abundance_vector
  ) %>% left_join(metadata_filtered %>% select(code, sp), by = c("sample" = "code"))
  
  # Compute mean abundance per host species
  mean_abundance_per_species <- df_abundance %>%
    group_by(sp) %>%
    summarise(mean_abundance = mean(abundance, na.rm = TRUE)) %>%
    filter(!is.na(sp))
  
  div <- setNames(mean_abundance_per_species$mean_abundance, as.character(mean_abundance_per_species$sp))
  
  # Filter tree for species present
  species_present <- intersect(tree$tip.label, names(div))
  if (length(species_present) < 4) {
    warning(paste("Too few species with data for", func, "- skipping"))
    next
  }
  
  pr_tree <- keep.tip(tree, species_present)
  
  # Ensure order matches
  div_filtered <- div[pr_tree$tip.label]
  
  # Fit Pagel's lambda
  tryCatch({
    fit_lambda <- fitContinuous(pr_tree, div_filtered, model = "lambda")
    results_lambda[[func]] <- fit_lambda$opt$lambda
  }, error = function(e) {
    warning(paste("Error fitting model for", func, ":", e$message))
    results_lambda[[func]] <- NA
  })
}

# === Create results data.frame and export ===
results_df <- data.frame(
  pathway = names(results_lambda),
  lambda = unlist(results_lambda)
)

results_df_ordered <- results_df %>% arrange(desc(lambda))

