# ========================================================
# Script: 08_phylogenetic_signal_gini.R
# Phylogenetic Signal Analysis (Pagel’s λ): Gini Indices per Function and Host Species
# --------------------------------------------------------
# This script is part of an ongoing project.
# The datasets required to run it are NOT publicly available.
# Example input files (not provided):
#   - gini_function_per_species_filtered_wide.xlsx
#   - tree.txt
# ========================================================

# === Required libraries ===
library(dplyr)
library(readxl)
library(ape)
library(geiger)
library(writexl)

# === Input files ===
gini_file <- "data/gini_function_per_species_filtered_wide.xlsx"  
tree_file <- "data/tree.txt"                                      

# === Load data ===
gini_data <- read_excel(gini_file)
tree      <- read.tree(tree_file)

# Rename first column to avoid conflicts
colnames(gini_data)[1] <- "pathway"
# Ensure species names are valid R names
colnames(gini_data)[-1] <- make.names(colnames(gini_data)[-1])

# === Functions to analyze ===
functions_to_analyze <- c(
  "1CMET2-PWY", "ARGSYN-PWY", ..., "URSIN-PWY", "VALSYN-PWY")  # Replace ... with full list

# === List to store results ===
results_lambda <- list()

# === Iterate over functions ===
for (func in functions_to_analyze) {
  cat("Processing:", func, "\n")
  
  # Extract row corresponding to the function
  row_func <- gini_data %>% filter(pathway == func)
  
  if (nrow(row_func) == 0) {
    warning(paste("Function not found:", func))
    next
  }
  
  # Create Gini vector per host species
  gini_vector <- as.numeric(row_func[ , -1])
  names(gini_vector) <- colnames(row_func)[-1]  # species names
  
  # Keep only species present in the tree
  species_present <- intersect(tree$tip.label, names(gini_vector))
  
  if (length(species_present) < 4) {
    warning(paste("Too few species with data for", func, "- skipping"))
    next
  }
  
  pr_tree <- keep.tip(tree, species_present)
  gini_vector_filtered <- gini_vector[pr_tree$tip.label]  # reorder vector to match tree
  
  # Validate names
  if (!all(names(gini_vector_filtered) == pr_tree$tip.label)) {
    warning(paste("Inconsistency between vector and tree for:", func))
    next
  }
  
  # Fit Pagel's lambda model
  tryCatch({
    fit_lambda <- fitContinuous(pr_tree, gini_vector_filtered, model = "lambda")
    results_lambda[[func]] <- fit_lambda$opt$lambda
  }, error = function(e) {
    warning(paste("Error fitting model for", func, ":", e$message))
    results_lambda[[func]] <- NA
  })
}

# === Create results data.frame ===
results_df <- data.frame(
  pathway = names(results_lambda),
  lambda  = unlist(results_lambda)
)

# Order by lambda descending
results_df_ordered <- results_df %>% arrange(desc(lambda))

