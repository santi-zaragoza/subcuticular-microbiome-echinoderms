# ========================================================
# Script: 02_otu_assignment_coverage.R
# Evaluates OTU assignment coverage in PICRUSt2
# --------------------------------------------------------
# This script is part of an ongoing project.
# The datasets required to run it are NOT publicly available.
# Example input files (not provided):
#   - OTU_table.tsv
#   - combined_marker_predicted_and_nsti.tsv
# ========================================================

# === Required libraries ===
# No external libraries are required for this basic calculation

# === Input files ===
# Place your files in the 'data' folder and modify the names if needed
otu_table <- "data/OTU_table.tsv"
nsti_table <- "data/combined_marker_predicted_and_nsti.tsv"

# === Read OTU table ===
otu_table <- read.table(otu_table, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
total_otus <- nrow(otu_table)

# === Read NSTI table (OTUs with functional predictions) ===
nsti_table <- read.table(nsti_table, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
predicted_otus <- nrow(nsti_table)

# === Calculate excluded OTUs ===
excluded_otus <- total_otus - predicted_otus
excluded_pct <- (excluded_otus / total_otus) * 100

# === Print summary ===
cat("=== PICRUSt2 OTU Assignment Coverage ===\n")
cat(sprintf("Total OTUs in original table:        %d\n", total_otus))
cat(sprintf("OTUs with functional prediction:     %d\n", predicted_otus))
cat(sprintf("OTUs excluded (no prediction):       %d (%.2f%%)\n", excluded_otus, excluded_pct))