# ========================================================
# Script: 04_differential_abundance_picrust2.R
# Differential Abundance Analysis (PICRUSt2 functions)
# --------------------------------------------------------
# This script is part of an ongoing project.
# The datasets required to run it are NOT publicly available.
# Example input files (not provided):
#   - path_abundance_data.tsv
#   - metadata.tsv
# ========================================================

# === Required libraries ===
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(gplots)
library(tools)

# Input files: place your data in the 'data' folder
abundance_file <- "data/path_abundance_data.tsv"
metadata_file  <- "data/metadata.tsv"

# Load data
abundance_data <- read_tsv(abundance_file, trim_ws = TRUE)
metadata       <- read_tsv(metadata_file, trim_ws = TRUE)

# Optional manual filtering of selected functions:
# If 'selected_functions' is defined, only the functions included
# among significant features (p_adjust < 0.05) will be shown.
# If it is not defined, all significant functions will be included.

#selected_functions <- c("1CMET2-PWY", "ARGSYN-PWY", "ARGSYNBSUB-PWY", ...) # example functions

# Parameters
sample_types <- c("sediment", "subcuticle")  # Types of Samples to Include

# Filter metadata by selected sample types
metadata_filtered <- metadata %>% filter(type %in% sample_types)

# Filter abundance data based on selected sample codes.
# Additional filtering can be added to separate by classes or sites.
abundance_filtered <- abundance_data %>%
  select(c("pathway", metadata_filtered$code))

# Differential analysis with ALDEx2
results <- pathway_daa(
  abundance = abundance_filtered %>% column_to_rownames("pathway"),
  metadata = metadata_filtered,
  group = "type",
  daa_method = "ALDEx2"
)

# Iterate over sub-methods
for (submethod in unique(results$method)) {
  result_sub <- results %>% filter(method == submethod)
  feature_with_p_0.05 <- result_sub %>% filter(p_adjust < 0.05)

  # Optional filtering by selected functions
  if (exists("selected_functions")) {
    feature_selected <- feature_with_p_0.05 %>%
      filter(feature %in% selected_functions)
  } else {
    feature_selected <- feature_with_p_0.05
  }

  if (nrow(feature_selected) > 0) {
    abundance_significant <- abundance_filtered %>%
      filter(pathway %in% feature_selected$feature)

    heatmap_matrix <- abundance_significant %>%
      select(c("pathway", metadata_filtered$code)) %>%
      column_to_rownames("pathway") %>%
      as.matrix()

    row_labels <- pathway_annotation(pathway = "MetaCyc", daa_results_df = result_sub, ko_to_kegg = FALSE) %>%
      filter(feature %in% rownames(heatmap_matrix)) %>%
      select(feature, description) %>%
      mutate(description = ifelse(is.na(description) | description == "", feature, description))

    rownames(heatmap_matrix) <- row_labels$description[match(rownames(heatmap_matrix), row_labels$feature)]

    metadata_filtered_ordered <- metadata_filtered %>% arrange(type)
    ordered_samples <- metadata_filtered_ordered$code
    heatmap_matrix <- heatmap_matrix[, ordered_samples]
    heatmap_matrix <- t(scale(t(heatmap_matrix)))

    group_colors <- c("#377eb8", "#00ced1")[as.numeric(as.factor(metadata_filtered_ordered$type))]
    sample_labels <- paste(metadata_filtered_ordered$code, metadata_filtered_ordered$site_loc, sep = "_")

    # Create the heatmap
    heatmap.2(
    heatmap_matrix,
    trace = "none",
    col = colorRampPalette(c("blue", "white", "red"))(100),
    dendrogram = "both",
    scale = "none",
    margins = c(1, 29),
    cexRow = 0.9,
    cexCol = 0.0001,
    key = FALSE,
    key.title = "Abundance",
    key.xlab = "Scaled values",
    density.info = "none",
    main = NULL,
    ColSideColors = group_colors,
    labRow = rownames(heatmap_matrix),
    labCol = NULL,
    keysize = 0.8,
    cex.main = 0.5, lwid = c(0.4, 4)
)
  }
}


