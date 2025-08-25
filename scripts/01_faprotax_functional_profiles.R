# ========================================================
# Script: 01_faprotax_functional_profiles.R
# Visualization of Functional Profiles (FAPROTAX)
# --------------------------------------------------------
# This script is part of an ongoing project.
# The datasets required to run it are NOT publicly available.
# Example input files (not provided):
#   - func_table_def.tsv
#   - metadata.xlsx
# ========================================================

# === Required libraries ===
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(scales)

# Input files: place your data in the 'data' folder
func_file <- "data/func_table_def.tsv"
metadata_file <- "data/metadata.xlsx"

# Read the data
func_data <- read.table(func_file, header = TRUE, row.names = 1, sep = "\t")
metadata <- read_excel(metadata_file)

# Filter functions with at least one abundance > 0
func_data <- func_data[rowSums(func_data) > 0, ]
func_data$group <- rownames(func_data)

# Filter samples of interest (e.g., water, sediment, subcuticle, or biofilm)
sample_type <- "sediment"  # Change here to the type of sample you want to analyze

code_samples <- subset(metadata, type == sample_type)
Funct_samples <- func_data[, colnames(func_data) %in% code_samples$code]
Funct_samples$group <- rownames(Funct_samples)

# ============================
# Transform data to long format and add metadata information
# ============================

# Metadata column to use for grouping (e.g., "site", "classes", "type", etc.)
metadata_column <- "site"  # Change according to what you want to analyze

# Convert to long format and add grouping information
func_long_samples <- Funct_samples %>%
  pivot_longer(-group, names_to = "Sample", values_to = "Abundance") %>%
  filter(Abundance > 0) %>%
  left_join(code_samples %>% select(code, all_of(metadata_column)), by = c("Sample" = "code")) %>%
  rename(Sample = !!metadata_column)  # 'Sample' indicates the selected metadata column (e.g., site, type, classes)


# Group by metadata and function, and calculate mean abundance
func_mean_group <- func_long_samples %>%
  group_by(Sample, group) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

# ============================
# Color palette and bubbleplot
# ============================

pal <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100)

# Values for the color scale
min_val <- min(func_mean_group$Abundance)
max_val <- max(func_mean_group$Abundance)
moderately_low_mid <- quantile(func_mean_group$Abundance, 0.9)

# Bubbleplot with averaged groups
ggplot(func_mean_group, aes(x = Sample, y = group, size = Abundance, fill = Abundance)) +
  geom_point(shape = 21, color = "black", alpha = 0.8) +
  scale_size_continuous(name = "Mean Abundance", range = c(1, 10)) +
  scale_fill_gradientn(
    name = "Mean Abundance",
    colours = pal,
    values = rescale(c(min_val, moderately_low_mid, max_val)),
    guide = "colorbar"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
    axis.text.y = element_text(size = 7)
  ) +
  labs(title = NULL, x = NULL, y = NULL)