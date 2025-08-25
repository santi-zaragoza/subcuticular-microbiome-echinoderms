# ========================================================
# Script: 03_weighted_nsti_analysis.R
# Nearest Sequenced Taxon Index (NSTI) Evaluation
# --------------------------------------------------------
# This script is part of an ongoing project.
# The datasets required to run it are NOT publicly available.
# Example input files (not provided):
#   - OTU_table.tsv
#   - combined_marker_predicted_and_nsti.tsv
#   - metadata.tsv
# ========================================================

# === Required libraries ===
library(tibble)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggpubr)

# Input files: place your data in the 'data' folder
otu_table <- "data/OTU_table.tsv"
nsti_table <- "data/combined_marker_predicted_and_nsti.tsv"
metadata <- "data/metadata.tsv"

# Load data
otu_table <- read.table(otu_table, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
nsti_table <- read.table(nsti_table, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
metadata <- read.table(metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Add OTU column for joining
otu_table2 <- otu_table %>% rownames_to_column("OTU")
nsti_df <- nsti_table %>% rownames_to_column("OTU") %>% select(OTU, metadata_NSTI)

# Convert OTU table to long format
long_otu <- otu_table2 %>%
  pivot_longer(cols = -OTU, names_to = "sample", values_to = "abundance")

# Join with NSTI values
merged <- long_otu %>%
  left_join(nsti_df, by = "OTU")

# Filter OTUs with available NSTI and abundance > 0
filtered <- merged %>%
  filter(!is.na(metadata_NSTI), abundance > 0)

# Calculate weighted NSTI per sample
weighted_nsti_df <- filtered %>%
  group_by(sample) %>%
  summarise(weighted_NSTI = sum(abundance * metadata_NSTI) / sum(abundance)) %>%
  ungroup()

# Join with metadata
merged_df <- weighted_nsti_df %>%
  left_join(metadata, by = c("sample" = "code"))

# Parameters
sample_type_filter <- "subcuticle"   # Change to any sample type, e.g., "water", "sediment"
class_column <- "classes"               # Use the metadata column you want to plot

# Filter by sample type
sample_df <- merged_df %>% filter(type == sample_type_filter)

# Convert class column to factor automatically
sample_df[[class_column]] <- as.factor(sample_df[[class_column]])

# Boxplot with colors and ordered classes
ggboxplot(sample_df, x = class_column, y = "weighted_NSTI",
          color = "black", fill = "white",
          ylab = "Weighted NSTI", xlab = class_column) +
    stat_compare_means(method = "kruskal.test") +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11,
                                   margin = margin(t = 10)),
        axis.title.x = element_text(margin = margin(t = 10), size = 13),
        axis.title.y = element_text(margin = margin(r = 10), size = 13)
    )
