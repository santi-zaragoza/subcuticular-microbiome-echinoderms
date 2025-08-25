# ========================================================
# Script: 06_gini_functional_contributions.R
# Evaluation of functional contribution inequality (Gini Index)
# --------------------------------------------------------
# This script is part of an ongoing project.
# The datasets required to run it are NOT publicly available.
# Example input files (not provided):
#   - path_abun_contrib.tsv
#   - metadata.tsv
# ========================================================

# === Required libraries ===
library(readr)
library(dplyr)
library(tidyr)
library(ineq)
library(writexl)
library(ggplot2)

# === Input files ===
contrib_file <- "data/path_abun_contrib.tsv"
metadata_file <- "data/metadata.tsv"

# === Load data ===
Contribution_data <- read_tsv(contrib_file, trim_ws = TRUE)
metadata          <- read_tsv(metadata_file, trim_ws = TRUE)


# === Functions of interest ===
selected_functions <- c(
  "1CMET2-PWY", "ARGSYN-PWY", ..., "URSIN-PWY", "VALSYN-PWY")

# === Filter metadata ===
sample_type <- "subcuticle"  # Replace with the sample type of interest
metadata_filtered <- metadata %>%
  filter(type == sample_type)
  

# === Filter contribution data (matching samples + functions) ===
Contribution_type <- Contribution_data %>%
  filter(sample %in% metadata_filtered$code) %>%
  filter(`function` %in% selected_functions) %>%
  left_join(metadata_filtered, by = c("sample" = "code"))


# === Compute Gini index per function × sample ===
# 'species' here refers to bacterial species contributing to each function
gini_fun_sample_species <- Contribution_type %>%
  group_by(type, `function`, sample) %>%
  summarise(
    gini = ineq(norm_taxon_function_contrib, na.rm = TRUE),
    .groups = "drop"
  )

# === Summarise Gini by function (mean, sd, sem, n) ===
gini_fun_species <- gini_fun_sample_species %>%
  group_by(type, `function`) %>%
  summarise(
    mean_gini = mean(gini, na.rm = TRUE),
    sd_gini = sd(gini, na.rm = TRUE),
    sem_gini = sd_gini / sqrt(n()),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(type, desc(mean_gini))

# === Wide format: functions as rows, types as columns ===
gini_fun_species_wide <- gini_fun_species %>%
  select(type, `function`, mean_gini) %>%
  pivot_wider(names_from = type, values_from = mean_gini)

# === Visualisation ===
ggplot(gini_fun_species, aes(x = mean_gini, y = `function`, color = type)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = paste("Mean Gini index per function in", sample_type, "samples"),
    x = "Mean Gini index",
    y = "Function"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Sample type"))