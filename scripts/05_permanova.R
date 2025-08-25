# ========================================================
# Script: 05_permanova.R
# PERMANOVA analysis in echinoderm microbiomes
# --------------------------------------------------------
# This script is part of an ongoing project.
# The datasets required to run it are NOT publicly available.
# Example input files (not provided):
#   - path_abundance_data.tsv
#   - metadata.tsv
# ========================================================

# === Required libraries ===
library(vegan)
library(tidyverse)
library(openxlsx) # For exporting Excel files

# === Input files ===
abundance_file <- "data/path_abundance_data.tsv"
metadata_file  <- "data/metadata.tsv"

# === Load abundance data ===
abundance_data <- read.table(abundance_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
abundance_data <- t(abundance_data)

# === Load metadata ===
metadata <- read.table(metadata_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
metadata <- metadata[match(rownames(abundance_data), rownames(metadata)), ]
stopifnot(all(rownames(abundance_data) == rownames(metadata)))

# === 2. Normalize abundances ===
abundances_relative <- sweep(abundance_data, 1, rowSums(abundance_data), FUN = "/")

# === 3. Bray-Curtis distance across all samples ===
dist_brays <- vegdist(abundances_relative, method = "bray")

# === 4. Global univariate PERMANOVA ===
adonis_type <- adonis2(dist_brays ~ type, data = metadata)

# === 5. Multivariate analysis for biological samples only (subcuticle and biofilm) ===
idx_bio <- metadata$type %in% c("subcuticle", "biofilm")
abundance_bio <- abundances_relative[idx_bio, ]
metadata_bio <- metadata[idx_bio, ]
dist_bio <- vegdist(abundance_bio, method = "bray")

adonis_bio_global <- adonis2(dist_bio ~ type + clases + sp + site_loc + site, data = metadata_bio)
adonis_clases     <- adonis2(dist_bio ~ clases, data = metadata_bio)
adonis_sp         <- adonis2(dist_bio ~ sp, data = metadata_bio)
adonis_site       <- adonis2(dist_bio ~ site, data = metadata_bio)
adonis_site_loc   <- adonis2(dist_bio ~ site_loc, data = metadata_bio)

# === 6. Subcuticle samples only ===
idx_sub <- metadata$type == "subcuticle"
abundance_sub <- abundances_relative[idx_sub, ]
metadata_sub <- metadata[idx_sub, ]
dist_sub <- vegdist(abundance_sub, method = "bray")

adonis_sub_global   <- adonis2(dist_sub ~ clases + sp + site_loc + site, data = metadata_sub)
adonis_sub_clases   <- adonis2(dist_sub ~ clases, data = metadata_sub)
adonis_sub_sp       <- adonis2(dist_sub ~ sp, data = metadata_sub)
adonis_sub_site     <- adonis2(dist_sub ~ site, data = metadata_sub)
adonis_sub_site_loc <- adonis2(dist_sub ~ site_loc, data = metadata_sub)

# === 7. PERMANOVA summary table ===
extract_results <- function(adonis_obj, model_name) {
  tibble(
    Model= model_name,
    R2 = adonis_obj["Model", "R2"],
    P_value = adonis_obj["Model", "Pr(>F)"]
  )
}

results <- bind_rows(
  extract_results(adonis_type,         "type (univariate, global with water/sediment)"),
  extract_results(adonis_bio_global,   "global multivariate (biofilm + subcuticle only)"),
  extract_results(adonis_clases,       "clases (biofilm + subcuticle)"),
  extract_results(adonis_sp,           "sp (biofilm + subcuticle)"),
  extract_results(adonis_site,         "site (biofilm + subcuticle)"),
  extract_results(adonis_site_loc,     "site_loc (biofilm + subcuticle)"),
  extract_results(adonis_sub_global,   "global multivariate (subcuticle only)"),
  extract_results(adonis_sub_clases,   "clases (subcuticle only)"),
  extract_results(adonis_sub_sp,       "sp (subcuticle only)"),
  extract_results(adonis_sub_site,     "site (subcuticle only)"),
  extract_results(adonis_sub_site_loc, "site_loc (subcuticle only)")
) %>% arrange(desc(R2))


