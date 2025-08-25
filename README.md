# subcuticular-microbiome-echinoderms

## 0. Note on Scripts and Data

All scripts described below are included in the scripts folder. The datasets required to run these analyses are not publicly available yet, as they are part of an ongoing project.


## 1. Project Objective
Scripts for functional and coevolutionary analysis of animal microbiome datasets, originally developed for the echinoderm subcuticular microbiome inferred with FAPROTAX and PICRUSt2. The focus is on differences in composition and functional structure between microbial communities, using analyses such as PERMANOVA, Gini index calculation, and Pagel’s λ.

## 2. Scripts Overview

Below is a list of all scripts included in this repository, with a brief description:

1. [01_faprotax_functional_profiles.R](scripts/01_faprotax_functional_profiles.R) – Visualizes functional profiles inferred with FAPROTAX using bubble plots.
2. [02_otu_assignment_coverage.R](scripts/02_otu_assignment_coverage.R) – Evaluates OTU assignment coverage and quality from PICRUSt2 predictions.
3. [03_weighted_nsti_analysis.R](scripts/03_weighted_nsti_analysis.R) – Computes abundance-weighted NSTI values per sample and compares across groups.
4. [04_differential_abundance_picrust2.R](scripts/04_differential_abundance_picrust2.R) – Performs differential abundance analysis of PICRUSt2-inferred functions and generates heatmaps.
5. [05_permanova.R](scripts/05_permanova.R) – Performs PERMANOVA to test for differences in community composition across groups.
6. [06_gini_functional_contributions.R](scripts/06_gini_functional_contributions.R) – Calculates Gini indices for functional contribution inequality per OTU.
7. [07_phylogenetic_signal_abundance.R](scripts/07_phylogenetic_signal_abundance.R) – Estimates phylogenetic signal (Pagel’s λ) for average functional abundances.
8. [08_phylogenetic_signal_gini.R](scripts/08_phylogenetic_signal_gini.R) – Estimates phylogenetic signal (Pagel’s λ) for Gini indices of functional contributions per host species.



## 3. Input File Requirements

All scripts require input files to follow a consistent format:

1. **OTU/Function Abundance Tables (.tsv)**
   - Rows: OTUs or functions
   - Columns/Header: Samples
   - Example:
     
     | OTU  | Sample1 | Sample2 | Sample3 |
     |:----:|:-------:|:-------:|:-------:|
     | OTU1 |   12    |    5    |    0    |
     | OTU2 |    0    |    8    |    3    |
    

 2. **Metadata Table (metadata.tsv)**
    - Rows: Samples
    - Columns/Header: Metadata variables (e.g., SampleID, HostSpecies, Location)
    - Example:

      | SampleID | HostClasses | Location  |
      |:--------:|:-----------:|:--------:|
      | Sample1  | Starfish    | Granada  |
      | Sample2  | Ophiura     | Ceuta    |
     
 
3. **Phylogenetic Tree (tree.txt)**
    - Format: Newick
   

4.  **NSTI values for each OTU (combined_marker_predicted_and_nsti.tsv)**
    - Rows: One OTU per row
    - Columns / Header:
      1. sequence – OTU identifier
      2. 16S_rRNA_Count – Number of 16S rRNA reads assigned to the OTU
      3. metadata_NSTI – NSTI value for the OTU
      4. best_domain – Taxonomic domain of closest reference genome (bac for bacteria)
      5. closest_reference_genome – Accession number of closest reference genome
    - Example:

      | sequence  | 16S_rRNA_Count | metadata_NSTI | best_domain | closest_reference_genome |
      |:---------:|:--------------:|:-------------:|:-----------:|:------------------------:|
      | OTU49955  | 1              | 0.129679      | bac         | GCA_018882785.1          |
      | OTU20398  | 4              | 0.026109      | bac         | GCF_014337155.1          |
      | OTU52804  | 1              | 0.147661      | bac         | GCF_001746755.1          |
      | OTU48340  | 1              | 0.044553      | bac         | GCA_023266095.1          |
     

 5. **Gini Indices per Function and Host Species (gini_function_per_species_filtered_wide.xlsx)**
    - Rows: Microbial functions (e.g., MetaCyc pathways)
    - Columns / Header:
      1. function – Identifier or name of the function
      2. Remaining columns – Gini indices for each host species (one column per species)
    - Example:
      
      | function  | *Antedon_mediterranea* | *Arbacia_lixula* | *Arbaciella_elegans* |
      |:---------:|:---------------------:|:----------------:|:------------------:|
      | PWY-6126  | 0.936761692           | 0.831952092      | 0.809817721        |
      | PWY-I9    | 0.935496167           | 0.704158826      | 0.748473582        |
     


## 4. Project Structure
### Script 01 – Visualization of Functional Profiles (FAPROTAX)
Processes and visualizes bacterial community functional profiles predicted with FAPROTAX. Functions with zero abundance are removed, and samples can be subset (e.g., subcuticle or sediment). The script is parameterized for use with other datasets.

Functional differences are visualized with bubble plots, where point size and color reflect mean relative abundances.

**Input**:
- func_table_def.tsv - matrix of functional abundances inferred with FAPROTAX
- metadata.tsv - sample metadata for grouping

**Output**:
- Bubble plots showing functional variation across samples

### Script 02 – OTU Assignment Coverage
Evaluates the coverage and quality of functional predictions from PICRUSt2. Uses two input files: the original OTU abundance table and the NSTI table (Nearest Sequenced Taxon Index). Calculates the proportion of OTUs with functional predictions and the percentage excluded.

**Input**:
- OTU_table.tsv - original OTU abundance matrix
- combined_marker_predicted_and_nsti.tsv - NSTI values for each OTU
  
**Output**:
- Summary table printed to console with total OTUs, predicted OTUs, and excluded OTUs.

### Script 03 - Nearest Sequenced Taxon Index (NSTI) Evaluation
Computes a weighted NSTI per sample, representing an abundance-weighted mean of NSTI values. Compares NSTI values across different sample groups using the non-parametric Kruskal-Wallis test.

**Input**:
- OTU_table.tsv
- combined_marker_predicted_and_nsti.tsv
- metadata.tsv
  
**Output**:
- Summary tables of weighted NSTI comparisons across groups.

### Script 04 – Differential Abundance Analysis (PICRUSt2 functions)
Performs differential abundance analysis (DAA) for functions inferred with PICRUSt2, comparing sample types. Implemented using the ggpicrust2 R package and the ALDEx2 method, which includes internal CLR transformation and Monte Carlo sampling. Both the Wilcoxon rank-sum test (non-parametric) and Welch’s t-test (parametric) are applied; conclusions prioritize Wilcoxon results for robustness.

For visualization, functions with adjusted p-values < 0.05 (Benjamini–Hochberg correction) were selected, annotated via pathway_annotation from ggpicrust2, and normalized by row (Z-score scaling). Heatmaps were then generated with heatmap.2 (gplots package), including hierarchical clustering of both functions and samples.

**Input**: 
- path_abundance_data.tsv - abundances inferred by PICRUSt2 corresponding to MetaCyc annotations
- metadata.tsv
  
**Output**:
- Heatmaps showing functional patterns across samples.

### Script 05 – PERMANOVA:
Performs multivariate analysis of variance with permutations using adonis2() from the vegan R package. Tests for significant differences in community composition between groups using Bray-Curtis dissimilarities. Can analyze the full dataset or subsets (sample types, taxonomic levels). Reports variance explained (R²) and p-values (999 permutations).

**Input**: 
- path_abundance_data.tsv
- metadata.tsv
  
**Output**:
- PERMANOVA summary tables for all tested factors.

### Script 06 – Functional Contribution Inequality (Gini Index)
Quantifies inequality in OTU functional contributions using the Gini index, implemented via the ineq R package. Operates on normalized OTU functional contributions. Higher Gini values indicate a few taxa dominate specific functions, while lower values indicate more even distribution.

**Input**: 
- path_abun_contrib.tsv - normalized functional contribution matrix
- metadata.tsv
  
**Output**:
- Gini index calculated for each function.

### Scripts 07 & 08 – Phylogenetic Signal Analysis (Pagel’s λ)
Estimate the phylogenetic signal of microbial functions or functional contribution inequality using the geiger R package.

#### Script 07 – Average Functional Abundances
**Input**:
- metadata.tsv
- path_abundance_data.tsv
- tree.txt – host phylogenetic tree
  
**Output**:
- Matrix of Pagel’s λ values for each function, where λ ~1 indicates strong phylogenetic structuring and λ ~0 indicates independence from phylogeny.

#### Script 08 – Gini Indices per Function and Host Species
**Input**:
- metadata.tsv
- gini_function_per_species_filtered_wide.xlsx - Gini indices per function and host species (from Script 6)
- tree.txt
  
**Output**:
- Matrix of Pagel’s λ values for Gini indices of functional contribution, interpreted as above.

