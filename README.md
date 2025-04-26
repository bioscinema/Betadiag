# Betadiag

Diagnostics and Remedy Tools for Beta Diversity Analysis

Betadiag is an R package designed to diagnose, evaluate, and correct issues in beta diversity analysis, particularly when working with dissimilarity matrices from ecological and microbiome studies. It provides tools to assess matrix properties, perform ordination and statistical evaluation, and apply remedies to ensure valid analysis.

## Installation

To install the development version of Betadiag from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")
library(devtools)
devtools::install_github("bioscinema/Betadiag")
library(Betadiag)
```

## Key Features

- Check whether dissimilarity matrices are metric and Euclidean
- Assess triangle inequality violations, collinearity, and nonlinearity
- Convert dissimilarity matrices into Gram matrices
- Perform Principal Coordinate Analysis (PCoA) and eigenvalue decomposition
- Conduct PERMANOVA and pseudo R² analysis for association testing
- Apply Higham or Tikhonov corrections to enforce positive semi-definiteness
- Compute weighted UniFrac distances directly from phyloseq objects

## End-to-end demo: weighted UniFrac on IBD 16S data

<details>
<summary>Click to expand full script</summary>

```r
# Install Betadiag (dev) and dependencies ---------------------------------
install.packages("devtools")
library(devtools)
devtools::install_github("bioscinema/Betadiag")

library(Betadiag)
library(phyloseq)
library(ggplot2)

# Load dataset ------------------------------------------------------------
load("RealData/IBD_16s_data_V4.RData")

# ---- Data cleaning ------------------------------------------------------
sample_data <- sample_data(phy1)
rows_with_na <- apply(
  sample_data[, c(1, 18, 50, 57, 115, 126, 127)],
  1,
  function(x) any(x %in% c("not providednot provided",
                           "-88",
                           "not provided",
                           "cd (from uc 7/17/2018)"))
)
physeq <- prune_samples(!rows_with_na, phy1)

metadata <- data.frame(
  diagnosis         = as.factor(physeq@sam_data$diagnosis),
  age_at_diagnosis  = as.numeric(physeq@sam_data$age_at__diagnosis),
  host_age          = as.numeric(physeq@sam_data$host_age),
  host_height       = as.numeric(physeq@sam_data$host_height),
  race              = physeq@sam_data$race,
  sex               = as.numeric(ifelse(physeq@sam_data$sex == "male", 1, 0)),
  smoking           = as.numeric(ifelse(physeq@sam_data$smoking == "n", 0, 1))
)

# ---- Diagnostics --------------------------------------------------------
wu_dist   <- phyloseq::distance(physeq, method = "wunifrac")
wu_matrix <- as.matrix(wu_dist)

wu.check <- check_distance(wu_matrix)

# Basic flags
wu.check$is.metric
wu.check$is.Euclidean
```

```text
[1] 0
[1] 0
```

```r
# ---- Baseline evaluation ------------------------------------------------
wu.evalution <- evaluate_beta(
  wu_matrix,
  as.data.frame(metadata[, 1]),  # Y
  as.data.frame(metadata[, -1]), # Z
  metadata
)

wu.evalution$pseudo_F
wu.evalution$permanova_p
wu.evalution$pseudo_R2
```

```text
[1] ...F...
[1] ...p-value...
[1] ...R²...
```

```r
# ---- Remedies -----------------------------------------------------------
# Higham projection
wu.higham           <- remedy_gram(wu_matrix, method = "Higham")
wu.higham.evalution <- evaluate_beta_gram(
  wu.higham,
  as.data.frame(metadata[, 1]),
  as.data.frame(metadata[, -1]),
  metadata
)
wu.higham.euclidean <- euclidean_check(wu.higham)

# Tikhonov ridge (ε = 0)
wu.tik           <- remedy_gram(wu_matrix, method = "Tikhonov", epsilon = 0)
wu.tik.evalution <- evaluate_beta_gram(
  wu.tik,
  as.data.frame(metadata[, 1]),
  as.data.frame(metadata[, -1]),
  metadata
)
wu.tik.euclidean <- euclidean_check(wu.tik)
```

```text
#> Higham Euclidean?  TRUE
#> Tikhonov Euclidean? TRUE
```

```r
# ---- Quick comparison table --------------------------------------------
data.frame(
  Matrix      = c("Raw", "Higham", "Tikhonov"),
  Euclidean   = c(wu.check$is.Euclidean,
                  wu.higham.euclidean$is.Euclidean,
                  wu.tik.euclidean$is.Euclidean),
  Pseudo_R2   = c(wu.evalution$pseudo_R2,
                  wu.higham.evalution$pseudo_R2,
                  wu.tik.evalution$pseudo_R2),
  PERM_pval   = c(wu.evalution$permanova_p,
                  wu.higham.evalution$permanova_p,
                  wu.tik.evalution$permanova_p)
)
```

```text
     Matrix Euclidean Pseudo_R2 PERM_pval
1       Raw     FALSE     ...     ...
2    Higham      TRUE     ...     ...
3 Tikhonov      TRUE     ...     ...
```

</details>
