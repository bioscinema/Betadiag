# Betadiag
Diagnostics and remedy tools for beta diversity analysis
This repository provides a suite of tools for diagnosing and addressing issues in beta diversity analysis. These tools help assess data quality, identify potential biases, and apply corrections to improve the robustness of diversity comparisons.

Features:
Diagnostic checks for beta diversity metrics;
Identification of potential biases and anomalies;
Remediation strategies for improving analysis accuracy;
Support for various ecological and microbiome datasets.

# Dissimilarity Analysis in R

## Overview
This R package provides tools for analyzing dissimilarity matrices, converting them into Gram matrices, performing diagnostic checks, and applying remedies for non-Euclidean or non-metric properties. It includes functions for Principal Coordinate Analysis (PCoA), PERMANOVA, and various remedial techniques to improve matrix properties.

## Installation

You can install the package by sourcing the provided R scripts manually:
```r
source("Dissimilarity_to_Gram.R")
source("Euclidean check.R")
source("pcoa_gower.R")
source("permanova_gower.R")
source("Remedy_Gram.R")
source("Triangle inequality check.R")
source("weighted_unifrac_calc.R")
```

## Functions

### 1. `Dissimilarity_to_Gram(D)`
Converts a dissimilarity matrix `D` into a Gram matrix.

### 2. `Euclidean.Check(G)`
Checks whether a Gram matrix `G` is Euclidean by analyzing its eigenvalues.

### 3. `pcoa_gower(G)`
Performs Principal Coordinate Analysis (PCoA) using a Gram matrix as input.

### 4. `permanova_gower(G, metadata, nperm=999)`
Performs PERMANOVA (Permutation Multivariate Analysis of Variance) on a Gram matrix `G` given metadata.

### 5. `Remedy_Gram(G, method, epsilon=NULL)`
Applies remedial techniques to correct a non-Euclidean Gram matrix. Methods include:
- `Higham`: Nearest Positive Semi-Definite (PSD) projection.
- `Tikhonov`: Relaxed Positive Definiteness adjustment.

### 6. `Triangle.Check(D)`
Checks if a dissimilarity matrix `D` satisfies the triangle inequality.

### 7. `weighted_unifrac_calc(physeq)`
Computes the weighted UniFrac distance matrix for a `phyloseq` object.

## Example Usage
A full pipeline demonstrating these functions is provided in `walk-through example.R`. Example steps include:

```r
# Load necessary libraries
library(phyloseq)

# Load and preprocess data
source("walk-through example.R")

# Compute weighted UniFrac distance
wu_dist <- weighted_unifrac_calc(physeq)
wu_matrix <- as.matrix(wu_dist)

# Check dissimilarity properties
wu_check <- Triangle.Check(wu_matrix)
wu_kernel <- Dissimilarity_to_Gram(wu_matrix)
wu_euc_check <- Euclidean.Check(wu_kernel)

# Perform PCoA and PERMANOVA
wu_pcoa <- pcoa_gower(wu_kernel)
wu_permanova <- permanova_gower(wu_kernel, metadata)

# Apply remedies
wu_higham <- Remedy_Gram(wu_kernel, method = "Higham")
wu_tik <- Remedy_Gram(wu_kernel, method = "Tikhonov", epsilon = 0)
```

## License
This package is distributed under the MIT License.

## Acknowledgments
This package makes use of `phyloseq` for microbial community analysis and various statistical techniques for matrix diagnostics and transformations.

