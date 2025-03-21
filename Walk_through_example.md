---
title: "Walk through Example"
author: ""
date: "`r format(Sys.Date(), '%Y-%m-%d')`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Load Libraries

```{r load-libraries}
library(phyloseq)
library(ggplot2)
source("Triangle inequality check.R")
source("Dissimilarity_to_Gram.R")
source("Euclidean check.R")
source("pcoa_gower.R")
source("permanova_gower.R")
source("Remedy_Gram.R")
source("MiRKAT_R.R")
```

## Load Data

```{r load-data}
load("IBD_16s_data_V4.RData")
```

## Data Cleaning

```{r data-cleaning}
sample_data <- sample_data(phy1)
rows_with_na <- apply(sample_data[, c(1,18, 50, 57, 115, 126, 127)], 1, function(x) any(x %in% c("not providednot provided", "-88", "not provided", "cd (from uc 7/17/2018)")))
physeq <- prune_samples(!rows_with_na, phy1)

metadata <- data.frame(
  diagnosis = as.factor(physeq@sam_data$diagnosis),
  age_at_diagnosis = as.numeric(physeq@sam_data$age_at__diagnosis),
  host_age = as.numeric(physeq@sam_data$host_age),
  host_height = as.numeric(physeq@sam_data$host_height),
  race = physeq@sam_data$race,
  sex = as.numeric(ifelse(physeq@sam_data$sex == "male", 1, 0)),
  smoking = as.numeric(ifelse(physeq@sam_data$smoking == "n", 0, 1))
)
```

## Proposed Pipeline with Diagnostics and Remedy

### Distance Calculation

```{r distance-calculation}
otu_table <- otu_table(physeq)
if (taxa_are_rows(physeq)) {
  otu_table <- t(otu_table)
}
otu_matrix <- as.matrix(otu_table)

wu_dist <- phyloseq::distance(physeq, method = "wunifrac")
wu_matrix = as.matrix(wu_dist)
```

### Step 1: Diagnostics

```{r diagnostics}
wu.check <- Triangle.Check(wu_matrix)
wu.check$is.metric # non-metric

plot(wu.check$collinearity.score, xlab = "Sample", ylab = "Collinearity Score") 
plot(wu.check$nonlinearity.score, xlab = "Sample", ylab = "Nonlinearity Score")

wu_kernel <- Dissimilarity_to_Gram(wu_matrix)
wu.euc.check <- Euclidean.Check(wu_kernel)
wu.euc.check

wu.pcoa <- pcoa_gower(wu_kernel)
wu.permanova <- permanova_gower(wu_kernel, metadata)
wu.permanova
```

### PCoA Visualization

```{r pcoa-visualization}
pcoa_points <- as.data.frame(wu.pcoa$scores[, c(1,2)])
colnames(pcoa_points) <- c("PC1", "PC2")
group_info <- metadata$diagnosis
pcoa_points$Group <- group_info

ggplot(pcoa_points, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Group), type = "norm", linetype = 2, size = 1) +  
  theme_minimal() +
  scale_color_manual(values = c("uc" = "#177cd5", "cd" = "#eb1f15", "ind" = "#2ca02c")) +
  labs(title = "PCoA based on Weighted Unifrac Distance", 
       x = paste0("PC1 (", round(wu.pcoa$PCo1.rate * 100, 2), "%)"),
       y = paste0("PC2 (", round(wu.pcoa$PCo2.rate * 100, 2), "%)"),
       color = "Group") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(wu.permanova$p_value, 3)), vjust = 2, hjust = 2, size = 5)
```

### Step 2: Remedy - Method 1 Higham Method

```{r remedy-higham}
wu.higham <- Remedy_Gram(wu_kernel, method = "Higham")
wu.higham.pcoa <- pcoa_gower(wu.higham)
wu.higham.permanova <- permanova_gower(wu.higham, metadata)
wu.higham.permanova
```

### Step 2: Remedy - Method 2 Tikhonov Method

```{r remedy-tikhonov}
wu.tik <- Remedy_Gram(wu_kernel, method = "Tikhonov", epsilon = 0)
wu.tik.pcoa <- pcoa_gower(wu.tik)
wu.tik.permanova <- permanova_gower(wu.tik, metadata)
wu.tik.permanova
```
