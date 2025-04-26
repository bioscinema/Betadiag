# Install devtools if not already installed
install.packages("devtools")
library(devtools)
devtools::install_github("bioscinema/Betadiag")
library(Betadiag)
library(phyloseq)
library(ggplot2)
load("RealData/IBD_16s_data_V4.RData")

################################################################
### data cleaning
################################################################
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

otu_table <- otu_table(physeq)
if (taxa_are_rows(physeq)) {
  otu_table <- t(otu_table)
}
otu_matrix <- as.matrix(otu_table)

################################################################
### Proposed pipeline with diagnostics and remedy
################################################################

### use weighted unifrac as an example
wu_dist <- phyloseq::distance(physeq, method = "wunifrac")
wu_matrix = as.matrix(wu_dist)

### step 1: diagnostics
wu.check <- check_distance(wu_matrix)
wu.check$is.metric # non-metric
wu.check$is.Euclidean # non-Euclidean(rate = 11.8%)

plot(wu.check$collinearity.score,xlab = "Sample",ylab = "Collinearity Score") 
plot(wu.check$nonlinearity.score,xlab = "Sample",ylab = "Nonlinearity Score") 
# This dataset suffer more from nonlinearity than collinearity using weighted Unifrac distance, hence tSNE or UMAP would be better

# understand the consequence

wu.evalution <- evaluate_beta(wu_matrix, as.data.frame(metadata[,1]), as.data.frame(metadata[,-1]),metadata)

### step 2: remedy
# Method 1: Higham
wu.higham <- remedy_gram(wu_matrix, method = "Higham")
wu.higham.evalution <- evaluate_beta_gram(wu.higham, as.data.frame(metadata[,1]), as.data.frame(metadata[,-1]),metadata)
wu.higham.euclidean <- euclidean_check(wu.higham)

# method 2: Tikhonov
wu.tik <- remedy_gram(wu_matrix, method = "Tikhonov",epsilon = 0)
wu.tik.evalution <- evaluate_beta_gram(wu.tik, as.data.frame(metadata[,1]), as.data.frame(metadata[,-1]),metadata)
wu.tik.euclidean <- euclidean_check(wu.tik)

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