library(phyloseq)

source("Triangle inequality check.R")
source("Dissimilarity_to_Gram.R")
source("Euclidean check.R")

source("pcoa_gower.R")
source("permanova_gower.R")

source("Remedy_Gram.R")

load("IBD_16s_data_V4.RData")

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
wu.check <- Triangle.Check(wu_matrix)
wu.check$is.metric # non-metric

plot(wu.check$collinearity.score,xlab = "Sample",ylab = "Collinearity Score") 
plot(wu.check$nonlinearity.score,xlab = "Sample",ylab = "Nonlinearity Score") 
# This dataset suffer more from nonlinearity than collinearity using weighted Unifrac distance, hence tSNE or UMAP would be better

wu_kernel <- Dissimilarity_to_Gram(wu_matrix)
wu.euc.check <- Euclidean.Check(wu_kernel) # rate = 11.8%, non-Euclidean

# understand the consequence
wu.pcoa <- pcoa_gower(wu_kernel)
wu.permanova <- permanova_gower(wu_kernel,metadata)

### step 2: remedy
# Method 1: Higham
wu.higham <- Remedy_Gram(wu_kernel, method = "Higham")
wu.higham.pcoa <- pcoa_gower(wu.higham)
wu.higham.permanova <- permanova_gower(wu.higham,metadata)

# method 2: Tikhonov
wu.tik <- Remedy_Gram(wu_kernel, method = "Tikhonov",epsilon = 0)
wu.tik.pcoa <- pcoa_gower(wu.tik)
wu.tik.permanova <- permanova_gower(wu.tik,metadata)

