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

pcoa_points <- as.data.frame(wu.evalution$scores[,c(1,2)])
colnames(pcoa_points) <- c("PC1", "PC2")
group_info <- metadata$diagnosis
pcoa_points$Group <- group_info

ggplot(pcoa_points, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Group), type = "norm", linetype = 2, size = 1) +  
  theme_minimal() +
  scale_color_manual(values = c("uc" = "#177cd5", "cd" = "#eb1f15", "ind" = "#2ca02c")) +
  labs(title = "PCoA based on Weighted Unifrac Distance", 
       x = paste0("PC1 (", round(wu.evalution$PCo1.rate * 100, 2), "%)"),
       y = paste0("PC2 (", round(wu.evalution$PCo2.rate * 100, 2), "%)"),
       color = "Group") +
  annotate("text", x = Inf, y = Inf, label = paste0("p = ", signif(wu.evalution$permanova_p, 3)), vjust = 2, hjust = 2, size = 5)

### step 2: remedy
# Method 1: Higham
wu.higham <- remedy_gram(wu_matrix, method = "Higham")
wu.higham.evalution <- evaluate_beta(wu.higham, as.data.frame(metadata[,1]), as.data.frame(metadata[,-1]),metadata)


# method 2: Tikhonov
wu.tik <- remedy_gram(wu_matrix, method = "Tikhonov",epsilon = 0)
wu.tik.evalution <- evaluate_beta(wu.tik, as.data.frame(metadata[,1]), as.data.frame(metadata[,-1]),metadata)

