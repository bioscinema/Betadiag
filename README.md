# Betadiag
Diagnostics and remedy tools for beta diversity analysis
This repository provides a suite of tools for diagnosing and addressing issues in beta diversity analysis. These tools help assess data quality, identify potential biases, and apply corrections to improve the robustness of diversity comparisons.

Features:
Diagnostic checks for beta diversity metrics;
Identification of potential biases and anomalies;
Remediation strategies for improving analysis accuracy;
Support for various ecological and microbiome datasets.

Walk-Through Example

In this example, we will use the IBD_16s_data_V4.RData dataset to demonstrate the diagnostic and remedy tools.

Data Cleaning

   1.Extract sample data and remove rows with missing or invalid values. 
   
   2. Create a metadata dataframe.
      
   3. Prepare the OTU table.

Proposed Pipeline with Diagnostics and Remedy

Step 1: Diagnostics
1.Calculate the weighted UniFrac distance:

wu_dist <- phyloseq::distance(physeq, method = "wunifrac")
wu_matrix <- as.matrix(wu_dist)

2. Perform diagnostic checks:

wu.check <- Triangle.Check(wu_matrix)
wu.check$is.metric # Check if the distance matrix is metric

plot(wu.check$collinearity.score, xlab = "Sample", ylab = "Collinearity Score")
plot(wu.check$nonlinearity.score, xlab = "Sample", ylab = "Nonlinearity Score")
 
3.Convert dissimilarity to Gram matrix and check Euclidean properties:

wu_kernel <- Dissimilarity_to_Gram(wu_matrix)
wu.euc.check <- Euclidean.Check(wu_kernel) # Check Euclidean properties

4.Understand the consequences:

wu.pcoa <- pcoa_gower(wu_kernel)
wu.permanova <- permanova_gower(wu_kernel, metadata)

Step 2: Remedy
 
 1.Apply the Higham method:

wu.higham <- Remedy_Gram(wu_kernel, method = "Higham")
wu.higham.pcoa <- pcoa_gower(wu.higham)
wu.higham.permanova <- permanova_gower(wu.higham, metadata)

 2.Apply the Tikhonov method:

wu.tik <- Remedy_Gram(wu_kernel, method = "Tikhonov", epsilon = 0)
wu.tik.pcoa <- pcoa_gower(wu.tik)
wu.tik.permanova <- permanova_gower(wu.tik, metadata)

Conclusion

By following these steps, you can diagnose and address issues in your beta diversity analysis, ensuring more reliable and robust results.
