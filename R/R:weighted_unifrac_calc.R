#' Calculate Weighted UniFrac Distance
#'
#' This function calculates the weighted UniFrac distance for a given phyloseq object.
#'
#' @param physeq A phyloseq object.
#' @return A dissimilarity matrix representing the weighted UniFrac distances.
#' @export
weighted_unifrac_calc <- function(physeq) {
  if (!inherits(physeq, "phyloseq")) {
    stop("Input must be a phyloseq object")
  }
  
  # Abundance table for all nodes
  otu <- as(otu_table(physeq), "matrix")
  otu <- sweep(otu, 2, colSums(otu), "/") 
  tree <- phy_tree(physeq)
  edge_matrix <- tree$edge
  edge_lengths <- tree$edge.length
  tip_labels <- tree$tip.label
  num_tips <- length(tip_labels) 
  num_nodes <- tree$Nnode 
  total_nodes <- num_tips + num_nodes 
  node_names <- c(tip_labels, paste0("N", (num_tips + 1):total_nodes))
  node_abundances <- matrix(0, nrow = total_nodes, ncol = ncol(otu))
  rownames(node_abundances) <- node_names
  colnames(node_abundances) <- colnames(otu)
  node_abundances[tip_labels, ] <- otu[tip_labels, ]  
 
  for (i in seq_len(nrow(edge_matrix))) {
    parent <- edge_matrix[i, 1]  
    child <- edge_matrix[i, 2]  
    node_abundances[node_names[parent], ] <- node_abundances[node_names[parent], ] + node_abundances[node_names[child], ]
  }
  
  # Calculate dissimilarity matrix
  samples <- colnames(otu)
  n_samples <- length(samples)
  dist_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
  rownames(dist_matrix) <- colnames(dist_matrix) <- samples
  for (i in 1:(n_samples - 1)) {
    for (j in (i + 1):n_samples) {
      sample1 <- samples[i]
      sample2 <- samples[j]
      diff_abundance <- abs(node_abundances[edge_matrix[,2], sample1] - node_abundances[edge_matrix[,2], sample2])
      weighted_diff <- sum(edge_lengths * diff_abundance, na.rm = TRUE)
      dist_matrix[i, j] <- weighted_diff
      dist_matrix[j, i] <- weighted_diff
    }
  }
  
  return(as.dist(dist_matrix))
}
