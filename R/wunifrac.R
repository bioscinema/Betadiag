#' Compute Weighted UniFrac distance matrix from a phyloseq object (metric version)
#'
#' This function computes the weighted UniFrac distance between samples using the phylogenetic tree 
#' and normalized abundances stored in a \code{phyloseq} object.
#'
#' @param physeq A \code{phyloseq} object with OTU table and rooted phylogenetic tree.
#'
#' @return A \code{dist} object representing the weighted UniFrac distance matrix between samples.
#'
#' @details This function manually computes the weighted UniFrac distances without using 
#' \code{phyloseq::UniFrac}. It supports custom workflows or further integration in pipelines.
#'
#' @importFrom stats as.dist
#' @importFrom phyloseq otu_table phy_tree
#' @export
wunifrac <- function(physeq) {
  if (!inherits(physeq, "phyloseq")) {
    stop("Input must be a phyloseq object")
  }
  
  # Extract and normalize OTU table
  otu <- as(otu_table(physeq), "matrix")
  otu <- sweep(otu, 2, colSums(otu), "/")  # column-wise relative abundance
  
  # Extract tree structure
  tree <- phy_tree(physeq)
  if (is.null(tree$edge.length)) {
    stop("The phylogenetic tree must have branch lengths (edge.length)")
  }
  
  edge_matrix <- tree$edge
  edge_lengths <- tree$edge.length
  tip_labels <- tree$tip.label
  num_tips <- length(tip_labels)
  num_nodes <- tree$Nnode
  total_nodes <- num_tips + num_nodes
  
  node_names <- c(tip_labels, paste0("N", (num_tips + 1):total_nodes))
  node_abundances <- matrix(0, nrow = total_nodes, ncol = ncol(otu),
                            dimnames = list(node_names, colnames(otu)))
  node_abundances[tip_labels, ] <- otu[tip_labels, , drop = FALSE]
  
  # Bottom-up tree traversal: accumulate child → parent
  for (i in seq_len(nrow(edge_matrix))) {
    parent <- edge_matrix[i, 1]
    child <- edge_matrix[i, 2]
    parent_name <- node_names[parent]
    child_name <- node_names[child]
    node_abundances[parent_name, ] <- node_abundances[parent_name, ] + node_abundances[child_name, ]
  }
  
  # Compute weighted UniFrac distances
  samples <- colnames(otu)
  n_samples <- length(samples)
  dist_matrix <- matrix(0, nrow = n_samples, ncol = n_samples,
                        dimnames = list(samples, samples))
  
  for (i in 1:(n_samples - 1)) {
    for (j in (i + 1):n_samples) {
      s1 <- samples[i]
      s2 <- samples[j]
      diff_abundance <- abs(node_abundances[node_names[edge_matrix[, 2]], s1] -
                              node_abundances[node_names[edge_matrix[, 2]], s2])
      weighted_diff <- sum(edge_lengths * diff_abundance, na.rm = TRUE)
      dist_matrix[i, j] <- dist_matrix[j, i] <- weighted_diff
    }
  }
  
  return(as.dist(dist_matrix))
}
