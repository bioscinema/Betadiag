#' Principal Coordinates Analysis (PCoA) with Gram Matrix
#'
#' This function performs Principal Coordinates Analysis (PCoA) using a Gram matrix as input.
#'
#' @param G A Gram matrix.
#' @return A list containing eigenvalues, scores, rates for the first two principal coordinates, and a boolean indicating if the matrix is valid.
#' @export
pcoa_gower <- function(G) {
  ### PCoA with Gram matrix as input
  eig <- eigen(G, symmetric = TRUE)
  values <- eig$values
  small.values = which(abs(values)<1e-16)
  values[small.values] = 0
  vectors <- eig$vectors
  pos_index <- values > 0
  positive_values <- values[pos_index]
  positive_vectors <- vectors[, pos_index, drop = FALSE]
  scores <- positive_vectors %*% diag(sqrt(positive_values))
  is_valid <- all(values >= 0)
  
  PCo1.rate = abs(values)[1]/sum(abs(values)[values>0])
  PCo2.rate = abs(values)[2]/sum(abs(values)[values>0])
  return(list(
    eigenvalues = values,
    scores = scores,
    PCo1.rate = PCo1.rate,
    PCo2.rate = PCo2.rate,
    is.valid = is_valid
  ))
}
