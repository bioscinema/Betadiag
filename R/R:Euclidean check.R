#' Check if Gram Matrix is Euclidean
#'
#' This function checks if a given Gram matrix is Euclidean.
#'
#' @param G A Gram matrix.
#' @return A list containing the ratio of negative eigenvalues to the total sum of eigenvalues and a boolean indicating if the matrix is Euclidean.
#' @export
Euclidean.Check <- function(G){
  # input Gram matrix
  eig <- eigen(G)
  lambda <- eig$values  
  small.values = which(abs(lambda)<1e-16)
  lambda[small.values] = 0
  neg_sum <- sum(abs(lambda[lambda < 0]))
  total_sum <- sum(abs(lambda))
  ratio <- neg_sum / total_sum
  
  # save results
  return(list(rate = ratio, is.Euclidean = ifelse(ratio>0,0,1)))
}
