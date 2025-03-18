#' Convert Dissimilarity Matrix to Gram Matrix
#'
#' This function converts a dissimilarity matrix to a Gram matrix.
#'
#' @param D A dissimilarity matrix.
#' @return A Gram matrix.
#' @export
Dissimilarity_to_Gram <- function(D){
  # convert dissimilarity matrix to Gram matrix
  R = ncol(D)
  J = diag(R) - (matrix(1,R,1) %*% matrix(1,1,R))/R
  Gram = -0.5 * J %*% D^2 %*% J
  return(Gram)
}
