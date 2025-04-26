#' Check Euclidean property of a Gram matrix
#'
#' This function determines whether a given Gram (inner-product) matrix can be
#' embedded in Euclidean space and computes its fraction of negative inertia
#' (FNI).
#'
#' @param G A symmetric Gram matrix (numeric matrix).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{FNI}{Fraction of negative inertia, calculated as the absolute
#'     sum of negative eigenvalues divided by the absolute sum of all
#'     eigenvalues.}
#'   \item{is.Euclidean}{Indicator of Euclidean embeddability
#'     (1 = Euclidean, 0 = non-Euclidean).}
#' }
#' @export
euclidean_check <- function(G) {
  eig <- eigen(G)
  lambda <- eig$values
  lambda[abs(lambda) < 1e-12] <- 0
  neg_sum <- sum(abs(lambda[lambda < 0]))
  total_sum <- sum(abs(lambda))
  ratio <- neg_sum / total_sum
  
  return(list(FNI = ratio, is.Euclidean = ifelse(ratio > 0, 0, 1)))
}