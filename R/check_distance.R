#' Check distance matrix properties
#'
#' This function checks whether a given distance matrix is metric and Euclidean, 
#' and computes collinearity and nonlinearity scores, fraction of negative inertia (FNI) for each sample.
#'
#' @param D A symmetric distance matrix (numeric matrix).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{is.metric}{Whether the distance matrix satisfies the triangle inequality (0/1).}
#'   \item{is.Euclidean}{Whether the distance matrix corresponds to a Euclidean space (0/1).}
#'   \item{collinearity.score}{A numeric vector of collinearity scores for each sample.}
#'   \item{nonlinearity.score}{A numeric vector of nonlinearity scores for each sample.}
#'   \item{FNI}{Fraction of negative eigenvalues (from Gram matrix), same as `rate` in `Euclidean.Check()`.}
#' }
#' @export
check_distance <- function(D) {
  triangle_res <- triangle_check(D)
  G <- dissimilarity_to_gram(D)
  euclidean_res <- euclidean_check(G)
  
  return(list(
    is.metric = triangle_res$is.metric,
    is.Euclidean = euclidean_res$is.Euclidean,
    collinearity.score = triangle_res$collinearity.score,
    nonlinearity.score = triangle_res$nonlinearity.score,
    FNI = euclidean_res$rate
  ))
}


# ---- Internal Functions ----

#' @noRd
triangle_check <- function(D) {
  R <- ncol(D)
  neg.counter <- 0
  zero.counter <- 0
  collinearity.score <- rep(0, R)
  nonlinearity.score <- rep(0, R)
  
  for (p in 1:R) {
    for (q in 1:R) {
      if (q != p) {
        r <- setdiff(1:R, c(p, q))
        numN <- sum(D[p, q] + D[q, r] - D[p, r] < 0)
        num0 <- sum(D[p, q] + D[q, r] - D[p, r] == 0)
        neg.counter <- neg.counter + numN / (R * (R - 1) * (R - 2))
        zero.counter <- zero.counter + num0 / (R * (R - 1) * (R - 2))
        collinearity.score[p] <- collinearity.score[p] + num0 / (R * (R - 1) * (R - 2))
        nonlinearity.score[p] <- nonlinearity.score[p] + numN / (R * (R - 1) * (R - 2))
      }
    }
  }
  
  return(list(
    is.metric = ifelse(neg.counter > 0, 0, 1),
    neg.counter = neg.counter,
    zero.counter = zero.counter,
    collinearity.score = collinearity.score,
    nonlinearity.score = nonlinearity.score
  ))
}

#' @noRd
dissimilarity_to_gram <- function(D) {
  R <- ncol(D)
  J <- diag(R) - matrix(1, R, R) / R
  Gram <- -0.5 * J %*% (D^2) %*% J
  return(Gram)
}

#' @noRd
euclidean_check <- function(G) {
  eig <- eigen(G)
  lambda <- eig$values
  lambda[abs(lambda) < 1e-12] <- 0
  neg_sum <- sum(abs(lambda[lambda < 0]))
  total_sum <- sum(abs(lambda))
  ratio <- neg_sum / total_sum
  
  return(list(rate = ratio, is.Euclidean = ifelse(ratio > 0, 0, 1)))
}
