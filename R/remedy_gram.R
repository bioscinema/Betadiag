#' Apply remedy methods to a dissimilarity matrix via Gram matrix correction
#'
#' This function applies remedial techniques to ensure the corresponding Gram matrix is 
#' positive semi-definite. It first converts the dissimilarity matrix to a Gram matrix, 
#' then applies either the \code{"Higham"} or \code{"Tikhonov"} method.
#'
#' @param D A symmetric dissimilarity matrix (distance matrix).
#' @param method Character string specifying the method. One of \code{"Higham"} or \code{"Tikhonov"}.
#' @param epsilon A small positive number used only for \code{"Tikhonov"} regularization.
#'
#' @return A modified Gram matrix that is positive semi-definite.
#' @export
remedy_gram <- function(D, method = c("Higham", "Tikhonov"), epsilon = NULL) {
  method <- match.arg(method)
  
  if (!is.matrix(D) || !is.numeric(D)) stop("D must be a numeric matrix.")
  if (!isSymmetric(D)) stop("D must be symmetric.")
  
  # Convert dissimilarity matrix to Gram matrix
  G <- dissimilarity_to_gram(D)
  
  # Eigen decomposition
  eig <- eigen(G)
  lambda <- eig$values
  vectors <- eig$vectors
  
  if (method == "Higham") {
    lambda0 <- pmax(lambda, 0)
    G_rem <- vectors %*% diag(lambda0) %*% t(vectors)
  }
  
  if (method == "Tikhonov") {
    if (is.null(epsilon)) stop("epsilon must be specified for method = 'Tikhonov'")
    
    lambda[abs(lambda) < 1e-16] <- 0
    lambda.min <- min(lambda)
    
    if (lambda.min == 0) {
      lambda0 <- lambda + epsilon
    } else if (lambda.min < 0) {
      lambda0 <- lambda + abs(lambda.min) + epsilon
    } else {
      lambda0 <- lambda
    }
    
    G_rem <- vectors %*% diag(lambda0) %*% t(vectors)
  }
  
  return(G_rem)
}

# ---- Internal Functions ----

#' @noRd
dissimilarity_to_gram <- function(D) {
  R <- ncol(D)
  J <- diag(R) - matrix(1, R, R) / R
  Gram <- -0.5 * J %*% (D^2) %*% J
  return(Gram)
}