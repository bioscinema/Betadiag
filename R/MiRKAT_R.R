#' MiRKAT_R Function
#'
#' Computes the R-squared statistic based on kernel matrices.
#'
#' @param G A square matrix representing the kernel matrix.
#' @param Y A response variable, either a vector or matrix.
#' @param metadata A data frame containing covariates.
#' @return A numeric value representing the R-squared statistic.
#' @export
MiRKAT_R <- function(G, Y, metadata) {
  # Ensure input validity
  if (!is.matrix(G) || !is.numeric(G)) stop("G must be a numeric matrix.")
  if (!is.data.frame(metadata)) stop("metadata must be a data frame.")
  
  n <- nrow(G)
  
  # Generate design matrices
  X <- model.matrix(~ ., data = metadata)
  Y <- model.matrix(~ ., data = Y)
  
  # Compute hat matrix and residuals
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  I_H <- diag(n) - H  
  R <- I_H %*% Y
  
  # Eigen decomposition of G
  eig <- eigen(G, symmetric = TRUE)
  Gd <- eig$vectors %*% diag(abs(eig$values)) %*% t(eig$vectors)
  
  # Compute R-squared statistic
  Rsquared <- sum(diag(t(R) %*% G %*% R)) / (sum(diag(Gd)) * sum(diag(t(R) %*% R)))
  
  return(Rsquared)
}
