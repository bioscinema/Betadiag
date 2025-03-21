#' Perform PERMANOVA Using Gower's Distance
#'
#' Computes a pseudo-F statistic and p-value for a PERMANOVA test using a given Gower dissimilarity matrix.
#'
#' @param G A square dissimilarity matrix.
#' @param metadata A data frame containing covariates.
#' @return A list containing the pseudo-F statistic and its associated p-value.
#' @export
permanova_gower <- function(G, metadata) {
  # Ensure G is a valid numeric matrix
  if (!is.matrix(G) || !is.numeric(G)) stop("G must be a numeric matrix.")
  if (!is.data.frame(metadata)) stop("metadata must be a data frame.")

  # Eigen decomposition
  eig <- eigen(G, symmetric = TRUE)
  lambda <- eig$values
  lambda[abs(lambda) < 1e-16] <- 0  # Set near-zero eigenvalues to zero
  
  # Check for non-Euclidean dissimilarity matrix
  if (any(lambda < 0)) {
    stop("The dissimilarity matrix is non-Euclidean!")
  }
  
  r <- sum(lambda == 0)  # Dimension loss
  n <- nrow(G)
  X <- model.matrix(~ ., data = metadata)
  p <- ncol(X)
  
  # Warn if between-sample degrees of freedom are low
  if (n - r < p) {
    warning("The between-sample degrees of freedom are low!")
  }
  
  # Compute hat matrix
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  
  # Compute sum of squares
  G_H <- H %*% G %*% H  
  SS_model <- sum(diag(G_H))
  
  I_H <- diag(n) - H  
  G_res <- I_H %*% G %*% I_H  
  SS_residual <- sum(diag(G_res))
  
  # Compute degrees of freedom
  df_model <- p - 1
  df_residual <- n - r - p
  
  # Compute pseudo-F statistic and p-value
  pseudo_F <- (SS_model / df_model) / (SS_residual / df_residual)
  p_value <- pf(pseudo_F, df_model, df_residual, lower.tail = FALSE)
  
  return(list(
    pseudo_F = pseudo_F,
    p_value = p_value
  ))
}
