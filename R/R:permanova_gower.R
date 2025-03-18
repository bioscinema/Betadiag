#' PERMANOVA with Gram Matrix
#'
#' This function performs PERMANOVA (Permutational Multivariate Analysis of Variance) using a Gram matrix and metadata.
#'
#' @param G A Gram matrix.
#' @param metadata A data frame containing metadata for the samples.
#' @param nperm The number of permutations to perform. Default is 999.
#' @return A list containing the pseudo-F statistic and the p-value.
#' @export
permanova_gower <- function(G, metadata, nperm = 999) {
  n <- nrow(G)
  
  X <- model.matrix(~ ., data = metadata)
  H <- X %*% solve(t(X) %*% X) %*% t(X)  
  G_H <- H %*% G %*% H  
  SS_model <- sum(diag(G_H))  
  
  I_H <- diag(n) - H  
  G_res <- I_H %*% G %*% I_H  
  SS_residual <- sum(diag(G_res)) 
  
  df_model <- ncol(X) - 1
  df_residual <- n - ncol(X)
  
  pseudo_F <- (SS_model / df_model) / (SS_residual / df_residual) 
  
  perm_F <- numeric(nperm)
  for (i in 1:nperm) {
    perm_index <- sample(1:n) 
    G_perm <- G[perm_index, perm_index]  
    G_H_perm <- H %*% G_perm %*% H  
    SS_model_perm <- sum(diag(G_H_perm))
    
    G_res_perm <- I_H %*% G_perm %*% I_H  
    SS_residual_perm <- sum(diag(G_res_perm))
    
    perm_F[i] <- (SS_model_perm / df_model) / (SS_residual_perm / df_residual)
  }
  
  p_value <- mean(perm_F >= pseudo_F)
  
  return(list(
    pseudo_F = pseudo_F,
    p_value = p_value
  ))
}
