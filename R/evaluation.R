#' Evaluate beta diversity association with outcome and confounders
#'
#' This function computes multiple evaluation metrics based on a dissimilarity matrix and metadata:
#' including PCoA eigen decomposition, PERMANOVA, and MiRKAT-style pseudo R².
#'
#' @param D A dissimilarity matrix (numeric and symmetric).
#' @param Y A data frame or one-column matrix representing the outcome variable(s).
#' @param Z A data frame of confounder variables.
#' @param metadata A data frame containing both outcome and confounders.
#' @param nperm Number of permutations for PERMANOVA (default = 999).
#'
#' @return A list containing:
#' \describe{
#'   \item{eigenvalues}{PCoA eigenvalues.}
#'   \item{scores}{PCoA scores from positive eigencomponents.}
#'   \item{PCo1.rate}{Proportion of variance explained by the first axis.}
#'   \item{PCo2.rate}{Proportion of variance explained by the second axis.}
#'   \item{PCoA.valid}{Whether all eigenvalues are non-negative.}
#'   \item{pseudo_F}{Pseudo F-statistic from PERMANOVA.}
#'   \item{pseudo_R2}{Pseudo R² from PERMANOVA.}
#'   \item{permanova_p}{Permutation p-value.}
#'   \item{MiRKAT_R2}{MiRKAT-style pseudo R² using residual projection.}
#' }
#' @importFrom stats model.matrix
#' @export
evaluate_beta <- function(D, Y, Z, metadata, nperm = 999) {
  # ---- Internal: Dissimilarity to Gram ----
  dissimilarity_to_gram <- function(D) {
    R <- ncol(D)
    J <- diag(R) - matrix(1, R, R) / R
    Gram <- -0.5 * J %*% (D^2) %*% J
    return(Gram)
  }
  
  # ---- Internal: PCoA ----
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
    
    PCo1.rate = abs(values)[1]/sum(abs(values))
    PCo2.rate = abs(values)[2]/sum(abs(values))
    return(list(
      eigenvalues = values,
      scores = scores,
      PCo1.rate = PCo1.rate,
      PCo2.rate = PCo2.rate,
      is.valid = is_valid
    ))
  }
  
  # ---- Internal: PERMANOVA ----
  permanova_gower <- function(G, metadata, nperm = 999) {
    n <- nrow(G)
    
    eig = eigen(G)
    values <- eig$values
    small.values = which(abs(values)<1e-16)
    values[small.values] = 0
    is.valid = all(values>=0)
    
    Gd = eig$vectors %*% diag(abs(values)) %*% t(eig$vectors)
    Gn = eig$vectors[,values>0] %*% diag(values[values>0]) %*% t(eig$vectors[,values>0])
    
    X <- model.matrix(~ ., data = metadata)
    H <- X %*% solve(t(X) %*% X) %*% t(X)  
    G_H <- H %*% Gn %*% H  
    SS_model <- sum(diag(G_H))  
    
    I_H <- diag(n) - H  
    G_res <- I_H %*% Gd %*% I_H  
    SS_residual <- sum(diag(G_res)) 
    
    df_model <- ncol(X) - 1
    df_residual <- n - ncol(X)
    
    pseudo_F <- (SS_model / df_model) / (SS_residual / df_residual) 
    
    pseudo_R2 = SS_model/(SS_model + SS_residual)
    
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
      pseudo_R2 = pseudo_R2,
      p_value = p_value,
      is.valid = is.valid
    ))
  }
  
  # ---- Internal: MiRKAT R2 ----
  mirkat_r <- function(G, Y, Z) {
    n <- nrow(G)
    X <- model.matrix(~ ., data = Z)
    Ymat <- model.matrix(~ ., data = Y)
    H <- X %*% solve(t(X) %*% X) %*% t(X)
    I_H <- diag(n) - H
    R <- I_H %*% Ymat
    eig <- eigen(G)
    Gd <- eig$vectors %*% diag(abs(eig$values)) %*% t(eig$vectors)
    Rsq <- sum(diag(t(R) %*% G %*% R)) / (sum(diag(Gd)) * sum(diag(t(R) %*% R)))
    return(Rsq)
  }
  
  # === Workflow ===
  G <- dissimilarity_to_gram(D)
  pcoa_res <- pcoa_gower(G)
  permanova_res <- permanova_gower(G, metadata, nperm = nperm)
  mirkat_r2 <- mirkat_r(G, Y, Z)
  
  return(list(
    eigenvalues = pcoa_res$eigenvalues,
    scores = pcoa_res$scores,
    PCo1.rate = pcoa_res$PCo1.rate,
    PCo2.rate = pcoa_res$PCo2.rate,
    PCoA.valid = pcoa_res$is.valid,
    pseudo_F = permanova_res$pseudo_F,
    pseudo_R2 = permanova_res$pseudo_R2,
    permanova_p = permanova_res$p_value,
    MiRKAT_R2 = mirkat_r2
  ))
}
