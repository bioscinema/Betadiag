permanova_gower <- function(G, metadata) {
  eig <- eigen(G)
  lambda <- eig$values
  lambda[abs(lambda)<1e-16] = 0
  if(sum(lambda<0)>0){
    stop("The dissimilarity matrix is non-Euclidean!")
  }
  else{
    r = sum(lambda==0) # dimension loss 
    n <- nrow(G)
    X <- model.matrix(~ ., data = metadata)
    p <- ncol(X)
    
    if(n-r<p){
      warning("The between-sample degree-of-freedom is low!")
    }
    
    H <- X %*% solve(t(X) %*% X) %*% t(X)  
    G_H <- H %*% G %*% H  
    SS_model <- sum(diag(G_H))  
      
    I_H <- diag(n) - H  
    G_res <- I_H %*% G %*% I_H  
    SS_residual <- sum(diag(G_res)) 
      
    df_model <- p - 1
    df_residual <- n - r - p
    pseudo_F <- (SS_model / df_model) / (SS_residual / df_residual) 
    p_value <- pf(pseudo_F, df_model, df_residual, 0, lower.tail = F)
  
    return(list(
      pseudo_F = pseudo_F,
      p_value = p_value
    ))
    }
}
  
  
  


  


  

}