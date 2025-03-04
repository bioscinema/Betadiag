Remedy_Gram <- function(G, method, epsilon = NULL){
  ## input G and implement remedial techniques on it
  ## specify epsilon if use Tikhonov
  if(method == "Higham"){
    ### Nearest PSD projection
    eig <- eigen(G)
    lambda <- eig$values  
    lambda0 <- pmax(lambda, 0)
    G_rem <- eig$vectors %*% diag(lambda0) %*% t(eig$vectors)
  }
  
  if(method == "Tikhonov"){
    ### Relaxed Positive Definiteness
    eig <- eigen(G)
    lambda <- eig$values  
    lambda[abs(lambda)<1e-16] = 0
    lambda.min = min(lambda)
    if(lambda.min==0){
      lambda = lambda + epsilon
    }
    if(lambda.min<0){
      lambda = lambda + abs(lambda.min) + epsilon
    }
    G_rem <- eig$vectors %*% diag(lambda) %*% t(eig$vectors)
  }
  
  return(G_rem)
}