MiRKAT_R <- function(G,Y,metadata){
  n <- nrow(G)
  X <- model.matrix(~ ., data = metadata)
  Y <- model.matrix(~ ., data = Y)
  H <- X %*% solve(t(X) %*% X) %*% t(X)  
  I_H <- diag(n) - H  
  R <- I_H %*% Y
  eig <- eigen(G)
  Gd <- eig$vectors %*% diag(abs(eig$values)) %*% t(eig$vectors)
  Rsquared <- sum(diag((t(R) %*% G %*% R)))/(sum(diag(Gd))*sum(diag(t(R) %*% R)))
  return(Rsquared)
}