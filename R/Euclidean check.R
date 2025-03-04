Euclidean.Check <- function(G){
  # input Gram matrix
  eig <- eigen(G)
  lambda <- eig$values  
  small.values = which(abs(lambda)<1e-16)
  lambda[small.values] = 0
  neg_sum <- sum(abs(lambda[lambda < 0]))
  total_sum <- sum(abs(lambda))
  ratio <- neg_sum / total_sum
  
  # save results
  return(list(rate = ratio, is.Euclidean = ifelse(ratio>0,0,1)))
}