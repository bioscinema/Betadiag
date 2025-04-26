euclidean_check <- function(G) {
  eig <- eigen(G)
  lambda <- eig$values
  lambda[abs(lambda) < 1e-12] <- 0
  neg_sum <- sum(abs(lambda[lambda < 0]))
  total_sum <- sum(abs(lambda))
  ratio <- neg_sum / total_sum
  
  return(list(rate = ratio, is.Euclidean = ifelse(ratio > 0, 0, 1)))
}