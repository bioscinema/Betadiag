#' Check Triangle Inequality for Dissimilarity Matrix
#'
#' This function checks the triangle inequality for a given dissimilarity matrix.
#'
#' @param D A dissimilarity matrix.
#' @return A list containing whether the matrix is metric, the count of negative and zero triangle inequalities, and collinearity and nonlinearity scores.
#' @export
Triangle.Check <- function(D){
  # input dissimilarity matrix
  R = ncol(D)
  neg.counter = 0
  zero.counter = 0
  collinearity.score = rep(0,R)
  nonlinearity.score = rep(0,R)
  
  for (p in 1:R) {
    for (q in 1:R){
      if(q!=p){
        r = seq(1,R)[-c(p,q)]
        numN = sum(D[p,q] + D[q,r] - D[p,r]<0)
        num0 = sum(D[p,q] + D[q,r] - D[p,r]==0)
        neg.counter = neg.counter + numN/(R*(R-1)*(R-2))
        zero.counter = zero.counter + num0/(R*(R-1)*(R-2))
        collinearity.score[p] = collinearity.score[p] + num0/(R*(R-1)*(R-2))
        nonlinearity.score[p] = nonlinearity.score[p] + numN/(R*(R-1)*(R-2))
        }
    }
  }
  
  return(list(is.metric = ifelse(neg.counter>0,0,1), 
              neg.counter = neg.counter,
              zero.counter = zero.counter,
              collinearity.score = collinearity.score,
              nonlinearity.score = nonlinearity.score))
}

