## Function to compute an independence matrix based on some regularization matrix and value lambda

getStructureLambda <- function(M_L, lambda) {
  d <- ncol(M_L)
  if (length(lambda) == 1) {
    lambda_Vector <- rep(lambda, d)
  } else {
    lambda_Vector <- c(lambda, -1, -1)
  }
  checkLambda <- function(column, M_L, lambda) {
    M_L[,column] > lambda[column]
  }
  out <- sapply(1:d, FUN = checkLambda, M_L = M_L, lambda = lambda_Vector)
  IndepMatrix <- matrix(as.numeric(out), ncol = d, nrow = d, byrow = FALSE)
  IndepMatrix[upper.tri(IndepMatrix, diag = TRUE)] <- 0
  return(IndepMatrix)
}
