## Function to compute an independence matrix based on some regularization matrix and value lambda

getStructureLambdaK <- function(M_L, k_percentage = NULL, k_bound = NULL, type = "floor") {

  getColContrib <- function(j, M_L, k_bound) {
    col_sum <- sum(M_L[,j])
    threshold <- k_bound * col_sum
    ## order in size
    set <- (j+1):nrow(M_L)
    ordering <- order(M_L[set,j], decreasing = TRUE)
    out <- M_L[set,j][ordering]
    ##
    calcToI <- function(entry, ordered_vector) {
      sum(ordered_vector[1:entry])
    }
    colsums_partial <- sapply(1:length(set), FUN = calcToI, ordered_vector = out)
    outset <- ordering[which(colsums_partial > threshold)]
    if (j == nrow(M_L)-1) {
      out <- c(rep(0,j), 1)
    } else {
      result <- rep(1, length(set))
      result[outset] <- 0
      out <- c(rep(0,j), result)
    }
    return(out)
  }

  getLowestColumn <- function(j, M_L, k_percentage, type = "floor") {
    set <- (j+1):nrow(M_L)
    evalBoundary <- function(x, type) {
      switch(type,
             floor = floor(x),
             ceiling = ceiling(x))
    }
    k <- evalBoundary(length(set) * k_percentage, type = type)
    result <- c(rep(0, (d - length(set))), rep(1, length(set)))
    if (k == 0) {
      return(result)
    } else {
      out <- j + order(M_L[set,j], decreasing = FALSE)[1:k]
      result[out] <- 0
    }
    return(result)
  }

  d <- nrow(M_L)
  if (!is.null(k_percentage)) {
    indep_matrix_incomplete <- sapply(1:(d-1), FUN = getLowestColumn, M_L = M_L, k_percentage = k_percentage, type = type)
  } else {
    indep_matrix_incomplete <- sapply(1:(d-1), FUN = getColContrib, M_L = M_L, k_bound = k_bound)
  }
  indep_matrix <- cbind(indep_matrix_incomplete, rep(0, d))
  return(indep_matrix)
}
