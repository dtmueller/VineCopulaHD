## Function calcStartOrdering20 for selection of an ordering for a R-vine structure using structural equation models and the Lasso

calcStartOrdering <- function(data, order_selection = "Lasso-z", lambda_first = NULL, alpha = 1, cores = 1) {

  if (cores != 1 | is.na(cores)) {
    if (is.na(cores))
      cores <- max(1, detectCores() - 1)
    if (cores > 1) {
      cl <- makeCluster(cores, outfile = "")
      registerDoParallel(cl)
      on.exit(try(stopCluster(), silent = TRUE))
      on.exit(try(closeAllConnections(), silent = TRUE), add = TRUE)
    }
  }

  d <- ncol(data)
  if(any(data < 0) | any(data > 1)) {
    data_u <- pobs(data)
    data_n <- as.data.frame(apply(data_u, 2, qnorm))
  } else {
    data_u <- data
    data_n <- as.data.frame(apply(data_u, 2, qnorm))
  }

  if (order_selection == "cor") {
    order_selection <- "cor-z"
  }

  if (order_selection == "cor-x" | order_selection == "cor-z") {
    if (order_selection == "cor-x") {
      data_n <- data
    } else {
      ## do nothing
    }
    weights <- apply(abs(cor(data_n)), 1, sum)
    infostats <- weights
    new_order <- order(weights, decreasing = TRUE)
    occurrence_table = NA
  } else {
    if (order_selection == "tau") {
      weights <- apply(abs(TauMatrix(data_u)), 1, sum)
      infostats <- weights
      new_order <- order(weights, decreasing = TRUE)
      occurrence_table = NA
    } else {
      if (order_selection == "Lasso-x" | order_selection == "Lasso" | order_selection == "Lasso-z") {
        ## define functions
        if (order_selection == "Lasso-x") {
          data_n <- data
        } else {
          ## do nothing
        }
        applyLassoOrdering <- function(i, dataset, lambda_first = NULL, alpha) {
          if (is.null(lambda_first)) {
            set.seed(501+i)
            model <- cv.glmnet(as.matrix(dataset[,-i]), matrix(dataset[,i]), family = "gaussian", intercept = FALSE, alpha = alpha)
            lambda_first <- model$lambda.1se
          } else {
            model <- glmnet(as.matrix(dataset[,-i]), matrix(dataset[,i]), family = "gaussian", intercept = FALSE, alpha = alpha)
          }
          out <- coef(model, lambda_first)[-1,] # it must be calculated but can not be reused
          result <- list(result = names(out)[which(out != 0)], infostats = lambda_first)
          return(result)
        }
        getPos <- function(x, names_dataset) {
          if(x %in% names_dataset) {
            which(x == names_dataset)
          } else {
            NULL
          }
        }
        getPosV <- Vectorize(getPos, vectorize.args = "x")
        ## define routine
        if (cores > 1) {
          clusterExport(cl, c("glmnet", "cv.glmnet"))
          occurrence <- parLapply(cl, as.list(1:d), fun = applyLassoOrdering, dataset = data_n, lambda_first = lambda_first, alpha = alpha)
          occurrences <- sapply(1:d, FUN = function(i, occurrence) {occurrence[[i]]$result}, occurrence = occurrence)
          occurrence_table <- table(unlist(occurrences))
          eval.list.fun <- function(i, outlist) { outlist[[i]]$infostats }
          infostats <- vapply(1:d, FUN = eval.list.fun, outlist = occurrence, FUN.VALUE = 1)
          weights_intermed <- names(occurrence_table[order(occurrence_table, decreasing = TRUE)])
          weights_intermed <- union(weights_intermed, setdiff(names(data_n), weights_intermed))
          new_order <- unlist(getPosV(weights_intermed, names_dataset = names(data_n)))
          stopCluster(cl)
        } else {
          occurrence <- lapply(as.list(1:d), FUN = applyLassoOrdering, dataset = data_n, lambda_first = lambda_first, alpha = alpha)
          occurrences <- sapply(1:d, FUN = function(i, occurrence) {occurrence[[i]]$result}, occurrence = occurrence)
          occurrence_table <- table(unlist(occurrences))
          eval.list.fun <- function(i, outlist) { outlist[[i]]$infostats }
          infostats <- vapply(1:d, FUN = eval.list.fun, outlist = occurrence, FUN.VALUE = 1)
          weights_intermed <- names(occurrence_table[order(occurrence_table, decreasing = TRUE)])
          weights_intermed <- union(weights_intermed, setdiff(names(data_n), weights_intermed))
          new_order <- unlist(getPosV(weights_intermed, names_dataset = names(data_n)))
        }
      }
    }
  }
  out <- list(new_order = new_order, infostats = infostats, occurrence_table = occurrence_table, occurrences = occurrences)
  return(out)
}
