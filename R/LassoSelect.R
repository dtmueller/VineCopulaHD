## Function LassoSelect for selection of an R-vine structure using structural equation models and the Lasso

LassoSelect <- function(data, order_selection = "Lasso", alpha = 1, cores = 1, logging = FALSE, filename) {
  d <- ncol(data)

  ## data conversion
  if(any(data < 0) | any(data > 1)) {
    data_u <- pobs(data)
    data_n <- as.data.frame(apply(data_u, 2, qnorm))
  } else {
    data_u <- data
    data_n <- as.data.frame(apply(data_u, 2, qnorm))
  }

  ## definition of internally used functions

  getFirstCoef <- function(path) {
    c(sum(path != 0), as.numeric(path[length(path) - sum(path != 0) + 1]))
  }
  getPos <- function(x, names_dataset) {
    if(x %in% names_dataset) {
      which(x == names_dataset)
    } else {
      NULL
    }
  }
  getPosV <- Vectorize(getPos, vectorize.args = "x")

  extractFirstCoef <- function(i, coord, coefs) {
    out <- coefs[coord[i,1],][coord[i,2]]
    if(is.na(out)) {
      0
    } else {
      out
    }
  }

  getRegressors1 <- function(j, RVM, dataset, alpha) {
    d <- ncol(RVM)
    regressors <- diag(RVM)[(j+1):d]

    model <- glmnet(as.matrix(dataset[,regressors]), matrix(dataset[,RVM[j,j]]),
                    family = "gaussian", intercept = FALSE, alpha = alpha)

    res <- model$beta != 0
    ncol.res <- ncol(res)
    start <- ncol.res - sapply(1:nrow(res), FUN = function(i, x) {sum(x[i,])}, x = res)

    coords <- cbind(1:nrow(res), start + 1) ## shift here
    first_coefs <- sapply(1:nrow(coords), FUN = extractFirstCoef, coord = coords, coefs = model$beta)

    temp <- start + first_coefs
    ordering <- order(temp, decreasing = FALSE)
    out <- start[ordering]

    lambda_out <- model$lambda[out]
    path_out <- getPosV(rownames(model$beta)[ordering], names_dataset = names(dataset))

    result <- list(regressor = path_out, lambda_out = lambda_out)
    return(result)
  }

  getRegressorsk <- function(j, possible_entries, dataset, alpha, RVM, i) {
    d <- ncol(RVM)
    whitelist <- RVM[(i+1):d,j]
    regressors <- vapply(possible_entries, FUN = function(x) {return(x$entry)}, FUN.VALUE = 1)
    columns  <- vapply(possible_entries, FUN = function(x) {return(x$column)}, FUN.VALUE = 1)
    regressors <- c(regressors, whitelist)
    penalty.factor <- rep(1, length(regressors))
    id_whitelist <- which(regressors %in% whitelist)
    penalty.factor[id_whitelist] <- 0

    model <- glmnet(as.matrix(dataset[,regressors]), matrix(dataset[,RVM[j,j]]),
                    family = "gaussian", intercept = FALSE, alpha = alpha,
                    penalty.factor = penalty.factor)

    res <- model$beta[-id_whitelist,] != 0
    ncol.res <- ncol(res)
    if (is.null(ncol.res)) {
      res <- matrix(res, nrow = 1)
    }
    start <- ncol(res) - sapply(1:nrow(res), FUN = function(i, x) {sum(x[i,])}, x = res)

    coords <- cbind(1:nrow(res), start + 1) ## shift here
    first_coefs <- sapply(1:nrow(coords), FUN = extractFirstCoef, coord = coords, coefs = model$beta)

    temp <- start + first_coefs
    ordering <- order(temp, decreasing = FALSE)
    out <- start[ordering]

    lambda_out <- model$lambda[out]
    path_out <- getPosV(rownames(model$beta)[ordering], names_dataset = names(dataset))

    result <- list(regressors = path_out, lambda_out = lambda_out, column = columns[which(regressors == path_out[[1]])])
    return(result)
  }

  cutResults1 <- function(results) {
    if (!is.null(results)) {
      return(results$regressor[-1])
    } else {
      return(integer(0))
    }
  }

  # cutResults2 <- function(results) {
  #   if (!is.null(results)) {
  #     return(results[-1])
  #   } else {
  #     return(integer(0))
  #   }
  # }
  #
  # calcFirstTreeApply <- function(j, results) {
  #   as.numeric(results[[j]][1])
  # }
  #
  # CheckFun <- function(m, RVS, y, k, d, col.x) {
  #   length(setdiff(union(y, as.vector(RVS[min(d, d-k+2):d, col.x])), union(as.vector(RVS[min(d, d-k+2):d, m]), RVS[m, m])))
  # }
  #
  # assignValues <- function(set, RVM, occurrence_table, dataset) {
  #   c(which(diag(RVM) == set[2,1]), getPos(names(which.max(occurrence_table[unlist(set[1, which(set[3,] == TRUE)])])), names_dataset = names(dataset)))
  # }
  #
  # sortValues <- function(X, RVM, i) {
  #   RVM[i,X[1]] <- X[2]
  # }


  # applyLassoCVEnd <- function(j, RVM, dataset, alpha) {
  #   set.seed(501)
  #   d <- ncol(RVM)
  #   regressors <- RVM[(j+1):d, j]
  #   response <- RVM[j,j]
  #   model <- cv.glmnet(as.matrix(dataset[,regressors]), matrix(dataset[,response]), family = "gaussian", intercept = FALSE, alpha = alpha)
  #   coefs <- apply(model$glmnet.fit$beta, 1, FUN = function(x) {sum(x == 0)})
  #
  #   lambdas <- model$lambda[coefs]
  #   names(lambdas) <- names(coefs)
  #   if (any(is.na(lambdas))) {
  #     lambdas[which(is.na(lambdas))] = 0
  #   }
  #   lambda_cv_1se <- model$lambda.1se
  #   lambda_cv_min <- model$lambda.min
  #   out <- list(lambdas = lambdas, lambda_cv_1se = lambda_cv_1se, lambda_cv_min = lambda_cv_min)
  #   return(out)
  # }

  ## outsourcing of function for determination of ordering of main diagonal
  out <- calcStartOrdering(data, order_selection = order_selection, alpha = alpha, cores = cores)
  new_order <- out$new_order
  occurrence_table <- out$occurrence_table
  occurrences <- out$occurrences
  if (length(occurrence_table) != d) {
    nameset <- setdiff(names(data), names(occurrence_table))
    add <- rep(0, length(data))
    names(add) <- nameset
    occurrence_table <- c(occurrence_table, add)
  }
  occurrence_table <- occurrence_table[unlist(getPosV(names(data), names_dataset = names(occurrence_table)))]
  if (logging == TRUE) {
    filename_save <- paste(paste("output/", paste(filename, "new_order", sep = "_"), sep = ""), ".rdata", sep = "")
    save(new_order, file = filename_save)
  }

  ## register cluster
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

  ## build initial objects
  RVM <- diag(rev(new_order))
  RVM[d, d-1] <- RVM[d,d]
  M_L <- matrix(0, d, d)

  col_matrix <- matrix(0, d, d)

  y <- data_n[,diag(RVM)[d-1]]
  z <- data_n[,diag(RVM)[d]]
  y_new <- y - mean(y)
  z_new <- y - mean(z)
  z_new <- z / (sum(z^2)/length(z))
  M_L[d, d-1] <- (1/length(z)) * abs((z_new %*% y_new))

  ## calculate first paths and tree T_1
  if (cores > 1) {
    clusterExport(cl, "glmnet")
    results <- parLapply(cl, as.list(1:(d-2)), fun = getRegressors1, RVM = RVM, dataset = data_n, alpha = alpha)
  } else {
    results <- lapply(as.list(1:(d-2)), FUN = getRegressors1, RVM = RVM, dataset = data_n, alpha = alpha)
  }

  if (logging == TRUE) {
    filename_save <- paste(paste("output/", paste(filename, "results", sep = "_"), sep = ""), ".rdata", sep = "")
    save(results, file = filename_save)
  }

  ## obtain lambda and parameter estimates
  regressors <- lapply(results, FUN = function(j) {j$regressor})
  lambdas <- sapply(results, FUN = function(j) {j$lambda_out})

  if (any(vapply(regressors, FUN = length, FUN.VALUE = 1) == 0)) {
    NULLs <- which(vapply(regressors, FUN = function(x) {length(x) == 0}, FUN.VALUE = 1) == TRUE)
    notNULLs <- setdiff(1:(d-2), NULLs)
    RVM[d, NULLs] <- diag(RVM)[NULLs + 1]
    RVM[d, notNULLs] <- vapply(regressors[notNULLs], FUN = function(x) {x[1]}, FUN.VALUE = 1)
    M_L[d, notNULLs] <- lambdas
  } else {
    RVM[d,1:(d-2)] <- vapply(regressors, FUN = function(x) {x[1]}, FUN.VALUE = 1)
    M_L[d,1:(d-2)] <- vapply(lambdas, FUN = function(x) {x[1]}, FUN.VALUE = 1)
  }

  col_matrix[d, 1:(d-1)] <- vapply(RVM[d,1:(d-1)], FUN = function(x, diagonal) {which(diagonal == x)}, diagonal = diag(RVM), FUN.VALUE = 1)

  proximity_fails <- matrix(1, d, d)
  proximity_fails[upper.tri(proximity_fails, diag = TRUE)] <- 0

  path_coef <- sapply(1:(d-2), FUN = function(j, results) {results[[j]]$regressor[-1]}, results = results)
  path_lambda <- sapply(1:(d-2), FUN = function(j, results) {results[[j]]$lambda_out[-1]}, results = results)

  ###################################
  symdiff <- function(x, y) { setdiff( c(x, y), intersect(x, y))}

  makeNodes <- function(j, RVM, i) {
    c(RVM[j,j], RVM[d:(i+1),j])
  }

  convertToList <- function(j, possible_entries, columns, y) {
    out <- list(entry = setdiff(possible_entries[[j]], y), column = columns[[j]])
    return(out)
  }

  checkProximityNovel <- function(j, nodes) {
    checked_list <- lapply(nodes[j:length(nodes)], FUN = symdiff, y = nodes[[j]])
    length_list <- unlist(vapply(checked_list, FUN = length, FUN.VALUE = 1))
    length_list_pos <- which(length_list == 2)
    possible_entries <- checked_list[length_list_pos]
    clean_entries <- lapply(as.list(1:length(length_list_pos)), FUN = convertToList, possible_entries = possible_entries, columns = length_list_pos + j-1, y = nodes[[j]])
    return(clean_entries)
  }

  evalPossibleEntries_s <- function(RVM, i) {
    nodes <- lapply(1:i, FUN = makeNodes, RVM = RVM, i = i)
    possible_entries <- lapply(1:(i-1), FUN = checkProximityNovel, nodes = nodes)
    return(possible_entries)
  }

  evalPossibleEntries_p <- function(RVM, i) {
    nodes <- lapply(1:i, FUN = makeNodes, RVM = RVM, i = i)
    possible_entries <- parLapply(cl, 1:(i-1), fun = checkProximityNovel, nodes = nodes)
    return(possible_entries)
  }

  # needForRecalc <- function(j, results, possible_entries) {
  #   if (length(results[[j]]) == 0) {
  #     return(FALSE)
  #   } else {
  #     if (results[[j]][1] %in% possible_entries[[j]]) {
  #       return(TRUE)
  #     } else {
  #       return(FALSE)
  #     }
  #   }
  # }
  #

  wrapperFunction <- function(j, possible_entries, RVM, dataset, alpha, i, path_coef, path_lambda) {
    if (length(path_coef[[j]]) == 0) {
      temp <- getRegressorsk(j, RVM = RVM, possible_entries = possible_entries[[j]], i = i, dataset = dataset, alpha = alpha)
      out <- list(reworked = TRUE, out = temp$regressors, lambda_out = temp$lambda_out, column = temp$column)
    } else {
      possible_entries_j <- sapply(possible_entries[[j]], FUN = function(item) {item$entry})
      path_j <- path_coef[[j]][1]
      out <- path_j %in% possible_entries_j
      if (out) {
        id <- which(path_j %in% possible_entries_j)
        column <- possible_entries[[j]][[id]]$column
        out <- list(reworked = FALSE, out = path_coef[[j]], lambda_out = path_lambda[[j]], column = column)
      } else {
        temp <- getRegressorsk(j, RVM = RVM, possible_entries = possible_entries[[j]], i = i, dataset = dataset, alpha = alpha)
        out <- list(reworked = TRUE, out = temp$regressors, lambda_out = temp$lambda_out, column = temp$column)
      }
    }
    return(out)
  }
  ####################################
  #   isolateResults <- function(j, out) {
  #   if (length(out[[j]]$results) == 0) {
  #     return(NA)
  #   } else {
  #     return(out[[j]]$results[1])
  #   }
  # }

  # loop over higher order trees)
  if (cores > 1) {
    clusterExport(cl, c("glmnet"))
    for(i in (d-1):2) {
      possible_entries <- evalPossibleEntries_p(RVM, i)
      temp_out <- parLapply(cl, as.list(1:(i-1)), fun = wrapperFunction, possible_entries = possible_entries, RVM = RVM, dataset = data_n, alpha = alpha, i = i, path_coef = path_coef, path_lambda = path_lambda)

      RVM[i, 1:(i-1)] <- vapply(1:(i-1), FUN = function (j, out) {out[[j]]$out[1]}, out = temp_out, FUN.VALUE = 1)
      path_coef <- sapply(1:(i-1), FUN = function (j, out) {out[[j]]$out[-1]}, out = temp_out)
      M_L[i, 1:(i-1)] <- vapply(1:(i-1), FUN = function (j, out) {out[[j]]$lambda_out[1]}, out = temp_out, FUN.VALUE = 1)
      path_lambda <- sapply(1:(i-1), FUN = function (j, out) {out[[j]]$lambda_out[-1]}, out = temp_out)

      proximity_fails[i, 1:(i-1)] <- as.numeric(FALSE == vapply(1:(i-1), FUN = function (j, out) {out[[j]]$reworked}, out = temp_out, FUN.VALUE = logical(1)))
      col_matrix[i, 1:(i-1)] <- vapply(1:(i-1), FUN = function (j, out) {out[[j]]$column}, out = temp_out, FUN.VALUE = 1)
      print(i)
    }
  } else {
    for(i in (d-1):2) {
      possible_entries <- evalPossibleEntries_s(RVM, i)
      temp_out <- lapply(as.list(1:(i-1)), FUN = wrapperFunction, possible_entries = possible_entries, RVM = RVM, dataset = data_n, alpha = alpha, i = i, path_coef = path_coef, path_lambda = path_lambda)

      RVM[i, 1:(i-1)] <- vapply(1:(i-1), FUN = function (j, out) {out[[j]]$out[1]}, out = temp_out, FUN.VALUE = 1)
      path_coef <- sapply(1:(i-1), FUN = function (j, out) {out[[j]]$out[-1]}, out = temp_out)
      M_L[i, 1:(i-1)] <- vapply(1:(i-1), FUN = function (j, out) {out[[j]]$lambda_out[1]}, out = temp_out, FUN.VALUE = 1)
      path_lambda <- sapply(1:(i-1), FUN = function (j, out) {out[[j]]$lambda_out[-1]}, out = temp_out)

      proximity_fails[i, 1:(i-1)] <- as.numeric(FALSE == vapply(1:(i-1), FUN = function (j, out) {out[[j]]$reworked}, out = temp_out, FUN.VALUE = logical(1)))
      col_matrix[i, 1:(i-1)] <- vapply(1:(i-1), FUN = function (j, out) {out[[j]]$column}, out = temp_out, FUN.VALUE = 1)
      print(i)
    }
  }

  ### check for remaining -1 entries in the M_L matrix - get coefficients for lambda
  # indexSet <- which(M_L == -1)
  # columns <- ceiling(indexSet / nrow(M_L))
  # rows <- indexSet %% nrow(M_L)
  # entries <- cbind(rows, columns)
  #
  # if (cores > 1) {
  #   clusterExport(cl, c("glmnet"))
  #   M_L[indexSet] <- parApply(cl, entries, 1, FUN = getLambdaValue, dataset = data_n, RVM = RVM, alpha = alpha)
  # } else {
  #   M_L[indexSet] <- apply(entries, 1, FUN = getLambdaValue, dataset = data_n, RVM = RVM, alpha = alpha)
  # }

  # if (cores > 1) {
  #   clusterExport(cl, c("glmnet", "cv.glmnet") )
  #   cv_lambda_end <- parLapply(cl, 1:(d-2), fun = applyLassoCVEnd, RVM = RVM, dataset = data_n, alpha = alpha)
  # } else {
  #   cv_lambda_end <- lapply(1:(d-2), FUN = applyLassoCVEnd, RVM = RVM, dataset = data_n, alpha = alpha)
  # }
  #
  # extractLambdas <- function(row, cv_lambda_end, d) {
  #   c(rep(0, row), cv_lambda_end[[row]]$lambdas)
  # }
  #
  # extractLambdaCV_1se <- function(row, cv_lambda_end) {
  #   cv_lambda_end[[row]]$lambda_cv_1se
  # }
  #
  # extractLambdaCV_min <- function(row, cv_lambda_end) {
  #   cv_lambda_end[[row]]$lambda_cv_min
  # }
  #
  # lambda_vectors <- sapply(1:(d-2), FUN = extractLambdas, cv_lambda_end = cv_lambda_end, d = d)
  # rownames(lambda_vectors) <- NULL
  # y <- data_n[,diag(RVM)[d-1]]
  # z <- data_n[,diag(RVM)[d]]
  # y_new <- y - mean(y)
  # z_new <- y - mean(z)
  # z_new <- z / (sum(z^2)/length(z))
  # lambda_21 <- (1/length(z)) * (z_new %*% y_new)
  # M_L <- cbind(lambda_vectors, c(rep(0, d-1), lambda_21), rep(0, d))
  #
  # lambda_cv_vector_1se <- sapply(1:(d-2), FUN = extractLambdaCV_1se, cv_lambda_end = cv_lambda_end)
  # lambda_cv_vector_min <- sapply(1:(d-2), FUN = extractLambdaCV_min, cv_lambda_end = cv_lambda_end)

  # out <- list(RVM = RVM, M_L = M_L, proximity_fails = proximity_fails, cv_lambda_ordering = out$infostats[rev(new_order)][1:(d-2)], cv_lambda_end_1se = lambda_cv_vector_1se,
  #             cv_lambda_end_min = lambda_cv_vector_min)

  .RVM <- RVM
  .M_L <- M_L
  .proximity_fails <- proximity_fails
  .occurrences <- occurrences[new_order]
  .proximity_matrix <- col_matrix

  remove(list = ls())
  gc()

  out <- list(RVM = .RVM, M_L = .M_L, proximity_fails = .proximity_fails, occurrences = .occurrences, proximity_matrix = .proximity_matrix)
  return(out)
}
