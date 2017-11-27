## Function to compute a clustered R-vine based on a graphical model

RVineClusterSelect <- function(data,
                               graph = NULL,
                               max_size = 50,
                               graph_select = NA,
                               nlambda = 30,
                               familyset = NA,
                               selectioncrit = "AIC",
                               indeptest = TRUE,
                               level = 0.05,
                               trunclevel = NA,
                               method = "mle",
                               fill_level = 2,
                               cores = 1,
                               batch_size = NA,
                               par_est = TRUE,
                               method_outside = "itau",
                               loop = TRUE) {

  ## data preprocessing
  if (any(data < 0) | any(data > 1)) {
    data_u <- pobs(data)
  } else {
    data_u <- data
  }
  # data_n <- as.data.frame(apply(data_u, 2, FUN = qnorm))

  print(paste("Graph  information", sum(graph), sep = " "))
  print(paste("Graph  information null", !is.null(graph), sep = " "))

  if (!is.null(graph) & sum(graph) == 0) {
    return("No dependence structure detected.")
  }

  names_dataset <- names(data_u) <- names(data)
  d <- ncol(data_u)
  n <- nrow(data_u)

  ## Definition of internal functions
  calcClusterSize <- function(adj_matrix) {
    max(components(graph_from_adjacency_matrix(adj_matrix, mode = "undirected"))$csize)
  }
  assignList <- function(x, vector) {
    which(vector == x)
  }
  getPos <- function(x, names_dataset) {
    if(x %in% names_dataset) {
      which(x == names_dataset)
    } else {
      NULL
    }
  }
  getPosV <- Vectorize(getPos, vectorize.args = "x")

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

  evalPossibleEntries_s_C <- function(RVM, i) {
    todo <- which(RVM[i,1:(i-1)] == "0")
    nodes <- lapply(1:i, FUN = makeNodes, RVM = RVM, i = i)
    out <- lapply(todo, FUN = checkProximityNovel, nodes = nodes)
    return(out)
  }

  evalPossibleEntries_p_C <- function(RVM, i) {
    todo <- which(RVM[i,1:(i-1)] == "0")
    nodes <- lapply(1:i, FUN = makeNodes, RVM = RVM, i = i)
    out <- parLapply(cl, todo, fun = checkProximityNovel, nodes = nodes)
    return(out)
  }

  evalPossibleEntries_s_Cc <- function(RVM, i) {
    todo <- which(RVM[i,1:(i-1)] == "0")
    nodes <- lapply(1:i, FUN = makeNodes, RVM = RVM, i = i)
    possible_entries <- lapply(todo, FUN = checkProximityNovel, nodes = nodes)
    return(possible_entries)
  }

  evalPossibleEntries_p_Cc <- function(RVM, i) {
    todo <- which(RVM[i,1:(i-1)] == "0")
    nodes <- lapply(1:i, FUN = makeNodes, RVM = RVM, i = i)
    possible_entries <- parLapply(cl, todo, fun = checkProximityNovel, nodes = nodes)
    return(possible_entries)
  }

  fitFunctionGraph <- function(fitSet,
                               data,
                               graph,
                               familyset = NA,
                               selectioncrit = "AIC",
                               indeptest = TRUE,
                               level = 0.05,
                               trunclevel = NA,
                               method = "itau") {
    graph <- graph[fitSet, fitSet] ## assign quadratic matrix

    out <- RVineGraphSelect(data[,fitSet],
                            graph = graph,
                            familyset = familyset,
                            selectioncrit = selectioncrit,
                            indeptest = indeptest,
                            level = level,
                            trunclevel = trunclevel,
                            method = method,
                            cores = 1)
    if (length(fitSet) == 2) {
      out$pair.AIC <- matrix(c(0,out$AIC,0,0), 2, 2)
      out$pair.BIC <- matrix(c(0,out$BIC,0,0), 2, 2)
    }
    out
  }

  getPossibleEntries <- function(j, Matrix) {
    d <- ncol(Matrix)
    columns <- (j+1):d
    nodes <- diag(Matrix)[columns]
    out <- lapply(1:length(columns), FUN = convertToList2, nodes = nodes, columns = columns)
    out
  }

  convertToList2 <- function(iterator, nodes, columns) {
    list(entry = nodes[[iterator]], columns = columns[[iterator]])
  }

  evalPossibleEntries_s_C_1 <- function(Matrix, i) {
    out <- as.list(rep(NA, i-1))
    out[which(Matrix[i,] == "0")] <- lapply(which(Matrix[i,] == "0"), FUN = getPossibleEntries, Matrix = Matrix)
    out
  }

  evalPossibleEntries_p_C_1 <- function(Matrix, i) {
    out <- as.list(rep(NA, i-1))
    out[which(Matrix[i,] == "0")] <- parLapply(cl, which(Matrix[i,] == "0"), fun = getPossibleEntries, Matrix = Matrix)
    out
  }

  # estPC is the main function called by parLapply over an entire tree level
  estPC_isolated <- function(j, i, V, diag_M,
                             familyset = NA,
                             selectioncrit = "AIC",
                             indeptest = TRUE,
                             level = 0.05,
                             method = "itau") {
    zr1 <- V$direct[[i]]

    if (j$entry == diag_M[j$column]) {
      # this is the direct case
      zr2 <- V$direct[[j$column]]
    } else {
      zr2 <- V$indirect[[j$column]]
    }
    cfit <- suppressWarnings(BiCopSelect(zr2,
                                         zr1,
                                         familyset,
                                         selectioncrit,
                                         indeptest,
                                         level,
                                         method = method))

    direct <- suppressWarnings(BiCopHfunc1(zr2,
                                           zr1,
                                           cfit,
                                           check.pars = FALSE))

    indirect <- suppressWarnings(BiCopHfunc2(zr2,
                                             zr1,
                                             cfit,
                                             check.pars = FALSE))
    out <- list(cfit = cfit,
                direct = direct,
                indirect = indirect)
    out
  }


  estPC <- function(i, open_entries, pos_entries, V, diag_M,
                    familyset = NA,
                    selectioncrit = "AIC",
                    indeptest = TRUE,
                    level = 0.05,
                    method = "itau") {

    temp <- lapply(pos_entries[[i]], FUN = estPC_isolated, open_entries[i], V, diag_M, familyset, selectioncrit, indeptest, level, method)
    opt <- which.min(sapply(temp, FUN = function(x) {x$cfit$AIC}))
    out <- list(id = pos_entries[[i]][[opt]]$entry,
                PC = temp[[opt]]$cfit,
                direct = temp[[opt]]$direct,
                indirect = temp[[opt]]$indirect)
    out
  }

  # works only for the first tree
  generatePseudoObs_1 <- function(j, V, Matrix, family, par, par2) {
    d <- ncol(Matrix)
    ell <- which(diag(Matrix) == Matrix[d,j])
    out.direct <- BiCopHfunc1(V$direct[[ell]],
                              V$direct[[j]],
                              family = family[d, j],
                              par = par[d, j],
                              par2 = par2[d, j])
    out.indirect <- BiCopHfunc2(V$direct[[ell]],
                                V$direct[[j]],
                                family = family[d, j],
                                par = par[d, j],
                                par2 = par2[d, j])
    out <- list(direct = out.direct,
                indirect = out.indirect)
    out
  }

  generatePseudoObs <- function(j, V, k, Matrix, MaxMat, family, par, par2) {
    ell <- which(diag(Matrix) == MaxMat[k,j])
    if (Matrix[k,j] == MaxMat[k,j]) {
      out.direct <- BiCopHfunc1(V$direct[[ell]],
                                V$direct[[j]],
                                family = family[k, j],
                                par = par[k, j],
                                par2 = par2[k, j])
      out.indirect <- BiCopHfunc2(V$direct[[ell]],
                                  V$direct[[j]],
                                  family = family[k, j],
                                  par = par[k, j],
                                  par2 = par2[k, j])
    } else {
      out.direct <- BiCopHfunc1(V$indirect[[ell]],
                                V$direct[[j]],
                                family = family[k, j],
                                par = par[k, j],
                                par2 = par2[k, j])
      out.indirect <- BiCopHfunc2(V$indirect[[ell]],
                                  V$direct[[j]],
                                  family = family[k, j],
                                  par = par[k, j],
                                  par2 = par2[k, j])
    }
    out <- list(direct = out.direct,
                indirect = out.indirect)
    out
  }

  print(paste("Function started at", Sys.time(), sep = " "))

  ## using the nonparanormal, we can operate on the u data
  if (is.null(graph)) {
    data_NPN <- huge.npn(data_u)
    out_path <- huge(data_NPN, nlambda = nlambda)
    if (!is.na(graph_select)) {
      graph <- huge.select(out_path, criterion = graph_select)
    } else {
      if (is.numeric(max_size)) {
        clusterSize <- sum(sapply(out_path$path, FUN = calcClusterSize) <= max_size)
        graph <- out_path$path[[clusterSize]]
      }
    }
  } else {
    ## graph is supplied
  }

  comp_info <- components(graph_from_adjacency_matrix(graph, mode = "undirected"))
  ## generate assignment list of
  component_list <- as.list(1:max(comp_info$membership))

  print(paste("Component list is of length", length(component_list), sep = " "))

  if (length(component_list) == 1) {
    out <- RVineGraphSelect(data = data_u, graph = graph, familyset = familyset, selectioncrit = selectioncrit, indeptest = indeptest, level = level, trunclevel = trunclevel,
                            method = method, cores = cores)
    ret <- list(RVM = out, graph = graph, time.taken = 0)
    return(ret)
  }

  component_assignment <- lapply(component_list, FUN = assignList, vector = comp_info$membership)
  sizes <- sapply(component_assignment, FUN = length)
  component_assignment <- component_assignment[order(sizes, decreasing = TRUE)]
  singleton_ids <- which(sapply(component_assignment, FUN = length) == 1)
  singletons <- unlist(component_assignment[singleton_ids])
  if (length(singleton_ids > 0)) {
    component_assignment <- component_assignment[-singleton_ids]
  }

  cl <- NULL
  if (cores > 1) {
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
  }

  if (!is.null(cl)) {
    fits <- parLapply(cl, component_assignment, fun = fitFunctionGraph, data = data_u,
                      graph = graph, familyset = familyset, selectioncrit = selectioncrit, indeptest = indeptest,
                      level = level, trunclevel = trunclevel, method = method)
  } else {
    fits <- lapply(component_assignment, FUN = fitFunctionGraph, data = data_u,
                   graph = graph, familyset = familyset, selectioncrit = selectioncrit, indeptest = indeptest,
                   level = level, trunclevel = trunclevel, method = method)
  }

  ## initialization of matrix objects
  Matrix <- matrix(0, d, d)
  MaxMat <- matrix(0, d, d)
  family <- matrix(0, d, d)
  par <- matrix(0, d, d)
  par2 <- matrix(0, d, d)
  pair.logLik <- matrix(0, d, d)
  pair.AIC <- matrix(0, d, d)
  pair.BIC <- matrix(0, d, d)
  emptau <- matrix(0, d, d)
  taildep.upper <- matrix(0, d, d)
  taildep.lower <- matrix(0, d, d)
  beta <- matrix(0, d, d)

  if (length(fits) > 0) {
    for(i in 1:length(fits)) {
      temp <- fits[[i]]
      size <- ncol(temp$Matrix)
      ## name exchange part
      RVM_vector <- as.vector(temp$Matrix)
      MaxMat_vector <- as.vector(temp$MaxMat)
      nnset <- which(RVM_vector > 0)
      RVM_vector[nnset] <- temp$names[RVM_vector[nnset]]
      MaxMat_vector[nnset] <- temp$names[MaxMat_vector[nnset]]
      RVM_matrix <- matrix(RVM_vector, length(temp$names))
      MaxMat_matrix <- matrix(MaxMat_vector, length(temp$names))
      ## finish name exchange part
      end <- sum(diag(Matrix) == 0)
      start <- end - size + 1
      diag(Matrix)[start:end] <- diag(RVM_matrix)
      lower_bound <- d
      upper_bound <- d - size + 2 # manual correction
      if (end < d) {
        diag(RVM_matrix) <- "0"
      }
      Matrix[upper_bound:lower_bound, start:(end-1)] <- RVM_matrix[-1, -size] # remove first row
      MaxMat[upper_bound:lower_bound, start:(end-1)] <- MaxMat_matrix[-1, -size] # remove first row
      family[upper_bound:lower_bound, start:end] <- temp$family[-1, ] # remove first row
      par[upper_bound:lower_bound, start:end] <- temp$par[-1, ] # remove first row
      par2[upper_bound:lower_bound, start:end] <- temp$par2[-1, ] # remove first row
      pair.logLik[upper_bound:lower_bound, start:end] <- temp$pair.logLik[-1, ] # remove first row
      ## new
      pair.AIC[upper_bound:lower_bound, start:end] <- temp$pair.AIC[-1, ] # remove first row
      pair.BIC[upper_bound:lower_bound, start:end] <- temp$pair.BIC[-1, ] # remove first row
      emptau[upper_bound:lower_bound, start:end] <- temp$emptau[-1, ] # remove first row
      taildep.upper[upper_bound:lower_bound, start:end] <- temp$taildep$upper[-1, ] # remove first row
      taildep.lower[upper_bound:lower_bound, start:end] <- temp$taildep$lower[-1, ] # remove first row
      beta[upper_bound:lower_bound, start:end] <- temp$beta[-1, ] # remove first row
    }
  }

  print(paste("Fits done at", Sys.time(), sep = " "))
  remove(fits)

  ## now, only singletons are remaining, check for them
  if (!is.null(singletons)) {
    diag(Matrix)[1:(start-1)] <- names_dataset[singletons]
  }

  ## whenever there are two clusters, the last entry of the main diagonal of the smaller cluster has no entry on the d-th row
  ## more sensible to parallelize longer tasks
  if (fill_level > 0) {
    diag_M <- diag(Matrix)
    V <- list()
    V$direct <- as.list(rep(NA, d))
    V$indirect <- as.list(rep(NA, d))
    V$direct <- lapply(1:d, FUN = function(i, data, diag_M) {data[,diag_M[i]]}, data = data_u, diag_M = diag_M)

    for (k in d:(d - fill_level + 1)) {
      open_entries <- which(Matrix[k,1:(k-1)] == "0")
      rem_entries <- setdiff(1:(k-1), open_entries)
      if (k == d) {
        if (!is.null(cl)) {
          pos_entries <- evalPossibleEntries_p_C_1(Matrix, d)[open_entries]
          if (!is.na(batch_size) & par_est == TRUE) {
            batch_assignment <- ceiling((1:length(open_entries))/batch_size)
            estimates <- list()
            print(paste("Open entries length is ", length(open_entries), sep = ""))
            for(i in unique(batch_assignment)) {
              worklist <- which(batch_assignment == i)
              estimates[worklist] <- parLapply(cl, (1:length(open_entries))[worklist], fun = estPC, open_entries, pos_entries, V, diag_M, familyset, selectioncrit, indeptest, level, method_outside)
              print(i)
            }
          } else {
            if (par_est == FALSE & loop == TRUE) {
              estimates <- list()
              print(paste("Open entries length is ", length(open_entries), sep = ""))
              for(i in 1:length(open_entries)) {
                estimates[[i]] <- estPC((1:length(open_entries))[[i]], open_entries, pos_entries, V, diag_M, familyset, selectioncrit, indeptest, level, method_outside)
                print(i)
              }
            } else {
              if (par_est == FALSE & loop == FALSE) {
                estimates <- lapply(1:length(open_entries), FUN = estPC, open_entries, pos_entries, V, diag_M, familyset, selectioncrit, indeptest, level, method_outside)
              } else {
                estimates <- parLapply(cl, 1:length(open_entries), fun = estPC, open_entries, pos_entries, V, diag_M, familyset, selectioncrit, indeptest, level, method_outside)
              }
            }
          }
          print("Estimates 1 okay")
          pseudoObs_rem <- lapply(rem_entries, FUN = generatePseudoObs_1, V = V, Matrix = Matrix, family = family, par = par, par2 = par2)
        } else {
          pos_entries <- evalPossibleEntries_s_C_1(Matrix, d)[open_entries]
          estimates <- lapply(1:length(open_entries), FUN = estPC, open_entries, pos_entries, V, diag_M, familyset, selectioncrit, indeptest, level, method_outside)
          pseudoObs_rem <- lapply(rem_entries, FUN = generatePseudoObs_1, V = V, Matrix = Matrix, family = family, par = par, par2 = par2)
        }
      } else {
        if (!is.null(cl)) {
          pos_entries <- evalPossibleEntries_p_C(Matrix, k)
          if (!is.na(batch_size) & (par_est == TRUE)) {
            batch_assignment <- ceiling((1:length(open_entries))/batch_size)
            estimates <- list()
            print(paste("Open entries length is", length(open_entries), sep = ""))
            for(i in unique(batch_assignment)) {
              worklist <- which(batch_assignment == i)
              estimates[worklist] <- parLapply(cl, (1:length(open_entries))[worklist], fun = estPC, open_entries, pos_entries, V, diag_M, familyset, selectioncrit, indeptest, level, method_outside)
              print(i)
            }
          } else {
            if (par_est == FALSE) {
              estimates <- lapply(1:length(open_entries), FUN = estPC, open_entries, pos_entries, V, diag_M, familyset, selectioncrit, indeptest, level, method_outside)
            } else {
              estimates <- parLapply(cl, 1:length(open_entries), fun = estPC, open_entries, pos_entries, V, diag_M, familyset, selectioncrit, indeptest, level, method_outside)
            }
          }
          print(paste("Estimates okay, level:", k, sep=""))
          pseudoObs_rem <- lapply(rem_entries, FUN = generatePseudoObs, k = k, V = V, Matrix = Matrix, MaxMat = MaxMat, family = family, par = par, par2 = par2)
        } else {
          pos_entries <- evalPossibleEntries_s_C(Matrix, k)
          estimates <- lapply(1:length(open_entries), FUN = estPC, open_entries, pos_entries, V, diag_M, familyset, selectioncrit, indeptest, level, method_outside)
          if (k > (d-fill_level+1)) {
            pseudoObs_rem <- lapply(rem_entries, FUN = generatePseudoObs, k = k, V = V, Matrix = Matrix, MaxMat = MaxMat, family = family, par = par, par2 = par2)
          }
        }
      }
      Matrix[k, open_entries] <- sapply(estimates, FUN = function(x) {x$id})
      family[k, open_entries] <- sapply(estimates, FUN = function(x) {x$PC$family})
      par[k, open_entries]  <- sapply(estimates, FUN = function(x) {x$PC$par})
      par2[k, open_entries] <- sapply(estimates, FUN = function(x) {x$PC$par2})
      pair.logLik[k, open_entries] <- sapply(estimates, FUN = function(x) {x$PC$logLik})
      pair.AIC[k, open_entries]  <- sapply(estimates, FUN = function(x) {x$PC$AIC})
      pair.BIC[k, open_entries] <- sapply(estimates, FUN = function(x) {x$PC$BIC})
      emptau[k, open_entries] <- sapply(estimates, FUN = function(x) {x$PC$emptau})
      taildep.upper[k, open_entries] <- sapply(estimates, FUN = function(x) {x$PC$taildep$upper})
      taildep.lower[k, open_entries] <- sapply(estimates, FUN = function(x) {x$PC$taildep$lower})
      print(paste("Fit done of fill_level", k, sep = ":"))
      ## generate pseudo observations

      if (k > (d-fill_level+1)) {
        V$direct <- as.list(rep(NA, k-1))
        V$indirect <- as.list(rep(NA, k-1))
        V$direct[open_entries] <- lapply(estimates, FUN = function(x) {x$direct})
        V$indirect[open_entries] <- lapply(estimates, FUN = function(x) {x$indirect})
        V$direct[rem_entries] <- lapply(pseudoObs_rem, FUN = function(x) {x$direct})
        V$indirect[rem_entries] <- lapply(pseudoObs_rem, FUN = function(x) {x$indirect})
        print(paste("V done of fill_level", k, sep = ":"))
      }
    }
  }

  Matrix_vector <- as.vector(Matrix)
  nnset <- which(Matrix_vector != "0")
  Matrix_vector[nnset] <- getPosV(Matrix_vector[nnset], names_dataset)
  Matrix <- matrix(as.numeric(Matrix_vector), d, d)

  print(paste("RVM conversion completed at", Sys.time(), sep = " "))

  dump <- list(Matrix = Matrix,
               family = family,
               par = par,
               par2 = par2,
               pair.logLik = pair.logLik,
               pair.AIC = pair.AIC,
               pair.BIC = pair.BIC,
               emptau = emptau,
               taildep.upper = taildep.upper,
               taildep.lower = taildep.lower,
               n = n,
               d = d)
  filename <- paste("output/",paste(sum(graph), ".rdata", sep = ""), sep ="")
  # save(dump, file = filename)
  # print(paste("Dump created at filename", filename, sep = ":"))

  ## complete Rvine matrix according to proximity condition and convert to sparse objects

  start.time <- Sys.time()
  # if (!is.null(cl)) {
  #   for(i in (d-fill_level):2) {
  #     completion_values <- evalPossibleEntries_p_Cc(Matrix, i)
  #     fill_in <- sapply(completion_values, FUN = function(x) {x[[1]]$entry})
  #     Matrix[i, which(Matrix[i, 1:(i-1)] == 0)] <- fill_in
  #     print(paste(paste("Row", i, sep = " "), Sys.time(), sep = " "))
  #   }
  # } else {
  #   for(i in (d-fill_level):2) {
  #     completion_values <- evalPossibleEntries_s_Cc(Matrix, i)
  #     fill_in <- sapply(completion_values, FUN = function(x) {x[[1]]$entry})
  #     Matrix[i, which(Matrix[i, 1:(i-1)] == 0)] <- fill_in
  #     print(paste(paste("Row", i, sep = " "), Sys.time(), sep = " "))
  #   }
  # }
  Matrix <- finalizeMatrixOuter(Matrix, d-fill_level)
  end.time <- Sys.time()
  time.taken <- as.numeric(difftime(end.time, start.time, units = "secs"))
  print(paste("RVM completion done at", Sys.time(), sep = " "))

  k <- sum(par != 0) + sum(par2 != 0)
  logLik <- sum(pair.logLik)
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + log(n) * k
  print(paste("Model stats done at", Sys.time(), sep = " "))
  # out_save <- list(Matrix = Matrix,
  #                  family = family,
  #                  par = par,
  #                  par2 = par2,
  #                  pair.logLik = pair.logLik,
  #                  logLik = logLik,
  #                  AIC = AIC,
  #                  BIC = BIC,
  #                  n = n,
  #                  names_dataset = names_dataset)
  # a <- runif(1,0,1)
  # print(a)
  # save(out_save, file = paste(a, ".rdata", sep = ""))

  out <- RVineMatrixS(Matrix = Matrix, family = family, par = par, par2 = par2, names = names_dataset, sparse = TRUE)
  out$pair.logLik <- pair.logLik
  out$logLik <- logLik
  out$AIC <- AIC
  out$BIC <- BIC
  out$nobs <- n
  out$pair.AIC <- pair.AIC
  out$pair.BIC <- pair.BIC
  out$emptau <- emptau
  out$taildep <- list(upper = taildep.upper, lower = taildep.lower)
  out$beta <- beta
  out$names <- names_dataset
  print(paste("Out object done at", Sys.time(), sep = " "))

  .graph <- graph
  .out <- out
  .time.taken <- time.taken

  remove(list = ls())
  ret <- list(RVM = .out, graph = .graph, time.taken = .time.taken)
  return(ret)
}
