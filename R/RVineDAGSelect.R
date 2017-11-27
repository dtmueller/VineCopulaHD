## Function taken from VineCopula package and adapted
## Function to estimate a sparse R-vine based on a DAG model

RVineDAGSelect <- function(data, familyset = NA, type = 0, selectioncrit = "AIC", indeptest = FALSE,
                           level = 0.05, trunclevel = NA, progress = FALSE,  weights = NA,
                           se = FALSE, rotations = TRUE, method = "mle",
                           maxparents = Inf, cores = 1) {
  ## preprocessing of arguments
  args <- preproc(c(as.list(environment()), call = match.call()),
                  check_data,
                  check_nobs,
                  check_if_01,
                  prep_familyset,
                  check_twoparams)
  list2env(args, environment())

  d <- ncol(data)
  n <- nrow(data)

  ## sanity checks
  if (d < 2)
    stop("Dimension has to be at least 2.")
  if (d == 2) {
    return(RVineCopSelect(data,
                          familyset = familyset,
                          Matrix = matrix(c(2, 1, 0, 1), 2, 2),
                          selectioncrit = selectioncrit,
                          indeptest = indeptest,
                          level = level,
                          trunclevel = trunclevel,
                          se = se,
                          rotations = rotations,
                          cores = cores))
  }
  if (!(selectioncrit %in% c("AIC", "BIC", "logLik")))
    stop("Selection criterion not implemented.")
  if (level < 0 & level > 1)
    stop("Significance level has to be between 0 and 1.")

  ## set defaults
  if (type == 0)
    type <- "RVine"
  if (type == 1)
    type <- "CVine"
  if (is.na(trunclevel))
    trunclevel <- d
  if (trunclevel == 0)
    familyset <- 0

  ## initialize object for results
  RVine <- list(Tree = NULL, Graph = NULL)
  warn <- NULL

  ## estimation in first tree ----------------------------
  # find optimal tree

  ## register parallel backend
  cl <- NULL
  if (cores != 1 | is.na(cores)) {
    if (is.na(cores))
      cores <- max(1, detectCores() - 1)
    if (cores > 1) {
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      on.exit(try(stopCluster(), silent = TRUE))
      on.exit(try(closeAllConnections(), silent = TRUE), add = TRUE)
    }
  }

  out <- initializeFirstGraph_DAG(data, maxparents, cl)
  .dags <- out$dags
  .h <- out$h
  g <- out$graph
  b <- min(setdiff(g$E$weights, 0)) / 2
  MST <- findMaxTree_DAG(g, mode = type, truncated = trunclevel < 1,  kdag = out$kdag, amat = out$amat, b)

  # estimate pair-copulas
  VineTree <- fit.FirstTreeCopulas_DAG(MST,
                                       data,
                                       type = familyset,
                                       selectioncrit,
                                       indeptest,
                                       level,
                                       se = se,
                                       weights = weights,
                                       method = method,
                                       cores = cores,
                                       DAG = out$kdag)
  # store results
  if (!is.null(VineTree$warn))
    warn <- VineTree$warn
  RVine$Tree[[1]] <- VineTree
  RVine$Graph[[1]] <- g
  oldVineGraph <- VineTree

  ## estimation in higher trees --------------------------
  for (tree in 2:(d - 1)) {
    ## old pseudo-observations are uneccessary in RVine object
    RVine$Tree[[tree - 1]]$E$Copula.CondData.1 <- NULL
    RVine$Tree[[tree - 1]]$E$Copula.CondData.2 <- NULL

    # only estimate pair-copulas if not truncated
    if (trunclevel == tree - 1)
      familyset <- 0
    # find optimal tree
    g <- buildNextGraph_DAG(VineTree, weights, b, cores > 1,
                            truncated = trunclevel < tree, out$kdag)
    MST <- findMaxTree_DAG(g, mode = type, truncated = trunclevel < tree, kdag = out$kdag, amat = out$amat, b)
    # estimate pair-copulas
    VineTree <- fit.TreeCopulas_DAG(MST,
                                    VineTree,
                                    type = familyset,
                                    selectioncrit,
                                    indeptest,
                                    level,
                                    se = se,
                                    progress,
                                    weights = weights,
                                    method = method,
                                    cores = cores,
                                    DAG = out$kdag)
    # store results
    if (!is.null(VineTree$warn))
      warn <- VineTree$warn
    RVine$Tree[[tree]] <- VineTree
    RVine$Graph[[tree]] <- g
  }
  if (any(is.na(data)))
    warning(" In ", args$call[1], ": ",
            "Some of the data are NA. ",
            "Only pairwise complete observations are used.",
            call. = FALSE)
  if (!is.null(warn))
    warning(" In ", args$call[1], ": ",
            warn,
            call. = FALSE)


  ## free memory and return results as 'RVineMatrix' object
  callexp <- match.call()
  .RVM <- as.RVM2(RVine, data, callexp)
  rm(list = ls())
  out <- list(RVM = .RVM, kdags = .dags, h = .h)
  out
}

wrapperDAG <- function(i,
                       indep,
                       pc.data,
                       type,
                       copulaSelectionBy,
                       testForIndependence,
                       testForIndependence.level,
                       se,
                       weights,
                       method) {
  if (!indep[i]) {
    pcSelect(pc.data[[i]],
             type,
             copulaSelectionBy,
             testForIndependence,
             testForIndependence.level,
             se,
             weights,
             method)
  } else {
    pcSelect(pc.data[[i]],
             type = 0,
             copulaSelectionBy,
             testForIndependence,
             testForIndependence.level,
             se,
             weights,
             method)
  }
}

# symdiff <- function(x, y) {setdiff(c(x, y), intersect(x, y))}

checkSepFirst_DAG <- function(pair, DAG) {
  if (dsep(DAG, names(DAG$nodes)[pair[1]], names(DAG$nodes)[pair[2]])) {
    TRUE
  } else {
    FALSE
  }
}

checkSepHigher_DAG <- function(j, k, D, DAG) {
  if (dsep(DAG, names(DAG$nodes)[j], names(DAG$nodes)[k], names(DAG$nodes)[D])) {
    TRUE
  } else {
    FALSE
  }
}

check.Indep.2 <- function(constraintSet, DAG){
  dsep(DAG, names(DAG$nodes)[constraintSet[1]], names(DAG$nodes)[constraintSet[2]], names(DAG$nodes)[constraintSet[3:length(constraintSet)]])
}

calc_amat <- function(dag, data.names){
  h.amat <- matrix(0, length(data.names), length(data.names))
  rownames(h.amat) <- colnames(h.amat) <- data.names
  h.amat[dag$arcs] <- 1
  h.amat
}

initializeFirstGraph_DAG <- function(data, maxparents, cl = NULL) {
  if (maxparents == Inf) {
    maxparents <- ncol(data) - 1
  }

  data <- as.data.frame(apply(data, 2, FUN = qnorm))
  if (!is.null(cl)){
    dags <- foreach(i = 1:maxparents, .packages = c("bnlearn")) %dopar% {
      hc(data, maxp = i)
    }
  } else {
    calc_hc <- function(i) {
      hc(data, maxp = i)
    }
    dags <- lapply(1:maxparents, calc_hc)
  }
  amats <- lapply(dags, calc_amat, data.names = names(data))
  h.amat <- Reduce('+', amats)
  h.amat <- h.amat + t(h.amat)
  h.amat <- h.amat / (max(h.amat) + 1) # rescale
  diag(h.amat) <- Inf
  rownames(h.amat) <- colnames(h.amat) <- NULL

  W <- diag(ncol(data))
  W[upper.tri(W, diag = FALSE)] <- h.amat[upper.tri(h.amat, diag = FALSE)]
  colnames(W) <- rownames(W) <- colnames(data)

  ## return full graph with edge weights
  graph <- graphFromWeightMatrix(W)
  out <- list(graph = graph, kdag = dags[[maxparents]], dags = dags, amat = W, h = h.amat)
  out
}

findMaxTree_DAG <- function(g, mode = "RVine", truncated = FALSE, kdag = NULL, amat = NULL, b) {
  A <- adjacencyMatrix(g)
  if (truncated == FALSE) {
    ## construct adjency matrix
    if (length(kdag$nodes) > length(g$V$names)) {
      ## set weights for matrix - mind the lower diagonal
      constraintSet <- mapply(c, g$E$conditionedSet, g$E$conditioningSet, SIMPLIFY = FALSE)
      dSeptest <- which(sapply(constraintSet, FUN = check.Indep.2, DAG = kdag))

      temp <- -amat[matrix(unlist(g$E$conditionedSet), ncol = 2, byrow = TRUE)] # switch the sign - here, the eligible edge weights are set according to the DAG
      temp[dSeptest] <- -b # manual change for d-separated edges

      # A <- matrix(Inf, length(g$V$names), length(g$V$names))
      A[matrix(unlist(g$E$nums)[,c(1,2)], ncol = 2)] <- temp
      A[matrix(unlist(g$E$nums)[,c(2,1)], ncol = 2)] <- temp
      # A <- -A
    }
    d <- ncol(A)

    if (mode == "RVine") {
      ## initialize
      tree <- NULL
      edges <- matrix(NA, d - 1, 2)
      w <- numeric(d - 1)
      i <- 1

      ## construct minimum spanning tree
      for (k in 1:(d - 1)) {
        # add selected edge to tree
        tree <- c(tree, i)

        m <- apply(as.matrix(A[, tree]), 2, min)
        a <- apply(as.matrix(A[, tree]), 2,
                   function(x) order(rank(x)))[1, ]
        b <- order(rank(m))[1]
        j <- tree[b]
        i <- a[b]

        # store edge and weight
        edges[k, ] <- c(j, i)
        w[k] <- A[i, j]

        ## adjust adjecency matrix to prevent loops
        for (t in tree)
          A[i, t] <- A[t, i] <- Inf
      }

      ## reorder edges for backwads compatibility with igraph output
      edges <- t(apply(edges, 1, function(x) sort(x)))
      edges <- edges[order(edges[, 2], edges[, 1]), ]

      ## delete unused edges from graph
      E <- g$E$nums
      in.tree <- apply(matrix(edges, ncol = 2), 1,
                       function(x) which((x[1] == E[, 1]) & (x[2] == E[, 2])))
      MST <- g
      g$E$todel <- rep(TRUE, nrow(E))
      if (any(g$E$todel)) {
        g$E$todel[in.tree] <- FALSE
        MST <- deleteEdges(g)
      }
    } else if (mode  == "CVine") {
      ## set root as vertex with minimal sum of weights
      A <- adjacencyMatrix(g)
      diag(A) <- 0
      sumtaus <- rowSums(A)
      root <- which.min(sumtaus)

      ## delete unused edges
      g$E$todel <- !((g$E$nums[, 2] == root) | (g$E$nums[, 1] == root))
      MST <- g
      if (any(g$E$todel ))
        MST <- deleteEdges(g)
    } else {
      stop("vine not implemented")
    }
  } else {
    MST <- g

    # get edge list
    edgesList <- g$E$nums
    uid <- sort(unique(as.vector(g$E$nums)))
    luid <- length(uid)

    if (luid > 2) {
      # transform to adjacency list
      adjacencyList <- lapply(uid, function(u)
        c(edgesList[edgesList[,1] == u,2],
          edgesList[edgesList[,2] == u,1]))

      # find a tree by traversing the graph
      dfsorder <- dfs(adjacencyList, 1)
      newEdgesList <- t(apply(dfsorder$E, 1, sort))

      ## delete unused edges
      MST$E$todel <- !duplicated(rbind(newEdgesList,
                                       edgesList))[-seq(1,luid-1)]
    }

    if (any(MST$E$todel))
      MST <- deleteEdges(MST)

  }

  ## return result
  MST
}

fit.FirstTreeCopulas_DAG <- function(MST, data.univ, type, copulaSelectionBy,
                                     testForIndependence, testForIndependence.level,
                                     se, weights = NA, method = "mle", cores = 1, DAG) {

  ## initialize estimation results with empty list
  d <- nrow(MST$E$nums)
  pc.data <- lapply(1:d, function(i) NULL)

  indep <- apply(MST$E$nums, 1, FUN = checkSepFirst_DAG, DAG = DAG)

  ## prepare for estimation and store names
  for (i in 1:d) {
    ## get edge and corresponding data
    a <- MST$E$nums[i, ]
    pc.data[[i]]$zr1 <- data.univ[, a[1]]
    pc.data[[i]]$zr2 <- data.univ[, a[2]]
    #         MST$E$Copula.Data.1[i] <- list(data.univ[, a[1]])
    #         MST$E$Copula.Data.2[i] <- list(data.univ[, a[2]])

    ## set names for this edge
    if (is.null(MST$V$names[a[1]])) {
      MST$E$Copula.CondName.1[i] <- a[1]
    } else {
      MST$E$Copula.CondName.1[i] <- MST$V$names[a[1]]
    }
    if (is.null(MST$V$names[a[2]])) {
      MST$E$Copula.CondName.2[i] <- a[2]
    } else {
      MST$E$Copula.CondName.2[i] <- MST$V$names[a[2]]
    }
    if (is.null(MST$V$names[a[1]]) || is.null(MST$V$names[a[2]])) {
      MST$E$Copula.Name[i] <- paste(a[1], a[2], sep = " , ")
    } else {
      MST$E$Copula.Name[i] <- paste(MST$V$names[a[1]],
                                    MST$V$names[a[2]],
                                    sep = " , ")
    }

  }

  ## estimate parameters and select family
  if (cores > 1) {
    pc.fits <- foreach(i = 1:length(pc.data),
                       .export = c("pcSelect",
                                   "fit.ACopula",
                                   "BiCopSelect",
                                   "wrapperDAG")) %dopar%
      wrapperDAG(i,
                 indep,
                 pc.data,
                 type,
                 copulaSelectionBy,
                 testForIndependence,
                 testForIndependence.level,
                 se,
                 weights,
                 method)

  } else {
    pc.fits <- lapply(as.list(1:length(pc.data)),
                      FUN = wrapperDAG,
                      indep = indep,
                      pc.data = pc.data,
                      type,
                      copulaSelectionBy,
                      testForIndependence,
                      testForIndependence.level,
                      se,
                      weights,
                      method)
  }

  ## store estimated model and pseudo-obversations for next tree
  for (i in 1:d) {
    MST$E$Copula.param[[i]] <- c(pc.fits[[i]]$par,
                                 pc.fits[[i]]$par2)
    MST$E$Copula.type[i] <- pc.fits[[i]]$family
    MST$E$fits[[i]] <- pc.fits[[i]]

    MST$E$Copula.CondData.1[i] <- list(pc.fits[[i]]$CondOn.1)
    MST$E$Copula.CondData.2[i] <- list(pc.fits[[i]]$CondOn.2)

    if (!is.null(pc.fits[[i]]$warn))
      MST$warn <- pc.fits[[i]]$warn
  }

  ## return results
  MST
}

## fit pair-copulas for vine trees 2,...
fit.TreeCopulas_DAG <- function(MST, oldVineGraph, type, copulaSelectionBy,
                                testForIndependence, testForIndependence.level,
                                se = se, progress, weights = NA, method = "mle",
                                cores = 1, DAG) {


  ## initialize estimation results with empty list
  d <- nrow(MST$E$nums)
  pc.data <- lapply(1:d, function(i) NULL)

  indep <- rep(NA, d)

  ## prepare for estimation
  for (i in 1:d) {
    ## get edge and corresponding data
    con <- MST$E$nums[i, ]
    temp <- oldVineGraph$E$nums[con, ]

    # condSet <- symdiff(temp[1,],temp[2,])
    condSet <- MST$E$conditionedSet[[i]]
    condingSet <- MST$E$conditioningSet[[i]]

    indep[i] <- checkSepHigher_DAG(condSet[1], condSet[2], condingSet, DAG = DAG)

    ## fetch corresponding data and names
    if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
      same <- temp[2, 1]
    } else {
      if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
        same <- temp[2, 2]
      }
    }

    if (temp[1, 1] == same) {
      zr1 <- oldVineGraph$E$Copula.CondData.2[con[1]]
      n1  <- oldVineGraph$E$Copula.CondName.2[con[1]]
    } else {
      zr1 <- oldVineGraph$E$Copula.CondData.1[con[1]]
      n1  <- oldVineGraph$E$Copula.CondName.1[con[1]]
    }
    if (temp[2, 1] == same) {
      zr2 <- oldVineGraph$E$Copula.CondData.2[con[2]]
      n2  <- oldVineGraph$E$Copula.CondName.2[con[2]]
    } else {
      zr2 <- oldVineGraph$E$Copula.CondData.1[con[2]]
      n2  <- oldVineGraph$E$Copula.CondName.1[con[2]]
    }

    if (is.list(zr1)) {
      zr1a <- as.vector(zr1[[1]])
      zr2a <- as.vector(zr2[[1]])
      n1a <- as.vector(n1[[1]])
      n2a <- as.vector(n2[[1]])
    } else {
      zr1a <- zr1
      zr2a <- zr2
      n1a <- n1
      n2a <- n2
    }

    if (progress == TRUE)
      message(n1a, " + ", n2a, " --> ", MST$E$names[i])

    pc.data[[i]]$zr1 <- zr1a
    pc.data[[i]]$zr2 <- zr2a

    #         MST$E$Copula.Data.1[i] <- list(zr1a)
    #         MST$E$Copula.Data.2[i] <- list(zr2a)

    MST$E$Copula.CondName.1[i] <- n1a
    MST$E$Copula.CondName.2[i] <- n2a
  }

  ## estimate parameters and select family
  if (cores > 1) {
    pc.fits <- foreach(i = 1:length(pc.data),
                       .export = c("pcSelect",
                                   "fit.ACopula",
                                   "BiCopSelect",
                                   "wrapperDAG")) %dopar%
      wrapperDAG(i,
                 indep,
                 pc.data,
                 type,
                 copulaSelectionBy,
                 testForIndependence,
                 testForIndependence.level,
                 se,
                 weights,
                 method)

  } else {
    pc.fits <- lapply(as.list(1:length(pc.data)),
                      FUN = wrapperDAG,
                      indep,
                      pc.data,
                      type,
                      copulaSelectionBy,
                      testForIndependence,
                      testForIndependence.level,
                      se,
                      weights,
                      method)
  }

  ## store estimated model and pseudo-obversations for next tree
  for (i in 1:d) {
    MST$E$Copula.param[[i]] <- c(pc.fits[[i]]$par,
                                 pc.fits[[i]]$par2)
    MST$E$Copula.type[i] <- pc.fits[[i]]$family
    MST$E$fits[[i]] <- pc.fits[[i]]
    MST$E$Copula.CondData.1[i] <- list(pc.fits[[i]]$CondOn.1)
    MST$E$Copula.CondData.2[i] <- list(pc.fits[[i]]$CondOn.2)
    if (!is.null(pc.fits[[i]]$warn))
      MST$warn <- pc.fits[[i]]$warn
  }

  ## return results
  MST
}

## initialize graph for next vine tree (possible edges)
buildNextGraph_DAG <- function(oldVineGraph, weights = NA, b, parallel,
                               truncated = FALSE, DAG) {

  d <- nrow(oldVineGraph$E$nums)

  ## initialize with full graph
  g <- makeFullGraph(d)
  g$V$names <- oldVineGraph$E$names
  g$V$conditionedSet <- oldVineGraph$E$conditionedSet
  g$V$conditioningSet <- oldVineGraph$E$conditioningSet

  ## get info for all edges
  if (parallel) {
    i <- NA  # dummy for CRAN checks
    out <- foreach(i = 1:nrow(g$E$nums)) %dopar%
      getEdgeInfo_DAG(i,
                      g = g,
                      oldVineGraph = oldVineGraph,
                      b = b,
                      weights = weights,
                      truncated = truncated,
                      DAG = DAG)
  } else {
    out <- lapply(1:nrow(g$E$nums),
                  getEdgeInfo_DAG,
                  g = g,
                  oldVineGraph = oldVineGraph,
                  b = b,
                  weights = weights,
                  truncated = truncated,
                  DAG = DAG)
  }

  ## annotate graph (same order as in old version of this function)
  g$E$weights         <- sapply(out, function(x) x$w)
  g$E$names           <- sapply(out, function(x) x$name)
  g$E$conditionedSet  <- lapply(out, function(x) x$nedSet)
  g$E$conditioningSet <- lapply(out, function(x) x$ningSet)
  g$E$todel           <- sapply(out, function(x) x$todel)

  ## delete edges that are prohibited by the proximity condition
  deleteEdges(g)
}

## function for obtaining edge information
getEdgeInfo_DAG <- function(i, g, oldVineGraph, b, weights,
                            truncated = FALSE, DAG) {

  ## get edge
  con <- g$E$nums[i, ]
  temp <- oldVineGraph$E$nums[con, ]


  ## check for proximity condition
  ok <- FALSE
  if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
    ok <- TRUE
    same <- temp[2, 1]
  } else {
    if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
      ok <- TRUE
      same <- temp[2, 2]
    }
  }

  ## dummy output
  w <- nedSet <- ningSet <- name <- NA
  todel <- TRUE

  # info if proximity condition is fulfilled ...
  if (ok) {
    ## infer conditioned set and conditioning set
    l1 <- c(g$V$conditionedSet[[con[1]]],
            g$V$conditioningSet[[con[1]]])
    l2 <- c(g$V$conditionedSet[[con[2]]],
            g$V$conditioningSet[[con[2]]])
    nedSet <- c(setdiff(l1, l2), setdiff(l2, l1))
    ningSet <- intersect(l1, l2)

    ## mark as ok
    todel <- FALSE

    if (truncated == FALSE) {
      ## get data
      if (temp[1, 1] == same) {
        zr1 <- oldVineGraph$E$Copula.CondData.2[con[1]]
      } else {
        zr1 <- oldVineGraph$E$Copula.CondData.1[con[1]]
      }
      if (temp[2, 1] == same) {
        zr2 <- oldVineGraph$E$Copula.CondData.2[con[2]]
      } else {
        zr2 <- oldVineGraph$E$Copula.CondData.1[con[2]]
      }
      if (is.list(zr1)) {
        zr1a <- as.vector(zr1[[1]])
        zr2a <- as.vector(zr2[[1]])
      } else {
        zr1a <- zr1
        zr2a <- zr2
      }

      ## calculate Kendall's tau
      keine_nas <- !(is.na(zr1a) | is.na(zr2a))
      if (checkSepHigher_DAG(nedSet[1], nedSet[2], ningSet, DAG)) {
        w <- b
      } else {
        w <- 0
      }

      ## get names
      name.node1 <- strsplit(g$V$names[con[1]], split = " *[,;] *")[[1]]
      name.node2 <- strsplit(g$V$names[con[2]], split = " *[,;] *")[[1]]

      ## set edge name
      nmdiff <- c(setdiff(name.node1, name.node2),
                  setdiff(name.node2, name.node1))
      nmsect <- intersect(name.node1, name.node2)
      name <- paste(paste(nmdiff, collapse = ","),
                    paste(nmsect, collapse = ","),
                    sep = " ; ")
    } else {
      w <- 1
    }
  }

  ## return edge information
  list(w = w,
       nedSet = nedSet,
       ningSet = ningSet,
       name = name,
       todel = todel)
}

# depth first search to build a tree without finding the MST
dfs <- function(adjacencyList, v, e = NULL, dfsorder = list()) {
  dfsorder$V <- c(dfsorder$V, v)
  dfsorder$E <- rbind(dfsorder$E, e)
  for (u in adjacencyList[[v]]) {
    if (!is.element(u, dfsorder$V)) {
      dfsorder <- dfs(adjacencyList, u, c(u,v), dfsorder)
    }
  }
  return(dfsorder)
}

pcSelect <- function(parameterForACopula, type, ...) {
  return(fit.ACopula(parameterForACopula$zr1,
                     parameterForACopula$zr2,
                     type,
                     ...))
}


## bivariate copula selection
fit.ACopula <- function(u1, u2, familyset = NA, selectioncrit = "AIC",
                        indeptest = FALSE, level = 0.05, se = FALSE,
                        weights = NA, method = "mle") {

  ## select family and estimate parameter(s) for the pair copula
  complete.i <- which(!is.na(u1 + u2))

  if (length(complete.i) < 10 || (length(familyset) == 1 && familyset == 0)) {
    out <- BiCop(0)
    ## add more information about the fit
    out$se  <- NA
    out$se2 <- NA
    out$nobs   <- 0
    out$logLik <- 0
    out$AIC    <- 0
    out$BIC    <- 0
    out$emptau <- NA
    out$p.value.indeptest <- NA

    if (length(complete.i) < 10) {
      out$warn <- paste("Insufficient data for at least one pair.",
                        "Independence has been selected automatically.")
    }
  } else {
    out <- suppressWarnings(BiCopSelect(u1[complete.i], u2[complete.i],
                                        familyset,
                                        selectioncrit,
                                        indeptest,
                                        level,
                                        weights = weights,
                                        rotations = FALSE,
                                        se = se,
                                        method = method))
    out$warn <- NULL
  }

  ## change rotation if family is not symmetric wrt the main diagonal
  if (out$family %in% c(23, 24, 26:30, 124, 224)) {
    out$family <- out$family + 10
  } else if (out$family %in% c(33, 34, 36:40, 134, 234)) {
    out$family <- out$family - 10
  }

  ## tawn copulas also have to change type
  if (out$family%/%100 == 1) {
    out$family <- out$family + 100
  } else if (out$family%/%200 == 1) {
    out$family <- out$family - 100
  }

  ## store pseudo-observations for estimation in next tree
  if (length(familyset) == 1 && familyset == 0) {
    out$CondOn.1 <- u1
    out$CondOn.2 <- u2
  } else {
    out$CondOn.1 <- suppressWarnings(BiCopHfunc1(u2, u1, out, check.pars = FALSE))
    out$CondOn.2 <- suppressWarnings(BiCopHfunc2(u2, u1, out, check.pars = FALSE))
  }

  ## return results
  out
}

## build R-Vine matrix object based on nested set of trees
as.RVM2 <- function(RVine, data, callexp) {

  ## initialize objects
  n <- length(RVine$Tree) + 1
  con <- list()
  nam <- RVine$Tree[[1]]$V$names
  nedSets <- list()
  crspParams <- list()
  crspTypes <- list()
  crspfits <- list()

  ## get selected pairs, families and estimated parameters
  for (k in 1:(n - 2)) {
    nedSets[[k]]    <- RVine$Tree[[k]]$E$conditionedSet
    crspParams[[k]] <- as.list(RVine$Tree[[k]]$E$Copula.param)
    crspTypes[[k]]  <- as.list(RVine$Tree[[k]]$E$Copula.type)
    crspfits[[k]]   <- as.list(RVine$Tree[[k]]$E$fits)

  }
  crspParams[[n - 1]] <- as.list(RVine$Tree[[n - 1]]$E$Copula.param)
  crspTypes[[n - 1]]  <- as.list(RVine$Tree[[n - 1]]$E$Copula.type)
  crspfits[[n - 1]]   <- as.list(RVine$Tree[[n - 1]]$E$fits)
  if (is.list(RVine$Tree[[1]]$E$conditionedSet)) {
    nedSets[[n - 1]] <- list(RVine$Tree[[n - 1]]$E$conditionedSet[[1]])
  } else {
    nedSets[[n - 1]] <- list(RVine$Tree[[n - 1]]$E$conditionedSet)
  }

  ## initialize matrices for RVineMatrix object
  Param <- array(dim = c(n, n))
  Params2 <- array(0, dim = c(n, n))
  Type <- array(dim = c(n, n))
  M <- matrix(NA, n, n)
  Ses     <- matrix(0, n, n)
  tmps    <- matrix(0, n, n)
  Se2s    <- matrix(0, n, n)
  emptaus <- matrix(0, n, n)
  pvals   <- matrix(0, n, n)
  logLiks <- matrix(0, n, n)

  ## store structure, families and parameters in matrices
  for (k in 1:(n - 1)) {
    w <- nedSets[[n - k]][[1]][1]

    M[k, k] <- w
    M[(k + 1), k] <- nedSets[[n - k]][[1]][2]

    Param[(k + 1), k]   <- crspParams[[n - k]][[1]][1]
    Params2[(k + 1), k] <- crspParams[[n - k]][[1]][2]
    Type[(k + 1), k]    <- crspTypes[[n - k]][[1]]
    tmpse               <- crspfits[[n - k]][[1]]$se
    Ses[(k + 1), k]     <- ifelse(is.null(tmpse), NA, tmpse)
    tmpse2              <- crspfits[[n - k]][[1]]$se2
    Se2s[(k + 1), k]    <- ifelse(is.null(tmpse2), NA, tmpse2)
    emptaus[(k + 1), k] <- crspfits[[n - k]][[1]]$emptau
    pvals[(k + 1), k]   <- crspfits[[n - k]][[1]]$p.value.indeptest
    logLiks[(k + 1), k] <- crspfits[[n - k]][[1]]$logLik

    if (k == (n - 1)) {
      M[(k + 1), (k + 1)] <- nedSets[[n - k]][[1]][2]
    } else {
      for (i in (k + 2):n) {
        for (j in 1:length(nedSets[[n - i + 1]])) {
          cs <- nedSets[[n - i + 1]][[j]]
          cty <- crspTypes[[n - i + 1]][[j]]
          if (cs[1] == w) {
            M[i, k] <- cs[2]
            Type[i, k] <- cty
            break
          } else if (cs[2] == w) {
            # correct family for rotations
            if (cty %in% c(23, 24, 26:30, 124, 224)) {
              cty <- cty + 10
            } else if (cty %in% c(33, 34, 36:40, 134, 234)) {
              cty <- cty - 10
            }
            # change type for Tawn
            if (cty%/%100 == 1) {
              cty <- cty + 100
            } else if (cty%/%200 == 1) {
              cty <- cty - 100
            }
            M[i, k] <- cs[1]
            Type[i, k] <- cty
            break
          }
        }
        Param[i, k]   <- crspParams[[n - i + 1]][[j]][1]
        Params2[i, k] <- crspParams[[n - i + 1]][[j]][2]
        tmpse         <- crspfits[[n - i + 1]][[j]]$se
        Ses[i, k]     <- ifelse(is.null(tmpse), NA, tmpse)
        tmpse2        <- crspfits[[n - i + 1]][[j]]$se2
        Se2s[i, k]    <- ifelse(is.null(tmpse2), NA, tmpse2)
        emptaus[i, k] <- crspfits[[n - i + 1]][[j]]$emptau
        pvals[i, k]   <- crspfits[[n - i + 1]][[j]]$p.value.indeptest
        logLiks[i, k] <- crspfits[[n - i + 1]][[j]]$logLik
        nedSets[[n - i + 1]][[j]]    <- NULL
        crspParams[[n - i + 1]][[j]] <- NULL
        crspTypes[[n - i + 1]][[j]]  <- NULL
        crspfits[[n - i + 1]][[j]] <- NULL
      }
    }
  }

  ## clean NAs
  M[is.na(M)] <- 0
  Type[is.na(Type)] <- 0

  ## create RVineMatrix object
  RVM <- RVineMatrix(M, family = Type, par = Param, par2 = Params2, names = nam)
  RVM$call <- callexp

  ## add information about pair-copulas
  if (!all(is.na(Ses[lower.tri(Ses)]))) {
    RVM$se <- Ses
    RVM$se2 <- Se2s
  }
  RVM$nobs <- crspfits[[1]][[1]]$nobs
  like <- suppressWarnings(RVineLogLik(data, RVM, calculate.V = FALSE))
  RVM$logLik <- like$loglik
  RVM$pair.logLik <- logLiks
  npar <- sum(RVM$family %in% allfams[onepar], na.rm = TRUE) +
    2 * sum(RVM$family %in% allfams[twopar], na.rm = TRUE)
  npar_pair <- matrix((RVM$family %in% allfams[onepar]) +
                        2 * (RVM$family %in% allfams[twopar]), nrow = n, ncol = n)
  N <- nrow(data)
  RVM$AIC <- -2 * like$loglik + 2 * npar
  RVM$pair.AIC <- -2 * logLiks + 2 * npar_pair
  RVM$BIC <- -2 * like$loglik + log(N) * npar
  RVM$pair.BIC <- -2 * logLiks + log(N) * npar_pair
  RVM$emptau <- emptaus

  ## return final object
  RVM
}


## functions for handling the tree structure -------------------------
graphFromWeightMatrix <- function(W) {
  d <- ncol(W)
  # get variable names
  nms <- colnames(W)
  if (is.null(nms))
    nms <- paste0("V", 1:d)
  # construct edge set
  E <- cbind(do.call(c, sapply(1:(d-1), function(i) seq.int(i))),
             do.call(c, sapply(1:(d-1), function(i) rep(i+1, i))))
  # add edge names
  E.names <- apply(E, 1, function(x) paste(nms[x[1]],  nms[x[2]], sep = ","))
  # set weights
  w <- W[upper.tri(W)]

  ## return results
  list(V = list(names = nms,
                conditionedSet = NULL,
                conditioningSet = NULL),
       E = list(nums = E,
                names = E.names,
                weights = w,
                conditionedSet = lapply(1:nrow(E), function(i) E[i, ]),
                conditioningSet = NULL))
}

makeFullGraph <- function(d) {
  ## create matrix of all combinations
  E <- cbind(do.call(c, lapply(1:(d-1), function(i) rep(i, d-i))),
             do.call(c, lapply(1:(d-1), function(i) (i+1):d)))
  E <- matrix(E, ncol = 2)

  ## output dummy list with edges set
  list(V = list(names = NULL,
                conditionedSet = NULL,
                conditioningSet = NULL),
       E = list(nums = E,
                names = NULL,
                weights = NULL,
                conditionedSet = E,
                conditioningSet = NULL))
}

adjacencyMatrix <- function(g) {
  ## create matrix of all combinations
  d <- length(g$V$names)
  v.all <- cbind(do.call(c, lapply(1:(d-1), function(i) seq.int(i))),
                 do.call(c, lapply(1:(d-1), function(i) rep(i+1, i))))

  ## fnd weight
  vals <- apply(v.all, 1, set_weight, E = g$E)

  ## create symmetric matrix of weights
  M <- matrix(0, d, d)
  M[upper.tri(M)] <- vals
  M <- M + t(M)
  diag(M) <- Inf

  ## return final matrix
  M
}

set_weight <- function(x, E) {
  ## convert weights so that minimum spanning tree algorithm can be applied
  is.edge <- (x[1] == E$nums[, 1]) & (x[2] == E$nums[, 2])
  if (!any(is.edge)) Inf else -E$weights[which(is.edge)]
}


deleteEdges <- function(g) {
  ## reduce edge list
  keep <- which(!g$E$todel)
  E <- list(nums            = matrix(g$E$nums[keep, ], ncol = 2),
            names           = g$E$names[keep],
            weights         = g$E$weights[keep],
            conditionedSet  = g$E$conditionedSet[keep],
            conditioningSet = g$E$conditioningSet[keep])

  ## return reduced graph
  list(V = g$V, E = E)
}

