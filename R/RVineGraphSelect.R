## Function taken from VineCopula package and adapted
## Function to estimate a sparse R-vine based on an undirected graphical model

RVineGraphSelect <- function(data, graph = NULL, familyset = NA, type = 0, selectioncrit = "AIC", indeptest = FALSE,
                              level = 0.05, trunclevel = NA, progress = FALSE,  weights = NA,
                              treecrit = "AIC", se = FALSE, rotations = TRUE, method = "itau", cores = 1) {

  if(is.null(graph))
    stop("No dependence structure supplied.")

  # if (any(is.na(familyset))) {
  #   familyset <- c(0,1,2,3,4,5,6,13,14,16,23,24,26,33,34,36)
  # } else {
  #   if (any(familyset %in% setdiff(VineCopula:::allfams, c(0,1,2,3,4,5,6,13,14,16,23,24,26,33,34,36)))) {
  #     familyset <- intersect(familyset, c(0,1,2,3,4,5,6,13,14,16,23,24,26,33,34,36))
  #   }
  # }

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
                          method = method,
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

  treecrit <- set_treecrit_graph6(treecrit, familyset, method)
  if (is.na(trunclevel))
    trunclevel <- d
  if (trunclevel == 0)
    familyset <- 0

  ## initialize object for results
  RVine <- list(Tree = NULL, Graph = NULL)
  warn <- NULL

  ## novel part
  graph_igraph <- graph_from_adjacency_matrix(graph, mode = "undirected")

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

  g <- initializeFirstGraph_graph6(data, treecrit, weights, graph, cl)
  MST <- findMaxTree_graph6(g, mode = type)

  # estimate pair-copulas
  VineTree <- fit.FirstTreeCopulas_graph6(MST,
                                          data,
                                          type = familyset,
                                          selectioncrit,
                                          indeptest,
                                          level,
                                          se = se,
                                          weights = weights,
                                          method = method,
                                          cores = cores)
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
    g <- buildNextGraph_graph6(VineTree, weights, treecrit = treecrit, cores > 1,
                               truncated = trunclevel < tree, graph_igraph)

    MST <- findMaxTree_graph6(g, mode = type, truncated = trunclevel < tree)
    # estimate pair-copulas
    VineTree <- fit.TreeCopulas_graph6(MST,
                                       VineTree,
                                       type = familyset,
                                       selectioncrit,
                                       indeptest,
                                       level,
                                       se = se,
                                       progress,
                                       weights = weights,
                                       method = method,
                                       cores = cores)
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
  .RVine <- RVine
  .data <- data
  .callexp <- match.call()
  rm(list = ls())
  as.RVM2(.RVine, .data, .callexp)
}

# bindMSTPC <- function(MST, edges = NULL, PC) {
#   if (edges)
#     id <- paste(MST$E$nums[,1], MST$E$nums[,2], sep = "a")
#   edges_id <- paste(edges[,1], edges[,2], sep = "a")
#   MST$E$PC <- lapply(as.list(1:nrow(MST$E$nums)), FUN = function(x) {BiCop(0)})
#   pos_ids <- which(edges_id %in% id)
#   ## find position of pos_ids in id
#   MST$E$PC[checkPosV(edges_id[pos_ids], id)] <- PC[pos_ids]
#   MST
# }
#
# checkPos <- function(item, checkList) {
#   which(checkList == item)
# }
#
# checkPosV <- Vectorize(checkPos, vectorize.args = "item")

# symdiff <- function(x, y) { setdiff( c(x, y), intersect(x, y))}

# checkSepFirst <- function(triplet, graph) {
#   if (are_adjacent(graph, triplet[1], triplet[2])) {
#     TRUE
#   } else {
#     FALSE
#   }
# }

checkSepHigher <- function(j, k, D, graph) {
  if (is_separator(graph, D)) {
    temp <- components(delete_vertices(graph, D))$membership
    j <- j - sum(D < j)
    k <- k - sum(D < k)
    if (temp[j] != temp[k]) {
      FALSE
    } else {
      TRUE
    }
  } else {
    TRUE
  }
}

set_treecrit_graph6 <- function(treecrit, famset, method) {
  ## check if function is appropriate or type is implemented
  if (is.function(treecrit)) {
    if (!all(names(formals(treecrit)) == c("u1", "u2", "weights")))
      stop("treecrit must be of the form 'function(u1, u2, weights)'")
    if (!is.numeric(treecrit(runif(10), runif(10), rep(1, 10))))
      stop("treecrit does not return a numeric value")
  } else if (all(treecrit == "tau")) {
    treecrit <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 10) {
        tau <- 0
      } else {
        complete.freq <- mean(!is.na(u1 + u2))
        tau <- abs(fasttau(u1[complete.i], u2[complete.i], weights))
        tau * sqrt(complete.freq)
      }
    }
  } else if (all(treecrit == "rho")) {
    treecrit <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 10) {
        tau <- 0
      } else {
        complete.freq <- mean(!is.na(u1 + u2))
        rho <- abs(cor(u1, u2, method = "spearman", use = "complete.obs"))
        rho * sqrt(complete.freq)
      }
    }
  } else if (all(treecrit == "AIC")) {
    treecrit <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 2) {
        AIC <- 0
      } else {
        temp <- suppressWarnings(BiCopSelect(u1[complete.i], u2[complete.i],
                                             familyset = famset,
                                             weights = weights,
                                             indeptest = TRUE,
                                             level = 0.05,
                                             method = method))
      }
      list(val = -temp$AIC, PC = temp)
    }
  } else if (all(treecrit == "BIC")) {
    treecrit <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 2) {
        BIC <- 0
      } else {
        temp <- suppressWarnings(BiCopSelect(u1[complete.i], u2[complete.i],
                                             familyset = famset,
                                             weights = weights,
                                             indeptest = TRUE,
                                             level = 0.05,
                                             method = method))
      }
      list(val = -temp$BIC, PC = temp)
    }
  } else if (all(treecrit == "cAIC")) {
    treecrit <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 2) {
        cAIC <- 0
      } else {
        fit <- suppressWarnings(BiCopSelect(u1[complete.i], u2[complete.i],
                                            familyset = famset,
                                            weights = weights,
                                            indeptest = TRUE,
                                            level = 0.05,
                                            method = method))
        n <- fit$nobs
        p <- fit$npars
        cAIC <- - (fit$AIC + 2 * p * (p + 1) / (n - p - 1))
      }
      list(val = cAIC, PC = fit)
    }
  } else {
    txt1 <- 'treecrit must be one of "tau", "rho", "AIC", "BIC", "cAIC"'
    txt2 <- 'or a function like function(u1, u2, weights) ... returning a numeric value.'
    stop(paste(txt1, txt2))
  }

  ## return treecrit function
  treecrit
}

calcTreeCrit <- function(pair, graph, data, treecrit, weights, cl = NULL) {
  if (graph[pair[1], pair[2]] != 0) {
    w <- treecrit(data[, pair[1]], data[, pair[2]], weights)
    PC <- w$PC
    w <- w$val
  } else {
    w <- 0
    PC <- BiCop(0)
    PC$nobs <- 0
    PC$logLik <- 0
    PC$AIC <- 0
    PC$BIC <- 0
    PC$emptau <- 0
    PC$p.value.indeptest <- 0
  }
  list(w = w, PC = PC)
}

initializeFirstGraph_graph6 <- function(data, treecrit, weights, graph, cl = NULL) {
  ## calculate edge weight for each possible pair
  # all.pairs <- combn(1:ncol(data), 2)
  all.pairs <- combn(ncol(data), 2)
  if (!is.null(cl)) {
    edge.ws <- parApply(cl, all.pairs, 2, calcTreeCrit, graph, data, treecrit, weights)
  } else {
    edge.ws <- apply(all.pairs, 2, calcTreeCrit, graph, data, treecrit, weights)
  }

  ws <- sapply(edge.ws, FUN = function(x) {x$w})
  PC <- lapply(edge.ws, FUN = function(x) {x$PC})

  # number of pairwise complete observations / all observations
  rel.nobs <- apply(all.pairs, 2,
                    function(ind)
                      mean(!is.na(data[, ind[1]] + data[, ind[2]])))

  ## store in symmetric matrix with appropriate names
  W <- diag(ncol(data))
  W[lower.tri(W)] <- ws
  # W <- t(W)
  colnames(W) <- rownames(W) <- colnames(data)

  # W <- diag(ncol(data))
  # W[which(graph != 0)] <- ws
  # colnames(W) <- rownames(W) <- colnames(data)

  ## return full graph with edge weights
  graphFromWeightMatrix_graph6(W, PC)
}

findMaxTree_graph6 <- function(g, mode = "RVine", truncated = FALSE) {

  if (truncated == FALSE) {
    ## construct adjency matrix
    A <- adjacencyMatrix_graph6(g)
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

        # find edge with minimal weight
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
        MST <- deleteEdges_graph6(g)
      }
    } else if (mode  == "CVine") {
      ## set root as vertex with minimal sum of weights
      A <- adjacencyMatrix_graph6(g)
      diag(A) <- 0
      sumtaus <- rowSums(A)
      root <- which.min(sumtaus)

      ## delete unused edges
      g$E$todel <- !((g$E$nums[, 2] == root) | (g$E$nums[, 1] == root))
      MST <- g
      if (any(g$E$todel ))
        MST <- deleteEdges_graph6(g)
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
      MST <- deleteEdges_graph6(MST)

  }

  ## return result
  MST
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

# not required any longer? Use TauMatrix instead
fasttau <- function(x, y, weights = NA) {
  if (any(is.na(weights))) {
    m <- length(x)
    n <- length(y)
    if (m == 0 || n == 0)
      stop("both 'x' and 'y' must be non-empty")
    if (m != n)
      stop("'x' and 'y' must have the same length")
    out <- .C("ktau",
              x = as.double(x),
              y = as.double(y),
              N = as.integer(n),
              tau = as.double(0),
              S = as.double(0),
              D = as.double(0),
              T = as.integer(0),
              U = as.integer(0),
              V = as.integer(0),
              PACKAGE = "VineCopula")
    ktau <- out$tau
  } else {
    ktau <- TauMatrix(matrix(c(x, y), length(x), 2), weights)[2, 1]
  }
  return(ktau)
}

## fit pair-copulas for the first vine tree
fit.FirstTreeCopulas_graph6 <- function(MST, data.univ, type, copulaSelectionBy,
                                        testForIndependence, testForIndependence.level,
                                        se, weights = NA, method = "mle", cores = 1) {

  ## initialize estimation results with empty list
  d <- nrow(MST$E$nums)
  pc.data <- lapply(1:d, function(i) NULL)

  ## prepare for estimation and store names
  for (i in 1:d) {
    ## get edge and corresponding data
    a <- MST$E$nums[i, ]
    pc.data[[i]]$zr1 <- data.univ[, a[1]]
    pc.data[[i]]$zr2 <- data.univ[, a[2]]
    pc.data[[i]]$PC <- MST$E$PC[[i]]
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
                       .export = c("pcSelect_graph6",
                                   "fit.ACopula_graph6",
                                   "BiCopSelect")) %dopar%
      pcSelect_graph6(pc.data[[i]],
                type,
                copulaSelectionBy,
                testForIndependence,
                testForIndependence.level,
                se,
                weights,
                method)

  } else {
    pc.fits <- lapply(X = pc.data,
                      FUN = pcSelect_graph6,
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
fit.TreeCopulas_graph6 <- function(MST, oldVineGraph, type, copulaSelectionBy,
                                   testForIndependence, testForIndependence.level,
                                   se = se, progress, weights = NA, method = "mle",
                                   cores = 1) {


  ## initialize estimation results with empty list
  d <- nrow(MST$E$nums)
  pc.data <- lapply(1:d, function(i) NULL)

  ## prepare for estimation
  for (i in 1:d) {
    ## get edge and corresponding data
    con <- MST$E$nums[i, ]
    temp <- oldVineGraph$E$nums[con, ]

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
    pc.data[[i]]$PC <- MST$E$PC[[i]]

    #         MST$E$Copula.Data.1[i] <- list(zr1a)
    #         MST$E$Copula.Data.2[i] <- list(zr2a)

    MST$E$Copula.CondName.1[i] <- n1a
    MST$E$Copula.CondName.2[i] <- n2a
  }

  if (d == 1) {
    pc.data[[1]]$PC <- MST$E$PC
  }

  ## estimate parameters and select family
  if (cores > 1) {
    pc.fits <- foreach(i = 1:length(pc.data),
                       .export = c("pcSelect_graph6",
                                   "fit.ACopula_graph6",
                                   "BiCopSelect")) %dopar%
      pcSelect_graph6(pc.data[[i]],
                type,
                copulaSelectionBy,
                testForIndependence,
                testForIndependence.level,
                se,
                weights,
                method)

  } else {
    pc.fits <- lapply(X = pc.data,
                      FUN = pcSelect_graph6,
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
buildNextGraph_graph6 <- function(oldVineGraph, treecrit, weights = NA, parallel,
                                  truncated = FALSE, graph) {

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
      getEdgeInfo_graph6(i,
                         g = g,
                         oldVineGraph = oldVineGraph,
                         treecrit = treecrit,
                         weights = weights,
                         truncated = truncated,
                         graph = graph)
  } else {
    out <- lapply(1:nrow(g$E$nums),
                  getEdgeInfo_graph6,
                  g = g,
                  oldVineGraph = oldVineGraph,
                  treecrit = treecrit,
                  weights = weights,
                  truncated = truncated,
                  graph = graph)
  }

  ## annotate graph (same order as in old version of this function)
  g$E$weights         <- sapply(out, function(x) x$w)
  g$E$names           <- sapply(out, function(x) x$name)
  g$E$conditionedSet  <- lapply(out, function(x) x$nedSet)
  g$E$conditioningSet <- lapply(out, function(x) x$ningSet)
  g$E$todel           <- sapply(out, function(x) x$todel)
  g$E$PC              <- lapply(out, function(x) x$PC)

  if (length(out) == 1) {
    g$E$PC <- out[[1]]$PC
  }

  ## delete edges that are prohibited by the proximity condition
  deleteEdges_graph6(g)
}

## function for obtaining edge information
getEdgeInfo_graph6 <- function(i, g, oldVineGraph, treecrit, weights,
                               truncated = FALSE, graph) {

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
  PC <- w <- nedSet <- ningSet <- name <- NA
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
      if (checkSepHigher(nedSet[1], nedSet[2], ningSet, graph)) {
        w <- treecrit(zr1a[keine_nas], zr2a[keine_nas], weights)
        PC <- w$PC
        w <- w$val
      } else {
        w <- 0
        PC <- BiCop(0)
        PC$nobs <- 0
        PC$logLik <- 0
        PC$AIC <- 0
        PC$BIC <- 0
        PC$emptau <- 0
        PC$p.value.indeptest <- 0
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
  } else {
  }

  ## return edge information
  list(w = w,
       nedSet = nedSet,
       ningSet = ningSet,
       name = name,
       todel = todel,
       PC = PC)
}


pcSelect_graph6 <- function(parameterForACopula, type, ...) {
  return(fit.ACopula_graph6(parameterForACopula$zr1,
                      parameterForACopula$zr2,
                      parameterForACopula$PC,
                      type,
                      ...))
}


## bivariate copula selection
fit.ACopula_graph6 <- function(u1, u2, PC, familyset = NA, selectioncrit = "AIC",
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
    # out <- suppressWarnings(BiCopSelect(u1[complete.i], u2[complete.i],
    #                                     familyset,
    #                                     selectioncrit,
    #                                     indeptest,
    #                                     level,
    #                                     weights = weights,
    #                                     rotations = FALSE,
    #                                     se = se,
    #                                     method = method))
    out <- PC
    # out$warn <- NULL
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
graphFromWeightMatrix_graph6 <- function(W, PC) {
  d <- ncol(W)
  # get variable names
  nms <- colnames(W)
  if (is.null(nms))
    nms <- paste0("V", 1:d)
  # construct edge set
  E <- t(combn(d,2))

  # E <- cbind(do.call(c, sapply(1:(d-1), function(i) seq.int(i))),
  #            do.call(c, sapply(1:(d-1), function(i) rep(i+1, i))))
  # add edge names
  E.names <- apply(E, 1, function(x) paste(nms[x[1]],  nms[x[2]], sep = ","))
  # set weights
  w <- W[lower.tri(W)]

  ## return results
  list(V = list(names = nms,
                conditionedSet = NULL,
                conditioningSet = NULL),
       E = list(nums = E,
                names = E.names,
                weights = w,
                conditionedSet = lapply(1:nrow(E), function(i) E[i, ]),
                conditioningSet = NULL,
                PC = PC))
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

adjacencyMatrix_graph6 <- function(g) {
  ## create matrix of all combinations
  d <- length(g$V$names)
  # v.all <- cbind(do.call(c, lapply(1:(d-1), function(i) seq.int(i))),
  #                do.call(c, lapply(1:(d-1), function(i) rep(i+1, i))))
  v.all <- t(combn(d,2))

  ## fnd weight
  vals <- apply(v.all, 1, set_weight, E = g$E)

  ## create symmetric matrix of weights
  M <- matrix(0, d, d)
  M[lower.tri(M)] <- vals
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


deleteEdges_graph6 <- function(g) {
  ## reduce edge list
  keep <- which(!g$E$todel)
  if (length(keep) == 1) {
    spec <- g$E$PC
  } else {
    spec <- g$E$PC[keep]
  }
  E <- list(nums            = matrix(g$E$nums[keep, ], ncol = 2),
            names           = g$E$names[keep],
            weights         = g$E$weights[keep],
            conditionedSet  = g$E$conditionedSet[keep],
            conditioningSet = g$E$conditioningSet[keep],
            PC = spec)

  ## return reduced graph
  list(V = g$V, E = E)
}

