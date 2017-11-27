## The Algorithm is an implementation of Soerensen & Janssens (2007) finding the best, i.e. weight-optimal k spanning trees of some given undirected graph.
## Input: some weighted graph
## Output: a list of graphs which can be analyzed if some favourable property exists
## Open: direct implementation of "favorable property" which shall be looked for - could speed up things a little
## PARALLEL PROGRAM

findKBestTreesP <- function(graph, k  = Inf, searchMode = "max", cores = 1) {
  library(igraph)
  library(optrees)
  library(doParallel)
  library(foreach)
  library(parallel)
  degree <- igraph::degree

  if (searchMode == "max") {
    E(graph)$weight <- -abs(E(graph)$weight)
    weight <- min(E(graph)$weight) * 2
    ## all edge weights are negative and "weight" is the smallest negative, so an edge with this weight is prefered
  } else {
    weight <- min(E(graph)$weight) / 2
    ## all edge weights are positive and "weight" is the smallest positive, so an edge with this weight is prefered
  }
  ## calculate number of spanning trees by Kirchhoff's formula
  graph.calc <- graph
  E(graph.calc)$weight <- 0
  number.of.trees <- det((diag(degree(graph)) - get.adjacency(graph, sparse = FALSE))[-1,-1]) # Kirchhoff's formula for number of spanning trees
  if (k > number.of.trees) {
    k <- number.of.trees - 1
  }

  ## register parallel backend
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

  ## aux functions
  identifyEdgeFromEdgelist <- function(rowcont, graph) {
    EL <- get.edgelist(graph)
    check1 <- which(rowcont[1] == EL[,1] & rowcont[2] == EL[,2])
    if(length(check1) == 0) {
      check2 <- which(rowcont[2] == EL[,1] & rowcont[1] == EL[,2])
      return(check2)
    } else {
      return(check1)
    }
  }
  assign.number.from.names <- function(arc, nodes) {
    which(nodes$name %in% arc)
  }
  assign.names.from.number <- function(arc, nodes) {
    nodes[arc[1:2]]$name
  }
  loop.check.fun <- function(id_internal, graph) { !E(graph)$status[id_internal] %in% c("included", "excluded") }
  loop.check.fun.v <- Vectorize(loop.check.fun, vectorize.args = "id_internal")

  assign.partition <- function(x, MST, P, loop.check.result) {
    if(loop.check.result[[x]] == FALSE) {
      out <- NULL
    } else {
      out <- P
      out$included <- rbind(out$included, get.edgelist(MST)[which(loop.check.result[0:(x-1)]),])
      out$excluded <- rbind(out$excluded, get.edgelist(MST)[x,])
    }
    return(out)
  }

  check.MST <- function(x, partition.list.internal.temp) {
    if(is.character(x$MST)) {
      if(x$MST == "not connected") {
        return(NULL)
      } else {
        out <- list(partition = partition.list.internal.temp[[x$id]], MST = x$MST)
        return(out)
      }
    } else {
      out <- list(partition = partition.list.internal.temp[[x$id]], MST = x$MST)
      return(out)
    }
  }
  ## end aux functions

  convertigraphOptTreeMST <- function(graph) {
    ## necessary, because in igraph only Prim's algorithm is implemented which might end up in not optimal results
    ## before converting, included edges have to be assigned a low weight
    ## before converting, excluded edges have to be deleted
    nodes.i <- V(graph)
    arcs.i <- get.edgelist(graph)
    weight.i <- E(graph)$weight
    nodes.o <- seq(1:length(nodes.i))
    arcs.i.l <- split(arcs.i, row(arcs.i))
    arcs.o.l <- lapply(arcs.i.l, assign.number.from.names, nodes = nodes.i)
    arcs.o <- matrix(unlist(arcs.o.l), ncol = 2, byrow = TRUE)
    arcs.all <- cbind(arcs.o, weight.i)
    MST <- getMinimumSpanningTree(nodes = nodes.o, arcs = arcs.all, algorithm = "Kruskal", show.data = FALSE, show.graph = FALSE)
    arcs.o <- MST$tree.arcs
    arcs.i <- t(apply(arcs.o, 1, assign.names.from.number, nodes = nodes.i))
    MST.out <- graph_from_edgelist(arcs.i, directed = FALSE) # since it is a spanning tree, all nodes are involved and do not need to be specified afterwards
    E(MST.out)$weight <- arcs.o[, 3] # hardcoding
    return(MST.out)
  }

  ## calculation of an MST given some included and excluded edges
  calculate.MST <- function(P_internal, graph, weight) {
    ## input: partition P
    ## full graph of which the MST has to be calculated
    ## weights for calculation of MST metric
    included.edges <- apply(P_internal$included, 1, FUN = identifyEdgeFromEdgelist, graph = graph)
    excluded.edges <- apply(P_internal$excluded, 1, FUN = identifyEdgeFromEdgelist, graph = graph)
    saveWeight <- E(graph)$weight[included.edges]
    E(graph)$weight[included.edges] <- rep(weight, length(included.edges)) # Assignment for included edges
    graph <- delete_edges(graph, E(graph)[excluded.edges])
    if (is_connected(graph)) {
      MST <- convertigraphOptTreeMST(graph)
      ## assign included/excluded status
      to.be.included.edges <- apply(P_internal$included, 1, FUN = identifyEdgeFromEdgelist, graph = MST) # mark edges which have been included by force
      E(MST)$weight[to.be.included.edges] <- saveWeight
    } else {
      MST <- "not connected"
    }
    return(MST)
  }

  ## partition function
  partition.MST <- function(P, MST, graph, weight, cores = 1) {
    partition.list.internal <- list()
    MST.list.internal <- list()
    #P1 <- P2 <- P
    included.edges <- apply(P$included, 1, FUN = identifyEdgeFromEdgelist, graph = MST)
    excluded.edges <- apply(P$excluded, 1, FUN = identifyEdgeFromEdgelist, graph = MST)
    E(MST)$status <- rep("null", length(E(MST)))
    E(MST)$status[included.edges] <- "included"
    E(MST)$status[excluded.edges] <- "excluded"

    loop.check.vector <- matrix(seq(1, length(E(MST)), 1), nrow = 1)
    loop.check.result <- loop.check.fun.v(loop.check.vector, graph = MST)

    partition.list.internal.temp <- as.list(seq(1, length(E(MST)), 1))
    partition.list.internal.temp <- lapply(partition.list.internal.temp, FUN = assign.partition, MST = MST, P = P, loop.check.result = loop.check.result)
    check.null <- lapply(partition.list.internal.temp, FUN = is.null) # NULL entrys can be removed because they are the same with the -1 entrys, i.e. one before
    partition.list.internal.temp[which(check.null == TRUE)] <- NULL

    if (cores > 1) {
      MST.out <- foreach(i = 1:length(partition.list.internal.temp),
                         .export = c("calculate.MST",
                                     "identifyEdgeFromEdgelist",
                                     "convertigraphOptTreeMST",
                                     "assign.number.from.names",
                                     "assign.names.from.number"),
                         .packages = c("optrees", "igraph")) %dopar% {
                           tempMST <- calculate.MST(partition.list.internal.temp[[i]], graph = graph, weight = weight)
                           out <- list(id = i, MST = tempMST)
                         }
    } else {
      MST.out <- list()
      tempMST <- lapply(partition.list.internal.temp, FUN = calculate.MST, graph = graph, weight = weight)
      MST.out <- list(id = i, MST = tempMST)
    }

    joint.partition.MST.list <- lapply(MST.out, check.MST, partition.list.internal.temp = partition.list.internal.temp)
    check.null <- lapply(joint.partition.MST.list, FUN = is.null)
    joint.partition.MST.list[which(check.null == TRUE)] <- NULL

    partition.list.internal <- lapply(joint.partition.MST.list, FUN = function(x) { x$partition })
    MST.list.internal <- lapply(joint.partition.MST.list, FUN = function(x) { x$MST })
    values.internal <- lapply(MST.list.internal, FUN = function(x) { sum(E(x)$weight)})
    out <- list(partition.list.internal = partition.list.internal, MST.list.internal = MST.list.internal, values.internal = values.internal)
    return(out)
  }

  MST <- minimum.spanning.tree(graph)
  E(MST)$status <- rep("null", length(E(MST)))

  ## initialize list for final storage
  result.list <- as.list(rep(NA, k))
  values.list <- as.list(rep(NA, k))
  ## initialize list for partition saving
  partition.list <- list()
  partition.list[[1]] <- list(included = matrix(0, nrow = 0, ncol = 2), excluded = matrix(0, nrow = 0, ncol = 2))
  MST.list <- list(MST)
  partition.values.list <- list(sum(E(MST)$weight))
  j <- 1

  while(j < k + 1) {
    ## step 1 of paper: get partition optimal
    id <- which.min(partition.values.list)
    ## step 2 of paper: save MST and value
    result.list[[j]] <- MST.list[[id]]
    values.list[[j]] <- partition.values.list[[id]]
    ## step 4 of paper: interchanged since otherwise we have to duplicate P_s
    temp.result <- partition.MST(partition.list[[id]], MST = MST.list[[id]], graph = graph, weight = weight, cores = cores)
    ## step 3 of paper: delete from list
    partition.list[[id]] <- NULL
    partition.values.list[[id]] <- NULL
    MST.list[[id]] <- NULL
    ## step 5: insert values from temp.result into partition.list, partition.values and MST.list
    if(length(temp.result$partition.list.internal) > 0) {
      partition.list[(length(partition.list) + 1) : (length(partition.list) + length(temp.result$partition.list.internal))] <- temp.result$partition.list.internal
      MST.list[(length(MST.list) + 1) : (length(MST.list) + length(temp.result$MST.list))] <- temp.result$MST.list.internal
      partition.values.list[(length(partition.values.list) + 1) : (length(partition.values.list) + length(temp.result$values))] <- temp.result$values.internal
    } else {
      # do nothing as no new MSTs have been found in the previous iteration
    }
    print(id)
    j <- j + 1
  }
  ## evaluation of corresponding weights
  return(list(MST = result.list, MST.weights = values.list, number.of.trees = number.of.trees))
}
