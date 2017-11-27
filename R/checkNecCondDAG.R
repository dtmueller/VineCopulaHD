checkNecCondDAG <- function(DAG) {
  degree <- igraph::degree
  dsep <- bnlearn::dsep
  calc.amat <- function(graph, data.names){
    h.amat <- matrix(0, length(data.names), length(data.names))
    colnames(h.amat) <- data.names
    rownames(h.amat) <- data.names
    h.amat[graph$arcs] <- 1
    h.amat
  }
  ## Input is either a bnlearn-DAG or an igraph-DAG
  if (is_igraph(DAG)) {
    DAG_bnlearn <- to_bn(DAG)
    DAG_igraph <- DAG
  } else {
    DAG_bnlearn <- DAG
    DAG_igraph <- graph_from_adjacency_matrix(calc.amat(DAG, names(DAG$nodes)), mode = "directed")
  }
  names <- V(DAG_igraph)$name
  V(DAG_igraph)$name <- seq(1:length(V(DAG_igraph)))
  DAG_bnlearn <- to_bn(DAG_igraph)

  warning_message <- "NULL"
  nodes <- V(DAG_igraph)$name
  parental_sets <- neighborhood(DAG_igraph, 1, mode = "in")
  l_parental_sets <- unlist(lapply(parental_sets, FUN = length))
  m_parental_sets <- max(l_parental_sets)
  selected_parental_sets <- parental_sets[which(l_parental_sets == m_parental_sets)]

  triplets <- combn(nodes, 3)
  out <- NULL

  checkPerm <- function(parental_sets, triplet) {
    trip1 <- triplet[c(1,2)]
    trip2 <- triplet[c(1,3)]
    trip3 <- triplet[c(2,3)]
    checkIn <- function(nbh_entry, trip) {
      all(trip %in% nbh_entry)
    }
    trip1_check <- which(unlist(lapply(parental_sets, FUN = checkIn, trip = trip1)))
    trip2_check <- which(unlist(lapply(parental_sets, FUN = checkIn, trip = trip2)))
    trip3_check <- which(unlist(lapply(parental_sets, FUN = checkIn, trip = trip3)))
    if (all(c(length(trip1_check) > 0, length(trip2_check) > 0, length(trip3_check) > 0))) {
      if (all(c(length(unique(c(trip1_check, trip2_check, trip3_check))) > length(unique(trip1_check)),
                length(unique(c(trip1_check, trip2_check, trip3_check))) > length(unique(trip2_check)),
                length(unique(c(trip1_check, trip2_check, trip3_check))) > length(unique(trip3_check))))) {
        return(triplet)
      } else {
        ## do nothing
      }
    } else {
      ## do nothing
    }
  }

  for (i in 1:ncol(triplets)){
    out <- checkPerm(selected_parental_sets, triplets[,i])
    if(!is.null(out)) {
      warning_message <- "cycles present => A1 violated"
      break
    }
  }

  degrees <- degree(DAG_igraph)
  if(any(degrees == 0)) {
    ## errorhandling in case of isolated vertices which must be independent of all
    isolated <- which(degrees == 0)
    skeleton <- as.undirected(DAG_igraph)
    DAG_igraph_wo_iso <- delete_vertices(DAG_igraph, isolated)
    top_sort <- rev(topo_sort(DAG_igraph_wo_iso))
    RVM <- diag(c(as.numeric(isolated),as.numeric(names(top_sort))))
  } else {
    top_sort <- rev(topo_sort(DAG_igraph))
    RVM <- diag(as.numeric(names(top_sort)))
  }
  d <- ncol(RVM)
  ## first tree
  applySelectFirstTree <- function(j, DAG_bnlearn, RVM) {
    parentset <- parents(DAG_bnlearn, as.character(RVM[j,j]))
    if(length(parentset) == 1) {
      RVM[d,j] <- as.numeric(parentset)
    } else {
      if(length(parentset) > 1) {
        ## attach nodes according to the negative top sort, i.e. from left to right in the RVM
        return(diag(RVM[(j+1):d,(j+1):d])[which(diag(RVM[(j+1):d,(j+1):d]) %in% parentset)][1]) ## take the first one
      } else {
        return(RVM[j+1,j+1]) ## dummy setting
      }
    }
  }
  RVM[d,1:(d-1)] <- apply(as.matrix(c(1:(d-1))), 1, FUN = applySelectFirstTree, DAG_bnlearn = DAG_bnlearn, RVM = RVM)

  ## first tree finished, complete according to proximity an going for the parents on a tree-to-tree basis
  ## Todo: hier könnte man ggf. noch einbauen, dass aufbauend auf T_1 gesucht werden soll welcher der nächste verfügbare knoten ist bzgl. Distanz

  applySelectHigherTree <- function(j, RVM, DAG_bnlearn, i) {
    possible_entries <- as.list(setdiff(diag(RVM[(j+1):d,(j+1):d]),RVM[(i+1):d,j]))
    allowed_entries <- unlist(possible_entries)[unlist(lapply(possible_entries, FUN = checkProximity, RVS = RVM, k = d-i+1, x = RVM[j,j]))]
    parentset <- parents(DAG_bnlearn, as.character(RVM[j,j]))
    parentset_open <- setdiff(as.numeric(parentset), RVM[(i+1):d,j])
    parentset_open_allowed <- intersect(allowed_entries, parentset_open)
    if(length(parentset_open_allowed) == 0) {
      ## complete only according to the proximity condition
      return(allowed_entries[1])
    } else {
      return(parentset_open_allowed[1]) ## muessten bereits korrekt sortiert sein
    }
  }

  for (i in (d-1):2) {
    RVM[i,1:(i-1)] <- apply(as.matrix(c(1:(i-1))), 1, FUN = applySelectHigherTree, RVM = RVM, DAG_bnlearn = DAG_bnlearn, i = i)
  }

  checkDsepFirst <- function(j, DAG_bnlearn, RVM) {
    d <- ncol(RVM)
    if (!dsep(DAG_bnlearn, as.character(RVM[j,j]), as.character(RVM[d,j]))) {
      return(1)
    } else {
      return(0)
    }
  }

  checkDsepHigher <- function(j, DAG_bnlearn, RVM, i) {
    d <- ncol(RVM)
    if (!dsep(DAG_bnlearn, as.character(RVM[j,j]), as.character(RVM[i,j]), as.character(RVM[(i+1):d,j]))) {
      return(1)
    } else {
      return(0)
    }
  }
  IndepMatrix <- matrix(0, d, d)
  IndepMatrix[d,1:(d-1)] <- apply(as.matrix(c(1:(d-1))), 1, FUN = checkDsepFirst, RVM = RVM, DAG_bnlearn = DAG_bnlearn)

  for (i in (d-1):2) {
    IndepMatrix[i,1:(i-1)] <- apply(as.matrix(c(1:(i-1))), 1, FUN = checkDsepHigher, RVM = RVM, DAG_bnlearn = DAG_bnlearn, i = i)
  }

  ratio <- sum(IndepMatrix)/choose(d,2)
  k <- max(degree(DAG_igraph, mode = "in"))
  k_prime <- (d + 1) - min(which(apply(IndepMatrix, 1, FUN = sum) > 0)) # actual truncation level of R-vine
  out <- list(warning_message = warning_message, RVM = list(RVM = RVM, IndepMatrix = IndepMatrix, names = names, k = k, k_prime = k_prime, ratio = ratio))
  return(out)
}
