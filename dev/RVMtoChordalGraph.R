RVMtoChordalGraph <- function(Matrix, family_Matrix, names) {
  require(igraph)
  d <- length(names)
  amat <- matrix(0, d, d)
  colnames(amat) <- rownames(amat) <- names
  for(i in d:2) {
    for(j in (i-1):1) {
      if (family_Matrix[i,j] != 0) {
        amat[Matrix[i,j], Matrix[j,j]] = d - i + 1 # weight will be tree number
      }
    }
  }
  graph <- graph_from_adjacency_matrix(amat, mode = "undirected", weighted = TRUE) # convert object to igraph
  E(graph)$label <- E(graph)$weight
  graph
}
