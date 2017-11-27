## Simple evaluation function for some cluster metrics

getClusterMetrics <- function(graph) {
  comp <- components(graph_from_adjacency_matrix(graph, mode = "undirected"))
  iterator <- which(comp$csize > 1)
  ## get subgraphs for adjacency
  getNumEdges <- function(cluster_id, comp, graph) {
    vertex_ids <- which(comp$membership == cluster_id)
    sum(graph[vertex_ids, vertex_ids]) / 2
  }
  n_vertices <- comp$csize[iterator]
  n_edges <- sapply(iterator, FUN = getNumEdges, comp = comp, graph = graph)
  factors_1 <- n_edges / n_vertices
  factors_2 <- n_edges / choose(n_vertices, 2)
  prod_1 <- prod(factors_1)
  prod_2 <- prod(factors_2)
  out <- list(prod_1 = prod_1, prod_2 = prod_2)
  return(out)
}

