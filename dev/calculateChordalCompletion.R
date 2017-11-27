## calculates a chordal completion and saves all eligible edge additions
## Output is a sequence of chordal graphs for each R-vine tree

calculateChordalCompletion <- function(start_graph, end_graph, data, method = "pearson") {
  require(igraph)
  require(ppcor)
  require(RColorBrewer)
  pcor <- ppcor::pcor

  if(!is.chordal(start_graph)$chordal | !is.chordal(end_graph)$chordal) {
    return("Start graph or end graph is not chordal.")
  }
  data_names <- V(start_graph)$name
  d <- length(data_names)

  full_graph <- graph.full(d, directed = FALSE)
  V(full_graph)$name <- data_names
  edgelist_for_vine <- as.list(rep(NA, choose(d, 2)))
  transform.Edgelist.To.List <- function(edgelist) {
    list(j = edgelist[1], k = edgelist[2], D = NULL)
  }
  edgelist_for_vine[1:length(E(start_graph))] <- apply(get.edgelist(start_graph), 1, FUN = transform.Edgelist.To.List)
  iterator <- length(E(start_graph)) + 1
  temp_graph <- start_graph
  E(temp_graph)$tree_level <- rep(1, length(E(temp_graph)))
  E(temp_graph)$label <- rep(1, length(E(temp_graph)))
  colorScheme <- brewer.pal(d-1, "Reds")[(d-1) : 1]
  E(temp_graph)$color <- rep(colorScheme[1], length(E(temp_graph)))

  # stages <- as.list(seq(length(E(start_graph)), 1))
  # cum_stages <- lapply(stages, FUN = function(x, stages) {i <- which(stages == x); x <- x + sum(unlist(stages[0:(i-1)]))}, stages = stages)

  ## DEFINITION OF SUPPORT FUNCTIONS

  count.Level <- function(edge) {
    if(length(edge) == 3){
      return(length(edge$D))
    } else {
      return(NULL)
    }
  }

  ## function for checking which edges can be added to obtain a chordal graph with an additional edge based on some graph
  check.Add.Edge.Chordal <- function(edge, graph) {
    is.chordal(add_edges(graph, edge))$chordal
  }

  eval.Pcor <- function(edge_name, data, method) {
    if(!is.null(edge_name)) {
      pcor.test(x = data[,edge_name$j], y = data[,edge_name$k], z = data[,edge_name$D], method = method)$estimate
    } else {
      return(0)
    }
  }

  check.Graph.Chordal <- function(graph, edgelist_difference, end_graph) {
    completion_list <- apply(edgelist_difference, 1, FUN = check.Add.Edge.Chordal, graph = graph)
    edgelist_difference_1 <- matrix(edgelist_difference[which(completion_list == TRUE),], ncol = 2, byrow = FALSE)
    edge_difference_1 <- graph_from_edgelist(edgelist_difference_1, directed = FALSE)
    edge_difference_2 <- graph.intersection(edge_difference_1, end_graph)
    completion_edges <- get.edgelist(edge_difference_2)
    return(completion_edges)
  }

  check.Graph.Chordal.Fallback <- function(graph, edgelist_difference, end_graph) {
    completion_list <- apply(edgelist_difference, 1, FUN = check.Add.Edge.Chordal, graph = graph)
    edgelist_difference_1 <- edgelist_difference[which(completion_list == TRUE),]
    return(edgelist_difference_1)
  }

  search.Cliques.aIo <- function(selected_edge, temp_graph) {
    ad_vertices <- adjacent_vertices(temp_graph, v = selected_edge, mode = "all")
    all_names <- union(selected_edge, union(names(unlist(ad_vertices[[1]])), names(unlist(ad_vertices[[2]]))))
    ## um hier etwas speed reinzubringen könnte man noch versuchen ein lower bound für die cliques zu finden, wenn das was hilft
    cut_graph <- induced_subgraph(temp_graph, vids = V(temp_graph)[which(names(V(temp_graph)) %in% all_names)])
    # upper_limit <- max(E(cut_graph)$treelevel) + 1
    cut_graph_added <- add_edges(cut_graph, edges = selected_edge)
    clique_added <- setdiff(cliques(cut_graph_added), cliques(cut_graph))
    selected_clique <- clique_added[which.max(unlist(lapply(clique_added, FUN = length)))]
    cross_check <- E(induced_subgraph(temp_graph, vids = V(temp_graph)[which(names(V(temp_graph)) %in% names(selected_clique[[1]]))]))$tree_level
    if (length(selected_clique[[1]]) != max(cross_check) + 2) {
      return(NULL)
      ## an edge is used which belongs to the considered tree or a higher tree
    } else {
      intersect_final <- setdiff(names(selected_clique[[1]]), selected_edge)
      out <- list(j = selected_edge[1], k = selected_edge[2], D = intersect_final)
      return(out)
    }
  }

  extract.Values <- function(x) {
    return(list(edge = c(x$j, x$k), level = length(x$D) + 1))
  }

  get.Edge <- function(list_element, node_name) {
    if (length(list_element) > 1) {
      if (node_name == list_element$j) {
        return(list_element$k)
      } else {
        if (node_name == list_element$k) {
          return(list_element$j)
        } else {
          return(NULL)
        }
      }
    } else {
      return(NULL)
    }
  }

  assignment.Function <- function(val, assignment_matrix) {
    if (val != 0) {
      assignment_matrix[which(val == assignment_matrix[,1]), 2]
    } else {
      return(0)
    }
  }

  # check.In <- function(entry, v) {
  #   (v[1] %in% names(entry) | v[2] %in% names(entry))
  # }
  #
  # identify.Cliques <- function(max_cliques_var, size, selected_edge) {
  #   cl1 <- which(unlist(lapply(max_cliques_var, FUN = length)) == size)
  #   cl2 <- which(unlist(lapply(max_cliques_var[cl1], FUN = check.In, selected_edge)))
  #   return(cl2)
  # }
  #
  # intersect.Cliques.Edge <- function(clique, edge) {
  #   intersect(names(clique), edge)
  # }
  # intersect.Cliques <- function(v, clique_list) {
  #   intersect(names(clique_list[[v[1]]]), names(clique_list[[v[2]]]))
  # }
  #
  # choose.Intersect <- function(intersect, selected_edge, i) {
  #   if(length(intersect) != i - 1) {
  #     return(NULL)
  #   } else {
  #     if (selected_edge[1] %in% intersect | selected_edge[2] %in% intersect) {
  #       return(NULL)
  #     } else {
  #       return(intersect)
  #     }
  #   }
  # }
  #
  # break.Temp.Graphs <- function(level, temp_graph) {
  #   edges_to_remove <- which(E(temp_graph)$tree_level > level)
  #   out <- delete_edges(temp_graph, E(temp_graph)[edges_to_remove])
  #   return(out)
  # }


  # ## Die funktion hier wird vermutlich noch angepasst werden müssen, dass es einfacher ist rauszufinden was genau die r-vine edge ist, d.h. mit conditioned und condinitiong set ist
  # search.Cliques.internal <- function(max_cliques_var, selected_edge) {
  #   ## only cliques of same size can be merged
  #   indices <- unlist(unique(lapply(max_cliques_var, FUN = length)))
  #   for(i in indices[order(indices, decreasing = TRUE)]) {
  #     relevant_cliques <- identify.Cliques(max_cliques_var, i, selected_edge)
  #     if (length(relevant_cliques) >= 2) {
  #       cl1 <- which(unlist(lapply(max_cliques_var, FUN = length)) == i)
  #       intersections <- lapply(max_cliques_var[cl1][relevant_cliques], FUN = intersect.Cliques.Edge, edge = selected_edge)
  #       length.intersections.larger0 <- which(lapply(intersections, FUN = length) > 0)
  #       if (length(length.intersections.larger0) >= 2) {
  #         intersects <- apply(combn(length(length.intersections.larger0), 2), 2, FUN = intersect.Cliques, clique_list = max_cliques_var[cl1][relevant_cliques][length.intersections.larger0])
  #         intersect_final <- unlist(lapply(intersects, FUN = choose.Intersect, selected_edge = selected_edge, i = i))
  #       } else {
  #         ## do nothing, next loop
  #       }
  #     } else {
  #       ## do nothing, next loop
  #     }
  #   }
  #   out <- list(j = selected_edge[1], k = selected_edge[2], D = intersect_final)
  #   return(out)
  # }

  # remove.Names <- function(x, y) { setdiff(names(x),y )}

  # check.Clique.Man <- function(clique, selected_edge) {
  #   ifelse(selected_edge[1] %in% names(clique) | selected_edge[2] %in% names(clique), TRUE, FALSE)
  # }

  # search.Cliques <- function(selected_edge, temp_graph, level) {
  #   ## identify relevant vertices
  #   ad_vertices <- adjacent_vertices(temp_graph, v = selected_edge, mode = "all")
  #   all_names <- union(selected_edge, union(names(unlist(ad_vertices[[1]])), names(unlist(ad_vertices[[2]]))))
  #   ## get cliques on this set
  #   relevant_cliques <- cliques(induced_subgraph(temp_graph, vids = V(temp_graph)[which(names(V(temp_graph)) %in% all_names)]))
  #   relevant_cliques_2 <- relevant_cliques[unlist(lapply(relevant_cliques, FUN = check.In, v = selected_edge))]
  #   joined_list_2 <- append(lapply(relevant_cliques_2, FUN = remove.Names, y = selected_edge[1]),lapply(relevant_cliques_2, FUN = remove.Names, y = selected_edge[2]))
  #   intersect <- joined_list_2[unlist(which(duplicated(joined_list_2)))]
  #   intersect_final <- intersect[which.max(unlist(lapply(intersect, FUN = length)))][[1]]
  #   out <- list(j = selected_edge[1], k = selected_edge[2], D = intersect_final)
  #   if(max(E(temp_graph)$tree_level) > length(intersect)) {
  #     return(NULL)
  #   } else {
  #     return(out)
  #   }
  # }

  # get.Poss.Triangulations <- function(graph, completion_edges) {
  #   return(apply(completion_edges, 1, FUN = search.Cliques, temp_graph = temp_graph))
  # }

  # join.Function <- function(index, edge_list, graph_list) {
  #   list(edge_list = edge_list[[index]], graph_list = graph_list[[index]])
  # }

  # eval.List.Triangulation <- function(list_element) {
  #   out <- apply(list_element$edge_list, 1, FUN = search.Cliques, temp_graph = list_element$graph_list, level = )
  #   out[lapply(out, FUN = is.null) == FALSE]
  # }
  #
  # eval.Max.Pcor.1 <- function(list_element, data, method) {
  #   if(!is.null(list_element)) {
  #     if(length(list_element) > 1) {
  #       eval_Pcor <- lapply(list_element, FUN = eval.Pcor, data = data, method = method)
  #       out <- c(which.max(abs(unlist(eval_Pcor))), max(abs(unlist(eval_Pcor))))
  #       return(out)
  #     } else {
  #       eval_Pcor <- eval.Pcor(list_element[[1]], data = data, method = method)
  #       out <- c(1, eval_Pcor)
  #       return(out)
  #     }
  #   } else {
  #     return(c(0,0))
  #   }
  # }
  ## END DEFINITION

  while (length(E(graph.difference(end_graph, temp_graph))) > 0) {
    ## all edges which can be added at that point to create higher order cliques and edgelist of previous graph, we consider all edges since sometimes no edge of the
    ## original graph can be added to find a chordal completion
    edge_difference <- graph.difference(full_graph, temp_graph)
    edgelist_difference <- get.edgelist(edge_difference)

    chordal_edges <- check.Graph.Chordal(temp_graph, edgelist_difference = edgelist_difference, end_graph = end_graph)

    if (nrow(chordal_edges) == 0) {
      ## case that no edge present in end_graph can be added for triangulation
      chordal_edges <- check.Graph.Chordal.Fallback(temp_graph, edgelist_difference = edgelist_difference, end_graph = end_graph)
    } else {
        ## do nothing
    }

    poss_triangulations <- apply(chordal_edges, 1, FUN = search.Cliques.aIo, temp_graph = temp_graph)
    if (length(poss_triangulations) == 0) {
      chordal_edges <- check.Graph.Chordal.Fallback(temp_graph, edgelist_difference = edgelist_difference, end_graph = end_graph)
      poss_triangulations <- apply(chordal_edges, 1, FUN = search.Cliques.aIo, temp_graph = temp_graph)
    }
    max_pcor <- which.max(abs(unlist(lapply(poss_triangulations, FUN = eval.Pcor, data  = data, method = method)))) # equivalent to all possible edges (i.e. "proximity condition")
    chosen_triangulation <- poss_triangulations[max_pcor][[1]]

    edgelist_for_vine[[iterator]] <- chosen_triangulation
    iterator <- iterator + 1
    selected_edge <- extract.Values(chosen_triangulation)
    temp_graph <- add_edges(temp_graph, selected_edge$edge, tree_level = selected_edge$level, label = selected_edge$level, color = colorScheme[selected_edge$level])
    ## alternativ über verschiedene edge Linienstaerken modelieren
  }

  while (length(E(graph.difference(full_graph, temp_graph))) > 0) {
    edge_difference <- graph.difference(full_graph, temp_graph)
    edgelist_difference <- get.edgelist(edge_difference)

    chordal_edges <- check.Graph.Chordal(temp_graph, edgelist_difference = edgelist_difference, end_graph = full_graph)
    poss_triangulations <- apply(chordal_edges, 1, FUN = search.Cliques.aIo, temp_graph = temp_graph)
    chosen_triangulation <- poss_triangulations[which(unlist(lapply(poss_triangulations, FUN = is.null)) == FALSE)][[1]]
    edgelist_for_vine[[iterator]] <- chosen_triangulation
    iterator <- iterator + 1
    selected_edge <- extract.Values(chosen_triangulation)
    temp_graph <- add_edges(temp_graph, selected_edge$edge, tree_level = selected_edge$level, label = selected_edge$level, color = colorScheme[selected_edge$level])
  }

  independence_graph <- graph.difference(temp_graph, end_graph)
  E(independence_graph)$independence <- rep(TRUE, length(E(independence_graph)))
  out_graph <- graph.union(temp_graph, independence_graph)

  tree_level <- E(temp_graph)$tree_level

  RVM <- matrix(NA, d, d)
  RVM[upper.tri(RVM)] <- 0
  for (i in 1:(d-1)) {
    k <- which.max(tree_level)
    RVM[i, i] <- edgelist_for_vine[[k]]$j
    entries_nonzero <- which(unlist(lapply(lapply(edgelist_for_vine, FUN = get.Edge, node_name = RVM[i, i]), FUN = length)) != 0)
    RVM[(i+1):d, i] <- unlist(lapply(edgelist_for_vine, FUN = get.Edge, node_name = RVM[i, i]))[length(entries_nonzero) : 1]
    edgelist_for_vine[entries_nonzero] <- 0
    tree_level[entries_nonzero] <- 0
  }
  RVM[d,d] <- RVM[d, d-1]

  indepMatrix <- matrix(0, d, d)
  indepMatrix[lower.tri(indepMatrix)] <- 1
  for(j in 1:(d-1)) {
    for(i in (j+1):d) {
      indepMatrix[i,j] <- ifelse(length(E(graph.intersection(graph_from_edgelist(matrix(c(RVM[j,j], RVM[i,j]), ncol = 2), directed = FALSE), independence_graph))) == 1, 0, -1)
    }
  }

  names_diag <- diag(RVM)
  diag_seq <- as.character(seq(1:length(names_diag)))
  assignment_matrix <- cbind(names_diag, diag_seq)
  colnames(assignment_matrix) <- NULL

  RVM <- apply(RVM, MARGIN = c(1,2), FUN = assignment.Function, assignment_matrix = assignment_matrix)
  RVM <- matrix(as.numeric(RVM), d, d, byrow = FALSE)

  out <- list(RVM = RVM, indepMatrix = indepMatrix, names = names_diag, end_graph = out_graph)
  .out <- out
  rm(list = ls())
  .out
}
