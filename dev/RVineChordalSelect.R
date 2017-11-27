## Prototype algorithm for estimation of chordal graph, finding the spanning trees and finding the best fitting R-vine tree sequences using these conditional independences

RVineChordalSelect <- function(data, num_trees, algo = "hc",
                               maxp = Inf,  restart = 0, perturb = 1,
                               tabu = 10, max.tabu = 10,
                               alpha = 0.05,
                               seed = 81542,
                               familyset = NA, selectioncrit = "AIC", indeptest = FALSE, level = 0.05, trunclevel = NA, rotations = TRUE,
                               cores = 1,
                               filename = NULL) {
  require(VCRef)
  require(bnlearn)
  require(igraph)
  require(gRbase)
  require(parallel)
  require(foreach)
  require(doParallel)

  d <- ncol(data)
  if (any(data < 0) | any(data > 1)) {
    data_u <- as.data.frame(pobs(data))
  } else {
    data_u <- data
  }
  data_normal <- as.data.frame(apply(data_u, 2, qnorm))
  if (restart > 0) {
    set.seed(seed)
  }
  dag_bnlearn <- NULL

  calc.Amat <- function(graph, data_names){
    h_amat <- matrix(0, length(data_names), length(data_names))
    colnames(h_amat) <- data_names
    rownames(h_amat) <- data_names
    h_amat[graph$arcs] <- 1
    h_amat
  }

  if (algo == "hc") {
    dag_bnlearn <- hc(data_normal, maxp = maxp, restart = restart, perturb = perturb)
    moral_graph_bnlearn <- moral(dag_bnlearn) # moral.graph is chordal
  } else {
    if (algo == "gs") {
      moral_graph_bnlearn <- gs(data_normal, alpha = alpha, undirected = TRUE)
    } else {
      if (algo == "iamb") {
        moral_graph_bnlearn <- iamb(data_normal, alpha = alpha, undirected = TRUE)
      } else {
        if (algo == "fast.iamb") {
          moral_graph_bnlearn <- fast.iamb(data_normal, alpha = alpha, undirected = TRUE)
        } else {
          if (algo == "inter.iamb") {
            moral_graph_bnlearn <- inter.iamb(data_normal, alpha = alpha, undirected = TRUE)
          } else {
            if (algo == "mmpc") {
              moral_graph_bnlearn <- mmpc(data_normal, alpha = alpha, undirected = TRUE)
            } else {
              if (algo == "si.hiton.pc") {
                moral_graph_bnlearn <- si.hiton.pc(data_normal, alpha = alpha, undirected = TRUE)
              } else {
                if (algo == "tabu") {
                  dag_bnlearn <- tabu(data_normal, tabu = tabu, max.tabu = max.tabu)
                  moral_graph_bnlearn <- moral(dag_bnlearn) # moral.graph is chordal
                } else {
                  if (algo == "rsmax2") {
                    dag_bnlearn <- rsmax2(data_normal, alpha = alpha)
                    moral_graph_bnlearn <- moral(dag_bnlearn) # moral.graph is chordal
                  } else {
                    if ( algo == "mmhc") {
                      dag_bnlearn <- mmhc(data_normal, alpha = alpha, restart = restart, perturb = perturb)
                      moral_graph_bnlearn <- moral(dag_bnlearn) # moral.graph is chordal
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  moral_graph_amat <- calc.Amat(moral_graph_bnlearn, names(data)) # calculate adjacency matrix
  moral_graph_igraph <- graph_from_adjacency_matrix(moral_graph_amat, mode = "undirected") # convert object to igraph
  chordality <- is.chordal(moral_graph_igraph, fillin = TRUE, newgraph = TRUE)
  if(!chordality$chordal) {
    moral_graph_igraph_tri <- igraph.from.graphNEL(minimalTriang(igraph.to.graphNEL(moral_graph_igraph)))
  } else {
    moral_graph_igraph_tri <- moral_graph_igraph
  }
  E(moral_graph_igraph_tri)$weight <- abs(cor(data_normal))[get.edgelist(moral_graph_igraph_tri)] # assign edge weights according to Kendall's tau
  if(d >= 15 & cores > 8) {
    MST_object <- findKBestTreesP(moral_graph_igraph_tri, k = num_trees, searchMode = "max", cores = cores) # calculate num.trees best MSTs

    filename_new <- paste(filename, "MST", sep = "")
    filename_save <- paste(paste("output/", filename_new, sep = ""), ".rdata", sep = "")
    save(MST_object, file = filename_save)

  } else {
    MST_object <- findKBestTrees(moral_graph_igraph_tri, k = num_trees, searchMode = "max") # calculate num.trees best MSTs
  }

  # min_max_pathlength <- maxPathOnTree(MST_object$MST, moral_graph_igraph_tri) ## diese Logik taugt evtl. nicht... zumindest nicht rein für den ersten Baum

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

  calc.Vine.unpar <- function(MST, end_graph, data, data_u, familyset, selectioncrit, indeptest, level, trunclevel, rotations) {
    graph_out <- calculateChordalCompletion(start_graph = MST,
                                            end_graph = moral_graph_igraph_tri,
                                            data = data)
    RVM <- RVineCopSelect2(data = data_u[graph_out$names],
                           familyset = NA,
                           Matrix = graph_out$RVM,
                           familyMatrix = graph_out$indepMatrix,
                           selectioncrit = selectioncrit,
                           indeptest = indeptest,
                           level = level,
                           trunclevel = trunclevel,
                           rotations = TRUE,
                           cores = cores)
    .out <- list(RVM = RVM, graph = graph_out$end_graph)
    rm(list = ls())
    .out
  }

  calc.Vine.par <- function(list_element, data_u, familyset, selectioncrit, indeptest, level, trunclevel, rotations, cores) {
    RVineCopSelect2(data = data_u[list_element$names],
                    familyset = NA,
                    Matrix = list_element$RVM,
                    familyMatrix = list_element$indepMatrix,
                    selectioncrit = selectioncrit,
                    indeptest = indeptest,
                    level = level,
                    trunclevel = trunclevel,
                    rotations = TRUE,
                    cores = cores)
  }

  if (cores > 1) {
    if (cores > 15) {
      graph_out <- foreach(i = 1:length(MST_object$MST),
                           .export = c("calculateChordalCompletion"),
                           .packages = c("igraph", "ppcor", "RColorBrewer", "VCRef")) %dopar% {
                             out <- calculateChordalCompletion(start_graph = MST_object$MST[[i]],
                                                        end_graph = moral_graph_igraph_tri,
                                                        data = data)
                             print(i)
                             out
                           }
      filename_new <- paste(filename, "RVM", sep = "")
      filename_save <- paste(paste("output/", filename_new, sep = ""), ".rdata", sep = "")
      save(graph_out, file = filename_save)
      vine_out <- lapply(graph_out, FUN = calc.Vine.par, data_u = data_u, familyset = familyset, selectioncrit = selectioncrit, indeptest = indeptest, level = level,
                         trunclevel = trunclevel, rotations = rotations, cores = cores)
    } else {
      vine_out <- foreach(i = 1:length(MST_object$MST),
                           .export = c("calculateChordalCompletion"),
                           .packages = c("igraph", "ppcor", "RColorBrewer", "VCRef")) %dopar% {
                             graph_out <- calculateChordalCompletion(start_graph = MST_object$MST[[i]],
                                                                     end_graph = moral_graph_igraph_tri,
                                                                     data = data)
                             RVM <- RVineCopSelect2(data = data_u[graph_out$names],
                                                    familyset = familyset,
                                                    Matrix = graph_out$RVM,
                                                    familyMatrix = graph_out$indepMatrix,
                                                    selectioncrit = selectioncrit,
                                                    indeptest = indeptest,
                                                    level = level,
                                                    trunclevel = trunclevel,
                                                    rotations = TRUE,
                                                    cores = 1)
                             out <- list(RVM = RVM, graph = graph_out$end_graph)
                           }
    }
  } else {
    vine_out <- lapply(MST_object$MST, FUN = calc.Vine.unpar, end_graph = moral_graph_igraph_tri, data = data, data_u = data_u, familyset = familyset,
                       selectioncrit = selectioncrit, indeptest = indeptest, level = level, trunclevel = trunclevel, rotations = rotations)
  }
  graph_models <- list(dag_bnlearn = dag_bnlearn,
                       moral_graph_bnlearn = moral_graph_bnlearn,
                       moral_graph_igraph = moral_graph_igraph,
                       moral_graph_igraph_tri = moral_graph_igraph_tri,
                       MST = MST_object$MST)
  .out <- list(vine_models = vine_out, graph_models = graph_models)
  rm(list = ls())
  .out
}
