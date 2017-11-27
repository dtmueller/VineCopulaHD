pRunChordal <- function(data, num_trees, algo, maxp = Inf, restart, perturb, tabu, max.tabu, alpha = 0.05, seed,
                        familyset = NA, selectioncrit = "AIC", indeptest = FALSE, indeptestlevel = 0.05, trunclevel = NA, rotations = TRUE,
                        filename, cores = 1) {

  VineTime <- system.time(model <- RVineChordalSelect(data = data,
                                                      num_trees = 1,
                                                      algo = "hc",
                                                      maxp = Inf,
                                                      restart = 0,
                                                      perturb = 1,
                                                      tabu = 10,
                                                      max.tabu = 10,
                                                      alpha = 0.05,
                                                      seed = 81542,
                                                      familyset = familyset,
                                                      selectioncrit = selectioncrit,
                                                      indeptest = indeptest,
                                                      level = indeptestlevel,
                                                      trunclevel = trunclevel,
                                                      cores = cores))[3]

  stats <- rep(NA, 17)
  d <- ncol(data)
  n.pc <- choose(d, 2)

  if(any(data < 0 | data > 1)) {
    data_u <- pobs(data)
  } else {
    data_u <- data
  }

  ## information with respect to number of parameters and model complexity
  stats[1] <- maxp # number of parents
  stats[2] <- trunclevel # trunclevel
  stats[3] <- length(which(model$vine_model[[1]]$par != 0)) + length(which(model$vine_model[[1]]$par2 != 0)) # number of parameters
  stats[4] <- length(which(model$vine_model[[1]]$family > 0)) # number of non independence copulas
  stats[5] <- n.pc - stats[4] # number of independence copulas
  stats[6] <- length(which(model$vine_model[[1]]$family == 1)) # Gaussian pair copulas
  stats[7] <- length(which(model$vine_model[[1]]$family > 1)) # non Gaussian pair copulas
  stats[8] <- NA # dummy setting for k_hat
  stats[9] <- (d + 1) - min(which(apply(model$vine_model[[1]]$family, 1, FUN = sum) > 0)) # actual truncation level of R-vine
  ## information with respect to goodness of fit measures
  stats[10] <- model$vine_model[[1]]$logLik # logLik
  stats[11] <- model$vine_model[[1]]$AIC # AIC
  stats[12] <- model$vine_model[[1]]$BIC # BIC
  ## information with respect to DAGs and graphical models
  stats[13] <- score(x = model$graph_model$dag_bnlearn, data = as.data.frame(apply(data_u, 2, qnorm)), type = "loglik-g")
  stats[14] <- nrow(model$graph_model$dag_bnlearn$arcs)
  stats[15] <- length(E(model$graph_model$moral_graph_igraph))
  stats[16] <- length(E(model$graph_model$moral_graph_igraph_tri))
  stats[17] <- VineTime

  name_model <- paste("model", ".rdata", sep = "")
  name_result <- paste("result", ".rdata", sep = "")
  filename_model <- paste("output/", paste(filename, name_model, sep = "_"), sep = "")
  filename_result <- paste("output/", paste(filename, name_result, sep = "_"), sep = "")
  save(model, file = filename_model)
  save(stats, file = filename_result)
  stats
}
