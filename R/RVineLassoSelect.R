## Wrapper function RVineLassoSelect for selection of an R-vine structure using structural equation models and the Lasso

RVineLassoSelect <- function(data, order_selection = "Lasso", alpha = 1, cores_structure = 1, logging = FALSE, filename = "", threshold_mode = "ST", lambda_T = 0.25^4, mu = 0.1,
                             familyset = NA, type = 0, selectioncrit = "AIC", indeptest = FALSE,
                             level = 0.05, trunclevel = NA, progress = FALSE,  weights = NA,
                             treecrit = "AIC", se = FALSE, rotations = TRUE, method = "itau", cores_vine = 1) {

  if (any(data < 0) | any(data > 1)) {
    data_u <- pobs(data)
    data_n <- as.data.frame(apply(data_u, 2, qnorm))
  } else {
    data_u <- data
    data_n <- as.data.frame(apply(data_u, 2, qnorm))
  }

  structure <- LassoSelect(data_n, order_selection = order_selection, alpha = alpha, cores = cores_structure, logging = logging, filename = filename)
  if (threshold_mode != "ST") {
    indep <- getStructureLambdaK(structure$M_L, k_percentage = mu, k_bound = k_bound, type = "floor")
  } else {
    indep <- getStructureLambda(structure$M_L, lambda = lambda_T)
  }

  model <- RVineCopSelectSparse(data_u,
                                familyset = familyset,
                                Matrix = structure$RVM,
                                familyMatrix = indep * (-1),
                                familyList = NULL,
                                selectioncrit = selectioncrit,
                                indeptest = indeptest,
                                level = level,
                                trunclevel = trunclevel,
                                method = method,
                                cores = cores_vine)

  out <- list()
  out$model <- model
  out$M_L <- structure$M_L
  out$indep <- indep

  out
}
