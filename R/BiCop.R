## Function taken from VineCopula package and adapted

# BiCop <- function(family, par, par2 = 0, tau = NULL, check.pars = TRUE) {
#     ## use tau to construct object (if provided)
#     if (!is.null(tau))
#         par <- BiCopTau2Par(family, tau)
#     if (missing(par) & (family == 0))
#         par <- 0
#     stopifnot(is.logical(check.pars))
#     if (length(c(family, par, par2)) != 3)
#         stop("family, par, and par2 have to be a single number.")
#
#     ## family/parameter consistency checks
#     if (check.pars) {
#         # check for consistency
#         BiCopCheck(family, par, par2, call = match.call()[1])
#         # warn if par2 is unused
#         if ((family %in% allfams[onepar]) && (par2 != 0)) {
#             warning("The ",
#                     BiCopName(family, short = FALSE),
#                     " copula has only one parameter; 'par2' is useless.")
#         }
#     }
#
#     # calculate dependence measures
#     if (family == 0) {
#         tau <- 0
#         beta <- 0
#         taildep <- list()
#         taildep$upper <- 0
#         taildep$lower <- 0
#     } else {
#         tau <- BiCopPar2Tau(family, par, par2, check.pars = FALSE)
#         taildep <- BiCopPar2TailDep(family, par, par2, check.pars = FALSE)
#         beta <- NA  # beta does not work for t copula (cdf disabled)
#         if (family != 2)
#             beta <- BiCopPar2Beta(family, par, par2, check.pars = FALSE)
#     }
#
#     ## get full family name and calculate number of parameters
#     familyname <- BiCopName(family, short = FALSE)
#     npars <- if (family == 0) 0 else ifelse(family %in% allfams[onepar], 1, 2)
#
#     ## return BiCop object
#     out <- list(family     = family,
#                 par        = par,
#                 par2       = par2,
#                 npars      = npars,
#                 familyname = familyname,
#                 tau        = tau,
#                 beta       = beta,
#                 taildep    = taildep,
#                 call       = match.call())
#     class(out) <- "BiCop"
#     out
# }

## sets of families
allfams <- c(0:10,
             13, 14, 16:20,
             23, 24, 26:30, 33, 34, 36:40,
             104, 114, 124, 134, 204, 214, 224, 234)
tawns <- which(allfams > 100)
onepar <- setdiff(which(allfams %% 10 %in% c(1, 3, 4, 5, 6)), tawns)
twopar <- seq_along(allfams)[-c(1, onepar)]
negfams <- c(1, 2, 5, 23, 24, 26:30, 33, 34, 36:40, 124, 134, 224, 234)
posfams <- c(1:10, 13, 14, 16:20, 104, 114, 204, 214)
# families with more dependence near (u, v) = (1, 1)
fams11 <- c(13, 4, 6, 7, 17, 8, 9, 19, 10, 104, 204)
# families with more dependence near (u, v) = (0, 0)
fams00 <- c(3, 14, 16, 7, 17, 18, 9, 19, 20, 114, 214)
# families with more dependence near (u, v) = (1, 0)
fams10 <- c(23, 34, 36, 27, 37, 38, 29, 39, 40, 134, 234)
# families with more dependence near (u, v) = (0, 1)
fams01 <- c(33, 24, 26, 27, 37, 28, 29, 39, 30, 124, 224)

print.BiCop <- function(x, ...) {
    cat("Bivariate copula: ")
    cat(x$familyname, " (par = ", round(x$par, 2), sep = "")
    if (x$family %in% allfams[twopar])
        cat(", par2 = ", round(x$par2, 2), sep = "")
    cat(", tau = ", round(x$tau, 2), sep = "")
    cat(") \n")

    ## return BiCop object invsibly
    invisible(x)
}


summary.BiCop <- function(object, ...) {
    ## print family name
    cat("Family\n")
    cat("------ \n")
    cat("No:   ", object$family)
    cat("\n")
    cat("Name: ", object$familyname)
    cat("\n")
    cat("\n")

    ## print parameters and standard errors
    cat("Parameter(s)\n")
    cat("------------\n")
    cat("par: ", as.character(round(object$par, 2)))
    if (!is.null(object$se))
        cat("  (SE = ", as.character(round(object$se[1], 2)), ")", sep = "")
    cat("\n")
    if (object$family %in% allfams[twopar]) {
        cat("par2:", as.character(round(object$par2, 2)))
        if (!is.null(object$se))
            cat("  (SE = ", as.character(round(object$se2, 2)), ")", sep = "")
    }
    cat("\n")

    ## show dependence measures
    #     object$rho <- BiCopPar2Rho(object)
    cat("Dependence measures\n")
    cat("-------------------\n")
    cat("Kendall's tau:   ", as.character(round(object$tau, 2)))
    if (!is.null(object$emptau)) {
        p <- object$p.value.indeptest
        cat(" (empirical = ",
            as.character(round(object$emptau, 2)),
            ", ",
            "p value ",
            ifelse(p < 0.01,
                   "< 0.01",
                   paste0("= ", as.character(round(p, 2)))),
            ")",
            sep = "")
    }
    cat("\n")
    cat("Upper TD:        ", as.character(round(object$taildep$upper, 2)), "\n")
    cat("Lower TD:        ", as.character(round(object$taildep$lower, 2)), "\n")
#     cat("Blomqvist's beta:", as.character(round(object$beta, 2)), "\n")
    cat("\n")

    ## print fit statistics if available
    if (!is.null(object$nobs)) {
        cat("Fit statistics\n")
        cat("--------------\n")
        cat("logLik: ", as.character(round(object$logLik, 2)), "\n")
        cat("AIC:   ", as.character(round(object$AIC, 2)), "\n")
        cat("BIC:   ", as.character(round(object$BIC, 2)), "\n")
        cat("\n")

    }

    ## return BiCop object invsibly
    invisible(object)
}
