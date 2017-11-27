\alias{RVineLassoSelect}
\name{RVineLassoSelect}
\title{Estimation of Vine Copulas based on structural equation models (SEMs) and the Lasso}
\usage{
  RVineLassoSelect(data, threshold_mode, threshold, ...)
}
\arguments{
  \item{data}{A \code{data.frame} or \code{matrix} containing copula data (i.e. data with uniform margins on [0,1]).}
  \item{threshold mode}{\code{ST} for single thresholding (default) or else for adaptive thresholding.}
  \item{lambda_T}{\code{lambda_T} for single thresholding (default \code{lambda_T = 0.25^4})}
  \item{mu}{\code{mu} for adapative thresholding (default \code{mu = 0.1})}
  \item{...}{Additional arguments as for families or estimation method of pair copula families, to be passed on to \code{VineCopula} package used for pair copula estimation.}

\description{
  The function \code{RVineLassoSelect} estimates a vine copula model based on a Lasso structure selection and thresholding.
}
\examples{
  library(VineCopulaHD)
  library(VineCopula)

  data(daxreturns) ## load some 15-dimensional stock return data attached to the VineCopula package
  head(daxreturns)
  model_Lasso <- RVineLassoSelect(data = daxreturns, threshold_model = graph = graph) ## fit an R-vine based on this undirected graphical model
  model_Lasso$model ## the R-vine object
  model_Lasso$M_L ## the regularization path matrix
  model_Lasso$indep ## the independence matrix
}
\author{
  Dominik Mueller
}
