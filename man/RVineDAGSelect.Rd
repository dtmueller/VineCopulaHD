\alias{RVineDAGSelect}
\name{RVineDAGSelect}
\title{Estimation of Vine Copulas based on DAGs}
\usage{
RVineDAGSelect(data, maxp = k, ...)
}
\arguments{
\item{data}{A \code{data.frame} or \code{matrix} containing copula data (i.e. data with uniform margins on [0,1]).}
\item{k}{Integer specifying the maximal number of parents for the DAG.}
\item{...}{Additional arguments as for families or estimation method of pair copula families, to be passed on to \code{VineCopula} package used for pair copula estimation.}
}
\description{
The function \code{RVineDAGSelect} estimates a vine copula model based on a Gaussian DAG (internally computed).
}
\examples{
library(VineCopulaHD)
library(VineCopula)

data(daxreturns) ## load some 15-dimensional stock return data attached to the VineCopula package
head(daxreturns)
daxreturns_1 <- RVineDAGSelect(data = daxreturns, maxp = 1) ## fit an R-vine based on a DAG with at most 1 parent
daxreturns_1$RVM ## the R-vine object
daxreturns_1$kdags ## the DAG object used for estimation (Markov tree)
daxreturns_5 <- RVineDAGSelect(data = daxreturns, maxp = 5, cores = 5) ## fit an R-vine based on a DAG with at most 5 parents and in parallel
daxreturns_5$kdags ## the DAGs used to find the R-vine structure
daxreturns_5$dags$h ## the matrix H as described in the paper for giving weights for the R-vine trees based on the union of the sekeletons

data(EuroStoxx50) ## load the EuroStoxx50 data as used in Section 6 of the paper
head(stoxx)
stoxx_1 <- RVineDAGSelect(data = stoxx, maxp = 1, cores = 5, method = "itau") ## fit an R-vine based on a DAG with at most 1 parent and in parallel on the EuroStoxx data
stoxx_1$RVM ## only the first tree contains non independence copulas since we represent a Markov Tree

}
\author{
Dominik Mueller
}
