\alias{checkNecCondDAG}
\name{checkNecCondDAG}
\title{Check whether on some given DAG the necessary conditions hold and compute an approximating R-vine}
\usage{checkNecCond3(DAG)}
\arguments{
  \item{DAG}{A DAG of class \code{bn}.}
}
\description{
  The function \code{checkNecCond3} determines whether a Gaussian DAG can be represented by a regular vine and computes an R-vine representation.
}
\examples{
  This function needs the graph package (only on Bioconductor), if you do not have it, run
  source("https://bioconductor.org/biocLite.R")
  biocLite("graph")

  library(bnlearn)
  library(graph)

  DAG <- hc(daxreturns, maxp = 2) ## fit a DAG with at most 2 parents on the daxreturns dataset
  out <- checkNecCondDAG(DAG) ## check necessary conditions for A1 and use presence of v-structure as proxy for A2
  out$warning_message ## status message of check for necessary conditions

  out$RVM$RVM ## R-vine matrix of the R-vine approximation
  out$RVM$IndepMatrix ## Independence matrix of the R-vine approximation
  out$RVM$k_prime ## empirical truncation level
  out$RVM$ratio ## percentage of independence pair copulas in the R-vine

}
\author{
  Dominik Mueller
}
