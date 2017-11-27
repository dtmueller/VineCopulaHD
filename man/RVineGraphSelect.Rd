\alias{RVineGraphSelect}
\name{RVineGraphSelect}
\title{Estimation of Vine Copulas based on an undirected graphical model}
\usage{
RVineGraphSelect(data, graph, ...)
}
\arguments{
\item{data}{A \code{data.frame} or \code{matrix} containing copula data (i.e. data with uniform margins on [0,1]).}
\item{graph}{An adjacency matrix of an undirected graph, i.e. a symmetric matrix.}
\item{...}{Additional arguments as for families or estimation method of pair copula families, to be passed on to \code{VineCopula} package used for pair copula estimation.}
}
\description{
The function \code{RVineGraphsSelect} estimates a vine copula model based on a undirected graphical model which needs to be provided.
}
\examples{
library(VineCopulaHD)
library(VineCopula)

data(daxreturns) ## load 15-dimensional stock return data attached to the VineCopula package
head(daxreturns)
graph <- matrix(0, ncol = ncol(daxreturns), nrow = ncol(daxreturns)) ## define raw matrix
graph[upper.tri(graph, diag = FALSE)] <- rbinom(choose(ncol(daxreturns),2), 1, 0.5)  ## assign random graph structure
graph <- graph + t(graph) ## generate symmetric matrix
model_graph <- RVineGraphSelect(data = daxreturns, graph = graph) ## fit an R-vine based on this undirected graphical model
model_graph ## the R-vine object
}
\author{
Dominik Mueller
}
