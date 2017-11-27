\alias{RVineClusterSelect}
\name{RVineClusterSelect}
\title{Estimation of Vine Copulas based on a clustered undirected graphical model}
\usage{
RVineClusterSelect(data, graph, max_size, fill_level, ...)
}
\arguments{
\item{data}{A \code{data.frame} or \code{matrix} containing copula data (i.e. data with uniform margins on [0,1]).}
\item{graph}{An adjacency matrix of an undirected graph, i.e. a symmetric matrix. If the graphical model has multiple connected components, the R-vine is estimated on the connected components. If it is a connected graph, the algorithm \code{RVineGraphSelect} is called.}
\item{max_size}{Threshold dimension for clustered graphical model.}
\item{fill_level}{If a graph with several connected components is employed, the first \code{fill_level} R-vine trees outside the connected components are estimated with pair copulas as well.}
\item{...}{Additional arguments as for families or estimation method of pair copula families, to be passed on to \code{VineCopula} package used for pair copula estimation.}
}
\description{
The function \code{RVineClusterSelect} estimates a vine copula model based on a clustered undirected graphical model which needs to be provided or is internally computed given some threshold dimension.
}
\examples{
library(VineCopulaHD)
library(VineCopula)

data(daxreturns) ## load 15-dimensional stock return data attached to the VineCopula package
head(daxreturns)
model_cluster <- RVineClusterSelect(data = daxreturns, max_size = 7, fill_level = 2) ## fit an R-vine based on threshold dimension 7 with fill level 2
model_cluster$RVM ## the R-vine object
model_cluster$graph ## the underlying graphical model
}
\author{
Dominik Mueller
}
