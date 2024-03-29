% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LatentDigraph.R
\name{LatentDigraph}
\alias{LatentDigraph}
\title{Construct a LatentDigraph object}
\usage{
LatentDigraph(L = matrix(0,1,1),
                     observedNodes = seq(1, length = nrow(L)),
                     latentNodes = integer(0))
}
\arguments{
\item{L}{see \code{\link{graphID}} for the appropriate form of L.}

\item{observedNodes}{a vector of positive integers representing
the vertex numbers of the observed nodes. These will correspond,
in order, to the first length(observedNodes) rows of L.}

\item{latentNodes}{a vector of positive integers representing
the vertex numbers of the latent nodes. These will correspond,
in order, to the last length(latentNodes) rows of L.}
}
\value{
An object representing the LatentDigraph
}
\description{
Creates an object representing a latent factor graph. The methods that are
currently available to be used on the latent factor graph include
\enumerate{
\item numObserved
\item numLatents
\item numNodes
\item toIn
\item toEx
\item L
\item observedNodes
\item latentNodes
\item parents
\item children
\item ancestors
\item descendants
\item trFrom
\item getTrekSystem
\item inducedSubgraph
\item stronglyConnectedComponent
\item plot
\item observedParents
\item getMixedGraph
}
see the individual function documentation for more information.
}
