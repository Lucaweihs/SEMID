% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotLatentDigraph}
\alias{plotLatentDigraph}
\title{Plot a latent factor graph}
\usage{
plotLatentDigraph(L, observedNodes, latentNodes, main = "")
}
\arguments{
\item{L}{see \code{\link{graphID}} for the appropriate form of L.}

\item{observedNodes}{a vector of positive integers representing
the vertex numbers of the observed nodes. These will correspond,
in order, to the first length(observedNodes) rows of L.}

\item{latentNodes}{a vector of positive integers representing
the vertex numbers of the latent nodes. These will correspond,
in order, to the last length(latentNodes) rows of L.}

\item{main}{the plot title.}
}
\value{
An object representing the LatentDigraph
}
\description{
Given an adjacency matrix representing the directed edges in a latent
factor graph, plots a representation of the graph. The latent nodes
should come last in L and the vertex labels should only be given for the
observed nodes.
}
