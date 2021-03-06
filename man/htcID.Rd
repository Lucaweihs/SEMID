% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/htcID.R
\name{htcID}
\alias{htcID}
\title{Determines which edges in a mixed graph are HTC-identifiable.}
\usage{
htcID(mixedGraph, tianDecompose = T)
}
\arguments{
\item{mixedGraph}{a \code{\link{MixedGraph}} object representing the L-SEM.}

\item{tianDecompose}{TRUE or FALSE determining whether or not the Tian
decomposition should be used before running the
current generic identification algorithm. In general
letting this be TRUE will make the algorithm faster and
more powerful.}
}
\value{
see the return value of \code{\link{generalGenericID}}.
}
\description{
Uses the half-trek criterion of Foygel, Draisma, and Drton (2012) determine
which edges in a mixed graph are generically identifiable. Depending on your
application it faster to use the \code{\link{graphID.htcID}} function
instead of this one, this function has the advantage of returning additional
information.
}
\references{
Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
generic identifiability of linear structural equation models.
\emph{Ann. Statist.} 40(3): 1682-1713.

Jin Tian. 2005. Identifying direct causal effects in linear models. In
\emph{Proceedings of the 20th national conference on Artificial intelligence
- Volume 1} (AAAI'05), Anthony Cohn (Ed.), Vol. 1. AAAI Press 346-352.
}
