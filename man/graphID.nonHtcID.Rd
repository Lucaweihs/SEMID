% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/htcID.R
\name{graphID.nonHtcID}
\alias{graphID.nonHtcID}
\title{Check for generic infinite-to-one via the half-trek criterion.}
\usage{
graphID.nonHtcID(L, O)
}
\arguments{
\item{L}{Adjacency matrix for the directed part of the path
diagram/mixed graph; an edge pointing from i to j is encoded as L[i,j]=1 and
the lack of an edge between i and j is encoded as L[i,j]=0. There should be
no directed self loops, i.e. no i such that L[i,i]=1.}

\item{O}{Adjacency matrix for the bidirected part of the path diagram/mixed
graph. Edges are encoded as for the L parameter. Again there should be no
self loops. Also this matrix will be coerced to be symmetric so it is only
necessary to specify an edge once, i.e. if O[i,j]=1 you may, but are not
required to, also have O[j,i]=1.}
}
\value{
TRUE if the graph could be determined to be generically
non-identifiable, FALSE if this test was inconclusive.
}
\description{
Checks if a mixed graph is infinite-to-one using the half-trek criterion
presented by Foygel, Draisma, and Drton (2012).
}
\references{
Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
generic identifiability of linear structural equation models.
\emph{Ann. Statist.} 40(3): 1682-1713.
}
