% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{validateMatrices}
\alias{validateMatrices}
\title{A helper function to validate input matrices.}
\usage{
validateMatrices(L, O)
}
\arguments{
\item{L}{See above description.}

\item{O}{See above description.}
}
\value{
This function has no return value.
}
\description{
This helper function validates that the two input matrices, L and O, are of
the appropriate form to be interpreted by the other functions. In particular
they should be square matrices of 1's and 0's with all 0's along their
diagonals. We do not require O to be symmetric here.
}
