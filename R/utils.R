#' A helper function to validate input matrices.
#'
#' This helper function validates that the two input matrices, L and O, are of
#' the appropriate form to be interpreted by the other functions. In particular
#' they should be square matrices of 1's and 0's with all 0's along their
#' diagonals. We do not require O to be symmetric here.
#'
#' @param L See above description.
#' @param O See above description.
#'
#' @return This function has no return value.
validateMatrices <- function(L, O) {
  if (!is.matrix(L) || !is.matrix(O)) {
    stop("L and O must be matrices.")
  } else if (length(unique(c(dim(L), dim(O)))) != 1) {
    stop("L and O must both be square matrices of the same dimensions.")
  }
  t1 = all(L %in% c(0, 1))
  t2 = all(O %in% c(0, 1))
  if (!t1 || !t2) {
    stop("L and O must contain only 1's and 0's.")
  } else if (any(diag(L) != 0) || any(diag(O) != 0)) {
    stop("L and O must have 0's along their diagonals.")
  }
}

#' A general identification algorithm template
#'
#' A function that encapsulates the general structure of our algorithms for
#' testing generic identifiability. Allows for various identification algorithms
#' to be used in concert.
#'
#' @export
#'
#' @inheritParams graphID
#' @param idStepFunctions a list of identification step functions
#'
#' @return a list
generalGenericID <- function(L, O, idStepFunctions) {
  validateMatrices(L, O)
  O = 1 * ((O + t(O)) != 0)
  diag(O) = 0
  m = nrow(L)

  mixedGraph = MixedGraph(L, O)
  unsolvedParents = lapply(1:m, function(node) { mixedGraph$allParents(node) })
  solvedParents = rep(list(numeric(0)), m)
  identifier = createIdentifierBaseCase(L, O)

  changeFlag = T
  while (changeFlag) {
    for (idStepFunction in idStepFunctions) {
      idResult = idStepFunction(mixedGraph, unsolvedParents,
                                solvedParents, identifier)
      changeFlag = length(idResult$identifiedEdges) != 0
      unsolvedParents = idResult$unsolvedParents
      solvedParents = idResult$solvedParents
      identifier = idResult$identifier
      if (changeFlag) {
        break
      }
    }
  }
  return(list(solvedParents = solvedParents,
              unsolvedParents = unsolvedParents,
              identifier = identifier))
}

#' Returns all subsets of a certain size
#'
#' For an input vector x, returns in a list, the collection of all subsets
#' of x of size k.
#'
#' @param x a vector from which to get subsets
#' @param k the size of the subsets returned
#'
#' @return a list of all subsets of x of a given size k
subsetsOfSize <- function(x, k) {
  if (k > length(x)) {
    return(list())
  }
  if (k == 0) {
    return(list(numeric(0)))
  }
  if (length(x) == k) {
    return(list(x))
  }
  return(combn(x, k, simplify = F))
}

#' Create an identifier base case
#'
#' Identifiers are functions that take as input a covariance matrix Sigma
#' corresponding to some mixed graph G and, from that covariance matrix,
#' identify some subset of the coefficients in the mixed graph G. This function
#' takes as input the matrices, L and O, defining G and creates an identifer
#' that does not identify any of the coefficients of G. This is useful as a
#' base case when building more complex identification functions.
#'
#' @inheritParams graphID
#'
#' @return a function that takes as input a covariance matrix compatible with
#'         the mixed graph defined by L/O and returns a list with two
#'         named components:
#'         Lambda - a matrix equal to L but with NA values instead of 1s,
#'         Omega - a matrix equal to O but with NA values instead of 1s.
#'         When building more complex identifiers these NAs will be replaced
#'         by the value that can be identified from Sigma.
createIdentifierBaseCase <- function(L, O) {
  # Redundant assignment puts L into the environment of the, below returned,
  # function
  validateMatrices(L, O)
  L <- L
  O <- O
  return(function(Sigma) {
    Lambda <- matrix(NA, nrow(L), nrow(L))
    Lambda[L == 0] <- 0
    Omega <- matrix(NA, nrow(O), nrow(O))
    Omega[(O == 0) && !(diag(nrow(O)))] <- 0
    return(list(Lambda = Lambda, Omega = Omega))
  })
}