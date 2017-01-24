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

createSimpleBiDirIdentifier <- function(idFunc) {
  idFunc <- idFunc
  return(
    function(Sigma) {
      m <- nrow(Sigma)
      identifiedParams <- idFunc(Sigma)
      Lambda <- identifiedParams$Lambda
      if (!any(is.na(Lambda))) {
        Omega = t(diag(m) - Lambda) %*% Sigma %*% (diag(m) - Lambda)
        return(list(Lambda = Lambda, Omega = Omega))
      }
      return(list(Lambda = Lambda, Omega = identifiedParams$Omega))
    }
  )
}

#' Globally identify the covariance matrix of a C-component
#'
#' The Tian decomposition of a mixed graph G allows one to globally identify
#' the covariance matrices Sigma' of special subgraphs of G called C-components.
#' This function takes the covariance matrix Sigma corresponding to G and
#' a collection of node sets which specify the C-component, and returns the
#' Sigma' corresponding to the C-component.
#'
#' @param Sigma the covariance matrix for the mixed graph G
#' @param internal an integer vector corresponding to the vertices of the
#'        C-component that are in the bidirected equivalence classes (if the
#'        graph is not-acyclic then these equivalence classes must be enlarged
#'        by combining two bidirected components if there are two vertices, one
#'        in each component, that are simultaneously on the same directed cycle).
#' @param incoming the parents of vertices in internal that are not in the set
#'        internal themselves
#' @param topOrder a topological ordering of c(internal, incoming) with respect
#'        to the graph G. For vertices in a strongly connected component the
#'        ordering is allowed to be arbitrary.
#'
#' @export
#'
#' @return This function has no return value.
tianSigmaForComponent <- function(Sigma, internal, incoming, topOrder) {
  if (length(incoming) == 0) {
    return(Sigma[topOrder, topOrder, drop = F])
  }
  newSigmaInv = matrix(0, length(topOrder), length(topOrder))
  for (j in 1:length(topOrder)) {
    node = topOrder[j]
    if (node %in% internal) {
      if (j == 1) {
        newSigmaInv[j, j] = newSigmaInv[j,j] + 1 / Sigma[node, node,drop = F]
      } else {
        inds = topOrder[1:(j - 1)]

        SigmaIndsInv = solve(Sigma[inds, inds, drop = F])
        schurInv = solve(Sigma[node,node,drop = F] -
                           Sigma[node, inds, drop = F] %*% SigmaIndsInv %*%
                           Sigma[inds, node, drop = F])
        newSigmaInv[j, j] = newSigmaInv[j,j] + schurInv
        meanMat = as.numeric(solve(Sigma[inds,inds]) %*% Sigma[inds, node] %*% schurInv)
        newSigmaInv[j, 1:(j-1)] = newSigmaInv[j, 1:(j-1)] - meanMat
        newSigmaInv[1:(j-1), j] = newSigmaInv[j, 1:(j-1)]
        newSigmaInv[1:(j-1), 1:(j-1)] = newSigmaInv[1:(j-1), 1:(j-1)] + SigmaIndsInv %*% Sigma[inds, node] %*% meanMat
      }
    }
  }
  newSigmaInv[topOrder %in% incoming, topOrder %in% incoming] =
    newSigmaInv[topOrder %in% incoming, topOrder %in% incoming] + diag(length(incoming))
  newSigma = solve(newSigmaInv)
  return(newSigma)
}

#' Identifies components in a tian decomposition
#'
#' @export
#' @return This function has no return value.
tianIdentifier <- function(idFuncs, cComponents) {
  idFuncs <- idFuncs
  cComponents <- cComponents
  return(
    function(Sigma) {
      Lambda = matrix(NA, ncol(Sigma), ncol(Sigma))
      Omega = matrix(NA, ncol(Sigma), ncol(Sigma))

      for (i in 1:length(cComponents)) {
        internal = cComponents[[i]]$internal
        incoming = cComponents[[i]]$incoming
        topOrder = cComponents[[i]]$topOrder

        newSigma = tianSigmaForComponent(Sigma, internal, incoming, topOrder)

        result = idFuncs[[i]](newSigma)
        internalInds = which(topOrder %in% internal)
        Lambda[topOrder, internal] = result$Lambda[,internalInds]
        Lambda[-topOrder, internal] = 0
        Omega[internal, internal] = result$Omega[internalInds, internalInds]
        Omega[-internal, internal] = 0
        Omega[internal, -internal] = 0
      }
      return(list(Lambda = Lambda, Omega = Omega))
    }
  )
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
#' @param tianDecompose True if the mixed graph be decomposed into Tian's
#'        C-components before running the identification algorithms.
#'
#' @return a list
generalGenericID <- function(L, O, idStepFunctions, tianDecompose = T) {
  O = 1 * ((O + t(O)) != 0)
  m = nrow(L)

  mixedGraph = MixedGraph(L, O)
  unsolvedParents = lapply(1:m, function(node) { mixedGraph$allParents(node) })
  solvedParents = rep(list(numeric(0)), m)

  if (!tianDecompose) {
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
    if (length(unlist(unsolvedParents)) == 0) {
      identifier = createSimpleBiDirIdentifier(identifier)
      solvedSiblings = lapply(1:m, FUN = function(x) { mixedGraph$allSiblings(x) })
      unsolvedSiblings = rep(list(integer(0)), m)
    } else {
      solvedSiblings = rep(list(integer(0)), m)
      unsolvedSiblings = lapply(1:m, FUN = function(x) { mixedGraph$allSiblings(x) })
    }
  } else {
    solvedSiblings = rep(list(integer(0)), m)
    cComps = mixedGraph$tianDecompose()

    compResults = vector("list", length(cComps))
    identifiers = vector("list", length(cComps))

    for (i in 1:length(cComps)) {
      result = generalGenericID(cComps[[i]]$L, cComps[[i]]$O, idStepFunctions, tianDecompose = F)
      topOrder = cComps[[i]]$topOrder
      compResults[[i]] = result
      for (j in 1:length(topOrder)) {
        solvedParents[[topOrder[j]]] = c(solvedParents[[topOrder[j]]], topOrder[result$solvedParents[[j]]])
        solvedSiblings[[topOrder[j]]] = c(solvedSiblings[[topOrder[j]]], topOrder[result$solvedSiblings[[j]]])
      }

      identifiers[[i]] = result$identifier
    }

    unsolvedParents = lapply(1:m, FUN = function(x) { setdiff(mixedGraph$allParents(x), solvedParents[[x]]) })
    unsolvedSiblings = lapply(1:m, FUN = function(x) { setdiff(mixedGraph$allSiblings(x), solvedSiblings[[x]]) })
    identifier = tianIdentifier(identifiers, cComps)
  }
  return(list(solvedParents = solvedParents,
              unsolvedParents = unsolvedParents,
              solvedSiblings = solvedSiblings,
              unsolvedSiblings = unsolvedSiblings,
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
