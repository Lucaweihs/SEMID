createIdentifierBaseCase <- function(L) {
  return(function(Sigma) {
    Lambda <- matrix(NA, nrow(L), nrow(L))
    Lambda[L == 0] = 0
    return(Lambda)
  })
}

createLinearIdentifierFunc <- function(idFunc, sources, targets, node, sourceHtrParents) {
  # These assignments may seem redundent but they are necessary as they
  # assign these variables to the local environment of the function call.
  # This allows them to persist and still be usable by the returned function.
  idFunc <- idFunc
  sources <- sources
  targets <- targets
  node <- node
  sourceHtrParents <- sourceHtrParents
  return(
    function(Sigma) {
      m <- nrow(Sigma)
      Lambda <- idFunc(Sigma)

      SigmaMinus = Sigma
      for (sourceInd in 1:length(sources)) {
        source = sources[sourceInd]
        parents = sourceHtrParents[[sourceInd]]
        if (length(parents) != 0) {
          SigmaMinus[source,] <- Sigma[source, , drop = F] -
            t(Lambda[parents, source, drop = F]) %*% Sigma[parents,, drop = F]
        }
      }

      if (abs(det(SigmaMinus[sources, targets, drop = F])) < 10^-10) {
        stop("In identification, found near-singular system. Is the input matrix generic?")
      }

      Lambda[targets, node] <-
        solve(SigmaMinus[sources, targets, drop = F],
              SigmaMinus[sources, node, drop = F])

      return(Lambda)
    }
  )
}


#' A helper function for edgewiseID that does one step through all the nodes
#' and tries to identify new edge coefficients.
#'
#' @return a list
#' @export
linearIdentifyStep = function(mixedGraph, unsolvedParents, solvedParents,
                              identifier) {
  changeFlag = F
  m = mixedGraph$numNodes()
  for (i in 1:m) {
    unsolved = unsolvedParents[[i]]
    if (length(unsolved) != 0) {
      allowedNodesTrueFalse = logical(m)
      for (j in 1:m) {
        if (i != j &&
            !mixedGraph$isSibling(i,j) &&
            length(intersect(mixedGraph$htrFrom(j), unsolved)) != 0 &&
            length(intersect(mixedGraph$htrFrom(i), unsolvedParents[[j]])) == 0) {
          allowedNodesTrueFalse[j] = TRUE
        }
      }
      if (all(!allowedNodesTrueFalse)) {
        next
      }
      allowedNodes = which(allowedNodesTrueFalse)

      subsetFound = F
      for (k in length(unsolved):1) {
        if (length(unsolved) == 1) {
          # Silly conditional required because the combn function
          # interprets inputs of size 1 as a length rather than an as
          # a vector from which to get subsets
          subsets = list(unsolved)
        } else {
          subsets = combn(unsolved, k, simplify = F)
        }
        for (subset in subsets) {
          allowedForSubsetTrueFalse = logical(length(allowedNodes))
          for (l in 1:length(allowedNodes)) {
            a = allowedNodes[l]
            allowedForSubsetTrueFalse[l] = all(intersect(mixedGraph$htrFrom(a),
                                                         unsolved) %in% subset)
          }
          allowedForSubset = allowedNodes[allowedForSubsetTrueFalse]
          if (length(allowedForSubset) == 0) {
            next
          }
          halfTrekSystemResult = mixedGraph$getHalfTrekSystem(allowedForSubset,
                                                              subset)
          if (halfTrekSystemResult$systemExists) {
            changeFlag = T
            subsetFound = T
            activeFrom = halfTrekSystemResult$activeFrom
            identifier = createLinearIdentifierFunc(identifier,
                                                    activeFrom,
                                                    subset,
                                                    i,
                                                    lapply(activeFrom, function(x) {
                                                      intersect(mixedGraph$htrFrom(i),
                                                                mixedGraph$allParents(x))
                                                    }))
            solvedParents[[i]] = sort(c(subset, solvedParents[[i]]))
            unsolvedParents[[i]] = setdiff(unsolved, subset)
            break
          }
        }
        if (subsetFound) {
          break
        }
      }
    }
  }
  return(list(changeFlag = changeFlag, unsolvedParents = unsolvedParents,
              solvedParents = solvedParents, identifier = identifier))
}

# trekSeparationIdentifyStep = function(env, maxSubsetSize) {
#   changeFlag = F
#   m = nrow(O)
#   for (i in 1:m) {
#     unsolved = env$unsolvedParents[[i]]
#     if (length(unsolved) != 0) {
#       for (j in unsolved) {
#         edgeIdentified = F
#         for (k in 1:min(length(unsolved), m - 1)) {
#           if (m == 2) {
#             subsets = list(3 - i)
#           } else {
#             subsets = combn(setdiff(1:m, i), k, simplify = F)
#           }
#
#           for (subset in subsets) {
#             if (T) {# TODO
#               edgeIdentified = T
#               env$identifier = createTrekSeparationIdentifier(
#                 env$identifier,
#                 # TODO)
#                 env$solvedParents[[i]] = sort(c(subset, env$solvedParents[[i]]))
#                 env$unsolvedParents[[i]] = setdiff(unsolved, subset)
#                 break
#             }
#           }
#         }
#       }
#     }
#   }
#   return(changeFlag)
# }

#' Edge-wise identification
#'
#' @inheritParams graphID
#'
#' @return a list
#' @export
edgewiseID = function(L, O) {
  validateMatrices(L, O)
  O = 1 * ((O + t(O)) != 0)
  diag(O) = 0
  m = nrow(L)

  mixedGraph = MixedGraph(L, O)

  unsolvedParents = lapply(1:m, function(node) { mixedGraph$allParents(node) })

  identifier = createIdentifierBaseCase(L)
  changeFlag = T
  solvedParents = rep(list(numeric(0)), m)
  while (changeFlag) {
    idResult = linearIdentifyStep(mixedGraph, unsolvedParents,
                                  solvedParents, identifier)
    changeFlag = idResult$changeFlag
    unsolvedParents = idResult$unsolvedParents
    solvedParents = idResult$solvedParents
    identifier = idResult$identifier
  }

  return(list(solvedParents = solvedParents,
              unsolvedParents = unsolvedParents,
              identifier = identifier))
}



