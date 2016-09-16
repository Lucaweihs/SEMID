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

createIdentifierBaseCase <- function(L) {
  return(function(Sigma) {
    Lambda <- matrix(NA, nrow(L), nrow(L))
    Lambda[L == 0] = 0
    return(Lambda)
  })
}

createLinearIdentifier <- function(idFunc, sources, targets, node, solvedNodeParents,
                                   sourceParentsToRemove) {
  # These assignments may seem redundent but they are necessary as they
  # assign these variables to the local environment of the function call.
  # This allows them to persist and still be usable by the returned function.
  idFunc <- idFunc
  sources <- sources
  targets <- targets
  node <- node
  solvedNodeParents <- solvedNodeParents
  sourceParentsToRemove <- sourceParentsToRemove
  return(
    function(Sigma) {
      m <- nrow(Sigma)
      Lambda <- idFunc(Sigma)

      SigmaMinus = Sigma
      for (sourceInd in 1:length(sources)) {
        source = sources[sourceInd]
        parentsToRemove = sourceParentsToRemove[[sourceInd]]
        if (length(parentsToRemove) != 0) {
          SigmaMinus[source,] <- Sigma[source, , drop = F] -
            t(Lambda[parentsToRemove, source, drop = F]) %*% Sigma[parentsToRemove,, drop = F]
        }
      }

      if (length(solvedNodeParents) != 0) {
        SigmaMinus[sources, node] = SigmaMinus[sources, node, drop = F] -
          SigmaMinus[sources, solvedNodeParents, drop = F] %*% Lambda[solvedNodeParents, node, drop = F]
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

createTrekSeparationIdentifier <- function(idFunc, sources, sinks, node, parent) {
  # These assignments may seem redundent but they are necessary as they
  # assign these variables to the local environment of the function call.
  # This allows them to persist and still be usable by the returned function.
  idFunc <- idFunc
  sources <- sources
  sinks <- sinks
  node <- node
  parent <- parent
  return(
    function(Sigma) {
      m <- nrow(Sigma)
      Lambda <- idFunc(Sigma)

      subSigmaNode = Sigma[sources, c(sinks, node), drop = F]
      subSigmaParent = Sigma[sources, c(sinks, parent), drop = F]

      Lambda[parent, node] = det(subSigmaNode) / det(subSigmaParent)
      return(Lambda)
    }
  )
}

#' A helper function for edgewiseID that does one step through all the nodes
#' and tries to identify new edge coefficients using half-treks.
#'
#' @return a list
#' @export
linearIdentifyStep = function(mixedGraph, unsolvedParents, solvedParents,
                              identifier) {
  changeFlag = F
  m = mixedGraph$numNodes()
  for (i in 1:m) {
    unsolved = unsolvedParents[[i]]
    solved = solvedParents[[i]]
    htrFromNode = mixedGraph$htrFrom(i)
    if (length(unsolved) != 0) {
      allowedNodesTrueFalse = logical(m)
      for (j in 1:m) {
        if (i != j &&
            !mixedGraph$isSibling(i,j) &&
            length(intersect(mixedGraph$htrFrom(j), unsolved)) != 0 &&
            length(intersect(htrFromNode, unsolvedParents[[j]])) == 0) {
          allowedNodesTrueFalse[j] = TRUE
        }
      }
      if (all(!allowedNodesTrueFalse)) {
        next
      }
      allowedNodes = which(allowedNodesTrueFalse)

      htrFromAllowedOrTrFromUnsolved = rep(list(c()), length(allowedNodes))
      for (j in 1:length(allowedNodes)) {
        a = allowedNodes[j]
        htrFromAllowedOrTrFromUnsolved[[j]] = mixedGraph$htrFrom(a)
        for (unsolvedForA in unsolvedParents[[a]]) {
          htrFromAllowedOrTrFromUnsolved[[j]] = c(htrFromAllowedOrTrFromUnsolved[[j]],
                                                  mixedGraph$trFrom(unsolvedForA))
        }
        htrFromAllowedOrTrFromUnsolved[[j]] = unique(htrFromAllowedOrTrFromUnsolved[[j]])
      }

      subsetFound = F
      for (k in length(unsolved):1) {
        subsets = subsetsOfSize(unsolved, k)
        for (subset in subsets) {
          allowedForSubsetTrueFalse = logical(length(allowedNodes))
          for (l in 1:length(allowedNodes)) {
            a = allowedNodes[l]
            allowedForSubsetTrueFalse[l] = all(intersect(htrFromAllowedOrTrFromUnsolved[[l]],
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
            identifier = createLinearIdentifier(identifier,
                                                activeFrom,
                                                subset,
                                                i,
                                                solvedParents[[i]],
                                                lapply(activeFrom, function(x) {
                                                  solvedParents[[x]]
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


#' A helper function for edgewiseID that does one step through all the nodes
#' and tries to identify new edge coefficients using trek separation.
#'
#' @return a list
#' @export
trekSeparationIdentifyStep = function(mixedGraph, unsolvedParents,
                                      solvedParents, identifier,
                                      maxSubsetSize = 3) {
  if (maxSubsetSize <= 0) {
    stop("Max subset size must be >= 1")
  }
  changeFlag = F
  m = mixedGraph$numNodes()
  for (i in 1:m) {
    unsolvedBefore = unsolvedParents[[i]]
    component = mixedGraph$stronglyConnectedComponent(i)
    if (length(unsolvedBefore) != 0 && length(component) == 1) {
      nonIDescendants = setdiff(1:m, c(mixedGraph$allDescendants(i), i))
      for (j in unsolvedBefore) {
        edgeIdentified = F
        for (k in 1:min(maxSubsetSize, m - 1)) {
          subsetsBig = subsetsOfSize(nonIDescendants, k)
          subsetsSmall = subsetsOfSize(setdiff(nonIDescendants, j), k - 1)

          for (subsetBig in subsetsBig) {
            for (subsetSmall in subsetsSmall) {
              systemWithJ = mixedGraph$getTrekSystem(subsetBig, c(subsetSmall, j))
              if (systemWithJ$systemExists) {
                systemWithoutJIEdge = mixedGraph$getTrekSystem(subsetBig,
                                                               c(subsetSmall, i),
                                                               c(j,i))
                if (!systemWithoutJIEdge$systemExists) {
                  print("HAPPENS")
                  changeFlag = T
                  edgeIdentified = T
                  identifier =
                    createTrekSeparationIdentifier(identifier, subsetBig,
                                                   subsetSmall, i, j)
                  solvedParents[[i]] = sort(c(j, solvedParents[[i]]))
                  unsolvedParents[[i]] = setdiff(unsolvedParents[[i]], j)
                  break
                }
              }
            }
            if (edgeIdentified) { break }
          }
          if (edgeIdentified) { break }
        }
      }
    }
  }
  return(list(changeFlag = changeFlag, unsolvedParents = unsolvedParents,
              solvedParents = solvedParents, identifier = identifier))
}

#' Edge-wise identification
#'
#' @inheritParams graphID
#'
#' @return a list
#' @export
edgewiseID = function(L, O, maxTrekSepSubsetSize = 3) {
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

    if (!changeFlag && maxTrekSepSubsetSize != 0) {
      idResult = trekSeparationIdentifyStep(mixedGraph, unsolvedParents,
                                            solvedParents, identifier,
                                            maxTrekSepSubsetSize)
      changeFlag = idResult$changeFlag
      unsolvedParents = idResult$unsolvedParents
      solvedParents = idResult$solvedParents
      identifier = idResult$identifier
    }
  }

  return(list(solvedParents = solvedParents,
              unsolvedParents = unsolvedParents,
              identifier = identifier))
}



