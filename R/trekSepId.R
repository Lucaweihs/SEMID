#' Create an trek separation identification function
#'
#' A helper function for \code{\link{trekSeparationIdentifyStep}}, creates an
#' identifier function based on its given parameters. This created identifier
#' function will identify the directed edge from 'parent' to 'node.'
#'
#' @inheritParams createHtcIdentifier
#' @param parent the parent of node for which the edge node -> parent should
#'               be generically identified.
#' @param solvedParents the parents of node that have been solved
#'
#' @return an identification function
createTrekSeparationIdentifier <- function(idFunc, sources, targets, node, parent,
    solvedParents) {
    # These assignments may seem redundent but they are necessary as they assign
    # these variables to the local environment of the function call.  This allows
    # them to persist and still be usable by the returned function.
    idFunc <- idFunc
    sources <- sources
    targets <- targets
    node <- node
    parent <- parent
    solvedParents <- solvedParents
    return(function(Sigma) {
        m <- nrow(Sigma)
        identifiedParams <- idFunc(Sigma)
        Lambda <- identifiedParams$Lambda

        SigmaMinus <- Sigma
        for (solvedParent in solvedParents) {
            SigmaMinus[sources, node] <- SigmaMinus[sources, node] - SigmaMinus[sources,
                solvedParent] * Lambda[solvedParent, node]
        }
        subSigmaNode <- SigmaMinus[sources, c(targets, node), drop = F]
        subSigmaParent <- SigmaMinus[sources, c(targets, parent), drop = F]

        Lambda[parent, node] <- det(subSigmaNode)/det(subSigmaParent)

        return(list(Lambda = Lambda, Omega = identifiedParams$Omega))
    })
}

#' Perform one iteration of trek separation identification.
#'
#' A function that does one step through all the nodes in a mixed graph
#' and tries to identify new edge coefficients using trek-separation as
#' described in Weihs, Robeva, Robinson, et al. (2017).
#'
#' @inheritParams htcIdentifyStep
#' @param maxSubsetSize a positive integer which controls the maximum subset
#'                      size considered in the trek-separation identification
#'                      algorithm. Making this parameter smaller means the
#'                      algorithm will be faster but less exhaustive (and hence
#'                      less powerful).
#'
#' @return see the return of \code{\link{htcIdentifyStep}}.
#' @export
trekSeparationIdentifyStep <- function(mixedGraph, unsolvedParents, solvedParents,
    identifier, maxSubsetSize = 3) {
    if (maxSubsetSize <= 0) {
        stop("Max subset size must be >= 1")
    }
    m <- mixedGraph$numNodes()
    identifiedEdges <- c()
    for (i in 1:m) {
        unsolvedBefore <- unsolvedParents[[i]]
        component <- mixedGraph$stronglyConnectedComponent(i)
        if (length(unsolvedBefore) != 0 && length(component) == 1) {
            allButI <- setdiff(1:m, i)
            iDescendants <- mixedGraph$descendants(i)
            nonIDescendants <- setdiff(allButI, iDescendants)
            for (j in unsolvedBefore) {
                edgeIdentified <- F
                for (k in 1:min(maxSubsetSize, m - 1)) {
                  sourceSets <- subsetsOfSize(1:m, k)
                  targetSets <- subsetsOfSize(setdiff(nonIDescendants, j), k - 1)

                  for (sources in sourceSets) {
                    for (targets in targetSets) {
                      systemWithJ <- mixedGraph$getTrekSystem(sources, c(targets,
                        j))
                      if (systemWithJ$systemExists) {
                        toRemoveOnRight <- as.integer(rbind(c(j, solvedParents[[i]]),
                          i))
                        systemWithoutEdges <- mixedGraph$getTrekSystem(sources, c(targets,
                          i), avoidRightEdges = toRemoveOnRight)
                        if (!systemWithoutEdges$systemExists) {
                          identifiedEdges <- c(identifiedEdges, i, j)
                          edgeIdentified <- T
                          identifier <- createTrekSeparationIdentifier(identifier,
                            sources, targets, i, j, solvedParents[[i]])
                          solvedParents[[i]] <- sort(c(j, solvedParents[[i]]))
                          unsolvedParents[[i]] <- setdiff(unsolvedParents[[i]], j)
                          break
                        }
                      }
                    }
                    if (edgeIdentified) {
                      break
                    }
                  }
                  if (edgeIdentified) {
                    break
                  }
                }
            }
        }
    }
    return(list(identifiedEdges = identifiedEdges, unsolvedParents = unsolvedParents,
        solvedParents = solvedParents, identifier = identifier))
}

#' Determines which edges in a mixed graph are edgewiseID+TS identifiable
#'
#' Uses the edgewise+TS identification criterion of Weihs, Robeva, Robinson, et
#' al. (2017) to determine which edges in a mixed graph are generically
#' identifiable. In particular this algorithm iterates between the half-trek,
#' edgewise, and trek-separation identification algorithms in an attempt to
#' identify as many edges as possible, this may be very slow.
#'
#' @export
#'
#' @inheritParams generalGenericID
#' @inheritParams semID
#' @inheritParams edgewiseIdentifyStep
#' @inheritParams trekSeparationIdentifyStep
#'
#' @return see the return of \code{\link{generalGenericID}}.
edgewiseTSID <- function(mixedGraph, tianDecompose = T, subsetSizeControl = 3, maxSubsetSize = 3) {
    eid <- function(mixedGraph, unsolvedParents, solvedParents, identifier) {
        return(edgewiseIdentifyStep(mixedGraph, unsolvedParents, solvedParents, identifier,
            subsetSizeControl = subsetSizeControl))
    }
    tsid <- function(mixedGraph, unsolvedParents, solvedParents, identifier) {
        return(trekSeparationIdentifyStep(mixedGraph, unsolvedParents, solvedParents,
            identifier, maxSubsetSize = maxSubsetSize))
    }
    result <- generalGenericID(mixedGraph, list(htcIdentifyStep, eid, tsid), tianDecompose = tianDecompose)
    result$call <- match.call()
    return(result)
}
