#' Create an edgewise identification function
#'
#' A helper function for \code{\link{edgewiseIdentifyStep}}, creates an
#' identifier function based on its given parameters. This created identifier
#' function will identify the directed edges from 'targets' to 'node.'
#'
#' @inheritParams createHtcIdentifier
#' @param solvedNodeParents the parents of node that have been solved
#' @param sourceParentsToRemove a list of the parents of the sources that should
#'        have their edge to their respect source removed.
#'
#' @return an identification function
createEdgewiseIdentifier <- function(idFunc, sources, targets, node, solvedNodeParents, 
    sourceParentsToRemove) {
    # These assignments may seem redundent but they are necessary as they assign
    # these variables to the local environment of the function call.  This allows
    # them to persist and still be usable by the returned function.
    idFunc <- idFunc
    sources <- sources
    targets <- targets
    node <- node
    solvedNodeParents <- solvedNodeParents
    sourceParentsToRemove <- sourceParentsToRemove
    return(function(Sigma) {
        m <- nrow(Sigma)
        identifiedParams <- idFunc(Sigma)
        Lambda <- identifiedParams$Lambda
        
        SigmaMinus <- Sigma
        for (sourceInd in 1:length(sources)) {
            source <- sources[sourceInd]
            parentsToRemove <- sourceParentsToRemove[[sourceInd]]
            if (length(parentsToRemove) != 0) {
                SigmaMinus[source, ] <- Sigma[source, , drop = F] - t(Lambda[parentsToRemove, 
                  source, drop = F]) %*% Sigma[parentsToRemove, , drop = F]
            }
        }
        
        if (length(solvedNodeParents) != 0) {
            SigmaMinus[sources, node] <- SigmaMinus[sources, node, drop = F] - SigmaMinus[sources, 
                solvedNodeParents, drop = F] %*% Lambda[solvedNodeParents, node, 
                drop = F]
        }
        
        if (abs(det(SigmaMinus[sources, targets, drop = F])) < 10^-10) {
            stop("In identification, found near-singular system. Is the input matrix generic?")
        }
        
        Lambda[targets, node] <- solve(SigmaMinus[sources, targets, drop = F], SigmaMinus[sources, 
            node, drop = F])
        
        return(list(Lambda = Lambda, Omega = identifiedParams$Omega))
    })
}

#' Perform one iteration of edgewise identification.
#'
#' A function that does one step through all the nodes in a mixed graph
#' and tries to identify new edge coefficients using the existence of
#' half-trek systems as described in Weihs, Robeva, Robinson, et al. (2017).
#'
#' @inheritParams htcIdentifyStep
#' @param subsetSizeControl a positive integer (Inf allowed) which controls
#'                          the size of edgesets searched in the edgewiseID
#'                          algorithm. Suppose, for example, this has value 3.
#'                          Then if a node i has n parents, this will restrict
#'                          the algorithm to only look at subsets of the parents
#'                          of size 1,2,3 and n-2, n-1, n. Making this
#'                          parameter smaller means the algorithm will be faster
#'                          but less exhaustive (and hence less powerful).
#'
#' @return see the return of \code{\link{htcIdentifyStep}}.
#' @export
edgewiseIdentifyStep <- function(mixedGraph, unsolvedParents, solvedParents, identifier, 
    subsetSizeControl = Inf) {
    if (subsetSizeControl <= 0) {
        stop("Invalid subset size control parameter.")
    }
    identifiedEdges <- c()
    m <- mixedGraph$numNodes()
    for (i in 1:m) {
        unsolved <- unsolvedParents[[i]]
        htrFromNode <- mixedGraph$htrFrom(i)
        if (length(unsolved) != 0) {
            allowedNodesTrueFalse <- logical(m)
            for (j in 1:m) {
                if (i != j && !mixedGraph$isSibling(i, j) && length(intersect(mixedGraph$htrFrom(j), 
                  unsolved)) != 0 && length(intersect(htrFromNode, unsolvedParents[[j]])) == 
                  0) {
                  allowedNodesTrueFalse[j] <- TRUE
                }
            }
            if (all(!allowedNodesTrueFalse)) {
                next
            }
            allowedNodes <- which(allowedNodesTrueFalse)
            
            htrFromAllowedOrTrFromUnsolved <- rep(list(c()), length(allowedNodes))
            for (j in 1:length(allowedNodes)) {
                a <- allowedNodes[j]
                htrFromAllowedOrTrFromUnsolved[[j]] <- mixedGraph$htrFrom(a)
                for (unsolvedForA in unsolvedParents[[a]]) {
                  htrFromAllowedOrTrFromUnsolved[[j]] <- c(htrFromAllowedOrTrFromUnsolved[[j]], 
                    mixedGraph$trFrom(unsolvedForA))
                }
                htrFromAllowedOrTrFromUnsolved[[j]] <- unique(htrFromAllowedOrTrFromUnsolved[[j]])
            }
            
            subsetFound <- F
            subsetSizes <- union(length(unsolved):min(max((length(unsolved) - subsetSizeControl + 
                1), 1), length(unsolved)), 1:min(length(unsolved), subsetSizeControl))
            
            for (k in subsetSizes) {
                subsets <- subsetsOfSize(unsolved, k)
                for (subset in subsets) {
                  allowedForSubsetTrueFalse <- logical(length(allowedNodes))
                  for (l in 1:length(allowedNodes)) {
                    a <- allowedNodes[l]
                    allowedForSubsetTrueFalse[l] <- all(intersect(htrFromAllowedOrTrFromUnsolved[[l]], 
                      unsolved) %in% subset)
                  }
                  allowedForSubset <- allowedNodes[allowedForSubsetTrueFalse]
                  if (length(allowedForSubset) == 0) {
                    next
                  }
                  halfTrekSystemResult <- mixedGraph$getHalfTrekSystem(allowedForSubset, 
                    subset)
                  
                  if (halfTrekSystemResult$systemExists) {
                    identifiedEdges <- c(identifiedEdges, as.integer(rbind(subset, 
                      i)))
                    subsetFound <- T
                    activeFrom <- halfTrekSystemResult$activeFrom
                    identifier <- createEdgewiseIdentifier(identifier, activeFrom, 
                      subset, i, solvedParents[[i]], lapply(activeFrom, function(x) {
                        solvedParents[[x]]
                      }))
                    solvedParents[[i]] <- sort(c(subset, solvedParents[[i]]))
                    unsolvedParents[[i]] <- setdiff(unsolved, subset)
                    break
                  }
                }
                if (subsetFound) {
                  break
                }
            }
        }
    }
    return(list(identifiedEdges = identifiedEdges, unsolvedParents = unsolvedParents, 
        solvedParents = solvedParents, identifier = identifier))
}

#' Determines which edges in a mixed graph are edgewiseID-identifiable
#'
#' Uses the edgewise identification criterion of Weihs, Robeva, Robinson, et al.
#' (2017) to determine which edges in a mixed graph are generically
#' identifiable.
#'
#' @export
#'
#' @inheritParams generalGenericID
#' @inheritParams semID
#' @inheritParams edgewiseIdentifyStep
#'
#' @return see the return of \code{\link{generalGenericID}}.
edgewiseID <- function(mixedGraph, tianDecompose = T, subsetSizeControl = 3) {
    eid <- function(mixedGraph, unsolvedParents, solvedParents, identifier) {
        return(edgewiseIdentifyStep(mixedGraph, unsolvedParents, solvedParents, identifier, 
            subsetSizeControl = subsetSizeControl))
    }
    result <- generalGenericID(mixedGraph, list(eid), tianDecompose = tianDecompose)
    result$call <- match.call()
    return(result)
}
