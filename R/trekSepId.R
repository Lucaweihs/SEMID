#' Create an trek seperation identification function
#'
#' TODO: Add details
#'
#' @return an identification function
createTrekSeparationIdentifier <- function(idFunc, sources, targets, node, parent,
                                           solvedParents) {
  # These assignments may seem redundent but they are necessary as they
  # assign these variables to the local environment of the function call.
  # This allows them to persist and still be usable by the returned function.
  idFunc <- idFunc
  sources <- sources
  targets <- targets
  node <- node
  parent <- parent
  solvedParents <- solvedParents
  return(
    function(Sigma) {
      m <- nrow(Sigma)
      identifiedParams <- idFunc(Sigma)
      Lambda <- identifiedParams$Lambda

      SigmaMinus = Sigma
      for (solvedParent in solvedParents) {
        SigmaMinus[sources, node] =
          SigmaMinus[sources, node] - SigmaMinus[sources, solvedParent] * Lambda[solvedParent, node]
      }
      subSigmaNode = SigmaMinus[sources, c(targets, node), drop = F]
      subSigmaParent = SigmaMinus[sources, c(targets, parent), drop = F]

      Lambda[parent, node] = det(subSigmaNode) / det(subSigmaParent)

      return(list(Lambda = Lambda, Omega = identifiedParams$Omega))
    }
  )
}

#' Perform one iteration of trek seperation identification.
#'
#' A function that does one step through all the nodes in a mixed graph
#' and tries to identify new edge coefficients using trek-separation.
#'
#' @return a list
#' @export
trekSeparationIdentifyStep = function(mixedGraph, unsolvedParents,
                                      solvedParents, identifier,
                                      maxSubsetSize = 3) {
  if (maxSubsetSize <= 0) {
    stop("Max subset size must be >= 1")
  }
  m = mixedGraph$numNodes()
  identifiedEdges = c()
 for (i in 1:m) {
    unsolvedBefore = unsolvedParents[[i]]
    component = mixedGraph$stronglyConnectedComponent(i)
    if (length(unsolvedBefore) != 0 && length(component) == 1) {
      allButI = setdiff(1:m, i)
      iDescendants = mixedGraph$allDescendants(i)
      nonIDescendants = setdiff(allButI, iDescendants)
      for (j in unsolvedBefore) {
        edgeIdentified = F
        for (k in 1:min(maxSubsetSize, m - 1)) {
          sourceSets = subsetsOfSize(1:m, k)
          targetSets = subsetsOfSize(setdiff(nonIDescendants, j), k - 1)

          for (sources in sourceSets) {
            for (targets in targetSets) {
              systemWithJ = mixedGraph$getTrekSystem(sources, c(targets, j))
              if (systemWithJ$systemExists) {
                toRemoveOnRight = as.numeric(rbind(c(j, solvedParents[[i]]), i))
                systemWithoutEdges = mixedGraph$getTrekSystem(sources,
                                                              c(targets, i),
                                                              toRemoveOnRight)
                if (!systemWithoutEdges$systemExists) {
                  identifiedEdges = c(identifiedEdges, i, j)
                  edgeIdentified = T
                  identifier =
                    createTrekSeparationIdentifier(identifier, sources,
                                                   targets, i, j, solvedParents[[i]])
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
  return(list(identifiedEdges = identifiedEdges, unsolvedParents = unsolvedParents,
              solvedParents = solvedParents, identifier = identifier))
}
