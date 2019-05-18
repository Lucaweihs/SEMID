#' Create an latent factor half-trek critierion identification function.
#'
#' A helper function for \code{\link{lfhtcIdentifyStep}}, creates an identifier
#' function based on its given parameters. This created identifier function will
#' identify the directed edges from 'targets' to 'node.'
#'
#' @param idFunc identification of edge coefficients often requires that other
#'        edge coefficients already be identified. This argument should be a
#'        function that produces all such identifications. The newly created
#'        identifier function will return these identifications along with its
#'        own.
#' @param sources the sources of the half-trek system.
#' @param targets the targets of the half-trek system (these should be the
#'        parents of node).
#' @param node the node for which all incoming edges are to be identified
#'        (the tails of which are targets).
#' @param htrSources the nodes in sources which are half-trek reachable from
#'        node. All incoming edges to these sources should be identified by
#'        idFunc for the newly created identification function to work.
#'
#' @return an identification function
#'
#' @references
#' TO BE WRITTEN
createLFHtcIdentifier <- function(idFunc, v, Y, Z, parents, reachableY) {
  # Necessary redundant assignments
  idFunc <- idFunc
  v <- v
  Y <- Y
  Z <- Z
  parents <- sort(parents)
  reachableY <- reachableY
  return(function(Sigma) {
    m <- nrow(Sigma)
    identifiedParams <- idFunc(Sigma)
    Lambda <- identifiedParams$Lambda

    targets <- c(parents, Z)
    SigmaMinus <- Sigma
    for (y in Y) {
      if (y %in% reachableY) {
        SigmaMinus[y, ] <- Sigma[y,] -
          (t(Lambda[, y, drop = F]) %*% Sigma)
      }
      SigmaMinus[y, Z] <- SigmaMinus[y,Z] -
        SigmaMinus[y, , drop = F] %*% Lambda[,Z]
    }

    if (abs(det(SigmaMinus[Y, targets, drop = F])) < 10^-10) {
      stop("In identification, found near-singular system. Is the input matrix generic?")
    }

    Lambda[parents, v] <-
      solve(SigmaMinus[Y, targets, drop = F],
            SigmaMinus[Y, v, drop = F])[seq(1, length = length(parents)), ]

    if (!any(is.na(Lambda))) {
      Omega <- t(diag(m) - Lambda) %*% Sigma %*% (diag(m) - Lambda)
      return(list(Lambda = Lambda, Omega = Omega))
    }
    return(list(Lambda = Lambda, Omega = identifiedParams$Omega))
  })
}

#' A helper function to validate that latent nodes in a LatentDigraph are sources.
#'
#' Produces an error not all latent nodes are sources.
#'
#' @param graph the LatentDigraph
validateLatentNodesAreSources <- function(graph) {
  for (node in graph$latentNodes()) {
    if (length(graph$parents(node)) != 0) {
      stop("Latent node in input graph is not a source node, i.e. it has a parent.")
    }
  }
}

createYZSets <- function(allowedYZPerLatent) {
  if (length(allowedYZPerLatent) == 0) {
    return(list(list(Y = integer(0), Z = integer(0))))
  }
  Y <- allowedYZPerLatent[[1]]$Y
  Z <- allowedYZPerLatent[[1]]$Z
  YZPairsRecurse <- createYZSets(allowedYZPerLatent[-1])
  YZPairsNew = vector("list", length(Y) * length(Z) * length(YZPairsRecurse))
  k = 1
  for (y in Y) {
    for (z in Z) {
      if (y != z) {
        for (YZPair in YZPairsRecurse) {
          YZPairsNew[[k]] = list(Y = c(y, YZPair$Y),
                              Z = c(z, YZPair$Z))

          k = k + 1
        }
      }
    }
  }
  return(YZPairsNew[seq(1, length = k - 1)])
}

getPossibleYZ <- function(graph, L, allowedForY, allowedForZ) {
  if (length(L) == 0) {
    return(list(list(Y = integer(0), Z = integer(0))))
  }
  allowedYZPerLatent = vector("list", length = length(L))
  for (i in 1:length(L)) {
    children = graph$children(L[i])
    allowedYZPerLatent[[i]] = list(Y = intersect(allowedForY, children),
                                   Z = intersect(allowedForZ, children))
  }
  return(createYZSets(allowedYZPerLatent))
}

#' Perform one iteration of latent factor HTC identification.
#'
#' A function that does one step through all the nodes in a latent factor graph
#' and tries to identify new edge coefficients using the existence of
#' half-trek systems as described in TO BE WRITTEN.
#'
#' @param graph a \code{\link{LatentDigraph}} object representing
#'         the latent factor graph. All latent nodes in this graph should be
#'         source nodes (i.e. have no parents).
#' @param unsolvedParents a list whose ith index is a vector of all the parents
#'        j of i in the graph which for which the edge j->i is not yet known to be
#'        generically identifiable.
#' @param solvedParents the complement of \code{unsolvedParents}, a list whose
#'        ith index is a vector of all parents j of i for which the edge i->j
#'        is known to be generically identifiable (perhaps by other algorithms).
#' @param identifier an identification function that must produce the
#'        identifications corresponding to those in solved parents. That is
#'        \code{identifier} should be a function taking a single argument Sigma
#'        (any generically generated covariance matrix corresponding
#'        to the latent factor graph) and returns a list with two named arguments
#' \describe{
#'   \item{Lambda}{denote the number of nodes in \code{graph} as n. Then
#'                 Lambda is an nxn matrix whose i,jth entry
#' \enumerate{
#'   \item equals 0 if i is not a parent of j,
#'   \item equals NA if i is a parent of j but \code{identifier} cannot
#'         identify it generically,
#'   \item equals the (generically) unique value corresponding to the weight
#'         along the edge i->j that was used to produce Sigma.
#' }}
#'   \item{Omega}{just as Lambda but for the error covariance matrix for the
#'                latent factor graph.}
#' }
#'        such that if j is in \code{solvedParents[[i]]} we must have that
#'        Lambda[j,i] is not NA.
#' @param subsetSizeControl the largest subset of latent nodes to consider.
#'
#' @return a list with four components:
#' \describe{
#'   \item{\code{identifiedEdges}}{a matrix rx2 matrix where r is the number
#'   of edges that where identified by this function call and
#'   \code{identifiedEdges[i,1] -> identifiedEdges[i,2]} was the ith edge
#'   identified}
#'   \item{\code{unsolvedParents}}{as the input argument but updated with
#'   any newly identified edges}
#'   \item{\code{solvedParents}}{as the input argument but updated with
#'   any newly identified edges}
#'   \item{\code{identifier}}{as the input argument but updated with
#'   any newly identified edges}
#' }
#'
#' @export
#'
#' @references
#' TO BE WRITTEN
lfhtcIdentifyStep <- function(graph, unsolvedParents, solvedParents, identifier,
                              subsetSizeControl = Inf) {
  validateLatentNodesAreSources(graph)
  identifiedEdges <- numeric(0)

  observedNodes <- graph$observedNodes()
  latentNodes <- graph$latentNodes()
  numObserved <- length(observedNodes)
  numLatents <- length(latentNodes)

  edgeMat <- graph$L()
  edgesBetweenObserved <- which(edgeMat[seq(1, length = numObserved),
                                        seq(1, length = numObserved), drop = F] == 1,
                                arr.ind = T)
  edgesBetweenObserved <- matrix(observedNodes[edgesBetweenObserved], ncol = 2)

  solvedNodes <- which(sapply(unsolvedParents, FUN = function(x) {
    length(x) == 0
  }))

  for (i in setdiff(observedNodes, solvedNodes)) {
    allParents <- graph$parents(i)
    latentParents <- intersect(latentNodes, allParents)
    observedParents <- intersect(observedNodes, allParents)

    childrenOfLatentParents <-
      lapply(latentParents, FUN = function(x) { graph$children(x) })
    latentParentHasGeq4Children <- sapply(childrenOfLatentParents,
                                     FUN = function(x) { length(x) >= 4 })
    latentsToControl <- latentParents[latentParentHasGeq4Children]

    for (k in seq(0, length = 1 + min(subsetSizeControl, length(latentsToControl)))) {
      for (L in subsetsOfSize(latentsToControl, k)) {
        Lcomp <- setdiff(latentParents, L)

        htrFromIAvoidingL <- intersect(observedNodes,
                                  union(graph$descendants(i),
                                        graph$trFrom(Lcomp)))
        siblingsOfIAvoidingL <- graph$children(Lcomp)
        maybeAllowedForY <- setdiff(observedNodes,
                               c(i, siblingsOfIAvoidingL,
                                 setdiff(htrFromIAvoidingL, solvedNodes)))
        allowedForZ <- setdiff(solvedNodes, allParents)
        possibleYZs <- getPossibleYZ(graph, L,
                                    allowedForY = maybeAllowedForY,
                                    allowedForZ = allowedForZ)
        for (YZ in possibleYZs) {
          Y = YZ$Y
          Z = YZ$Z
          if (length(unique(Y)) == length(Y) &&
              length(unique(Z)) == length(Z) &&
              length(intersect(Y, Z)) == 0) {
            latentParentsOfZandI <- intersect(latentNodes, graph$parents(c(i,Z)))
            latentParentsOfZandINotInL <- setdiff(latentParentsOfZandI, L)
            htrFromZandIAvoidingL <- union(graph$descendants(c(i,Z)),
                                           graph$trFrom(latentParentsOfZandINotInL))
            allowed <- setdiff(maybeAllowedForY,
                               setdiff(htrFromZandIAvoidingL, solvedNodes))
            allowed <- setdiff(allowed,
                               c(graph$children(latentParentsOfZandINotInL),
                                 Z, i))
            if (any(!(Y %in% allowed))) {
              next
            }

            trekSystemResults <- graph$getTrekSystem(
              allowed, observedParents,
              avoidLeftNodes = c(Y, L), avoidRightNodes = c(Z, L),
              avoidLeftEdges = edgesBetweenObserved)

            if (trekSystemResults$systemExists) {
              identifiedEdges <- c(identifiedEdges, as.integer(rbind(observedParents, i)))
              activeFrom <- trekSystemResults$activeFrom
              identifier <- createLFHtcIdentifier(identifier, v = i,
                                                  Y = c(activeFrom, Y),
                                                  Z = Z,
                                                  parents = observedParents,
                                                  reachableY = htrFromZandIAvoidingL)
              solvedParents[[i]] <- observedParents
              unsolvedParents[[i]] <- integer(0)
              solvedNodes <- c(i, solvedNodes)
              break
            }
          }
        }
        if (i %in% solvedNodes) {
          break
        }
      }
      if (i %in% solvedNodes) {
        break
      }
    }
  }
  return(list(identifiedEdges = matrix(identifiedEdges, byrow = T, ncol = 2),
              unsolvedParents = unsolvedParents,
              solvedParents = solvedParents,
              identifier = identifier))
}