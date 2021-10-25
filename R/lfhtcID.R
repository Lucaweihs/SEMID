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

# Helper function for getPossibleYZ()
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

# Returns a list of all possible YZ-pairs with |Y| = |Z| <= |L|
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
lfhtcIdentifyStep <- function(graph, unsolvedParents, solvedParents, activeFroms, Zs, Ls, identifier,
                              subsetSizeControl = Inf) {
  # Sanity check
  validateLatentNodesAreSources(graph)

  # Variable to store all newly identified edges
  identifiedEdges <- numeric(0)

  # Collect basic infos from graph
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

  # Only latent nodes with >=  children may be possibly in L
  childrenOfLatentNodes <- lapply(latentNodes, FUN = function(x) { graph$children(x) })
  latentNodeHasGeq4Children <- vapply(childrenOfLatentNodes,
                                        FUN = function(x) { length(x) >= 4 },
                                        logical(1))
  latentsToControl <- latentNodes[latentNodeHasGeq4Children]

  # Loop over all unsolved nodes
  for (i in setdiff(observedNodes, solvedNodes)) {

    # Collect basic info of unsolved node i
    allParents <- graph$parents(i)
    latentParents <- intersect(latentNodes, allParents)
    observedParents <- intersect(observedNodes, allParents)

    # # Only latent parents of i with >= 4 children may possibly be in L
    # childrenOfLatentParents <-
    #   lapply(latentParents, FUN = function(x) { graph$children(x) })
    # latentParentHasGeq4Children <- vapply(childrenOfLatentParents,
    #                                  FUN = function(x) { length(x) >= 4 },
    #                                  logical(1))
    # latentsToControl <- latentParents[latentParentHasGeq4Children]

    # Loop over possible cardinalities of the L
    for (k in seq(0, length = 1 + min(subsetSizeControl, length(latentsToControl)))) {
      # Loop over all subsets L in latentsToControl with cardinality k, i.e. |L| = k
      for (L in subsetsOfSize(latentsToControl, k)) {
        Lcomp <- setdiff(latentParents, L)  # All latent nodes of i not in L

        # Find set maybeAllowdForY. These are all nodes that are
        # - not i
        # - no siblings of i avoiding L (all latent parents of Y that are at the same
        #                                time latent parents of i have to be in L!)
        # - solved if half-trek reachable from i avoiding L
        htrFromIAvoidingL <- intersect(observedNodes,
                                  union(graph$descendants(i),
                                        graph$trFrom(Lcomp)))  # trFrom(l) finds all half-treks
                                                               # starting with i<-l-> ...
        siblingsOfIAvoidingL <- graph$children(Lcomp)
        maybeAllowedForY <- setdiff(observedNodes,
                               c(i, siblingsOfIAvoidingL,
                                 setdiff(htrFromIAvoidingL, solvedNodes)))

        # Each element in Z has to be solved and has to be no parent of i
        # (Since in this case we would have sided intersection in the system of half-treks)
        allowedForZ <- setdiff(solvedNodes, allParents)

        # All possible YZ-pairs with |Y| = |Z| <= |L|
        # NOTE: The definition of Y is DIFFERENT than the one in the paper!!!
        # If k = 0 (i.e. L is empty set), then we only have one YZ-pair: Y and Z both empty
        possibleYZs <- getPossibleYZ(graph, L,
                                    allowedForY = maybeAllowedForY,
                                    allowedForZ = allowedForZ)

        # Loop over all pairs (Y,Z) with |Y| = |Z| = |L| = k
        for (YZ in possibleYZs) {
          Y = YZ$Y
          Z = YZ$Z
          # Sanity check of Y and Z
          if (length(unique(Y)) == length(Y) &&
              length(unique(Z)) == length(Z) &&
              length(intersect(Y, Z)) == 0) {

            # Find allowed set for Y. This is now possibly since Z (and L) is known.
            # Consists of all nodes in maybeAllowedForY that are
            # - not in (Z or i)
            # - no siblings of (Z or i) avoiding L (all latent parents of Y that are at the same
            #                                time latent parents of (Z or i) have to be in L!)
            # - solved if half-trek reachable from (Z or i) avoiding L
            latentParentsOfZandI <- intersect(latentNodes, graph$parents(c(i,Z)))
            latentParentsOfZandINotInL <- setdiff(latentParentsOfZandI, L)
            htrFromZandIAvoidingL <- union(graph$descendants(c(i,Z)),
                                           graph$trFrom(latentParentsOfZandINotInL))
            allowed <- setdiff(maybeAllowedForY,
                               setdiff(htrFromZandIAvoidingL, solvedNodes))
            allowed <- setdiff(allowed,
                               c(graph$children(latentParentsOfZandINotInL),
                                 Z, i))

            # Sanity check: Maybe we have nodes in Y that are not allowed any more
            if (any(!(Y %in% allowed))) {
              next
            }

            # Check if there is a half-trek system from allowed nodes to pa(i)
            # - avoid (Y,L) on LHS
            # - avoid (Z,L) on RHS
            #   (These two conditions enforce no sided intersection of half-trek system!)
            # - Avoid starting with "<-" vetween obsrved nodes.
            #   This forces trek system to be half-trek system.
            trekSystemResults <- graph$getTrekSystem(
              allowed, observedParents,
              avoidLeftNodes = c(Y, L), avoidRightNodes = c(Z, L),
              avoidLeftEdges = edgesBetweenObserved)

            # If half-trek system exists we know that all edges between pa(i) and i are identified
            if (trekSystemResults$systemExists) {
              identifiedEdges <- c(identifiedEdges, as.integer(rbind(observedParents, i)))
              activeFrom <- trekSystemResults$activeFrom
              identifier <- createLFHtcIdentifier(identifier, v = i,
                                                  Y = c(activeFrom, Y),  # This recovers Y from paper
                                                  Z = Z,
                                                  parents = observedParents,
                                                  reachableY = htrFromZandIAvoidingL) # half-trek reachable elements of Y are known
              solvedParents[[i]] <- observedParents
              unsolvedParents[[i]] <- integer(0)
              activeFroms[[i]] <- c(activeFrom, Y)
              Zs[[i]] <- Z
              Ls[[i]] <- L
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
              activeFroms = activeFroms,
              Zs = Zs,
              Ls = Ls,
              identifier = identifier))
}




lfhtcIdentifyStep2 <- function(graph, unsolvedParents, solvedParents, activeFroms, Zs, Ls, identifier,
                              subsetSizeControl = Inf) {
  # Sanity check
  validateLatentNodesAreSources(graph)

  # Variable to store all newly identified edges
  identifiedEdges <- numeric(0)

  # Collect basic infos from graph
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

  # Only latent nodes with >=  children may be possibly in L
  childrenOfLatentNodes <- lapply(latentNodes, FUN = function(x) { graph$children(x) })
  latentNodeHasGeq4Children <- vapply(childrenOfLatentNodes,
                                      FUN = function(x) { length(x) >= 4 },
                                      logical(1))
  latentsToControl <- latentNodes[latentNodeHasGeq4Children]

  # Loop over all unsolved nodes
  for (i in setdiff(observedNodes, solvedNodes)) {

    # Collect basic info of unsolved node i
    allParents <- graph$parents(i)
    latentParents <- intersect(latentNodes, allParents)
    observedParents <- intersect(observedNodes, allParents)

    # Loop over possible cardinalities of the L
    for (k in seq(0, length = 1 + min(subsetSizeControl, length(latentsToControl)))) {

      # Loop over all subsets L in latentsToControl with cardinality k, i.e. |L| = k
      for (L in subsetsOfSize(latentsToControl, k)) {

        # Allowed nodes for Z
        childrenOfL = graph$children(L)
        allowedForZ <- setdiff(intersect(solvedNodes,childrenOfL), c(i, observedParents))

        for (Z in subsetsOfSize(allowedForZ, k)) {

          # Allowed set for Y. This is now possibly since Z and L is known.
          # Consists of all nodes in maybeAllowedForY that are
          # - not in (Z or i)
          # - no siblings of (Z or i) avoiding L (all latent parents of Y that are at the same
          #                                time latent parents of (Z or i) have to be in L!)
          # - solved if half-trek reachable from (Z or i) avoiding L
          latentParentsOfZandI <- intersect(latentNodes, graph$parents(c(i,Z)))
          latentParentsOfZandINotInL <- setdiff(latentParentsOfZandI, L)
          htrFromZandIAvoidingL <- union(graph$descendants(c(i,Z)),
                                         graph$trFrom(latentParentsOfZandINotInL))
          allowed <- setdiff(observedNodes,
                             setdiff(htrFromZandIAvoidingL, solvedNodes))
          allowed <- setdiff(allowed,
                             c(graph$children(latentParentsOfZandINotInL),
                               Z, i))

          # Check if there is a half-trek system from allowed nodes to the observed parents pa(i) and Z
          # Avoid starting with edge between two observed nodes.
          # Avoid ending a half-trek in Z with an observed edge.
          # This forces trek system to be half-trek system of correct form.
          avoidRightEdges <- matrix(edgesBetweenObserved[is.element(edgesBetweenObserved[,2], Z),], ncol=2)
          trekSystemResults <- graph$getTrekSystem(
            allowed, c(observedParents, Z),
            avoidLeftEdges = edgesBetweenObserved, avoidRightEdges = avoidRightEdges)

          # If half-trek system exists we know that all edges between pa(i) and i are identified
          if (trekSystemResults$systemExists) {
            identifiedEdges <- c(identifiedEdges, as.integer(rbind(observedParents, i)))
            Y <- trekSystemResults$activeFrom
            identifier <- createLFHtcIdentifier(identifier, v = i,
                                                Y = Y,
                                                Z = Z,
                                                parents = observedParents,
                                                reachableY = htrFromZandIAvoidingL) # half-trek reachable elements of Y are known
            solvedParents[[i]] <- observedParents
            unsolvedParents[[i]] <- integer(0)
            activeFroms[[i]] <- Y
            Zs[[i]] <- Z
            Ls[[i]] <- L
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
  return(list(identifiedEdges = matrix(identifiedEdges, byrow = T, ncol = 2),
              unsolvedParents = unsolvedParents,
              solvedParents = solvedParents,
              activeFroms = activeFroms,
              Zs = Zs,
              Ls = Ls,
              identifier = identifier))
}




#' Create an latent identifier base case
#'
#' Identifiers are functions that take as input a covariance matrix Sigma
#' corresponding to some latent digraph G and, from that covariance matrix,
#' identify some subset of the coefficients in the latent digraph G. This function
#' takes as input the matrix L defining G and creates an identifier
#' that does not identify any of the coefficients of G. This is useful as a
#' base case when building more complex identification functions.
#'
#' @export
#'
#' @param graph a \code{\link{LatentDigraph}} object representing
#'         the latent factor graph. All latent nodes in this graph should be
#'         source nodes (i.e. have no parents).
#'
#' @return a function that takes as input a covariance matrix compatible with
#'         the latent digraph defined by L and returns a list with two
#'         named components:
#'         Lambda - a matrix equal to L but with NA values instead of 1s,
#'         Omega - a matrix equal to O but with NA values for coefficients not equal to zero.
#'         When building more complex identifiers these NAs will be replaced
#'         by the value that can be identified from Sigma.
createLFIdentifierBaseCase <- function(graph) {

  L <- graph$L()
  observedNodes <- graph$observedNodes()
  latentNodes <- graph$latentNodes()
  numObservedNodes <- length(observedNodes)
  childrenOfLatentNodes <- lapply(latentNodes, graph$children)
  indicesOmega <- sapply(childrenOfLatentNodes, combn, 2, simplify=FALSE)
  indicesOmega <- t(do.call(cbind, indicesOmega))
  indicesOmega <- indicesOmega[!duplicated(indicesOmega),]

  return(function(Sigma) {
    Lambda <- matrix(NA, numObservedNodes, numObservedNodes)
    Lambda[L[observedNodes, observedNodes] == 0] <- 0
    Omega <- matrix(0, numObservedNodes, numObservedNodes)
    Omega[indicesOmega] <- NA
    Omega[indicesOmega[,c(2,1)]] <- NA
    diag(Omega) <- NA
    return(list(Lambda = Lambda, Omega=Omega))
  })
}


#' Checks that a LatentDigraph has appropriate node numbering
#'
#' Checks that the input latent digraph has nodes numbered from 1
#' to latentDigraph$numObserved()+latentDigraph$numLatents(). The first latentDigraph$numObserved()
#' nodes correspond to the observed nodes in the graph, all other nodes are considered unobserved.
#' Throws an error if this is not true.
#'
#' @param graph a \code{\link{LatentDigraph}} object representing
#'         the latent factor graph. All latent nodes in this graph should be
#'         source nodes (i.e. have no parents).
latentDigraphHasSimpleNumbering <- function(graph) {
  observedNodes <- graph$observedNodes()
  latentNodes <- graph$latentNodes()
  if ((any(observedNodes != 1:graph$numObserved())) ||
      (any(latentNodes != (graph$numObserved()+1):(graph$numObserved()+graph$numLatents())))) {
    stop(paste("Currently only latent graphs whose vertices are numbered from 1",
               "to graph$numObserved()+graph$numLatents() in order are supported."))
  }
}



#' Determines which edges in a latent digraph are LF-HTC-identifiable.
#'
#' Uses the latent factor half-trek criterion to determine
#' which edges in a latent digraph are generically identifiable.
#'
#' @export
#'
#' @param graph a \code{\link{LatentDigraph}} object representing
#'         the latent factor graph. All latent nodes in this graph should be
#'         source nodes (i.e. have no parents).
#'
#' @return returns a list with 5 components:
#' \describe{
#'   \item{\code{solvedParents}}{a list whose ith element contains a vector
#'   containing the subsets of parents of node i for which the edge j->i could
#'   be shown to be generically identifiable.}
#'   \item{\code{unsolvedParents}}{as for \code{solvedParents} but for the
#'   unsolved parents.}
#'   \item{\code{identifier}}{a function that takes a (generic) covariance
#'   matrix corresponding to the graph and identifies the edges parameters
#'   from solvedParents and solvedSiblings. See \code{\link{htcIdentifyStep}}
#'   for a more in-depth discussion of identifier functions.}
#'   \item{\code{graph}}{a latent digraph object of the graph.}
#'   \item{\code{call}}{the call made to this function.}
#' }
#'
#' @references
#' TO BE WRITTEN
lfhtcID <- function(graph, version=1){

  # Check the graph
  latentDigraphHasSimpleNumbering(graph)

  unsolvedParents <- lapply(graph$observedNodes(), graph$observedParents)
  solvedParents <- rep(list(numeric(0)), length(graph$observedNodes()))
  activeFroms <- rep(list(numeric(0)), length(graph$observedNodes()))
  Zs <- rep(list(numeric(0)), length(graph$observedNodes()))
  Ls <- rep(list(numeric(0)), length(graph$observedNodes()))
  identifier <- createLFIdentifierBaseCase(graph)

  changeFlag <- T
  while (changeFlag) {
    if (version==1){
      idResult <- lfhtcIdentifyStep(graph, unsolvedParents, solvedParents, activeFroms, Zs, Ls, identifier)
    }
    else {
      idResult <- lfhtcIdentifyStep2(graph, unsolvedParents, solvedParents, activeFroms, Zs, Ls, identifier)
    }
    changeFlag <- (nrow(idResult$identifiedEdges) != 0)
    unsolvedParents <- idResult$unsolvedParents
    solvedParents <- idResult$solvedParents
    activeFroms <- idResult$activeFroms
    Zs <- idResult$Zs
    Ls <- idResult$Ls
    identifier <- idResult$identifier
  }

  result <- list()
  class(result) <- "LfhtcIDResult"
  result$solvedParents <- solvedParents
  result$unsolvedParents <- unsolvedParents
  result$identifier <- identifier
  result$graph <- graph
  result$activeFroms <- activeFroms
  result$Zs <- Zs
  result$Ls <- Ls
  result$call <- match.call()
  result$LatentDigraph <- graph
  return(result)
}






#' Prints a LfhtcIDResult object
#'
#' Prints a LfhtcIDResult object as returned by
#' \code{\link{lfhtcID}}. Invisibly returns its argument via
#' \code{\link{invisible}(x)} as most print functions do.
#'
#' @export
#'
#' @param x the LfhtcIDResult object
#' @param ... optional parameters, currently unused.
print.LfhtcIDResult <- function(x, ...) {
  cat("Call: ")
  print(x$call)

  observedNodes <- x$LatentDigraph$observedNodes()
  nObserved <- length(observedNodes)
  latentNodes <- x$LatentDigraph$latentNodes()
  nLatent <- length(latentNodes)
  solvedParents <- x$solvedParents

  cat(paste("\nLatent Digraph Info\n"))
  cat(paste("# observed nodes:", nObserved, "\n"))
  cat(paste("# latent nodes:", nLatent, "\n"))
  cat(paste("# total nr. of edges between observed nodes:", sum(x$LatentDigraph$L()[observedNodes,]), "\n"))

  cat(paste("\nGeneric Identifiability Summary\n"))
  cat(paste("# nr. of edges between observed nodes shown gen. identifiable:", length(unlist(solvedParents)),
            "\n"))

  cat("# gen. identifiable edges: ")
  edges <- character(min(length(unlist(solvedParents))/2, 11))
  k <- 0
  for (i in 1:nObserved) {
    if (length(solvedParents[[i]]) != 0) {
      for (j in solvedParents[[i]]) {
        k <- k + 1
        edges[k] <- paste(j, "->", i, sep = "")

        if (k == 10) {
          edges[11] <- "..."
          break
        }
      }
      if (k == 10) {
        break
      }
    }
  }
  if (length(edges) == 0) {
    cat("None\n")
  } else {
    cat(paste(paste(edges, collapse = ", "), "\n"))
  }

  invisible(x)
}




















