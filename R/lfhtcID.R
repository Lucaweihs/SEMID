#' Create a latent-factor half-trek critierion identification function.
#'
#' A helper function for \code{\link{lfhtcIdentifyStep}}, creates an identifier
#' function based on its given parameters. This created identifier function will
#' identify the directed edges from 'targets' to 'node.'
#'
#' @export
#'
#' @param idFunc identification of edge coefficients often requires that other
#'        edge coefficients already be identified. This argument should be a
#'        function that produces all such identifications. The newly created
#'        identifier function will return these identifications along with its
#'        own.
#' @param v the node for which all incoming edges are to be identified
#'        (the tails of which are targets).
#' @param Y the sources of the latent-factor half-trek system.
#' @param Z the nodes that are reached from Y via an latent-factor half-trek of the form
#'        \code{y <- h -> z} where \code{h} is an element of L.
#' @param parents the parents of node v.
#' @param reachableY the nodes in Y which are latent-factor half-trek reachable
#'        from Z or v  by avoiding the nodes in L. All incoming edges to these
#'        nodes should be identified by idFunc the newly created identification function to work.
#'
#' @return an identification function
#'
#' @references
#' Barber, R. F., Drton, M., Sturma, N., and Weihs L. (2022).
#' Half-Trek Criterion for Identifiability of Latent Variable Models.
#' \emph{arXiv preprint} arXiv:2201.04457
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
#' Produces an error if not all latent nodes are sources.
#'
#' @export
#'
#' @param graph the LatentDigraph
validateLatentNodesAreSources <- function(graph) {
  for (node in graph$latentNodes()) {
    if (length(graph$parents(node)) != 0) {
      stop("Latent node in input graph is not a source node, i.e. it has a parent.")
    }
  }
}



#' Perform one iteration of latent-factor HTC identification.
#'
#' A function that does one step through all the nodes in a latent-factor graph
#' and tries to identify new edge coefficients using the existence of
#' latent-factor half-trek systems.
#'
#' @param graph a \code{\link{LatentDigraph}} object representing
#'         the latent-factor graph. All latent nodes in this graph should be
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
#'        to the latent-factor graph) and returns a list with two named arguments
#' @param activeFroms list. If node i is solved then the ith index is a vector
#'        containing the nodes Y otherwise it is empty.
#' @param Zs list. If node i is solved then the ith index is a vector
#'        containing the nodes Z otherwise it is empty.
#' @param Ls list. If node i is solved then the ith index is a vector
#'        containing the nodes Z otherwise it is empty.
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
#'   \item{\code{activeFroms}}{as the input argument but updated with
#'   any newly solved node}
#'   \item{\code{Zs}}{as the input argument but updated with
#'   any newly solved node}
#'   \item{\code{Ls}}{as the input argument but updated with
#'   any newly solved node}
#' }
#'
#' @export
#'
#' @references
#' Barber, R. F., Drton, M., Sturma, N., and Weihs L. (2022).
#' Half-Trek Criterion for Identifiability of Latent Variable Models.
#' \emph{arXiv preprint} arXiv:2201.04457
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
#' corresponding to some latent digraph \code{G} and, from that covariance matrix,
#' identify some subset of the coefficients coresponding to the direct causal effects
#' in the latent digraph \code{G}. This function
#' takes as input the digraph \code{G} and creates an identifier
#' that does not identify any of the direct causal effects. This is useful as a
#' base case when building more complex identification functions.
#'
#' @export
#'
#' @param graph a \code{\link{LatentDigraph}} object representing
#'         the latent-factor graph. All latent nodes in this graph should be
#'         source nodes (i.e. have no parents).
#'
#' @return a function that takes as input a covariance matrix compatible with
#'         the latent digraph defined by \code{L} and returns a list with two
#'         named components:
#'         \describe{
#'         \item{\code{Lambda}}{a matrix equal to the observed part of \code{graph$L()} but with \code{NA} values
#'         instead of 1s}
#'         \item{\code{Omega}}{a matrix equal to \code{graph$O()} but with \code{NA} values for coefficients
#'         not equal to zero.}
#'         }
#'         When building more complex identifiers these NAs will be replaced
#'         by the value that can be identified from the covariance matrix corresponding to \code{G}.
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
#' @export
#'
#' @param graph a \code{\link{LatentDigraph}} object representing
#'         the latent-factor graph. All latent nodes in this graph should be
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
#' Uses the latent-factor half-trek criterion to determine
#' which edges in a latent digraph are generically identifiable.
#'
#' @export
#'
#' @param graph a \code{\link{LatentDigraph}} object representing
#'         the latent-factor graph. All latent nodes in this graph should be
#'         source nodes (i.e. have no parents).
#'
#' @return returns a list with 8 components:
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
#'   \item{\code{activeFroms}}{list. If node i is solved then the ith index
#'   is a vector containing the nodes Y otherwise it is empty.}
#'   \item{\code{Zs}}{list. If node i is solved then the ith index is a
#'   vector containing the nodes Z otherwise it is empty.}
#'   \item{\code{Ls}}{list. If node i is solved then the ith index is a
#'   vector containing the nodes L otherwise it is empty.}
#' }
#'
#' @references
#' Barber, R. F., Drton, M., Sturma, N., and Weihs L. (2022).
#' Half-Trek Criterion for Identifiability of Latent Variable Models.
#' Ann. Statist. 50(6):3174--3196. <doi:10.1214/22-AOS2221>.
lfhtcID <- function(graph){

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
    idResult <- lfhtcIdentifyStep(graph, unsolvedParents, solvedParents, activeFroms, Zs, Ls, identifier)
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

  observedNodes <- x$graph$observedNodes()
  nObserved <- length(observedNodes)
  latentNodes <- x$graph$latentNodes()
  nLatent <- length(latentNodes)
  solvedParents <- x$solvedParents

  cat(paste("\nLatent Digraph Info\n"))
  cat(paste("# observed nodes:", nObserved, "\n"))
  cat(paste("# latent nodes:", nLatent, "\n"))
  cat(paste("# total nr. of edges between observed nodes:", sum(x$graph$L()[observedNodes,]), "\n"))

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
