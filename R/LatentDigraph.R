#' A helper function to validate if input nodes are valid.
#'
#' Produces an error if outside bounds.
#'
#' @param nodes the input nodes, expected to be from the collection
#'        1:(number of nodes in the graph)
#' @param numNodes the number of observed nodes in the graph.
validateNodes <- function(nodes, numNodes) {
  if (any(is.na(nodes)) || (length(nodes) > 1 &&
                            (min(nodes) < 1 || max(nodes) > numNodes))) {
    stop(paste(
      "nodes given to function",
      deparse(sys.calls()[[sys.nframe() - 1]]),
      "must be nodes in the graph."
    ))
  }
}

#' A helper function to validate that there are no variable arguments
#'
#' Produces an error if there are variable arguments.
#'
#' @param ... the variable arguments
validateVarArgsEmpty <- function(...) {
  if (length(list(...)) != 0) {

    args = list(...)
    a = "\nThe arguments\n"
    for (i in 1:length(args)) {
      name = names(args)[i]
      if (length(name) == 0 || name == "") {
        a = paste(a, "  [[", i, "]]  =  ", toString(args[[i]]), "\n", sep = "")
      } else {
        a = paste(a, "  ", name, "  =  ", toString(args[[i]]), "\n", sep = "")
      }
    }
    a = paste(a, "are unused. Check for spelling errors.")
    stop(a)
  }
}


#' Construct LatentDigraphFixedOrder object
#'
#' Creates an object representing a directed graph with some number of
#' nodes which are latent (unobserved).
#'
#' @name LatentDigraphFixedOrder
#' @usage LatentDigraphFixedOrder(L = matrix(0,1,1), numObserved = nrow(L))
#' @export
#'
#' @param L see \code{\link{graphID}} for the appropriate form of L. The
#'        first numObserved rows of L correspond to the observed nodes in the
#'        graph, all other nodes are considered unobserved.
#' @param numObserved a non-negative integer representing the number of observed
#'        nodes in the graph.
#'
#' @return An object representing the LatentDigraphFixedOrder
setConstructorS3("LatentDigraphFixedOrder", function(L = matrix(0, 1, 1),
                                                     numObserved = nrow(L)) {
  validateMatrix(L)
  numObserved = as.integer(numObserved)

  # Sanity check
  if (length(numObserved) != 1 || numObserved < 0 || numObserved > nrow(L)) {
    stop("numObserved must be a nonnegative integer of size no more than nrow(L).")
  }

  digraph <- igraph::graph.adjacency(L, mode = "directed")

  R.oo::extend(
    R.oo::Object(),
    "LatentDigraphFixedOrder",
    .L = L,
    .numObserved = numObserved,
    .numLatents = nrow(L) - numObserved,
    .digraph = digraph
  )
})

#' Get directed adjacency matrix.
#'
#' @name L
#' @export L
#'
#' @param this the graph object
L <- function(this, ...) {
  UseMethod("L")
}

#' @rdname   L
#' @name     L.LatentDigraphFixedOrder
#' @export
setMethodS3("L", "LatentDigraphFixedOrder", function(this, ...) {
  validateVarArgsEmpty(...)
  return(this$.L)
}, appendVarArgs = F)

#' Number of nodes in the graph.
#'
#' @name numNodes
#' @export numNodes
#'
#' @param this the graph object
#'
numNodes <- function(this, ...) {
  UseMethod("numNodes")
}

#' @rdname   numNodes
#' @name     numNodes.LatentDigraphFixedOrder
#' @export
setMethodS3("numNodes", "LatentDigraphFixedOrder", function(this, ...) {
  validateVarArgsEmpty(...)
  return(nrow(this$.L))
}, appendVarArgs = F)

#' Number of observed nodes in the graph.
#'
#' @name numObserved
#' @export numObserved
#'
#' @param this the graph object
#' @param ... ignored
#'
numObserved <- function(this, ...) {
  UseMethod("numObserved")
}

#' @rdname   numObserved
#' @name     numObserved.LatentDigraphFixedOrder
#' @export
setMethodS3("numObserved", "LatentDigraphFixedOrder", function(this, ...) {
  validateVarArgsEmpty(...)
  return(this$.numObserved)
}, appendVarArgs = F)

#' Number of latent nodes in the graph.
#'
#' @name numLatents
#' @export numLatents
#'
#' @param this the graph object
#' @param ... ignored
#'
numLatents <- function(this, ...) {
  UseMethod("numLatents")
}

#' @rdname   numLatents
#' @name     numLatents.LatentDigraphFixedOrder
#' @export
setMethodS3("numLatents", "LatentDigraphFixedOrder", function(this, ...) {
  validateVarArgsEmpty(...)
  return(this$.numLatents)
}, appendVarArgs = F)

#' All parents of a collection of nodes.
#'
#' Returns all parents of the collection (does not necessarily include the
#' input nodes themselves unless they are parents of one another).
#'
#' @name parents
#' @export parents
#'
#' @param this the graph object.
#' @param nodes nodes the nodes of which to find the parents.
#'
#' @return the observed parents.
parents <- function(this, nodes, ...) {
  UseMethod("parents")
}

#' @param includeObserved if TRUE includes observed nodes in the returned set.
#' @param includeLatents if TRUE includes latent nodes in the returned set.
#'
#' @rdname   parents
#' @name     parents.LatentDigraphFixedOrder
#' @export
setMethodS3("parents", "LatentDigraphFixedOrder", function(this, nodes,
                                                           includeObserved = T,
                                                           includeLatents = T,
                                                           ...) {
  validateVarArgsEmpty(...)
  validateNodes(nodes, this$numNodes())
  nobs <- this$numObserved()
  nlat <- this$numLatents()
  includedNodes = c(
    if (includeObserved) { seq(1, length = nobs) } else { integer(0) },
    if (includeLatents) { seq(nobs + 1, length = nlat) } else { integer(0) })
  return(which(rowSums(this$.L[includedNodes, nodes, drop = F]) != 0))
}, appendVarArgs = F)

#' All children of a collection of nodes.
#'
#' Returns all children of the collection (does not necessarily include the
#' input nodes themselves unless they are parents of one another).
#'
#' @name children
#' @export children
#'
#' @param this the graph object.
#' @param nodes nodes the nodes of which to find the children
#'
#' @return the observed children
children <- function(this, nodes, ...) {
  UseMethod("children")
}

#' @inheritParams parents.LatentDigraphFixedOrder
#' @rdname   children
#' @name     children.LatentDigraphFixedOrder
#' @export
setMethodS3("children", "LatentDigraphFixedOrder", function(this, nodes,
                                                            includeObserved = T,
                                                            includeLatents = T,
                                                            ...) {
  validateVarArgsEmpty(...)
  validateNodes(nodes, this$numNodes())
  nobs <- this$numObserved()
  nlat <- this$numLatents()
  includedNodes = c(
    if (includeObserved) { seq(1, length = nobs) } else { integer(0) },
    if (includeLatents) { seq(nobs + 1, length = nlat) } else { integer(0) })
  return(which(colSums(this$.L[nodes, includedNodes, drop = F]) != 0))
}, appendVarArgs = F)

#' All ancestors of a collection of nodes
#'
#' Finds all the ancestors of a collection of nodes. These ancestors DO
#' include the nodes themselves (every node is considered an ancestor of itself).
#'
#' @name ancestors
#' @export ancestors
#'
#' @param this the graph object
#' @param nodes the nodes from which to find all ancestors
#'
#' @return the ancestors of the nodes in the observed part of the graph.
ancestors <- function(this, nodes, ...) {
  UseMethod("ancestors")
}

#' @inheritParams parents.LatentDigraphFixedOrder
#' @rdname   ancestors
#' @name     ancestors.LatentDigraphFixedOrder
#' @export
setMethodS3("ancestors", "LatentDigraphFixedOrder", function(this, nodes,
                                                             includeObserved = T,
                                                             includeLatents = T,
                                                             ...) {
  validateVarArgsEmpty(...)
  validateNodes(nodes, this$numNodes())
  ancestors = as.integer(unique(unlist(
    igraph::neighborhood(
      this$.digraph,
      nodes = nodes,
      order = this$numNodes(),
      mode = "in"
    )
  )))

  nobs <- this$numObserved()
  return(ancestors[(includeObserved | (ancestors > nobs)) &
                     (includeLatents | (ancestors <= nobs))])
}, appendVarArgs = F)

#' Get descendants of a collection of observed nodes
#'
#' Finds all descendants of a collection of nodes, this DOES include the nodes
#' themselves (every node is considered a descendant of itself).
#'
#' @name descendants
#' @export descendants
#'
#' @param this the graph object
#' @param nodes the nodes from which to get the descendants.
descendants <- function(this, nodes, ...) {
  UseMethod("descendants")
}

#' @inheritParams parents.LatentDigraphFixedOrder
#' @rdname   descendants
#' @name     descendants.LatentDigraphFixedOrder
#' @export
setMethodS3("descendants", "LatentDigraphFixedOrder", function(this, nodes,
                                                               includeObserved = T,
                                                               includeLatents = T,
                                                               ...) {
  validateVarArgsEmpty(...)
  validateNodes(nodes, this$numNodes())
  descendants <- unique(as.integer(unlist(
    igraph::neighborhood(
      this$.digraph,
      order = this$numNodes(),
      nodes = nodes,
      mode = "out"
    )
  )))
  nobs <- this$numObserved()
  return(descendants[(includeObserved | (descendants > nobs)) &
                     (includeLatents | (descendants <= nobs))])
}, appendVarArgs = F)

#' Helper function to create a graph encoding trek reachable relationships.
#' @name createTrGraph
#' @export createTrGraph
#'
#' @param this the graph object
#' @param ... ignored
createTrGraph <- function(this, ...) {
  UseMethod("createTrGraph")
}

#' @rdname   createTrGraph
#' @name     createTrGraph.LatentDigraphFixedOrder
#' @export
setMethodS3("createTrGraph", "LatentDigraphFixedOrder", function(this, ...) {
  validateVarArgsEmpty(...)
  numObs <- this$numObserved()
  numLats <- this$numLatents()
  m <- numObs + numLats
  adjMat <- matrix(0, 2 * m, 2 * m)

  # Left nodes point to each other in the opposite direction of L
  adjMat[1:m, 1:m] <- t(this$.L)

  # Left nodes point to their corresponding right nodes
  adjMat[cbind(1:m, m + 1:m)] <- 1

  # If i -> j then right i point to right j
  adjMat[m + 1:m, m + 1:m] <- this$.L

  return(igraph::graph.adjacency(adjMat, mode = "directed"))
}, appendVarArgs = F, private = T)

#' Trek reachable nodes.
#'
#' Gets all nodes that are trek reachable from a collection of nodes.
#'
#' @name trFrom
#' @export trFrom
#'
#' @param this the graph object
#' @param nodes the nodes from which to find trek-reachable nodes.
trFrom <- function(this, nodes, ...) {
  UseMethod("trFrom")
}

#' @inheritParams getTrekSystem
#' @inheritParams parents.LatentDigraphFixedOrder
#'
#' @rdname trFrom
#' @name   trFrom.LatentDigraphFixedOrder
#' @export
setMethodS3("trFrom", "LatentDigraphFixedOrder", function(this, nodes,
                                                          avoidLeftNodes = integer(0),
                                                          avoidRightNodes = integer(0),
                                                          includeObserved = T,
                                                          includeLatents = T,
                                                          ...) {
  validateVarArgsEmpty(...)
  numObs <- this$numObserved()
  numLats <- this$numLatents()
  m <- numObs + numLats
  validateNodes(c(nodes, avoidLeftNodes, avoidRightNodes), m)

  if (is.null(this$.trGraph)) {
    this$.trGraph <- this$createTrGraph()
  }

  allowed <- setdiff(1:(2 * m), c(avoidLeftNodes, avoidRightNodes + m))
  trFrom <-
    unique(as.integer(
      igraph::bfs(
        this$.trGraph,
        root = nodes,
        mode = "out",
        unreachable = F,
        restricted = allowed
      )$order
    ))
  trFrom <- trFrom[!is.na(trFrom)]
  trFrom[trFrom > m] <- trFrom[trFrom > m] - m
  trFrom <- unique(trFrom)
  return(trFrom[(includeObserved | (trFrom > numObs)) &
                (includeLatents | (trFrom <= numObs))])
}, appendVarArgs = F)

#' Helper function to create a flow graph.
#'
#' @name createTrekFlowGraph
#' @export createTrekFlowGraph
#'
#' @param this the graph object
#' @param ... ignored
createTrekFlowGraph <- function(this, ...) {
  UseMethod("createTrekFlowGraph")
}

#' @rdname   createTrekFlowGraph
#' @name     createTrekFlowGraph.LatentDigraphFixedOrder
#' @export
setMethodS3("createTrekFlowGraph", "LatentDigraphFixedOrder", function(this, ...) {
  validateVarArgsEmpty(...)
  if (is.null(this$.trGraph)) {
    this$.trGraph <- this$createTrGraph()
  }
  adjMat <- as.matrix(igraph::get.adjacency(this$.trGraph))

  # Create the flow graph from adjMat. All vertices and edges have capacity 1
  flowGraph <- FlowGraph(adjMat, 1, adjMat)

  return(flowGraph)
}, appendVarArgs = F, private = T)

#' Determines if a trek system exists in the mixed graph.
#'
#' @name getTrekSystem
#' @export getTrekSystem
#'
#' @param this the graph object
#' @param fromNodes the start nodes
#' @param toNodes the end nodes
#' @param avoidLeftNodes a collection of nodes to avoid on the left
#' @param avoidRightNodes a collection of nodes to avoid on the right
#' @param avoidLeftEdges a collection of edges between observed nodes
#'                         in the graph that should not be used on any right
#'                         hand side of any trek in the trek system.
#' @param avoidRightEdges a collection of edges between observed noes
#'                          in the graph that should not be used on any right
#'                          hand side of any trek in the trek system.
#' @param ... ignored
getTrekSystem <-
  function(this, fromNodes, toNodes,
           avoidLeftNodes = integer(0),
           avoidRightNodes = integer(0),
           avoidLeftEdges = integer(0),
           avoidRightEdges = integer(0), ...) {
    UseMethod("getTrekSystem")
  }

#' @rdname   getTrekSystem
#' @name     getTrekSystem.LatentDigraphFixedOrder
#' @export
setMethodS3("getTrekSystem", "LatentDigraphFixedOrder",
            function(this, fromNodes, toNodes,
                     avoidLeftNodes = integer(0),
                     avoidRightNodes = integer(0),
                     avoidLeftEdges = integer(0),
                     avoidRightEdges = integer(0),
                     ...) {
  validateVarArgsEmpty(...)
  numObs <- this$numObserved()
  numLats <- this$numLatents()
  m <- numObs + numLats
  validateNodes(c(fromNodes, toNodes, avoidLeftNodes,
                  avoidRightNodes, as.integer(avoidLeftEdges),
                  as.integer(avoidRightEdges)), m)

  if (is.null(this$.trekFlowGraph)) {
    this$.trekFlowGraph <- this$createTrekFlowGraph()
  }

  # Left edges
  if (is.vector(avoidLeftEdges)) {
    avoidLeftEdges <- matrix(avoidLeftEdges, byrow = T, ncol = 2)
  }
  if (nrow(avoidLeftEdges) != 0) {
    this$.trekFlowGraph$updateEdgeCapacities(t(avoidLeftEdges)[2:1,], 0)
  }

  # Right edges
  if (is.vector(avoidRightEdges)) {
    avoidRightEdges <- matrix(avoidRightEdges, byrow = T, ncol = 2)
  }
  if (nrow(avoidRightEdges) != 0) {
    this$.trekFlowGraph$updateEdgeCapacities(t(avoidRightEdges + m), 0)
  }

  # Left nodes
  if (length(avoidLeftNodes) != 0) {
    this$.trekFlowGraph$updateVertexCapacities(avoidLeftNodes, 0)
  }

  # Right nodes
  if (length(avoidRightNodes) != 0) {
    this$.trekFlowGraph$updateVertexCapacities(avoidRightNodes + m, 0)
  }

  flowResult <- this$.trekFlowGraph$flowBetween(fromNodes, m + toNodes)

  # Left edges
  if (nrow(avoidLeftEdges) != 0) {
    this$.trekFlowGraph$updateEdgeCapacities(t(avoidLeftEdges)[2:1,], 1)
  }

  # Right edges
  if (nrow(avoidRightEdges) != 0) {
    this$.trekFlowGraph$updateEdgeCapacities(t(avoidRightEdges + m), 1)
  }

  # Left nodes
  if (length(avoidLeftNodes) != 0) {
    this$.trekFlowGraph$updateVertexCapacities(avoidLeftNodes, 1)
  }

  # Right nodes
  if (length(avoidRightNodes) != 0) {
    this$.trekFlowGraph$updateVertexCapacities(avoidRightNodes + m, 1)
  }

  return(list(
    systemExists = (flowResult$value == length(toNodes)),
    activeFrom = flowResult$activeSources
  ))
}, appendVarArgs = F)

#' Strongly connected component
#'
#' Get the strongly connected component for a node i in the graph the graph.
#'
#' @name stronglyConnectedComponent
#' @export stronglyConnectedComponent
#'
#' @param this the graph object
#' @param node the node for which to get the strongly connected component.
stronglyConnectedComponent <- function(this, node, ...) {
  UseMethod("stronglyConnectedComponent")
}

#' @rdname   stronglyConnectedComponent
#' @name     stronglyConnectedComponent.LatentDigraphFixedOrder
#' @export
setMethodS3("stronglyConnectedComponent", "LatentDigraphFixedOrder",
            function(this, node, ...) {
  if (is.null(this$.stronglyConnectedMembership)) {
    this$.stronglyConnectedMembership <- igraph::components(this$.digraph,
                                                            "strong")$membership
  }
  return(which(this$.stronglyConnectedMembership[node] == this$.stronglyConnectedMembership))
}, appendVarArgs = F)


#######################################
### The LatentDigraph wrapper class ###
#######################################

#' Construct a LatentDigraph object
#'
#' Creates an object representing a latent factor graph. The methods that are
#' currently available to be used on the latent factor graph include
#' \enumerate{
#' \item numObserved
#' \item numLatents
#' \item numNodes
#' \item toIn
#' \item toEx
#' \item L
#' \item observedNodes
#' \item latentNodes
#' \item parents
#' \item children
#' \item ancestors
#' \item descendants
#' \item trFrom
#' \item getTrekSystem
#' \item inducedSubgraph
#' \item stronglyConnectedComponent
#' \item plot
#' \item observedParents
#' \item getMixedGraph
#' }
#' see the individual function documentation for more information.
#'
#' @name LatentDigraph
#' @usage LatentDigraph(L = matrix(0,1,1),
#'                      observedNodes = seq(1, length = nrow(L)),
#'                      latentNodes = integer(0))
#' @export LatentDigraph
#'
#' @param L see \code{\link{graphID}} for the appropriate form of L.
#' @param observedNodes a vector of positive integers representing
#'        the vertex numbers of the observed nodes. These will correspond,
#'        in order, to the first length(observedNodes) rows of L.
#' @param latentNodes a vector of positive integers representing
#'        the vertex numbers of the latent nodes. These will correspond,
#'        in order, to the last length(latentNodes) rows of L.
#'
#' @return An object representing the LatentDigraph
setConstructorS3("LatentDigraph", function(L = matrix(0, 1, 1),
                                     observedNodes = seq(1, length = nrow(L)),
                                     latentNodes = integer(0)) {
  # Sanity check
  vertexNums = c(observedNodes, latentNodes)
  if (nrow(L) == 0) {
    vertexNums <- as.integer(NA)
    vertexNumsToInternal <- as.integer(NA)
  } else if (any(vertexNums %% 1 != 0) || any(vertexNums < 1) ||
             length(unique(vertexNums)) != length(vertexNums) ||
             length(vertexNums) != nrow(L)) {
    stop(
      paste(
        "observedNodes and latentNodes must be all unique positive",
        "integers and must have total length == nrow(L)."
      )
    )
  } else {
    vertexNumsToInternal <- rep(NA, max(vertexNums))
    vertexNumsToInternal[vertexNums] <- 1:nrow(L)
  }

  internalGraph <- LatentDigraphFixedOrder(L, length(observedNodes))
  # in internalGraph, nodes 1:length(observedNodes) are always the observed nodes

  R.oo::extend(
    R.oo::Object(),
    "LatentDigraph",
    .L = L,
    .observedNodes = observedNodes,
    .latentNodes = latentNodes,
    .vertexNums = vertexNums,
    .vertexNumsToInternal = vertexNumsToInternal,
    .internalGraph = internalGraph
  )
})

#' @rdname   numObserved
#' @name     numObserved.LatentDigraph
#' @export
setMethodS3("numObserved", "LatentDigraph", function(this, ...) {
  validateVarArgsEmpty(...)
  return(length(this$.observedNodes))
}, appendVarArgs = F)

#' @rdname   numLatents
#' @name     numLatents.LatentDigraph
#' @export
setMethodS3("numLatents", "LatentDigraph", function(this, ...) {
  validateVarArgsEmpty(...)
  return(length(this$.latentNodes))
}, appendVarArgs = F)

#' @rdname   numNodes
#' @name     numNodes.LatentDigraph
#' @export
setMethodS3("numNodes", "LatentDigraph", function(this, ...) {
  validateVarArgsEmpty(...)
  return(nrow(this$.L))
}, appendVarArgs = F)

#' Transforms a vector of given node indices into their internal numbering
#'
#' @name toIn
#' @export toIn
#'
#' @param this the graph object
#' @param nodes the nodes to transform
#' @param ... ignored
toIn <- function(this, nodes, ...) {
  UseMethod("toIn")
}

#' @rdname   toIn
#' @name     toIn.LatentDigraph
#' @export
setMethodS3("toIn", "LatentDigraph", function(this, nodes, ...) {
  validateVarArgsEmpty(...)
  if (any(nodes <= 0)) {
    stop("All nodes must be positive integers")
  }
  if (is.vector(nodes) || is.null(nodes)) {
    return(this$.vertexNumsToInternal[nodes])
  } else if (is.matrix(nodes)) {
    return(matrix(this$.vertexNumsToInternal[nodes], nrow = nrow(nodes)))
  } else {
    stop("Unsupported input nodes.")
  }

}, appendVarArgs = F, private = T)

#' Transforms a vector of node indices in the internal rep. into external numbering
#'
#' @name toEx
#' @export toEx
#'
#' @param this the graph object
#' @param nodes the nodes to transform
#' @param ... ignored
toEx <- function(this, nodes, ...) {
  UseMethod("toEx")
}

#' @rdname   toEx
#' @name     toEx.LatentDigraph
#' @export
setMethodS3("toEx", "LatentDigraph", function(this, nodes, ...) {
  validateVarArgsEmpty(...)
  if (any(nodes <= 0)) {
    stop("All nodes must be positive integers")
  }
  if (is.vector(nodes) || is.null(nodes)) {
    return(this$.vertexNums[nodes])
  } else if (is.matrix(nodes)) {
    return(matrix(this$.vertexNums[nodes], nrow = nrow(nodes)))
  } else {
    stop("Unsupported input nodes.")
  }
}, appendVarArgs = F, private = T)

#' @rdname   L
#' @name     L.LatentDigraph
#' @export
setMethodS3("L", "LatentDigraph", function(this, ...) {
  validateVarArgsEmpty(...)
  return(this$.L)
}, appendVarArgs = F)

#' Get all observed nodes in the graph.
#'
#' @name observedNodes
#' @export observedNodes
#'
#' @param this the graph object
#' @param ... ignored
observedNodes <- function(this, ...) {
  UseMethod("observedNodes")
}

#' @rdname   observedNodes
#' @name     observedNodes.LatentDigraph
#' @export
setMethodS3("observedNodes", "LatentDigraph", function(this, ...) {
  validateVarArgsEmpty(...)
  return(this$.observedNodes)
}, appendVarArgs = F)

#' Get all latent nodes in the graph.
#'
#' @name latentNodes
#' @export latentNodes
#'
#' @param this the graph object
#' @param ... ignored
latentNodes <- function(this, ...) {
  UseMethod("latentNodes")
}

#' @rdname   latentNodes
#' @name     latentNodes.LatentDigraph
#' @export
setMethodS3("latentNodes", "LatentDigraph", function(this, ...) {
  validateVarArgsEmpty(...)
  return(this$.latentNodes)
}, appendVarArgs = F)

#' @rdname   parents
#' @name     parents.LatentDigraph
#' @export
setMethodS3("parents", "LatentDigraph", function(this, nodes,
                                                 includeObserved = T,
                                                 includeLatents = T,
                                                 ...) {
  return(this$toEx(this$.internalGraph$parents(this$toIn(nodes),
                                               includeObserved = includeObserved,
                                               includeLatents = includeLatents,
                                               ...)))
}, appendVarArgs = F)

#' @inheritParams parents.LatentDigraphFixedOrder
#' @rdname   children
#' @name     children.LatentDigraph
#' @export
setMethodS3("children", "LatentDigraph", function(this, nodes,
                                                  includeObserved = T,
                                                  includeLatents = T,
                                                  ...) {
  return(this$toEx(this$.internalGraph$children(this$toIn(nodes),
                                                includeObserved = includeObserved,
                                                includeLatents = includeLatents,
                                                ...)))
}, appendVarArgs = F)

#' @inheritParams parents.LatentDigraphFixedOrder
#' @rdname   ancestors
#' @name     ancestors.LatentDigraph
#' @export
setMethodS3("ancestors", "LatentDigraph", function(this, nodes,
                                                   includeObserved = T,
                                                   includeLatents = T,
                                                   ...) {
  return(this$toEx(this$.internalGraph$ancestors(this$toIn(nodes),
                                                 includeObserved = includeObserved,
                                                 includeLatents = includeLatents,
                                                 ...)))
}, appendVarArgs = F)

#' @rdname   descendants
#' @name     descendants.LatentDigraph
#' @export
setMethodS3("descendants", "LatentDigraph", function(this, nodes,
                                                     includeObserved = T,
                                                     includeLatents = T,
                                                     ...) {
  return(this$toEx(this$.internalGraph$descendants(this$toIn(nodes),
                                                   includeObserved = includeObserved,
                                                   includeLatents = includeLatents,
                                                   ...)))
}, appendVarArgs = F)

#' @inheritParams parents.LatentDigraphFixedOrder
#' @rdname   trFrom
#' @name     trFrom.LatentDigraph
#' @export
setMethodS3("trFrom", "LatentDigraph", function(this, nodes,
                                                avoidLeftNodes = integer(0),
                                                avoidRightNodes = integer(0),
                                                includeObserved = T,
                                                includeLatents = T,
                                                ...) {
  return(this$toEx(
    this$.internalGraph$trFrom(this$toIn(nodes),
                               avoidLeftNodes = this$toIn(avoidLeftNodes),
                               avoidRightNodes = this$toIn(avoidRightNodes),
                               includeObserved = includeObserved,
                               includeLatents = includeLatents,
                               ...)
  ))
}, appendVarArgs = F)

#' @rdname   getTrekSystem
#' @name     getTrekSystem.LatentDigraph
#' @export
setMethodS3("getTrekSystem", "LatentDigraph", function(this,
                                                 fromNodes,
                                                 toNodes,
                                                 avoidLeftNodes = integer(0),
                                                 avoidRightNodes = integer(0),
                                                 avoidLeftEdges = integer(0),
                                                 avoidRightEdges = integer(0),
                                                 ...) {
  l <-
    this$.internalGraph$getTrekSystem(
      this$toIn(fromNodes),
      this$toIn(toNodes),
      avoidLeftNodes = this$toIn(avoidLeftNodes),
      avoidRightNodes = this$toIn(avoidRightNodes),
      avoidLeftEdges = this$toIn(avoidLeftEdges),
      avoidRightEdges = this$toIn(avoidRightEdges),
      ...
    )
  return(list(
    systemExists = l$systemExists,
    activeFrom = this$toEx(l$activeFrom)
  ))
}, appendVarArgs = F)

#' Get the induced subgraph on a collection of nodes
#'
#' @name inducedSubgraph
#' @export inducedSubgraph
#'
#' @param this the graph object
#' @param nodes the nodes on which to create the induced subgraph.
inducedSubgraph <- function(this, nodes, ...) {
  UseMethod("inducedSubgraph")
}

#' @rdname   inducedSubgraph
#' @name     inducedSubgraph.LatentDigraph
#' @export
setMethodS3("inducedSubgraph", "LatentDigraph", function(this, nodes, ...) {
  validateVarArgsEmpty(...)
  observed <- intersect(nodes, this$observedNodes())
  latents <- intersect(nodes, this$latentNodes())
  if (length(unique(nodes)) != length(observed) + length(latents)) {
    stop("Input nodes to inducedSubgraph are not a part of the graph")
  }
  observedIn <- this$toIn(observed)
  latentsIn <- this$toIn(latents)
  nodesIn <- c(observedIn, latentsIn)
  newL <- this$.internalGraph$.L[nodesIn, nodesIn, drop = F]
  return(LatentDigraph(newL, observed, latents))
}, appendVarArgs = F)

#' @rdname   stronglyConnectedComponent
#' @name     stronglyConnectedComponent.LatentDigraph
#' @export
setMethodS3("stronglyConnectedComponent", "LatentDigraph", function(this, node, ...) {
  return(this$toEx(this$.internalGraph$stronglyConnectedComponent(this$toIn(node))))
}, appendVarArgs = F)

#' Plots the latent digraph
#'
#' @param x the LatentDigraph object
#' @param ... additional plotting arguments. Currently ignored.
#'
#' @rdname   plot
#' @name     plot.LatentDigraph
#' @export
setMethodS3("plot", "LatentDigraph", function(x, ...) {
  plotLatentDigraph(x$L(), x$observedNodes(), x$latentNodes())
}, appendVarArgs = F)

#' Get the observed parents on a collection of nodes
#'
#' @name observedParents
#' @export observedParents
#'
#' @param this the graph object
#' @param nodes the nodes on which to get the observed parents
#' @param ... ignored
observedParents <- function(this, nodes, ...) {
  UseMethod("observedParents")
}

#' @rdname   observedParents
#' @name     observedParents.LatentDigraph
#' @export
setMethodS3("observedParents", "LatentDigraph", function(this, nodes, ...) {
  parents = this$parents(nodes)
  observedParents <- unique(parents[! parents %in% this$latentNodes()])
  return(observedParents)
}, appendVarArgs = F)




#' Get the corresponding mixed graph
#'
#' Only works for graphs where the latent nodes are source nodes
#'
#' @name getMixedGraph
#' @export getMixedGraph
#'
#' @param this the LatentDigraph object
#' @param ... ignored
getMixedGraph <- function(this, ...) {
  UseMethod("getMixedGraph")
}

#' @rdname   getMixedGraph
#' @name     getMixedGraph.LatentDigraph
#' @export
setMethodS3("getMixedGraph", "LatentDigraph", function(this, ...) {

  if (length(this$parents(this$latentNodes())) != 0){
    stop("In LatentDigraph: mixed graphs are only available when latent nodes have no parents.")
  }

  # Remove observed edges
  latentL <- this$.L

  observedNodes <- this$toIn(this$.observedNodes)
  latentNodes <- this$toIn(this$.latentNodes)

  latentL[observedNodes, observedNodes] <- 0

  latentGraph <-  LatentDigraph(latentL, observedNodes, latentNodes)
  reachableByLatentTrek = lapply(observedNodes, latentGraph$trFrom, includeLatents = FALSE)

  O = matrix(0, length(observedNodes), length(observedNodes))
  if (length(observedNodes) >= 1){
    for (i in 1:length(observedNodes)){
      if (length(reachableByLatentTrek[[i]]) >= 1){
        for (j in 1:length(reachableByLatentTrek[[i]])){
          O[observedNodes[i],reachableByLatentTrek[[i]][j]] <- 1
        }
      }
    }
  }

  diag(O) <- 0


  return(MixedGraph(as.matrix(this$.L[observedNodes, observedNodes]), O, vertexNums=this$.observedNodes))
}, appendVarArgs = F)
