#' Construct MixedGraph object
#'
#' Creates an object representing a mixed graph. The methods that are currently
#' available to be used on the mixed graph include
#' \enumerate{
#' \item ancestors
#' \item descendants
#' \item parents
#' \item siblings
#' \item isSibling
#' \item htrFrom
#' \item trFrom
#' \item getHalfTrekSystem
#' \item getTrekSystem
#' \item inducedSubgraph
#' \item L
#' \item O
#' \item nodes
#' \item numNodes
#' \item stronglyConnectedComponent
#' \item tianComponent
#' \item tianDecompose
#' }
#' see the individual function documentation for more information.
#'
#' @name MixedGraph
#'
#' @export MixedGraph
#'
#' @param L see \code{\link{graphID}} for the appropriate form of L.
#' @param O as for L.
#' @param vertexNums the labeling of the vertices in the graph in the order
#'        of the rows of L and O. Labels must be positive integers.
#'
#' @return An object representing the MixedGraph
setConstructorS3("MixedGraph", function(L = matrix(0, 1, 1),
                                        O = matrix(0, 1, 1),
                                        vertexNums = seq(1, length = nrow(L))) {
  validateMatrices(L, O)
  numObserved = nrow(L)
  numLatents = sum(O) / 2
  numNodes = numObserved + numLatents
  LwithLatents = matrix(0, nrow = numNodes, numNodes)
  LwithLatents[seq(1, length = numObserved), seq(1, length = numObserved)] = L

  siblings = which((upper.tri(O) * O) == 1, arr.ind = T)
  for (i in seq(1, length = nrow(siblings))) {
    LwithLatents[numObserved + i, siblings[i,]] = 1
  }
  if (length(vertexNums)==0){
    latentNodes = integer(0)
  } else {
    latentNodes = seq(max(vertexNums) + 1, length = numLatents)
  }
  R.oo::extend(R.oo::Object(),
               "MixedGraph",
               .internalGraph = LatentDigraph(LwithLatents,
                                              observedNodes = vertexNums,
                                              latentNodes = latentNodes),
               .vertexNums = vertexNums,
               .L = L,
               .O = O)
})

#' @rdname   toIn
#' @name     toIn.MixedGraph
#' @export
setMethodS3("toIn", "MixedGraph", function(this, nodes, ...) {
  return(this$.internalGraph$toIn(nodes, ...))
}, appendVarArgs = F, private = T)

#' @rdname   toEx
#' @name     toEx.LatentDigraph
#' @export
setMethodS3("toEx", "MixedGraph", function(this, nodes, ...) {
  return(this$.internalGraph$toEx(nodes, ...))
}, appendVarArgs = F, private = T)

#' @rdname   L
#' @name     L.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("L", "MixedGraph", function(this, ...) {
  return(this$.L)
}, appendVarArgs = F)

#' Get adjacency matrix for bidirected part.
#'
#' @name O
#' @export O
#' @param ... ignored.
#' @param this the mixed graph object
O <- function(this, ...) {
    UseMethod("O")
}

#' @rdname   O
#' @name     O.MixedGraph
#' @export
setMethodS3("O", "MixedGraph", function(this, ...) {
    return(this$.O)
}, appendVarArgs = F)

#' Get all nodes in the graph.
#'
#' @name nodes
#' @export nodes
#'
#' @param this the mixed graph object
#' @param ... ignored.
nodes <- function(this, ...) {
    UseMethod("nodes")
}

#' @rdname   nodes
#' @name     nodes.MixedGraph
#' @export
setMethodS3("nodes", "MixedGraph", function(this, ...) {
    return(this$.vertexNums)
}, appendVarArgs = F)

#' @rdname   numNodes
#' @name     numNodes.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("numNodes", "MixedGraph", function(this, ...) {
    return(nrow(this$.L))
}, appendVarArgs = F)

#' All siblings of a collection of nodes
#'
#' @name siblings
#' @export siblings
#'
#' @param this the mixed graph object
#' @param nodes a vector of nodes of which to find the siblings.
#' @param ... ignored.
#' @return a vector of all of the siblings.
siblings <- function(this, nodes, ...) {
  UseMethod("siblings")
}

#' @rdname   siblings
#' @name     siblings.MixedGraph
#' @export
setMethodS3("siblings", "MixedGraph", function(this, nodes, ...) {
  return(this$toEx(which(rowSums(this$.O[, this$toIn(nodes), drop = F]) != 0)))
}, appendVarArgs = F)

#' Are two nodes siblings?
#'
#' @name isSibling
#' @export isSibling
#'
#' @param this the mixed graph object
#' @param node1 a node
#' @param node2 a second node
#' @param ... ignored.
#' @return TRUE if the nodes are siblings in the graph, FALSE otherwise
isSibling <- function(this, node1, node2, ...) {
  UseMethod("isSibling")
}

#' @rdname   isSibling
#' @name     isSibling.MixedGraph
#' @export
setMethodS3("isSibling", "MixedGraph", function(this, node1, node2, ...) {
    return(this$.O[this$toIn(node1), this$toIn(node2)] != 0)
}, appendVarArgs = F)

#' @rdname   parents
#' @name     parents.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("parents", "MixedGraph", function(this, nodes, ...) {
    return(this$.internalGraph$parents(nodes, includeLatents = F, ...))
}, appendVarArgs = F)

#' @rdname   children
#' @name     children.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("children", "MixedGraph", function(this, nodes, ...) {
  return(this$.internalGraph$children(nodes, includeLatents = F, ...))
}, appendVarArgs = F)

#' @rdname   ancestors
#' @name     ancestors.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("ancestors", "MixedGraph", function(this, nodes, ...) {
  return(this$.internalGraph$ancestors(nodes, includeLatents = F, ...))
}, appendVarArgs = F)

#' @rdname   descendants
#' @name     descendants.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("descendants", "MixedGraph", function(this, nodes, ...) {
  return(this$.internalGraph$descendants(nodes, includeLatents = F, ...))
}, appendVarArgs = F)

#' Half trek reachable nodes.
#'
#' @name htrFrom
#' @export htrFrom
#' @param this the mixed graph object
#' @param nodes the nodes from which to get all half-trek reachable nodes.
#' @param avoidLeftNodes a collection of nodes to avoid on the left
#' @param avoidRightNodes a collection of nodes to avoid on the right
#' @param ... ignored.
#' @return a vector of all nodes half-trek reachable from node.
htrFrom <- function(this, nodes, ...) {
  UseMethod("htrFrom")
}

#' @rdname   htrFrom
#' @name     htrFrom.MixedGraph
#' @export
setMethodS3("htrFrom", "MixedGraph", function(this, nodes,
                                              avoidLeftNodes = integer(0),
                                              avoidRightNodes = integer(0),
                                              ...) {
  avoidLeftNodes = c(avoidLeftNodes, setdiff(this$nodes(), nodes))
  return(this$trFrom(nodes,
                     avoidLeftNodes = avoidLeftNodes,
                     avoidRightNodes = avoidRightNodes))
}, appendVarArgs = F)

#' @rdname   trFrom
#' @name     trFrom.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("trFrom", "MixedGraph", function(this, nodes,
                                             avoidLeftNodes = integer(0),
                                             avoidRightNodes = integer(0),
                                             ...) {
  return(this$.internalGraph$trFrom(nodes,
                                    avoidLeftNodes = avoidLeftNodes,
                                    avoidRightNodes = avoidRightNodes,
                                    includeLatents = F))
}, appendVarArgs = F)

#' Determines if a half-trek system exists in the mixed graph.
#'
#' @name getHalfTrekSystem
#' @export getHalfTrekSystem
#'
#' @param this the mixed graph object
#' @param fromNodes the nodes from which the half-trek system should start.
#'        If length(fromNodes) > length(toNodes) will find if there exists
#'        any half-trek system from any subset of fromNodes of size
#'        length(toNodes) to toNodes.
#' @param toNodes the nodes where the half-trek system should end.
#' @param ... ignored.
#'
#' @return a list with two named components, \code{systemExists} (TRUE if a
#'         system exists, FALSE otherwise) and \code{activeFrom} (the subset
#'         of fromNodes from which the maximal half-trek system was started).
getHalfTrekSystem <- function(this, fromNodes, toNodes, ...) {
  UseMethod("getHalfTrekSystem")
}

#' @inheritParams getTrekSystem.LatentDigraph
#' @rdname   getHalfTrekSystem
#' @name     getHalfTrekSystem.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("getHalfTrekSystem", "MixedGraph", function(this,
                                                        fromNodes,
                                                        toNodes,
                                                        avoidLeftNodes = integer(0),
                                                        avoidRightNodes = integer(0),
                                                        avoidRightEdges = integer(0),
                                                        ...) {
  if (is.null(this$.edgesBetweenObserved)) {
    observedNodes <- this$nodes()
    edgesBetweenObserved <- which(this$L() == 1, arr.ind = T)
    this$.edgesBetweenObserved <- matrix(observedNodes[edgesBetweenObserved], ncol = 2)
  }
  return(
    this$getTrekSystem(fromNodes, toNodes,
                       avoidLeftNodes = avoidLeftNodes,
                       avoidRightNodes = avoidRightNodes,
                       avoidLeftEdges = this$.edgesBetweenObserved,
                       avoidRightEdges = avoidRightEdges)
  )
}, appendVarArgs = F)

#' @inheritParams getTrekSystem.LatentDigraph
#' @rdname   getHalfTrekSystem
#' @name     getHalfTrekSystem.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("getTrekSystem", "MixedGraph", function(this,
                                                    fromNodes,
                                                    toNodes,
                                                    avoidLeftNodes = integer(0),
                                                    avoidRightNodes = integer(0),
                                                    avoidLeftEdges = integer(0),
                                                    avoidRightEdges = integer(0),
                                                    ...) {

  return(
    this$.internalGraph$getTrekSystem(fromNodes, toNodes,
                                      avoidLeftNodes = avoidLeftNodes,
                                      avoidRightNodes = avoidRightNodes,
                                      avoidLeftEdges = avoidLeftEdges,
                                      avoidRightEdges = avoidRightEdges)
  )
}, appendVarArgs = F)

#' @rdname   stronglyConnectedComponent
#' @name     stronglyConnectedComponent.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("stronglyConnectedComponent", "MixedGraph", function(this, node, ...) {
  return(this$.internalGraph$stronglyConnectedComponent(node))
}, appendVarArgs = F)

#' @rdname   inducedSubgraph
#' @name     inducedSubgraph.MixedGraph
#' @param ... ignored.
#' @export
setMethodS3("inducedSubgraph", "MixedGraph", function(this, nodes, ...) {
    nodesIn <- this$toIn(nodes)
    newL <- this$L()[nodesIn, nodesIn]
    newO <- this$O()[nodesIn, nodesIn]

    return(MixedGraph(newL, newO, vertexNums = nodes))
}, appendVarArgs = F)

#' Performs the tian decomposition on the mixed graph
#'
#' Uses the Tian decomposition to break the mixed graph into c-components.
#' These c-components are slightly different than those from Tian (2005)
#' in that if they graph is not acyclic the bidirected components are
#' combined whenever they are connected by a directed loop.
#'
#' @name tianDecompose
#' @export tianDecompose
#'
#' @param this the mixed graph object
#'
#' @references
#' Jin Tian. 2005. Identifying direct causal effects in linear models. In
#' \emph{Proceedings of the 20th national conference on Artificial intelligence
#' - Volume 1} (AAAI'05), Anthony Cohn (Ed.), Vol. 1. AAAI Press 346-352.
tianDecompose <- function(this) {
  UseMethod("tianDecompose")
}

#' @rdname   tianDecompose
#' @name     tianDecompose.MixedGraph
#' @export
setMethodS3("tianDecompose", "MixedGraph", function(this) {
  if (!is.null(this$.cComponents)) {
      return(this$.cComponents)
  }
  if (is.null(this$.strongComponentsInternalOrdering)) {
    nodes <- this$nodes()
    components <- list()
    while (length(nodes) != 0) {
      node <- nodes[1]
      component <- this$stronglyConnectedComponent(node)
      components <- c(components, list(this$toIn(component)))
      nodes <- setdiff(nodes, component)
    }
    this$.strongComponentsInternalOrdering <- components
  }
  components <- this$.strongComponentsInternalOrdering
  numComponents <- length(components)

  O = this$O()
  L = this$L()
  shrunkO <- matrix(0, numComponents, numComponents)
  shrunkL <- matrix(0, numComponents, numComponents)
  for (i in seq(1, length = numComponents - 1)) {
      for (j in seq(i + 1, length = numComponents - i)) {
          shrunkO[i, j] <- 1 * any(O[components[[i]], components[[j]]] != 0)
          shrunkL[i, j] <- 1 * any(L[components[[i]], components[[j]]] != 0)
          shrunkL[j, i] <- 1 * any(L[components[[j]], components[[i]]] != 0)
      }
  }
  shrunkO <- shrunkO + t(shrunkO)

  biGraph <- igraph::graph.adjacency(shrunkO, mode = "undirected")
  biComponents <- igraph::components(biGraph)$membership
  shrunkTopOrder <- as.integer(igraph::topological.sort(
    igraph::graph.adjacency(shrunkL, mode = "directed"), mode = "out"))
  topOrder <- unlist(components[shrunkTopOrder])

  cComponents <- rep(list(list()), max(biComponents))

  for (i in seq(1, length = max(biComponents))) {
      superNodes <- which(biComponents == i)
      internal <- topOrder[topOrder %in% unlist(components[superNodes])]
      incoming <- topOrder[topOrder %in%
                             setdiff(which(rowSums(L[,internal,drop = F]) != 0), internal)]
      allOrdered <- topOrder[topOrder %in% c(internal, incoming)]

      indsInt <- which(allOrdered %in% internal)
      indsInc <- which(allOrdered %in% incoming)

      newL <- matrix(0, length(allOrdered), length(allOrdered))
      newO <- newL
      newL[c(indsInt, indsInc), indsInt] <- L[c(internal, incoming), internal]
      newO[indsInt, indsInt] <- O[internal, internal]

      cComponents[[i]]$internal <- this$toEx(internal)
      cComponents[[i]]$incoming <- this$toEx(incoming)
      cComponents[[i]]$topOrder <- this$toEx(allOrdered)
      cComponents[[i]]$L <- newL
      cComponents[[i]]$O <- newO
  }

  this$.cComponents <- cComponents
  return(cComponents)
}, appendVarArgs = F)


#' Returns the Tian c-component of a node
#'
#' @name tianComponent
#' @export tianComponent
#'
#' @param this the mixed graph object
#' @param node the node for which to return its c-component
tianComponent <- function(this, node) {
    UseMethod("tianComponent")
}

#' @rdname   tianComponent
#' @name     tianComponent.MixedGraph
#' @export
setMethodS3("tianComponent", "MixedGraph", function(this, node) {
    cComponents <- this$tianDecompose()
    for (i in seq(1, length = length(cComponents))) {
        if (node %in% cComponents[[i]]$internal) {
            return(cComponents[[i]])
        }
    }
    stop("No tian component for node, was the node mispecified?")
}, appendVarArgs = F)

#' Plots the mixed graph
#'
#' @param x the mixed graph object
#' @param ... additional plotting arguments. Currently ignored.
#'
#' @rdname   plot
#' @name     plot.MixedGraph
#' @export
setMethodS3("plot", "MixedGraph", function(x, ...) {
    plotMixedGraph(x$L(), x$O(), vertexLabels = x$nodes())
}, appendVarArgs = F)

