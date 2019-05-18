### First we define a class, MixedGraphFixedOrder, representing mixed graphs whose
### vertex labeling is fixed. We then write a wrapper class, MixedGraph, which
### allows for different vertex labelings.

#' Construct MixedGraphFixedOrder object
#'
#' Creates an object representing a mixed graph.
#'
#' @name MixedGraphFixedOrder
#' @usage MixedGraphFixedOrder(L = matrix(0,1,1), O = matrix(0,1,1))
#' @export
#'
#' @param L see \code{\link{graphID}} for the appropriate form of L.
#' @param O as for L.
#'
#' @return An object representing the MixedGraphFixedOrder
setConstructorS3("MixedGraphFixedOrder", function(L = matrix(0, 1, 1), O = matrix(0,
    1, 1)) {
    validateMatrices(L, O)
    O <- 1 * ((O + t(O)) != 0)

    dirGraph <- igraph::graph.adjacency(L, mode = "directed")

    R.oo::extend(R.oo::Object(), "MixedGraphFixedOrder", .L = L, .O = O, .dirGraph = dirGraph)
})


#' Number of nodes in the graph.
#'
#' @name numNodes
#' @export numNodes
#'
#' @param this the mixed graph object
#'
numNodes <- function(this) {
    UseMethod("numNodes")
}

#' @rdname   numNodes
#' @name     numNodes.MixedGraphFixedOrder
#' @export
setMethodS3("numNodes", "MixedGraphFixedOrder", function(this) {
    return(nrow(this$.L))
}, appendVarArgs = F)

#' All siblings of a collection of nodes
#'
#' @name siblings
#' @export siblings
#'
#' @param this the mixed graph object
#' @param nodes a vector of nodes of which to find the siblings.
#'
#' @return a vector of all of the siblings.
siblings <- function(this, nodes) {
    UseMethod("siblings")
}

#' @rdname   siblings
#' @name     siblings.MixedGraphFixedOrder
#' @export
setMethodS3("siblings", "MixedGraphFixedOrder", function(this, nodes) {
    return(which(rowSums(this$.O[, nodes, drop = F]) != 0))
}, appendVarArgs = F)

#' Are two nodes siblings?
#'
#' @name isSibling
#' @export isSibling
#'
#' @param this the mixed graph object
#' @param node1 a node
#' @param node2 a second node
#'
#' @return TRUE if the nodes are siblings in the graph, FALSE otherwise
isSibling <- function(this, node1, node2) {
    UseMethod("isSibling")
}

#' @rdname   isSibling
#' @name     isSibling.MixedGraphFixedOrder
#' @export
setMethodS3("isSibling", "MixedGraphFixedOrder", function(this, node1,
    node2) {
    return(this$.O[node1, node2] != 0)
}, appendVarArgs = F)

#' All parents a collection of nodes.
#'
#' @name parents
#' @export parents
#'
#' @param this the mixed graph object.
#' @param nodes nodes the nodes of which to find the parents.
#'
#' @return a vector of parents of the nodes.
parents <- function(this, nodes) {
    UseMethod("parents")
}

#' @rdname   parents
#' @name     parents.MixedGraphFixedOrder
#' @export
setMethodS3("parents", "MixedGraphFixedOrder", function(this, nodes) {
    return(which(rowSums(this$.L[, nodes, drop = F]) != 0))
}, appendVarArgs = F)

#' All ancestors of a collection of nodes
#'
#' Finds all the ancestors of a collection of nodes. These ancestors DO include
#' the nodes themselves (every node is considered an ancestor of itself).
#'
#' @name ancestors
#' @export ancestors
#'
#' @param this the mixed graph object
#' @param nodes the nodes from which to find all ancestors
ancestors <- function(this, nodes) {
    UseMethod("ancestors")
}

#' @rdname   ancestors
#' @name     ancestors.MixedGraphFixedOrder
#' @export
setMethodS3("ancestors", "MixedGraphFixedOrder", function(this, nodes) {
    return(as.integer(unique(unlist(igraph::neighborhood(this$.dirGraph, nodes = nodes,
        order = length(igraph::V(this$.dirGraph)), mode = "in")))))
}, appendVarArgs = F)

#' Helper function to create a graph encoding htr relationships.
#'
#' @name createHtrGraph
#' @export createHtrGraph
#'
#' @param this the mixed graph object
createHtrGraph <- function(this) {
    UseMethod("createHtrGraph")
}

#' @rdname   createHtrGraph
#' @name     createHtrGraph.MixedGraphFixedOrder
#' @export
setMethodS3("createHtrGraph", "MixedGraphFixedOrder", function(this) {
    m <- this$numNodes()
    adjMat <- matrix(0, 2 * m, 2 * m)

    # Left nodes point to right nodes
    adjMat[cbind(1:m, m + 1:m)] <- 1

    # If i <-> j, then left i should point to right j and left j to right i
    adjMat[1:m, m + 1:m] <- this$.O + adjMat[1:m, m + 1:m]

    # If i -> j then right i should point to right j
    adjMat[m + 1:m, m + 1:m] <- this$.L

    return(igraph::graph.adjacency(adjMat, mode = "directed"))
}, appendVarArgs = F, private = TRUE)

#' Half trek reachable nodes.
#'
#' @name htrFrom
#' @export htrFrom
#' @param this the mixed graph object
#' @param node the node from which to get all half-trek reachable nodes.
#'
#' @return a vector of all nodes half-trek reachable from node.
htrFrom <- function(this, node) {
    UseMethod("htrFrom")
}

#' @rdname   htrFrom
#' @name     htrFrom.MixedGraphFixedOrder
#' @export
setMethodS3("htrFrom", "MixedGraphFixedOrder", function(this, node) {
    if (!is.null(this$.htrFrom)) {
        return(this$.htrFrom[[node]])
    }
    if (is.null(this$.htrGraph)) {
        this$.htrGraph <- this$createHtrGraph()
    }
    m <- this$numNodes()
    htrFrom <- igraph::neighborhood(this$.htrGraph, order = 2 * m, nodes = 1:m, mode = "out",
        mindist = 1)
    for (i in 1:m) {
        htrFrom[[i]] <- as.integer(htrFrom[[i]]) - m
    }
    this$.htrFrom <- htrFrom
    return(this$.htrFrom[[node]])
}, appendVarArgs = F)

#' Helper function to create a flow graph.
#'
#' @name createHalfTrekFlowGraph
#' @export createHalfTrekFlowGraph
#'
#' @param this the mixed graph object
createHalfTrekFlowGraph <- function(this) {
    UseMethod("createHalfTrekFlowGraph")
}

#' @rdname   createHalfTrekFlowGraph
#' @name     createHalfTrekFlowGraph.MixedGraphFixedOrder
#' @export
setMethodS3("createHalfTrekFlowGraph", "MixedGraphFixedOrder", function(this) {
    if (is.null(this$.htrGraph)) {
        this$.htrGraph <- this$createHtrGraph()
    }
    adjMat <- as.matrix(igraph::get.adjacency(this$.htrGraph))
    # Create the flow graph from adjMat with all vertices and edges having a capacity
    # of 1
    flowGraph <- FlowGraph(adjMat, rep(1, 2 * this$numNodes()), adjMat)
    return(flowGraph)
}, appendVarArgs = F, private = TRUE)

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
#'
#' @return a list with two named components, \code{systemExists} (TRUE if a
#'         system exists, FALSE otherwise) and \code{activeFrom} (the subset
#'         of fromNodes from which the maximal half-trek system was started).
getHalfTrekSystem <- function(this, fromNodes, toNodes) {
    UseMethod("getHalfTrekSystem")
}

#' @rdname   getHalfTrekSystem
#' @name     getHalfTrekSystem.MixedGraphFixedOrder
#' @export
setMethodS3("getHalfTrekSystem", "MixedGraphFixedOrder", function(this,
    fromNodes, toNodes) {
    if (is.null(this$.halfTrekFlowGraph)) {
        this$.halfTrekFlowGraph <- this$createHalfTrekFlowGraph()
    }
    flowResult <- this$.halfTrekFlowGraph$flowBetween(fromNodes, this$numNodes() +
        toNodes)
    return(list(systemExists = (flowResult$value == length(toNodes)), activeFrom = flowResult$activeSources))
}, appendVarArgs = F)


#' Helper function to create a graph encoding trek reachable relationships.
#' @name createTrGraph
#' @export createTrGraph
#'
#' @param this the mixed graph object
createTrGraph <- function(this) {
    UseMethod("createTrGraph")
}

#' @rdname   createTrGraph
#' @name     createTrGraph.MixedGraphFixedOrder
#' @export
setMethodS3("createTrGraph", "MixedGraphFixedOrder", function(this) {
    m <- this$numNodes()
    adjMat <- matrix(0, 2 * m, 2 * m)

    # Left nodes point to each other in the opposite direction of L
    adjMat[1:m, 1:m] <- t(this$.L)

    # Left nodes point to their corresponding right nodes
    adjMat[cbind(1:m, m + 1:m)] <- 1

    # If i <-> j, then left i should point to right j and left j to right i
    adjMat[1:m, m + 1:m] <- 1 * ((this$.O + adjMat[1:m, m + 1:m]) != 0)

    # If i -> j then right i point to right j
    adjMat[m + 1:m, m + 1:m] <- this$.L

    return(igraph::graph.adjacency(adjMat, mode = "directed"))
}, appendVarArgs = F, private = T)

#' Trek reachable nodes.
#'
#' Like \code{\link{htrFrom}} but for the nodes that are trek-reachable from a
#' node
#'
#' @name trFrom
#' @export trFrom
#'
#' @param this the mixed graph object
#' @param node the node from which to find trek-reachable nodes.
trFrom <- function(this, node) {
    UseMethod("trFrom")
}

#' @rdname   trFrom
#' @name     trFrom.MixedGraphFixedOrder
#' @export
setMethodS3("trFrom", "MixedGraphFixedOrder", function(this, node) {
    if (is.null(this$.trGraph)) {
        this$.trGraph <- this$createTrGraph()
    }
    m <- this$numNodes()
    trFrom <- as.integer(igraph::neighborhood(this$.trGraph, order = 2 * m, nodes = node,
        mode = "out")[[1]])
    trFrom[trFrom > m] <- trFrom[trFrom > m] - m
    return(unique(trFrom))
}, appendVarArgs = F)

#' Helper function to create a flow graph.
#'
#' @name createTrekFlowGraph
#' @export createTrekFlowGraph
#'
#' @param this the mixed graph object
createTrekFlowGraph <- function(this) {
    UseMethod("createTrekFlowGraph")
}

#' @rdname   createTrekFlowGraph
#' @name     createTrekFlowGraph.MixedGraphFixedOrder
#' @export
setMethodS3("createTrekFlowGraph", "MixedGraphFixedOrder", function(this) {
    if (is.null(this$.trGraph)) {
        this$.trGraph <- this$createTrGraph()
    }
    adjMat <- as.matrix(igraph::get.adjacency(this$.trGraph))

    # Create the flow graph from adjMat. All vertices and edges have capacity 1
    flowGraph <- FlowGraph(adjMat, rep(1, 2 * this$numNodes()), adjMat)

    return(flowGraph)
}, appendVarArgs = F, private = T)

#' Determines if a trek system exists in the mixed graph.
#'
#' @name getTrekSystem
#' @export getTrekSystem
#'
#' @param this the mixed graph object
#' @param fromNodes the start nodes
#' @param toNodes the end nodes
#' @param avoidEdgesOnRight a collection of edges in the graph that should not
#'                          be used on any right hand side of any trek in the
#'                          trek system.
getTrekSystem <- function(this, fromNodes, toNodes, avoidEdgesOnRight) {
    UseMethod("getTrekSystem")
}

#' @rdname   getTrekSystem
#' @name     getTrekSystem.MixedGraphFixedOrder
#' @export
setMethodS3("getTrekSystem", "MixedGraphFixedOrder", function(this,
    fromNodes, toNodes, avoidEdgesOnRight = NULL) {
    if (is.null(this$.trekFlowGraph)) {
        this$.trekFlowGraph <- this$createTrekFlowGraph()
    }

    if (length(avoidEdgesOnRight) != 0) {
        if (is.vector(avoidEdgesOnRight)) {
            avoidEdgesOnRight <- matrix(avoidEdgesOnRight, byrow = T, ncol = 2)
        }
        if (any(this$.L[avoidEdgesOnRight] == 0)) {
            stop("Some edge in avoidEdgesOnRight is not an edge in the graph")
        }
        this$.trekFlowGraph$updateEdgeCapacities(t(avoidEdgesOnRight + this$numNodes()),
            0)
    }
    flowResult <- this$.trekFlowGraph$flowBetween(fromNodes, this$numNodes() + toNodes)
    if (length(avoidEdgesOnRight) != 0) {
        this$.trekFlowGraph$updateEdgeCapacities(t(avoidEdgesOnRight + this$numNodes()),
            1)
    }
    return(list(systemExists = (flowResult$value == length(toNodes)), activeFrom = flowResult$activeSources))
}, appendVarArgs = F)

#' Strongly connected components
#'
#' Get the strongly connected components of a graph
#'
#' @name stronglyConnectedComponents
#' @export stronglyConnectedComponents
#'
#' @param this the mixed graph object
stronglyConnectedComponents <- function(this) {
    UseMethod("stronglyConnectedComponents")
}

#' @rdname   stronglyConnectedComponents
#' @name     stronglyConnectedComponents.MixedGraphFixedOrder
#' @export
setMethodS3("stronglyConnectedComponents", "MixedGraphFixedOrder", function(this) {
    if (is.null(this$.stronglyConnectedComponents)) {
        this$.stronglyConnectedComponents <- igraph::components(this$.dirGraph, "strong")$membership
    }
    numComponents <- max(this$.stronglyConnectedComponents)
    components <- vector("list", numComponents)
    for (i in 1:numComponents) {
        components[[i]] <- which(this$.stronglyConnectedComponents == i)
    }
    return(components)
}, appendVarArgs = F)


#' Strongly connected component
#'
#' Get the strongly connected component for a node i in the directed part of
#' the graph.
#'
#' @name stronglyConnectedComponent
#' @export stronglyConnectedComponent
#'
#' @param this the mixed graph object
#' @param node the node for which to get the strongly connected component.
stronglyConnectedComponent <- function(this, node) {
    UseMethod("stronglyConnectedComponent")
}

#' @rdname   stronglyConnectedComponent
#' @name     stronglyConnectedComponent.MixedGraphFixedOrder
#' @export
setMethodS3("stronglyConnectedComponent", "MixedGraphFixedOrder", function(this,
    node) {
    if (is.null(this$.stronglyConnectedComponents)) {
        this$.stronglyConnectedComponents <- igraph::components(this$.dirGraph, "strong")$membership
    }
    return(which(this$.stronglyConnectedComponents[node] == this$.stronglyConnectedComponents))
}, appendVarArgs = F)

#' Get descendants of a node
#'
#' Finds all descendants of a node, this DOES include the node itself (every
#' node is considered a descendant of itself).
#'
#' @name descendants
#' @export descendants
#'
#' @param this the mixed graph object
#' @param node the node from which to get the descendants.
descendants <- function(this, node) {
    UseMethod("descendants")
}

#' @rdname   descendants
#' @name     descendants.MixedGraphFixedOrder
#' @export
setMethodS3("descendants", "MixedGraphFixedOrder", function(this, node) {
    return(as.integer(igraph::neighborhood(this$.dirGraph, order = this$numNodes(),
        nodes = node, mode = "out")[[1]]))
}, appendVarArgs = F)


### The MixedGraph wrapper class

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
#' @usage MixedGraph(L = matrix(0,1,1), O = matrix(0,1,1),
#'                   vertexNums = 1:nrow(L))
#' @export MixedGraph
#'
#' @param L see \code{\link{graphID}} for the appropriate form of L.
#' @param O as for L.
#' @param vertexNums the labeling of the vertices in the graph in the order
#'        of the rows of L and O. Labels must be positive integers.
#'
#' @return An object representing the MixedGraph
setConstructorS3("MixedGraph", function(L = matrix(0, 1, 1), O = matrix(0,
    1, 1), vertexNums = 1:nrow(L)) {

    internalGraph <- MixedGraphFixedOrder(L, O)

    if (nrow(L) == 0) {
        vertexNums <- c()
        vertexNumsToInternal <- c()
    } else if (any(vertexNums%%1 != 0) || any(vertexNums < 1) || length(unique(vertexNums)) !=
        length(vertexNums) || length(vertexNums) != nrow(L)) {
        stop(paste("vertexNums must be all unique positive", "integers and must have length == nrow(L)"))
    } else {
        vertexNumsToInternal <- rep(NA, max(vertexNums))
        vertexNumsToInternal[vertexNums] <- 1:nrow(L)
    }

    R.oo::extend(R.oo::Object(), "MixedGraph", .L = L, .O = O, .internalGraph = internalGraph,
        .vertexNums = vertexNums, .vertexNumsToInternal = vertexNumsToInternal)
})

#' Transforms a vector of given node indices into their internal numbering
#'
#' @name toIn
#' @export toIn
#'
#' @param this the mixed graph object
#' @param nodes the nodes to transform
toIn <- function(this, nodes) {
    UseMethod("toIn")
}

#' @rdname   toIn
#' @name     toIn.MixedGraph
#' @export
setMethodS3("toIn", "MixedGraph", function(this, nodes) {
    return(this$.vertexNumsToInternal[nodes])
}, appendVarArgs = F, private = T)

#' Transforms a vector of node indices in the internal rep. into external numbering
#'
#' @name toEx
#' @export toEx
#'
#' @param this the mixed graph object
#' @param nodes the nodes to transform
toEx <- function(this, nodes) {
    UseMethod("toEx")
}

#' @rdname   toEx
#' @name     toEx.MixedGraph
#' @export
setMethodS3("toEx", "MixedGraph", function(this, nodes) {
    return(this$.vertexNums[nodes])
}, appendVarArgs = F, private = T)

#' Get adjacency matrix for directed part.
#'
#' @name L
#' @export L
#'
#' @param this the mixed graph object
L <- function(this) {
    UseMethod("L")
}

#' @rdname   L
#' @name     L.MixedGraph
#' @export
setMethodS3("L", "MixedGraph", function(this) {
    return(this$.L)
}, appendVarArgs = F)

#' Get adjacency matrix for bidirected part.
#'
#' @name O
#' @export O
#'
#' @param this the mixed graph object
O <- function(this) {
    UseMethod("O")
}

#' @rdname   O
#' @name     O.MixedGraph
#' @export
setMethodS3("O", "MixedGraph", function(this) {
    return(this$.O)
}, appendVarArgs = F)

#' Get all nodes in the graph.
#'
#' @name nodes
#' @export nodes
#'
#' @param this the mixed graph object
nodes <- function(this) {
    UseMethod("nodes")
}

#' @rdname   nodes
#' @name     nodes.MixedGraph
#' @export
setMethodS3("nodes", "MixedGraph", function(this) {
    return(this$.vertexNums)
}, appendVarArgs = F)

#' @rdname   numNodes
#' @name     numNodes.MixedGraph
#' @export
setMethodS3("numNodes", "MixedGraph", function(this) {
    return(this$.internalGraph$numNodes())
}, appendVarArgs = F)

#' @rdname   siblings
#' @name     siblings.MixedGraph
#' @export
setMethodS3("siblings", "MixedGraph", function(this, nodes) {
    return(this$toEx(this$.internalGraph$siblings(this$toIn(nodes))))
}, appendVarArgs = F)

#' @rdname   isSibling
#' @name     isSibling.MixedGraph
#' @export
setMethodS3("isSibling", "MixedGraph", function(this, node1, node2) {
    return(this$.internalGraph$isSibling(this$toIn(node1), this$toIn(node2)))
}, appendVarArgs = F)

#' @rdname   parents
#' @name     parents.MixedGraph
#' @export
setMethodS3("parents", "MixedGraph", function(this, nodes) {
    return(this$toEx(this$.internalGraph$parents(this$toIn(nodes))))
}, appendVarArgs = F)

#' @rdname   ancestors
#' @name     ancestors.MixedGraph
#' @export
setMethodS3("ancestors", "MixedGraph", function(this, nodes) {
    return(this$toEx(this$.internalGraph$ancestors(this$toIn(nodes))))
}, appendVarArgs = F)

#' @rdname   htrFrom
#' @name     htrFrom.MixedGraph
#' @export
setMethodS3("htrFrom", "MixedGraph", function(this, node) {
    return(this$toEx(this$.internalGraph$htrFrom(this$toIn(node))))
}, appendVarArgs = F)

#' @rdname   getHalfTrekSystem
#' @name     getHalfTrekSystem.MixedGraph
#' @export
setMethodS3("getHalfTrekSystem", "MixedGraph", function(this, fromNodes,
    toNodes) {
    l <- this$.internalGraph$getHalfTrekSystem(this$toIn(fromNodes), this$toIn(toNodes))

    return(list(systemExists = l$systemExists, activeFrom = this$toEx(l$activeFrom)))
}, appendVarArgs = F)

#' @rdname   trFrom
#' @name     trFrom.MixedGraph
#' @export
setMethodS3("trFrom", "MixedGraph", function(this, node) {
    return(this$toEx(this$.internalGraph$trFrom(this$toIn(node))))
}, appendVarArgs = F)

#' @rdname   getTrekSystem
#' @name     getTrekSystem.MixedGraph
#' @export
setMethodS3("getTrekSystem", "MixedGraph", function(this, fromNodes,
    toNodes, avoidEdgesOnRight = NULL) {
    l <- this$.internalGraph$getTrekSystem(this$toIn(fromNodes), this$toIn(toNodes),
        this$toIn(avoidEdgesOnRight))
    return(list(systemExists = l$systemExists, activeFrom = this$toEx(l$activeFrom)))
}, appendVarArgs = F)

#' @rdname   stronglyConnectedComponent
#' @name     stronglyConnectedComponent.MixedGraph
#' @export
setMethodS3("stronglyConnectedComponent", "MixedGraph", function(this,
    node) {
    return(this$toEx(this$.internalGraph$stronglyConnectedComponent(this$toIn(node))))
}, appendVarArgs = F)

#' @rdname   descendants
#' @name     descendants.MixedGraph
#' @export
setMethodS3("descendants", "MixedGraph", function(this, node) {
    return(this$toEx(this$.internalGraph$descendants(this$toIn(node))))
}, appendVarArgs = F)

#' Get the induced subgraph on a collection of nodes
#'
#' @name inducedSubgraph
#' @export inducedSubgraph
#'
#' @param this the mixed graph object
#' @param nodes the nodes on which to create the induced subgraph.
inducedSubgraph <- function(this, nodes) {
    UseMethod("inducedSubgraph")
}

#' @rdname   inducedSubgraph
#' @name     inducedSubgraph.MixedGraph
#' @export
setMethodS3("inducedSubgraph", "MixedGraph", function(this, nodes) {
    nodesIn <- this$toIn(nodes)
    newL <- this$.internalGraph$.L[nodesIn, nodesIn]
    newO <- this$.internalGraph$.O[nodesIn, nodesIn]

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
    components <- this$.internalGraph$stronglyConnectedComponents()
    numComponents <- length(components)

    shrunkO <- matrix(0, numComponents, numComponents)
    shrunkL <- matrix(0, numComponents, numComponents)

    if (numComponents > 1) {
        for (i in 1:(numComponents - 1)) {
            for (j in (i + 1):numComponents) {
                shrunkO[i, j] <- 1 * any(this$.internalGraph$.O[components[[i]],
                  components[[j]]] != 0)
                shrunkL[i, j] <- 1 * any(this$.internalGraph$.L[components[[i]],
                  components[[j]]] != 0)
            }
        }
        shrunkO <- shrunkO + t(shrunkO)
    }

    biGraph <- igraph::graph.adjacency(shrunkO, mode = "undirected")
    biComponents <- igraph::components(biGraph)$membership
    shrunkTopOrder <- as.integer(igraph::topological.sort(igraph::graph.adjacency(shrunkL,
        mode = "directed"), mode = "out"))
    topOrder <- unlist(components[shrunkTopOrder])

    cComponents <- rep(list(list()), max(biComponents))

    for (i in 1:max(biComponents)) {
        superNodes <- which(biComponents == i)
        internal <- topOrder[topOrder %in% unlist(components[superNodes])]
        incoming <- topOrder[topOrder %in% setdiff(this$.internalGraph$parents(internal),
            internal)]
        allOrdered <- topOrder[topOrder %in% c(internal, incoming)]

        indsInt <- which(allOrdered %in% internal)
        indsInc <- if (length(incoming) != 0) {
            which(allOrdered %in% incoming)
        } else {
            c()
        }
        newL <- matrix(0, length(internal) + length(incoming), length(internal) +
            length(incoming))
        newO <- newL
        newL[c(indsInt, indsInc), indsInt] <- this$.internalGraph$.L[c(internal,
            incoming), internal]
        newO[indsInt, indsInt] <- this$.internalGraph$.O[internal, internal]

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
    for (i in 1:length(cComponents)) {
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
    plotMixedGraph(x$L(), x$O(), vertexLabels = x$.vertexNums)
}, appendVarArgs = F)

