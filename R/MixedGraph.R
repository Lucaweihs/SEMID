###
# First we define a class, MixedGraphFixedOrder, representing mixed graphs whose
# vertex labeling is fixed. We then write a wrapper class, MixedGraph, which
# allows for different vertex labelings.
###

#' Construct MixedGraphFixedOrder object
#'
#' Creates an object representing a mixed graph.
#'
#' @name MixedGraphFixedOrder
#' @usage MixedGraphFixedOrder(L = matrix(0,1,1), O = matrix(0,1,1))
#' @export MixedGraphFixedOrder
#'
#' @inheritParams graphID
#'
#' @return An object representing the MixedGraphFixedOrder
NULL
setConstructorS3("MixedGraphFixedOrder",
                 function(L = matrix(0,1,1), O = matrix(0,1,1)) {
                   validateMatrices(L, O)
                   O = 1 * ((O + t(O)) != 0)

                   dirGraph = igraph::graph.adjacency(L, mode = "directed")

                   extend(Object(),
                          "MixedGraphFixedOrder",
                          .L = L,
                          .O = O,
                          .dirGraph = dirGraph)
                 })


#' Number of nodes in the graph.
#'
#' @name     numNodes
#' @export   numNodes
#' @export   numNodes.MixedGraphFixedOrder
#'
NULL
setMethodS3("numNodes", "MixedGraphFixedOrder", function(this) {
  return(nrow(this$.L))
}, appendVarArgs = F)

#' All siblings of a collection of nodes
#'
#' @name     allSiblings
#' @export   allSiblings
#' @export   allSiblings.MixedGraphFixedOrder
NULL
setMethodS3("allSiblings", "MixedGraphFixedOrder", function(this, nodes) {
  return(which(rowSums(this$.O[,nodes,drop = F]) != 0))
}, appendVarArgs = F)

#' Are two nodes siblings?
#'
#' @name     isSibling
#' @export   isSibling
#' @export   isSibling.MixedGraphFixedOrder
NULL
setMethodS3("isSibling", "MixedGraphFixedOrder", function(this, node1, node2) {
  return(this$.O[node1, node2] != 0)
}, appendVarArgs = F)

#' All parents a collection of nodes.
#'
#' @name     allParents
#' @export   allParents
#' @export   allParents.MixedGraphFixedOrder
NULL
setMethodS3("allParents", "MixedGraphFixedOrder", function(this, nodes) {
  return(which(rowSums(this$.L[,nodes,drop=F]) != 0))
}, appendVarArgs = F)

#' Helper function to create a graph encoding htr relationships.
#'
#' @name     createHtrGraph
#' @export   createHtrGraph
#' @export   createHtrGraph.MixedGraphFixedOrder
NULL
setMethodS3("createHtrGraph", "MixedGraphFixedOrder", function(this) {
  m = this$numNodes()
  adjMat = matrix(0, 2 * m, 2 * m)

  # Left nodes point to right nodes
  adjMat[cbind(1:m, m + 1:m)] = 1

  # If i <-> j, then left i should point to right j and left j to right i
  adjMat[1:m, m + 1:m] = this$.O + adjMat[1:m, m + 1:m]

  # If i -> j then right i should point to right j
  adjMat[m + 1:m, m + 1:m] = this$.L

  return(igraph::graph.adjacency(adjMat, mode = "directed"))
}, appendVarArgs = F)

#' Half trek reachable nodes.
#'
#' @name     htrFrom
#' @export   htrFrom
#' @export   htrFrom.MixedGraphFixedOrder
NULL
setMethodS3("htrFrom", "MixedGraphFixedOrder", function(this, node) {
  if (!is.null(this$.htrFrom)) {
    return(this$.htrFrom[[node]])
  }
  if (is.null(this$.htrGraph)) {
    this$.htrGraph = this$createHtrGraph()
  }
  m = this$numNodes()
  htrFrom = igraph::neighborhood(this$.htrGraph, order = 2 * m, nodes = 1:m, mode = "out", mindist = 1)
  for (i in 1:m) {
    htrFrom[[i]] = as.numeric(htrFrom[[i]]) - m
  }
  this$.htrFrom = htrFrom
  return(this$.htrFrom[[node]])
}, appendVarArgs = F)

#' Helper function to create a flow graph.
#'
#' @name     createHalfTrekFlowGraph
#' @export   createHalfTrekFlowGraph
#' @export   createHalfTrekFlowGraph.MixedGraphFixedOrder
NULL
setMethodS3("createHalfTrekFlowGraph", "MixedGraphFixedOrder", function(this) {
  if (is.null(this$.htrGraph)) {
    this$.htrGraph = this$createHtrGraph()
  }
  adjMat = as.matrix(igraph::get.adjacency(this$.htrGraph))
  # Create the flow graph from adjMat with all vertices and edges having a
  # capacity of 1
  flowGraph = FlowGraph(adjMat, rep(1, 2 * this$numNodes()), adjMat)
  return(flowGraph)
}, appendVarArgs = F)

#' Determines if a half trek system exists in the mixed graph.
#'
#' @name     getHalfTrekSystem
#' @export   getHalfTrekSystem
#' @export   getHalfTrekSystem.MixedGraphFixedOrder
NULL
setMethodS3("getHalfTrekSystem", "MixedGraphFixedOrder", function(this, fromNodes,
                                                              toNodes) {
  if (is.null(this$.halfTrekFlowGraph)) {
    this$.halfTrekFlowGraph = this$createHalfTrekFlowGraph()
  }
  flowResult = this$.halfTrekFlowGraph$flowBetween(fromNodes,
                                                   this$numNodes() + toNodes)
  return(list(systemExists = (flowResult$value == length(toNodes)),
              activeFrom = flowResult$activeSources))
}, appendVarArgs = F)


#' Helper function to create a graph encoding trek reachable relationships.
#'
#' @name     createTrGraph
#' @export   createTrGraph
#' @export   createTrGraph.MixedGraphFixedOrder
NULL
setMethodS3("createTrGraph", "MixedGraphFixedOrder", function(this) {
  m = this$numNodes()
  adjMat = matrix(0, 2 * m, 2 * m)

  # Left nodes point to each other in the opposite direction of L
  adjMat[1:m, 1:m] = t(this$.L)

  # Left nodes point to their corresponding right nodes
  adjMat[cbind(1:m, m + 1:m)] = 1

  # If i <-> j, then left i should point to right j and left j to right i
  adjMat[1:m, m + 1:m] = 1 * ((this$.O + adjMat[1:m, m + 1:m]) != 0)

  # If i -> j then right i point to right j
  adjMat[m + 1:m, m + 1:m] = this$.L

  return(igraph::graph.adjacency(adjMat, mode = "directed"))
}, appendVarArgs = F)

#' Trek reachable nodes.
#'
#' @name     trFrom
#' @export   trFrom
#' @export   trFrom.MixedGraphFixedOrder
NULL
setMethodS3("trFrom", "MixedGraphFixedOrder", function(this, node) {
  if (is.null(this$.trGraph)) {
    this$.trGraph = this$createTrGraph()
  }
  m = this$numNodes()
  trFrom = as.numeric(
    igraph::neighborhood(this$.trGraph, order = 2 * m, nodes = node,
                          mode = "out")[[1]])
  trFrom[trFrom > m] = trFrom[trFrom > m] - m
  return(unique(trFrom))
}, appendVarArgs = F)

#' Helper function to create a flow graph.
#'
#' @name     createTrekFlowGraph
#' @export   createTrekFlowGraph
#' @export   createTrekFlowGraph.MixedGraphFixedOrder
NULL
setMethodS3("createTrekFlowGraph", "MixedGraphFixedOrder", function(this) {
  if (is.null(this$.trGraph)) {
    this$.trGraph = this$createTrGraph()
  }
  adjMat = as.matrix(igraph::get.adjacency(this$.trGraph))

  # Create the flow graph from adjMat. All vertices and edges have capacity 1
  flowGraph = FlowGraph(adjMat, rep(1, 2 * this$numNodes()), adjMat)

  return(flowGraph)
}, appendVarArgs = F)

#' Determines if a trek system exists in the mixed graph.
#'
#' @name     getTrekSystem
#' @export   getTrekSystem
#' @export   getTrekSystem.MixedGraphFixedOrder
NULL
setMethodS3("getTrekSystem", "MixedGraphFixedOrder", function(this, fromNodes, toNodes,
                                                    avoidEdgesOnRight = NULL) {
  if (is.null(this$.trekFlowGraph)) {
    this$.trekFlowGraph = this$createTrekFlowGraph()
  }

  if (length(avoidEdgesOnRight) !=0 ) {
    if (is.vector(avoidEdgesOnRight)) {
      avoidEdgesOnRight = matrix(avoidEdgesOnRight, byrow = T, ncol = 2)
    }
    if (any(this$.L[avoidEdgesOnRight]  == 0)) {
      stop("Some edge in avoidEdgesOnRight is not an edge in the graph")
    }
    this$.trekFlowGraph$updateEdgeCapacities(t(avoidEdgesOnRight + this$numNodes()),
                                             0)
  }
  flowResult = this$.trekFlowGraph$flowBetween(fromNodes,
                                               this$numNodes() + toNodes)
  if (length(avoidEdgesOnRight) != 0) {
    this$.trekFlowGraph$updateEdgeCapacities(t(avoidEdgesOnRight + this$numNodes()),
                                             1)
  }
  return(list(systemExists = (flowResult$value == length(toNodes)),
              activeFrom = flowResult$activeSources))
}, appendVarArgs = F)

#' Strongly connected components
#'
#' Get the strongly connected components of a graph
#'
#' @name     stronglyConnectedComponents
#' @export   stronglyConnectedComponents
#' @export   stronglyConnectedComponents.MixedGraphFixedOrder
NULL
setMethodS3("stronglyConnectedComponents", "MixedGraphFixedOrder", function(this) {
  if (is.null(this$.stronglyConnectedComponents)) {
    this$.stronglyConnectedComponents = igraph::components(this$.dirGraph,
                                                           "strong")$membership
  }
  numComponents = max(this$.stronglyConnectedComponents)
  components = vector("list", numComponents)
  for (i in 1:numComponents) {
    components[[i]] = which(this$.stronglyConnectedComponents == i)
  }
  return(components)
}, appendVarArgs = F)


#' Strongly connected component
#'
#' Get the strongly connected component for a node i in the directed part of
#' the graph.
#'
#' @name     stronglyConnectedComponent
#' @export   stronglyConnectedComponent
#' @export   stronglyConnectedComponent.MixedGraphFixedOrder
NULL
setMethodS3("stronglyConnectedComponent", "MixedGraphFixedOrder", function(this, node) {
  if (is.null(this$.stronglyConnectedComponents)) {
    this$.stronglyConnectedComponents = igraph::components(this$.dirGraph,
                                                           "strong")$membership
  }
  return(which(this$.stronglyConnectedComponents[node] == this$.stronglyConnectedComponents))
}, appendVarArgs = F)

#' Get descendents of a node
#'
#' @name     allDescendants
#' @export   allDescendants
#' @export   allDescendants.MixedGraphFixedOrder
NULL
setMethodS3("allDescendants", "MixedGraphFixedOrder", function(this, node) {
  return(as.numeric(
    igraph::neighborhood(this$.dirGraph, order = this$numNodes(),
                         nodes = node, mode = "out", mindist = 1)[[1]]))
}, appendVarArgs = F)



###
# The MixedGraph wrapper class
###

#' Construct MixedGraph object
#'
#' Creates an object representing a mixed graph.
#'
#' @name MixedGraph
#' @usage MixedGraph(L = matrix(0,1,1), O = matrix(0,1,1))
#' @export MixedGraph
#'
#' @inheritParams graphID
#'
#' @return An object representing the MixedGraph
NULL
setConstructorS3("MixedGraph",
                 function(L = matrix(0,1,1), O = matrix(0,1,1), vertexNums = 1:nrow(L)) {

                   internalGraph = MixedGraphFixedOrder(L, O)

                   if (nrow(L) == 0) {
                     vertexNums = c()
                     vertexNumsToInternal = c()
                   } else if (vertexNums %% 1 != 0 || any(vertexNums < 1) ||
                              length(unique(vertexNums)) != length(vertexNums) ||
                              length(vertexNums) != nrow(L)) {
                     stop(paste("vertexNums must be all unique positive",
                                "integers and must have length == nrow(L)"))
                   } else {
                     vertexNumsToInternal = rep(NA, max(vertexNums))
                     vertexNumsToInternal[vertexNums] = 1:nrow(L)
                   }

                   extend(Object(),
                          "MixedGraph",
                          .L = L,
                          .O = O,
                          .internalGraph = internalGraph,
                          .vertexNums = vertexNums,
                          .vertexNumsToInternal = vertexNumsToInternal)
                 })

#' Transforms a vector of given node indices into their internal numbering
#'
#' @name     toIn
#' @export   toIn
#' @export   toIn.MixedGraph
#'
NULL
setMethodS3("toIn", "MixedGraph", function(this, nodes) {
  return(this$.vertexNumsToInternal[nodes])
}, appendVarArgs = F)

#' Transforms a vector of node indices in the internal rep. into external numbering
#'
#' @name     toEx
#' @export   toEx
#' @export   toEx.MixedGraph
#'
NULL
setMethodS3("toEx", "MixedGraph", function(this, nodes) {
  return(this$.vertexNums[nodes])
}, appendVarArgs = F)

#' Number of nodes in the graph.
#'
#' @name     numNodes
#' @export   numNodes
#' @export   numNodes.MixedGraph
#'
NULL
setMethodS3("numNodes", "MixedGraph", function(this) {
  return(this$.internalGraph$numNodes())
}, appendVarArgs = F)

#' All siblings of a collection of nodes.
#'
#' @name     allSiblings
#' @export   allSiblings
#' @export   allSiblings.MixedGraph
NULL
setMethodS3("allSiblings", "MixedGraph", function(this, nodes) {
  return(this$toEx(this$.internalGraph$allSiblings(this$toIn(nodes))))
}, appendVarArgs = F)

#' Are two nodes siblings?
#'
#' @name     isSibling
#' @export   isSibling
#' @export   isSibling.MixedGraph
NULL
setMethodS3("isSibling", "MixedGraph", function(this, node1, node2) {
  return(this$.internalGraph$isSibling(this$toIn(node1), this$toIn(node2)))
}, appendVarArgs = F)

#' All parents of a collection of nodes.
#'
#' @name     allParents
#' @export   allParents
#' @export   allParents.MixedGraph
NULL
setMethodS3("allParents", "MixedGraph", function(this, nodes) {
  return(this$toEx(this$.internalGraph$allParents(this$toIn(nodes))))
}, appendVarArgs = F)

#' Half trek reachable nodes.
#'
#' @name     htrFrom
#' @export   htrFrom
#' @export   htrFrom.MixedGraph
NULL
setMethodS3("htrFrom", "MixedGraph", function(this, node) {
  return(this$toEx(this$.internalGraph$htrFrom(this$toIn(node))))
}, appendVarArgs = F)

#' Determines if a half trek system exists in the mixed graph.
#'
#' @name     getHalfTrekSystem
#' @export   getHalfTrekSystem
#' @export   getHalfTrekSystem.MixedGraph
NULL
setMethodS3("getHalfTrekSystem", "MixedGraph", function(this, fromNodes,
                                                        toNodes) {
  l = this$.internalGraph$getHalfTrekSystem(this$toIn(fromNodes), this$toEx(toNodes))

  return(list(systemExists = l$systemExists,
              activeFrom = this$toEx(l$activeFrom)))
}, appendVarArgs = F)

#' Trek reachable nodes.
#'
#' @name     trFrom
#' @export   trFrom
#' @export   trFrom.MixedGraph
NULL
setMethodS3("trFrom", "MixedGraph", function(this, node) {
  return(this$toEx(this$.internalGraph$trFrom(this$toIn(node))))
}, appendVarArgs = F)

#' Determines if a trek system exists in the mixed graph.
#'
#' @name     getTrekSystem
#' @export   getTrekSystem
#' @export   getTrekSystem.MixedGraph
NULL
setMethodS3("getTrekSystem", "MixedGraph", function(this, fromNodes, toNodes,
                                                    avoidEdgesOnRight = NULL) {
  l = this$.internalGraph$getTrekSystem(this$toIn(fromNodes),
                                        this$toIn(toNodes),
                                        this$toIn(avoidEdgesOnRight))
  return(list(systemExists = l$systemExists,
              activeFrom = this$toEx(l$activeFrom)))
}, appendVarArgs = F)

#' Strongly connected component
#'
#' Get the strongly connected component for a node i in the directed part of
#' the graph.
#'
#' @name     stronglyConnectedComponent
#' @export   stronglyConnectedComponent
#' @export   stronglyConnectedComponent.MixedGraph
NULL
setMethodS3("stronglyConnectedComponent", "MixedGraph", function(this, node) {
  return(this$toEx(this$.internalGraph$stronglyConnectedComponent(this$toIn(node))))
}, appendVarArgs = F)

#' Get descendents of a node
#'
#' @name     allDescendants
#' @export   allDescendants
#' @export   allDescendants.MixedGraph
NULL
setMethodS3("allDescendants", "MixedGraph", function(this, node) {
  return(this$toEx(this$.internalGraph$allDescendants(this$toIn(node))))
}, appendVarArgs = F)

#' Get the induced subgraph on a collection of nodes
#'
#' @name     inducedSubgraph
#' @export   inducedSubgraph
#' @export   inducedSubgraph.MixedGraph
NULL
setMethodS3("inducedSubgraph", "MixedGraph", function(this, nodes) {
  nodesIn = this$toIn(nodes)
  newL = this$.internalGraph$.L[nodesIn, nodesIn]
  newO = this$.internalGraph$.O[nodesIn, nodesIn]

  return(MixedGraph(newL, newO, vertexNums = nodes))
}, appendVarArgs = F)

#' Performs the tian decomposition
#'
#' @name     tianDecompose
#' @export   tianDecompose
#' @export   tianDecompose.MixedGraph
NULL
setMethodS3("tianDecompose", "MixedGraph", function(this) {

  components = this$.internalGraph$stronglyConnectedComponents()
  numComponents = length(components)

  shrunkO = matrix(0, numComponents, numComponents)
  shrunkL = matrix(0, numComponents, numComponents)

  if (numComponents > 1) {
    for (i in 1:(numComponents - 1)) {
      for (j in (i+1):numComponents) {
        shrunkO[i,j] = 1 * any(this$.internalGraph$.O[components[[i]], components[[j]]] != 0)
        shrunkL[i,j] = 1 * any(this$.internalGraph$.L[components[[i]], components[[j]]] != 0)
      }
    }
    shrunkO = shrunkO + t(shrunkO)
  }

  biGraph = igraph::graph.adjacency(shrunkO, mode = "undirected")
  biComponents = igraph::components(biGraph)$membership
  shrunkTopOrder = as.numeric(igraph::topological.sort(igraph::graph.adjacency(shrunkL, mode = "directed"), mode = "out"))
  topOrder = unlist(components[shrunkTopOrder])

  cComponents = rep(list(list()), max(biComponents))

  for (i in 1:max(biComponents)) {
    superNodes = which(biComponents == i)
    internal = topOrder[topOrder %in% unlist(components[superNodes])]
    incoming = topOrder[topOrder %in% setdiff(this$.internalGraph$allParents(internal), internal)]
    allOrdered = topOrder[topOrder %in% c(internal, incoming)]

    indsInt = which(allOrdered %in% internal)
    indsInc = if (length(incoming) != 0) { which(allOrdered %in% incoming) } else { c() }
    newL = matrix(0, length(internal) + length(incoming), length(internal) + length(incoming))
    newO = newL
    newL[c(indsInt, indsInc), indsInt] =
      this$.internalGraph$.L[c(internal, incoming), internal]
    newO[indsInt, indsInt] = this$.internalGraph$.O[internal, internal]

    cComponents[[i]]$internal = this$toEx(internal)
    cComponents[[i]]$incoming = this$toEx(incoming)
    cComponents[[i]]$topOrder = allOrdered
    cComponents[[i]]$L = newL
    cComponents[[i]]$O = newO
  }

  return(cComponents)
}, appendVarArgs = F)

#' Plot the mixed graph
#'
#' @name     plot.MixedGraph
#' @export   plot.MixedGraph
NULL
setMethodS3("plot", "MixedGraph", function(this) {
  plotMixedGraph(this$.L, this$.O, vertexLabels = this$.vertexNums)
}, appendVarArgs = F)

