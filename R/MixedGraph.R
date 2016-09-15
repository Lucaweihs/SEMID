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
                 function(L = matrix(0,1,1), O = matrix(0,1,1)) {
                   validateMatrices(L, O)
                   O = 1 * ((O + t(O)) != 0)

                   dirGraph = igraph::graph.adjacency(L, mode = "directed")

                   extend(Object(),
                          "MixedGraph",
                          .L = L,
                          .O = O,
                          .dirGraph = dirGraph)
                 })


#' Number of nodes in the graph.
#'
#' @name     numNodes
#' @export   numNodes
#' @export   numNodes.MixedGraph
#'
NULL
setMethodS3("numNodes", "MixedGraph", function(this) {
  return(nrow(this$.L))
}, appendVarArgs = F)

#' All siblings of a node.
#'
#' @name     allSiblings
#' @export   allSiblings
#' @export   allSiblings.MixedGraph
NULL
setMethodS3("allSiblings", "MixedGraph", function(this, node) {
  return(which(this$.O[node,] != 0))
}, appendVarArgs = F)

#' Are two nodes siblings?
#'
#' @name     isSibling
#' @export   isSibling
#' @export   isSibling.MixedGraph
NULL
setMethodS3("isSibling", "MixedGraph", function(this, node1, node2) {
  return(this$.O[node1, node2] != 0)
}, appendVarArgs = F)

#' All parents of a node.
#'
#' @name     allParents
#' @export   allParents
#' @export   allParents.MixedGraph
NULL
setMethodS3("allParents", "MixedGraph", function(this, node) {
  return(which(this$.L[,node] != 0))
}, appendVarArgs = F)

#' Helper function to create a graph encoding htr relationships.
#'
#' @name     createHtrGraph
#' @export   createHtrGraph
#' @export   createHtrGraph.MixedGraph
NULL
setMethodS3("createHtrGraph", "MixedGraph", function(this) {
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
#' @export   htrFrom.MixedGraph
NULL
setMethodS3("htrFrom", "MixedGraph", function(this, node) {
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
#' @export   createHalfTrekFlowGraph.MixedGraph
NULL
setMethodS3("createHalfTrekFlowGraph", "MixedGraph", function(this) {
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
#' @export   getHalfTrekSystem.MixedGraph
NULL
setMethodS3("getHalfTrekSystem", "MixedGraph", function(this, fromNodes,
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
#' @export   createTrGraph.MixedGraph
NULL
setMethodS3("createTrGraph", "MixedGraph", function(this) {
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
#' @export   trFrom.MixedGraph
NULL
setMethodS3("trFrom", "MixedGraph", function(this, node) {
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
#' @export   createTrekFlowGraph.MixedGraph
NULL
setMethodS3("createTrekFlowGraph", "MixedGraph", function(this) {
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
#' @export   getTrekSystem.MixedGraph
NULL
setMethodS3("getTrekSystem", "MixedGraph", function(this, fromNodes, toNodes,
                                                    avoidEdgeOnRight = NULL) {
  if (is.null(this$.trekFlowGraph)) {
    this$.trekFlowGraph = this$createTrekFlowGraph()
  }

  if (!is.null(avoidEdgeOnRight)) {
    if (this$.L[avoidEdgeOnRight[1], avoidEdgeOnRight[2]]  == 0) {
      stop("avoidEdgeOnRight is not an edge in the graph")
    }
    this$.trekFlowGraph$updateEdgeCapacities(avoidEdgeOnRight + this$numNodes(),
                                             0)
  }
  flowResult = this$.trekFlowGraph$flowBetween(fromNodes,
                                               this$numNodes() + toNodes)
  if (!is.null(avoidEdgeOnRight)) {
    this$.trekFlowGraph$updateEdgeCapacities(avoidEdgeOnRight + this$numNodes(),
                                             1)
  }
  return(list(systemExists = (flowResult$value == length(toNodes)),
              activeFrom = flowResult$activeSources))
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
#' @export   allDescendants.MixedGraph
NULL
setMethodS3("allDescendants", "MixedGraph", function(this, node) {
  return(as.numeric(
    igraph::neighborhood(this$.dirGraph, order = this$numNodes(),
                         nodes = node, mode = "out", mindist = 1)[[1]]))
}, appendVarArgs = F)
