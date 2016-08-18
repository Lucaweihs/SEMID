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
  L = this$.L
  O = this$.O
  m = this$numNodes()
  htrGraphAdjMatrix = matrix(0, 2 * m, 2 * m)
  for (i in 1:m) {
    # Left i points to right i
    htrGraphAdjMatrix[i, i + m] = 1
    for (j in 1:m) {
      if (O[i,j] == 1) {
        # If bidirected edge from i to j then
        # left i should point to right j
        htrGraphAdjMatrix[i, j + m] = 1
      }

      if (L[i,j] == 1) {
        # If directed edge from i to j then
        # right i should point to right j
        htrGraphAdjMatrix[i + m, j + m] = 1
      }
    }
  }
  return(igraph::graph.adjacency(htrGraphAdjMatrix, mode = "directed"))
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
  m = this$numNodes()
  htrGraph = this$createHtrGraph()
  htrFrom = igraph::neighborhood(htrGraph, order = 2 * m, nodes = 1:m, mode = "out", mindist = 1)
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
  m = this$numNodes()
  adjMat = matrix(0, 3 * m + 2, 3 * m + 2)
  for (i in 1:m) {
    # Left-i points to right-i-in, right-i-in points to right-i-out
    adjMat[i, i + m] = 1 # Left-i -> right-i-in
    adjMat[i + m, i + 2 * m] = 1 # Right-i-in -> right-i-out

    for (j in 1:m) {
      if (i != j) {
        if (this$.O[i,j] == 1) {
          # If bidirected edge from i to j then
          # left-i should point to right-j-in
          adjMat[i, j + m] = 1
        }

        if (this$.L[i,j] == 1) {
          # If directed edge from i to j then
          # right-i-out should point to right-j-in
          adjMat[i + 2 * m, j + m] = 1
        }
      }
    }
  }
  adjMat[3 * m + 1, 1:m] = 1 # Source points to all lefts
  adjMat[(2 * m + 1):(3 * m), 3 * m + 2] = 1 # All right-outs point to target
  g = igraph::graph.adjacency(adjMat, mode = "directed")
  igraph::E(g)$capacity = 1 # All edges (and thus vertices) have capacity 1
  sOutEdgeIds = get.edge.ids(g, as.numeric(rbind(rep(3 * m + 1, m), 1:m)))
  tInEdgeIds = get.edge.ids(g, as.numeric(rbind((2 * m + 1):(3 * m), rep(3 * m + 2, m))))
  return(list(flowGraph = g, s = 3 * m + 1, t = 3 * m + 2,
              sOutEdgeIds = sOutEdgeIds, tInEdgeIds = tInEdgeIds))
}, appendVarArgs = F)

#' Determines if a trek system exists in the mixed graph.
#'
#' @name     getHalfTrekSystem
#' @export   getHalfTrekSystem
#' @export   getHalfTrekSystem.MixedGraph
NULL
setMethodS3("getHalfTrekSystem", "MixedGraph", function(this, fromNodes,
                                                              toNodes) {
  if (is.null(this$.halfTrekFlowGraph)) {
    flowGraphList = this$createHalfTrekFlowGraph()
    this$.halfTrekFlowGraph = flowGraphList$flowGraph
    this$.s = flowGraphList$s
    this$.t = flowGraphList$t
    this$.sOutEdgeIds = flowGraphList$sOutEdgeIds
    this$.tInEdgeIds = flowGraphList$tInEdgeIds
  }

  igraph::E(this$.halfTrekFlowGraph)$capacity[this$.sOutEdgeIds] = 0
  igraph::E(this$.halfTrekFlowGraph)$capacity[this$.sOutEdgeIds[fromNodes]] = 1
  igraph::E(this$.halfTrekFlowGraph)$capacity[this$.tInEdgeIds] = 0
  igraph::E(this$.halfTrekFlowGraph)$capacity[this$.tInEdgeIds[toNodes]] = 1

  flowResult = igraph::graph.maxflow(this$.halfTrekFlowGraph, this$.s, this$.t)
  return(list(systemExists = (flowResult$value == length(toNodes)),
              activeFrom = which(flowResult$flow[this$.sOutEdgeIds] == 1)))
}, appendVarArgs = F)
