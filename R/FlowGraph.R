#' Construct FlowGraph object
#'
#' Creates an object representing a flow graph.
#'
#' @name FlowGraph
#' @usage FlowGraph(L = matrix(0,1,1), vertexCaps = 1, edgeCaps = matrix(1,1,1))
#' @export FlowGraph
#'
#' @inheritParams graphID
#'
#' @return An object representing the FlowGraph
NULL
setConstructorS3("FlowGraph",
                 function(L = matrix(0,1,1), vertexCaps = 1,
                          edgeCaps = matrix(1,1,1)) {
                   validateMatrices(L, L)
                   m = nrow(L)
                   adjMat = matrix(0, 2 * m + 2, 2 * m + 2)

                   # Each in-part of a vertex should point to its corresponding
                   # out part, i.e. i should point to i + m
                   adjMat[cbind(1:m, 1:m + m)] = 1

                   # All of the out parts of the vertices should point to the
                   # in parts of the vertices they used to point to
                   adjMat[1:m + m, 1:m] = L

                   # Source vertex points to all in-parts
                   s = 2 * m + 1
                   adjMat[s, 1:m] = 1

                   # Terminal vertex has all out parts as parents
                   t = 2 * m + 2
                   adjMat[m + 1:m, t] = 1

                   # Create flow graph and set all edge capacities to 0 initially
                   flowGraph = igraph::graph.adjacency(adjMat)
                   igraph::E(flowGraph)$capacity = 0

                   # Set vertex capacities
                   vertexCapEdges = igraph::get.edge.ids(flowGraph,
                                                         rbind(1:m, 1:m + m))
                   igraph::E(flowGraph)$capacity[vertexCapEdges] = vertexCaps

                   # Set edge capacities
                   edgeArrInds = which(L != 0, arr.ind = T)
                   numEdges = nrow(edgeArrInds)
                   edgeInds = igraph::get.edge.ids(flowGraph,
                                                   t(edgeArrInds) +
                                                     rbind(rep(m, numEdges),
                                                           rep(0, numEdges))
                                                   )
                   igraph::E(flowGraph)$capacity[edgeInds] = edgeCaps[edgeArrInds]

                   # Edge indices for source to in-nodes and forom out-nodes
                   # to terminal node.
                   sOutIndices = igraph::get.edge.ids(flowGraph, rbind(rep(s, m), 1:m))
                   tInIndices = igraph::get.edge.ids(flowGraph, rbind(m + 1:m, rep(t, m)))

                   extend(Object(),
                          "FlowGraph",
                          .adjMat = adjMat,
                          .m = m,
                          .s = s,
                          .t = t,
                          .sOutIndices = sOutIndices,
                          .tInIndices = tInIndices,
                          .flowGraph = flowGraph,
                          .vertexCapEdges = vertexCapEdges)
                 })


#' Flow from one set of nodes to another.
#'
#' @name     flowBetween
#' @export   flowBetween
#' @export   flowBetween.FlowGraph
#'
NULL
setMethodS3("flowBetween", "FlowGraph", function(this, sources, sinks) {
  totalCap = sum(igraph::E(this$.flowGraph)$capacity)
  igraph::E(this$.flowGraph)$capacity[this$.sOutIndices[sources]] = totalCap
  igraph::E(this$.flowGraph)$capacity[this$.tInIndices[sinks]] = totalCap
  flowResult = igraph::max_flow(this$.flowGraph, this$.s, this$.t)
  igraph::E(this$.flowGraph)$capacity[this$.sOutIndices[sources]] = 0
  igraph::E(this$.flowGraph)$capacity[this$.tInIndices[sinks]] = 0
  return(list(value = flowResult$value,
              activeSources = which(flowResult$flow[this$.sOutIndices] != 0)))
}, appendVarArgs = F)

#' Update vertex capacities
#'
#' @name     updateVertexCapacities
#' @export   updateVertexCapacities
#' @export   updateVertexCapacities.FlowGraph
#'
NULL
setMethodS3("updateVertexCapacities", "FlowGraph", function(this, vertices, newCaps) {
  igraph::E(this$.flowGraph)$capacity[this$.vertexCapEdges[vertices]] = newCaps
}, appendVarArgs = F)

#' Update edge capacities
#'
#' @name     updateEdgeCapacities
#' @export   updateEdgeCapacities
#' @export   updateEdgeCapacities.FlowGraph
NULL
setMethodS3("updateEdgeCapacities", "FlowGraph", function(this, edges, newCaps) {
  if (is.vector(edges)) {
    edges = matrix(edges, nrow = 2)
  }
  edgeIds = igraph::get.edge.ids(this$.flowGraph,
                                 edges + rbind(rep(this$.m, ncol(edges)), 0))
  igraph::E(this$.flowGraph)$capacity[edgeIds] = newCaps
}, appendVarArgs = F)