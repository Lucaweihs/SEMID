#' Construct FlowGraph object
#'
#' Creates an object representing a flow graph.
#'
#' @name FlowGraph
#' @usage FlowGraph(L = matrix(0,1,1), vertexCaps = 1, edgeCaps = matrix(1,1,1))
#' @export
#'
#' @param L the adjacency matrix for the flow graph. The (i,j)th of L should be
#'        a 1 if there is an edge from i to j and 0 otherwise.
#' @param vertexCaps the capacity of the vertices in the flow graph, should
#'        either be a single number or a vector whose ith entry is the capacity
#'        of vertex i.
#' @param edgeCaps the capacities of the edges in the the flow graph, should
#'        be a matrix of the same dimensions as L with (i,j)th entry the
#'        capacity of the i->j edge.
#'
#' @return An object representing the FlowGraph
setConstructorS3("FlowGraph", function(L = matrix(0, 1, 1), vertexCaps = 1,
    edgeCaps = matrix(1, 1, 1)) {
    validateMatrices(L, L)
    m <- nrow(L)
    adjMat <- matrix(0, 2 * m + 2, 2 * m + 2)

    # Each in-part of a vertex should point to its corresponding out part, i.e.  i
    # should point to i + m
    adjMat[cbind(1:m, 1:m + m)] <- 1

    # All of the out parts of the vertices should point to the in parts of the
    # vertices they used to point to
    adjMat[1:m + m, 1:m] <- L

    # Source vertex points to all in-parts
    s <- 2 * m + 1
    adjMat[s, 1:m] <- 1

    # Terminal vertex has all out parts as parents
    t <- 2 * m + 2
    adjMat[m + 1:m, t] <- 1

    # Create flow graph and set all edge capacities to 0 initially
    flowGraph <- igraph::graph.adjacency(adjMat)
    igraph::E(flowGraph)$capacity <- 0

    # Set vertex capacities
    vertexCapEdges <- igraph::get.edge.ids(flowGraph, rbind(1:m, 1:m + m))
    igraph::E(flowGraph)$capacity[vertexCapEdges] <- vertexCaps

    # Set edge capacities
    edgeArrInds <- which(L != 0, arr.ind = T)
    numEdges <- nrow(edgeArrInds)
    edgeInds <- igraph::get.edge.ids(flowGraph, t(edgeArrInds) + rbind(rep(m, numEdges),
        rep(0, numEdges)))
    igraph::E(flowGraph)$capacity[edgeInds] <- edgeCaps[edgeArrInds]

    # Edge indices for source to in-nodes and forom out-nodes to terminal node.
    sOutIndices <- igraph::get.edge.ids(flowGraph, rbind(rep(s, m), 1:m))
    tInIndices <- igraph::get.edge.ids(flowGraph, rbind(m + 1:m, rep(t, m)))

    R.oo::extend(R.oo::Object(), "FlowGraph", .adjMat = adjMat, .m = m, .s = s, .t = t,
        .sOutIndices = sOutIndices, .tInIndices = tInIndices, .flowGraph = flowGraph,
        .vertexCapEdges = vertexCapEdges)
})


#' Flow from one set of nodes to another.
#'
#' @name flowBetween
#' @export flowBetween
#'
#' @param this the flow graph object
#' @param sources the nodes from which flow should start.
#' @param sinks the nodes at which the flow should end.
#'
#' @return a list with two named components, \code{value} (the size of the
#'         computed flow) and \code{activeSources} (a vector representing the
#'         subset of sources which have non-zero flow out of them for the found
#'         max-flow).
flowBetween <- function(this, sources, sinks) {
    UseMethod("flowBetween")
}

#' @rdname   flowBetween
#' @name     flowBetween.FlowGraph
#' @export
setMethodS3("flowBetween", "FlowGraph", function(this, sources, sinks) {
    igraph::E(this$.flowGraph)$capacity[this$.sOutIndices[sources]] <- 1
    igraph::E(this$.flowGraph)$capacity[this$.tInIndices[sinks]] <- 1
    flowResult <- igraph::max_flow(this$.flowGraph, this$.s, this$.t)
    igraph::E(this$.flowGraph)$capacity[this$.sOutIndices[sources]] <- 0
    igraph::E(this$.flowGraph)$capacity[this$.tInIndices[sinks]] <- 0
    return(list(value = flowResult$value, activeSources = which(flowResult$flow[this$.sOutIndices] !=
        0)))
}, appendVarArgs = F)

#' Update vertex capacities.
#'
#' @name updateVertexCapacities
#' @export updateVertexCapacities
#'
#' @param this the flow graph object
#' @param vertices the vertices to update.
#' @param newCaps the new capacities for the vertices.
updateVertexCapacities <- function(this, vertices, newCaps) {
    UseMethod("updateVertexCapacities")
}

#' @rdname   updateVertexCapacities
#' @name     updateVertexCapacities.FlowGraph
#' @export
setMethodS3("updateVertexCapacities", "FlowGraph", function(this, vertices,
    newCaps) {
    igraph::E(this$.flowGraph)$capacity[this$.vertexCapEdges[vertices]] <- newCaps
}, appendVarArgs = F)

#' Update edge capacities.
#'
#' @name updateEdgeCapacities
#' @export updateEdgeCapacities
#'
#' @param this the flow graph object
#' @param edges the vertices to update (as a 2xr matrix with ith row
#'        corresponding to the edge edges[i,1]->edges[i,2].
#' @param newCaps the new capacities for the edges
updateEdgeCapacities <- function(this, edges, newCaps) {
    UseMethod("updateEdgeCapacities")
}

#' @rdname   updateEdgeCapacities
#' @name     updateEdgeCapacities.FlowGraph
#' @export
setMethodS3("updateEdgeCapacities", "FlowGraph", function(this, edges,
    newCaps) {
    if (is.vector(edges)) {
        edges <- matrix(edges, nrow = 2)
    }
    edgeIds <- igraph::get.edge.ids(this$.flowGraph, edges + rbind(rep(this$.m, ncol(edges)),
        0))
    igraph::E(this$.flowGraph)$capacity[edgeIds] <- newCaps
}, appendVarArgs = F)
