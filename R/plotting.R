#' Plot a mixed graph
#'
#' Given adjacency matrices representing the directed and bidirected portions
#' of a mixed graph, plots a representation of the graph.
#'
#' @inheritParams graphID
#' @param main the plot title.
#' @param vertexLabels labels to use for the vertices.
#'
#' @export
plotMixedGraph <- function(L, O, main = "", vertexLabels = 1:nrow(L)) {
    if (nrow(L) == 0) {
        stop("Can only plot graphs with >= 1 vertex")
    }
    R.utils::withSeed({
        dirGraph <- igraph::graph.adjacency(L)
        biGraph <- igraph::graph.adjacency(O, mode = "undirected")

        layout <- igraph::layout_in_circle(dirGraph)

        dirEdgeList <- igraph::get.edgelist(dirGraph)
        dirCurveMat <- (L * t(L)) != 0
        curvedOrNot <- dirCurveMat[dirEdgeList]

        igraph::V(dirGraph)$name <- vertexLabels
        plot(dirGraph, edge.color = "blue", layout = layout, edge.curved = curvedOrNot,
            vertex.size = 30, vertex.color = "lightgrey", main = main)

        biEdgeList <- igraph::get.edgelist(biGraph)
        biCurveMat <- ((L + t(L)) != 0) * O * (1 - dirCurveMat)
        curvedOrNot <- biCurveMat[biEdgeList] != 0

        igraph::V(biGraph)$name <- vertexLabels
        plot(biGraph, add = T, edge.color = "red", layout = layout, edge.curved = curvedOrNot,
            edge.arrow.mode = 3, vertex.size = 30, vertex.color = "lightgrey")
    }, 123)
}

#' Plot a latent factor graph
#'
#' Given an adjacency matrix representing the directed edges in a latent
#' factor graph, plots a representation of the graph. The latent nodes
#' should come last in L and the vertex labels should only be given for the
#' observed nodes.
#'
#' @inheritParams graphID
#' @inheritParams LatentDigraph
#' @param main the plot title.
#' @inherit LatentDigraph
#' @importFrom grDevices rgb
plotLatentDigraph <- function(L, observedNodes, latentNodes, main = "") {
  if (nrow(L) == 0) {
    stop("Can only plot graphs with >= 1 vertex")
  }
  numObs = length(observedNodes)
  numLats = length(latentNodes)

  R.utils::withSeed({
    g <- igraph::graph.adjacency(L)
    layout <- igraph::layout_nicely(g)

    edgeList <- igraph::get.edgelist(g)
    curveMat <- (L * t(L)) != 0
    curvedOrNot <- curveMat[edgeList]

    isLatentEdge <- edgeList[,1] > numObs
    edgeColors <- rep("blue", length(isLatentEdge))
    edgeColors[isLatentEdge] <- "red"

    igraph::V(g)$name <- c(observedNodes, latentNodes)

    plot(g, edge.color = edgeColors, layout = layout, edge.curved = curvedOrNot,
         vertex.size = 30,
         vertex.color = c(rep("lightblue", numObs), rep(rgb(1, 0, 0, 0.5), numLats)),
         main = main)
  }, 123)
}