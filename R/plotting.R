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
        igraph::V(dirGraph)$name <- vertexLabels
        biGraph <- igraph::graph.adjacency(O, mode = "undirected")
        igraph::V(biGraph)$name <- vertexLabels
        
        layout <- igraph::layout_in_circle(dirGraph)
        
        dirEdgeList <- igraph::get.edgelist(dirGraph)
        dirCurveMat <- (L * t(L)) != 0
        curvedOrNot <- dirCurveMat[dirEdgeList]
        
        plot(dirGraph, edge.color = "blue", layout = layout, edge.curved = curvedOrNot, 
            vertex.size = 30, vertex.color = "lightgrey", main = main)
        
        biEdgeList <- igraph::get.edgelist(biGraph)
        biCurveMat <- ((L + t(L)) != 0) * O * (1 - dirCurveMat)
        curvedOrNot <- biCurveMat[biEdgeList] != 0
        
        plot(biGraph, add = T, edge.color = "red", layout = layout, edge.curved = curvedOrNot, 
            edge.arrow.mode = 3, vertex.size = 30, vertex.color = "lightgrey")
    }, 123)
}
