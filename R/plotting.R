#' Plot a mixed graph
#'
#' Given adjacency matrices representing the directed and bidirected portions
#' of a mixed graph, plots a representation of the graph.
#'
#' @inheritParams graphID
#' @param vertexLabels labels to use for the vertices.
#'
#' @export
plotMixedGraph <- function(L, O, main = "", vertexLabels = 1:nrow(L)) {
  R.utils::withSeed({
  m = ncol(L)
  dirGraph = igraph::graph.adjacency(L)
  igraph::V(dirGraph)$name = vertexLabels
  biGraph = igraph::graph.adjacency(O, mode = "undirected")
  igraph::V(biGraph)$name = vertexLabels

  layout = igraph::layout_in_circle(dirGraph)

  plot(dirGraph, edge.color = "blue", layout = layout, vertex.size = 30,
       vertex.color = "lightgrey", main = main)

  biEdgeList = igraph::get.edgelist(biGraph)
  curvedOrNot = logical(nrow(biEdgeList))
  if (nrow(biEdgeList) != 0) {
    for (i in 1:nrow(biEdgeList)) {
      v1 = biEdgeList[i,1]
      v2 = biEdgeList[i,2]
      curvedOrNot[i] = (L[v1,v2] || L[v2,v1])
    }
  }
  plot(biGraph, add = T, edge.color = "red", layout = layout,
       edge.curved = curvedOrNot, edge.arrow.mode = 3, vertex.size = 30,
       vertex.color = "lightgrey")
  }, 123)
}
