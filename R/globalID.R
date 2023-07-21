#' Get bidirected components of a mixed graph
#'
#' Returns induced subgraphs of connected bidirected components with more than 1 node.
#'
#' @export
#'
#' @param graph a \code{\link{MixedGraph}} object representing
#'         the mixed graph.
#'
#' @return list, where each object is a MixedGraph with at least two nodes.
bidirectedComponents <- function(graph){
  biGraph <- igraph::graph.adjacency(graph$O(), mode = "undirected")
  igraph::V(biGraph)$name <- graph$nodes()
  comps = igraph::components(biGraph)
  res= list()
  count = 1
  for (i in 1:comps$no){
    if (comps$csize[i] > 1){
      res[[count]] = graph$inducedSubgraph(as.integer(names(comps$membership[comps$membership==i])))
      count = count + 1
    }
  }
  return(res)
}

#' Determines whether a mixed graph is globally identifiable.
#'
#' Uses the criterion in Theorem 2 of the paper by Drton, Foygel and Sullivant (2011) to determine
#' whether a mixed graph is globally identifiable.
#'
#' @export
#'
#' @param graph a \code{\link{MixedGraph}} object representing
#'         the mixed graph.
#'
#' @return TRUE if the graph is globally identifiable, FALSE otherwise.
#'
#' @references
#' Drton, M., Barber, R. F., and Sullivant S. (2011).
#' Half-Trek Criterion for Identifiability of Latent Variable Models.
#' Ann. Statist. 39 (2011), no. 2, 865--886 <doi:10.1214/10-AOS859>.
globalID <- function(graph){

  # Check whether graph is acyclic (Theorem 1)
  dirGraph <- igraph::graph.adjacency(graph$L())
  if (!igraph::is_dag(dirGraph)){
    return(FALSE)
  }

  # Otherwise check criterion (Theorem 2)
  comps = bidirectedComponents(graph)
  while (length(comps)!=0){

    comp = comps[[1]]
    comps = comps[-1]

    dirGraph <- igraph::graph.adjacency(comp$L())
    igraph::V(dirGraph)$name <- comp$nodes()
    sinks <- as.integer(names(which(igraph::degree(dirGraph, mode = "out")==0)))
    if (length(sinks)==1){
      return(FALSE)
    } else {
      for (s in sinks){
        ancestors = comp$ancestors(s)
        if (length(ancestors) > 1){
          ancestralGraph = comp$inducedSubgraph(ancestors)
          new_comps = bidirectedComponents(ancestralGraph)
          comps = c(comps, new_comps)
        }
      }
    }
  }
  return(TRUE)
}