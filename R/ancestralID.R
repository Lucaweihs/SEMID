#' Get ancestors of nodes in a graph.
#'
#' Get the ancestors of a collection of nodes in a graph g, the ancestors DO
#' include the the nodes themselves.
#'
#' @param g the graph (as an igraph).
#' @param nodes the nodes in the graph of which to get the ancestors.
#'
#' @return a sorted vector of all ancestor nodes.
ancestors <- function(g, nodes) {
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  as.numeric(sort(graph.bfs(g, nodes, neimode = "in", unreachable = F)$order,
                  na.last = NA))
  #sort(unique(unlist(neighborhood(g, vcount(g), nodes=nodes, mode="in"))))
}

#' Get parents of nodes in a graph.
#'
#' Get the parents of a collection of nodes in a graph g, the parents DO include
#' the input nodes themselves.
#'
#' @param nodes the nodes in the graph of which to get the parents.
#' @inheritParams ancestors
#'
#' @return a sorted vector of all parent nodes.
parents <- function(g, nodes) {
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  sort(unique(unlist(neighborhood(g, 1, nodes = nodes, mode = "in"))))
}

#' Get siblings of nodes in a graph.
#'
#' Get the siblings of a collection of nodes in a graph g, the siblings DO
#' include the input nodes themselves.
#'
#' @param nodes the nodes in the graph of which to get the siblings.
#' @inheritParams ancestors
#'
#' @return a sorted vector of all siblings of nodes.
siblings <- function(g, nodes) {
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  sort(unique(unlist(neighborhood(g, 1, nodes = nodes, mode = "all"))))
}

#' Getdescendants of nodes in a graph.
#'
#' Gets the descendants of a collection of nodes in a graph (all nodes that can
#' be reached by following directed edges from those nodes). Descendants DO
#' include the nodes themselves.
#'
#' @param nodes the nodes in the graph of which to get the descendants.
#' @inheritParams ancestors
#'
#' @return a sorted vector of all descendants of nodes.
descendants <- function(g, nodes) {
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  as.numeric(sort(graph.bfs(g, nodes, neimode = "out", unreachable = F)$order,
                  na.last = NA))
  #sort(unique(unlist(neighborhood(g, vcount(g), nodes=nodes, mode="out"))))
}

#' Get all HTR nodes from a set of nodes in a graph.
#'
#' Gets all vertices in a graph that are half-trek reachable from a set of
#' nodes.
#' WARNING: Often the half-trek reachable nodes from a vertex v are defined to
#' not include the vertex v or its siblings. We DO NOT follow this convention,
#' the returned set will include input nodes and their siblings.
#'
#' @inheritParams getMixedCompForNode
#' @param nodes the nodes in the graph of which to get the HTR nodes.
#'
#' @return a sorted list of all half-trek reachable nodes.
htr <- function(dG, bG, nodes) {
  if (!is.directed(dG) || is.directed(bG) || vcount(dG) != vcount(bG)) {
    stop("dG is undirected or bG is directed.")
  }
  if (vcount(dG) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  return(descendants(dG, siblings(bG, nodes)))
}

#' Get the mixed component of a node in a mixed subgraph.
#'
#' For an input mixed graph H and set of nodes A, let GA be the subgraph of
#' H on the nodes A. This function returns the mixed component of GA containing
#' a specified node.
#'
#' @param dG a directed graph representing the directed part of the mixed graph.
#' @param bG an undirected graph representing the undirected part of the mixed
#'        graph.
#' @param subNodes an ancestral set of nodes in the mixed graph, this set should
#'        include the node for which the mixed component sould be found.
#' @param node the node for which the mixed component is found.
#' @return a list with two named elements:
#'          biNodes - the nodes of the mixed graph in the biDirected component
#'                    containing nodeName w.r.t the ancestral set of nodes
#'          inNodes - the nodes in the graph which are not part of biNodes
#'                    but which are a parent of some node in biNodes.
getMixedCompForNode <- function(dG, bG, subNodes, node) {
  VdG = V(dG)
  VbG = V(bG)
  m = vcount(dG)
  if (is.null(VdG$names) || is.null(VbG$names) ||
      any(VdG$names != 1:m) || any(VbG$names != 1:m)) {
    stop(paste("Input graphs to getMixedCompForNode must have vertices named",
               "1:m in order."))
  }

  bidirectedComp =
    as.numeric(sort(graph.bfs(bG, root = node, restricted = subNodes - 1,
                              neimode = "total", unreachable = F)$order))
  incomingNodes =
    intersect(setdiff(parents(dG, bidirectedComp), bidirectedComp),
              subNodes)
  return(list(biNodes = bidirectedComp, inNodes = incomingNodes))
}


#' Size of largest HT system Y satisfying the HTC for a node v except perhaps
#' having |parents(v)| < |Y|.
#'
#' For an input mixed graph H, constructs the Gflow graph as described in Foygel
#' et al. (2012) for a subgraph G of H. A max flow algorithm is then run on
#' Gflow to determine the largest half-trek system in G to a particular node's
#' parents given a set of allowed nodes. Here G should consist of a bidirected
#' part and nodes which are not in the bidirected part but are a parent of some
#' node in the bidirected part. G should contain the node for which to compute
#' the max flow.
#'
#' @param allowedNodes the set of allowed nodes.
#' @param biNodes the set of nodes in the subgraph G which are part of the
#'        bidirected part.
#' @param inNodes the nodes of the subgraph G which are not in the bidirected
#'        part but are a parent of some node in the bidirected component.
#' @param node the node (as an integer) for which the maxflow the largest half
#' trek system
#' @inheritParams graphID
#'
#' @return See title.
#'
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
getMaxFlow <- function(L, O, allowedNodes, biNodes, inNodes, node) {
  if (!(node %in% biNodes) || any(O[biNodes, inNodes] != 0)) {
    stop(paste("When getting max flow either some in-nodes were connected to",
               "some bi-nodes or node was not in biNodes."))
  }
  if (length(intersect(allowedNodes, c(node, which(O[node,] == 1)))) != 0) {
    stop("Allowed nodes contained siblings of input node or the node itself.")
  }
  allowedNodes = union(allowedNodes, inNodes)
  m = length(biNodes) + length(inNodes)

  if (m == 1) {
    return(0)
  }

  oldNumToNewNum = numeric(m)
  oldNumToNewNum[biNodes] = 1:length(biNodes)
  oldNumToNewNum[inNodes] = (length(biNodes) + 1):m

  nodesToUse = c(biNodes, inNodes)
  allowedNodes = intersect(allowedNodes, nodesToUse)
  L[c(inNodes, biNodes), inNodes] = 0
  O[inNodes, inNodes] = 0
  L = L[nodesToUse, nodesToUse]
  O = O[nodesToUse, nodesToUse]

  # 1 & 2 = source & target
  # 2 + {1,...,m} = L(i) for i=1,...,m
  # 2 + m + {1,...,m} = R(i)-in for i=1,...,m
  # 2 + 2*m + {1,...,m} = R(i)-out for i=1,...,m

  Cap.matrix <- matrix(0, 2 + 3*m, 2 + 3*m)

  for (i in 1:m) {
    # edge from L(i) to R(i)-in, and to R(j)-in for all siblings j of i
    Cap.matrix[2 + i, 2 + m + c(i, which(O[i,] == 1))] <- 1
    # edge from R(i)-in to R(i)-out
    Cap.matrix[2 + m + i, 2 + 2*m + i] <- 1
    # edge from R(i)-out to R(j)-in for all directed edges i->j
    Cap.matrix[2 + 2*m + i, 2 + m + which(L[i,] == 1)] <- 1
  }

  allowedNodes = oldNumToNewNum[allowedNodes]
  node = oldNumToNewNum[node]
  Cap.matrix[1, 2 + allowedNodes] = 1
  Cap.matrix[2 + 2*m + which(L[,node] == 1), 2] = 1

  return(igraph::graph.maxflow(igraph::graph.adjacency(Cap.matrix),
                               source = 1, target = 2)$value)
}

#' Determine generic identifiability of an acyclic mixed graph using ancestral
#' decomposition.
#'
#' For an input, acyclic, mixed graph attempts to determine if the graph is
#' generically identifiable using decomposition by ancestral subsets. See
#' algorithm 1 of Drton and Weihs (2015).
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return The vector of nodes that could be determined to be generically
#' identifiable using the above algorithm.
#'
#' @references
#' {Drton}, M. and {Weihs}, L. (2015) Generic Identifiability of Linear
#' Structural Equation Models by Ancestor Decomposition. arXiv 1504.02992
graphID.ancestralID <- function(L, O) {
  m <- nrow(L)
  validateMatrices(L, O)
  O <- 1 * ((O + t(O)) != 0)

  dG = igraph::graph.adjacency(L)
  newOrder = as.numeric(topological.sort(dG))
  L = L[newOrder, newOrder]
  O = O[newOrder, newOrder]

  dG = igraph::graph.adjacency(L)
  bG = igraph::graph.adjacency(O)
  V(dG)$names = 1:m
  V(bG)$names = 1:m

  # Generates a list where, for each node v, we have a vector
  # corresponding to all the nodes that could ever be in a
  # half-trek system for v
  halfTrekSources = vector("list", length = m)
  for (i in 1:m) {
    halfTrekSources[[i]] = siblings(bG, ancestors(dG, i))
  }

  # A matrix determining which nodes are half-trek reachable from each node
  Dependence.matrix <- O + diag(m)
  for (i in 1:m) {
    Dependence.matrix <- ((Dependence.matrix + Dependence.matrix %*% L) > 0)
  }

  Solved.nodes <- rep(0, m)
  Solved.nodes[which(colSums(L) == 0)] <- 1 # nodes with no parents
  change <- 1
  count <- 1
  while (change == 1) {
    change <- 0

    Unsolved.nodes <- which(Solved.nodes == 0)
    for (i in Unsolved.nodes) {
      # A <- (Solved Nodes 'union' nodes not htr from i) \ ({i} 'union' sibs(i))
      A = setdiff(c(which(Solved.nodes > 0), which(Dependence.matrix[i,] == 0)),
                  c(i, which(O[i,] == 1)))
      # A <- A intersect (nodes that can ever be in a HT system for i)
      A = intersect(A, halfTrekSources[[i]])

      mixedCompList = getMixedCompForNode(dG, bG, ancestors(dG, c(i,A)), i)
      flow = getMaxFlow(L, O, A,
                        mixedCompList$biNodes, mixedCompList$inNodes, i)
      if (flow == sum(L[,i])) {
        change <- 1
        count <- count + 1
        Solved.nodes[i] <- count
        next
      }

      mixedCompList = getMixedCompForNode(dG, bG, ancestors(dG, i), i)
      A = intersect(A, unlist(mixedCompList))
      flow = getMaxFlow(L, O, A,
                        mixedCompList$biNodes, mixedCompList$inNodes, i)
      if (flow == sum(L[,i])) {
        change <- 1
        count <- count + 1
        Solved.nodes[i] <- count
        next
      }
    }
  }

  if (all(Solved.nodes == 0)) {
    Solved.nodes <- NULL
  } else {
    Solved.nodes <- order(Solved.nodes)[(1 + sum(Solved.nodes == 0)):m]
  }

  return(newOrder[Solved.nodes])
}