#' Create an ancestral identification function.
#'
#' A helper function for \code{\link{ancestralIdentifyStep}}, creates an identifier function
#' based on its given parameters. This created identifier function will
#' identify the directed edges from 'targets' to 'node.'
#'
#' @inheritParams createHtcIdentifier
#' @param ancestralSubset an ancestral subset of the graph containing \code{node}.
#' @param cComponent a list corresponding to the connected component containing \code{node}
#'        in the subgraph induced by \code{ancestralSubset}. See
#'        \code{\link{tianDecompose}} for how such connected component lists are formed.
#'
#' @return an identification function
createAncestralIdentifier <- function(idFunc, sources, targets, node, htrSources,
    ancestralSubset, cComponent) {
    # Necessary 'redundant' assignments
    idFunc <- idFunc
    sources <- sources
    targets <- targets
    node <- node
    htrSources <- htrSources
    ancestralSubset <- sort(ancestralSubset)
    cComponent <- cComponent
    return(function(Sigma) {
        m <- nrow(Sigma)
        identifiedParams <- idFunc(Sigma)
        Lambda <- identifiedParams$Lambda

        topOrder <- cComponent$topOrder

        afterAnc <- function(x) {
            which(x == ancestralSubset)
        }
        topOrderAfterAnc <- sapply(topOrder, afterAnc)
        sourcesAfterAnc <- sapply(sources, afterAnc)
        targetsAfterAnc <- sapply(targets, afterAnc)
        htrSourcesAfterAnc <- sapply(htrSources, afterAnc)
        internalAfterAnc <- sapply(cComponent$internal, afterAnc)
        incomingAfterAnc <- sapply(cComponent$incoming, afterAnc)
        nodeAfterAnc <- afterAnc(node)

        SigmaAnc <- Sigma[ancestralSubset, ancestralSubset]
        SigmaAncTian <- tianSigmaForComponent(SigmaAnc, internalAfterAnc, incomingAfterAnc,
            topOrderAfterAnc)

        afterAncTian <- function(x) {
            which(x == topOrderAfterAnc)
        }
        sourcesAfterAncTian <- sapply(sourcesAfterAnc, afterAncTian)
        targetsAfterAncTian <- sapply(targetsAfterAnc, afterAncTian)
        htrSourcesAfterAncTian <- sapply(htrSourcesAfterAnc, afterAncTian)
        nodeAfterAncTian <- afterAncTian(nodeAfterAnc)

        LambdaAfterAncTian <- Lambda[ancestralSubset, ancestralSubset, drop = F]
        LambdaAfterAncTian <- LambdaAfterAncTian[topOrderAfterAnc, topOrderAfterAnc,
            drop = F]

        SigmaMinus <- SigmaAncTian
        for (source in htrSourcesAfterAncTian) {
            SigmaMinus[source, ] <- SigmaAncTian[source, ] - t(LambdaAfterAncTian[,
                source, drop = F]) %*% SigmaAncTian
        }

        if (abs(det(SigmaMinus[sourcesAfterAncTian, targetsAfterAncTian, drop = F])) <
            10^-10) {
            stop("In identification, found near-singular system. Is the input matrix generic?")
        }

        sols <- solve(SigmaMinus[sourcesAfterAncTian, targetsAfterAncTian, drop = F],
            SigmaMinus[sourcesAfterAncTian, nodeAfterAncTian, drop = F])

        undoAncTian <- function(x) {
            ancestralSubset[topOrderAfterAnc[x]]
        }
        Lambda[undoAncTian(targetsAfterAncTian), node] <- sols

        return(list(Lambda = Lambda, Omega = identifiedParams$Omega))
    })
}

#' Perform one iteration of ancestral identification.
#'
#' A function that does one step through all the nodes in a mixed graph
#' and tries to determine if directed edge coefficients are generically
#' identifiable by leveraging decomposition by ancestral subsets. See
#' Algorithm 1 of Drton and Weihs (2015); this version of the algorithm
#' is somewhat different from Drton and Weihs (2015) in that it also works
#' on cyclic graphs.
#'
#' @inheritParams htcIdentifyStep
#'
#' @return a list with four components:
#' \describe{
#'   \item{\code{identifiedEdges}}{a matrix rx2 matrix where r is the number
#'   of edges that where identified by this function call and
#'   \code{identifiedEdges[i,1] -> identifiedEdges[i,2]} was the ith edge
#'   identified}
#'   \item{\code{unsolvedParents}}{as the input argument but updated with
#'   any newly identified edges}
#'   \item{\code{solvedParents}}{as the input argument but updated with
#'   any newly identified edges}
#'   \item{\code{identifier}}{as the input argument but updated with
#'   any newly identified edges}
#' }
#'
#' @export
#'
#' @references
#' {Drton}, M. and {Weihs}, L. (2015) Generic Identifiability of Linear
#' Structural Equation Models by Ancestor Decomposition. arXiv 1504.02992
ancestralIdentifyStep <- function(mixedGraph, unsolvedParents, solvedParents, identifier) {
    identifiedEdges <- c()
    m <- mixedGraph$numNodes()
    ancestralComps <- rep(list(list()), m)
    solvedNodes <- which(sapply(unsolvedParents, FUN = function(x) {
        length(x) == 0
    }))
    for (i in 1:m) {
        unsolved <- unsolvedParents[[i]]
        if (length(unsolved) != 0) {

            if (length(ancestralComps[[i]]) == 0) {
                nodeAncestors <- mixedGraph$ancestors(i)
                ancGraph <- mixedGraph$inducedSubgraph(nodeAncestors)
                tianComp <- ancGraph$tianComponent(i)
                tianComp$mixedGraph <- MixedGraph(tianComp$L, tianComp$O, vertexNums = tianComp$topOrder)
                ancestralComps[[i]] <- tianComp
                ancestralComps[[i]]$ancestors <- nodeAncestors
            }

            nodeParents <- ancGraph$parents(i)

            # Using the first ancestor graph
            ancGraph <- ancestralComps[[i]]$mixedGraph
            htrFromNode <- ancGraph$htrFrom(i)
            allowedNodes <- setdiff(solvedNodes, ancGraph$siblings(i))
            allowedNodes <- union(setdiff(ancGraph$nodes(), htrFromNode), allowedNodes)
            allowedNodes <- intersect(ancGraph$nodes(), allowedNodes)
            if (length(allowedNodes) >= length(nodeParents)) {
                halfTrekSystemResult <- ancGraph$getHalfTrekSystem(allowedNodes,
                  nodeParents)
                if (halfTrekSystemResult$systemExists) {
                  identifiedEdges <- c(identifiedEdges, as.integer(rbind(nodeParents,
                    i)))
                  activeFrom <- halfTrekSystemResult$activeFrom
                  identifier <- createAncestralIdentifier(identifier, activeFrom,
                    nodeParents, i, intersect(activeFrom, htrFromNode), ancestralComps[[i]]$ancestors,
                    ancestralComps[[i]])
                  solvedParents[[i]] <- nodeParents
                  unsolvedParents[[i]] <- integer(0)
                  solvedNodes <- c(i, solvedNodes)
                }
                next
            }

            # Using the second ancestor graph
            nodeAncestors <- mixedGraph$ancestors(i)
            nodeAncSibs <- union(nodeAncestors, mixedGraph$siblings(nodeAncestors))
            allowedAncSibs <- intersect(nodeAncSibs, setdiff(union(solvedNodes, setdiff(1:m,
                mixedGraph$htrFrom(i))), c(i, mixedGraph$siblings(i))))
            ancGraph <- mixedGraph$inducedSubgraph(mixedGraph$ancestors(c(i, allowedAncSibs)))
            tianComp <- ancGraph$tianComponent(i)
            allowedNodes <- union(intersect(allowedAncSibs, tianComp$internal), tianComp$incoming)

            htrFromNode <- ancGraph$htrFrom(i)
            allowedNodes <- setdiff(solvedNodes, ancGraph$siblings(i))
            allowedNodes <- union(setdiff(ancGraph$nodes(), htrFromNode), allowedNodes)
            allowedNodes <- intersect(ancGraph$nodes(), allowedNodes)

            tianAncGraph <- MixedGraph(tianComp$L, tianComp$O, tianComp$topOrder)
            if (length(allowedNodes) >= length(nodeParents)) {
                halfTrekSystemResult <- tianAncGraph$getHalfTrekSystem(allowedNodes,
                  nodeParents)
                if (halfTrekSystemResult$systemExists) {
                  identifiedEdges <- c(identifiedEdges, as.integer(rbind(nodeParents,
                    i)))
                  activeFrom <- halfTrekSystemResult$activeFrom
                  identifier <- createAncestralIdentifier(identifier, activeFrom,
                    nodeParents, i, intersect(activeFrom, htrFromNode), ancGraph$nodes(),
                    tianComp)
                  solvedParents[[i]] <- nodeParents
                  unsolvedParents[[i]] <- integer(0)
                  solvedNodes <- c(i, solvedNodes)
                }
            }
        }
    }
    return(list(identifiedEdges = identifiedEdges, unsolvedParents = unsolvedParents,
        solvedParents = solvedParents, identifier = identifier))
}

#' Determines which edges in a mixed graph are ancestralID-identifiable
#'
#' Uses the an identification criterion of Drton and Weihs (2015); this version
#' of the algorithm is somewhat different from Drton and Weihs (2015) in that it
#' also works on cyclic graphs. The original version of the algorithm can be
#' found in the function \code{\link{graphID.ancestralID}}.
#'
#' @export
#'
#' @inheritParams generalGenericID
#' @inheritParams semID
#' @inheritParams ancestralIdentifyStep
#'
#' @return see the return of \code{\link{generalGenericID}}.
ancestralID <- function(mixedGraph, tianDecompose = T) {
    result <- generalGenericID(mixedGraph, list(ancestralIdentifyStep), tianDecompose = tianDecompose)
    result$call <- match.call()
    return(result)
}

#' Get getAncestors of nodes in a graph.
#'
#' Get the getAncestors of a collection of nodes in a graph g, the getAncestors DO
#' include the the nodes themselves.
#'
#' @param g the graph (as an igraph).
#' @param nodes the nodes in the graph of which to get the getAncestors.
#'
#' @return a sorted vector of all ancestor nodes.
getAncestors <- function(g, nodes) {
    if (vcount(g) == 0 || length(nodes) == 0) {
        return(numeric(0))
    }
    as.integer(sort(graph.bfs(g, nodes, mode = "in", unreachable = F)$order, na.last = NA))
    # sort(unique(unlist(neighborhood(g, vcount(g), nodes=nodes, mode='in'))))
}

#' Get getParents of nodes in a graph.
#'
#' Get the getParents of a collection of nodes in a graph g, the getParents DO include
#' the input nodes themselves.
#'
#' @param nodes the nodes in the graph of which to get the getParents.
#' @inheritParams getAncestors
#'
#' @return a sorted vector of all parent nodes.
getParents <- function(g, nodes) {
    if (vcount(g) == 0 || length(nodes) == 0) {
        return(numeric(0))
    }
    sort(unique(unlist(neighborhood(g, 1, nodes = nodes, mode = "in"))))
}

#' Get getSiblings of nodes in a graph.
#'
#' Get the getSiblings of a collection of nodes in a graph g, the getSiblings DO
#' include the input nodes themselves.
#'
#' @param nodes the nodes in the graph of which to get the getSiblings.
#' @inheritParams getAncestors
#'
#' @return a sorted vector of all getSiblings of nodes.
getSiblings <- function(g, nodes) {
    if (vcount(g) == 0 || length(nodes) == 0) {
        return(numeric(0))
    }
    sort(unique(unlist(neighborhood(g, 1, nodes = nodes, mode = "all"))))
}

#' Get descendants of nodes in a graph.
#'
#' Gets the descendants of a collection of nodes in a graph (all nodes that can
#' be reached by following directed edges from those nodes). Descendants DO
#' include the nodes themselves.
#'
#' @param nodes the nodes in the graph of which to get the descendants.
#' @inheritParams getAncestors
#'
#' @return a sorted vector of all descendants of nodes.
getDescendants <- function(g, nodes) {
    if (vcount(g) == 0 || length(nodes) == 0) {
        return(numeric(0))
    }
    as.integer(sort(graph.bfs(g, nodes, mode = "out", unreachable = F)$order,
        na.last = NA))
    # sort(unique(unlist(neighborhood(g, vcount(g), nodes=nodes, mode='out'))))
}

#' Get all HTR nodes from a set of nodes in a graph.
#'
#' Gets all vertices in a graph that are half-trek reachable from a set of
#' nodes.
#' WARNING: Often the half-trek reachable nodes from a vertex v are defined to
#' not include the vertex v or its getSiblings. We DO NOT follow this convention,
#' the returned set will include input nodes and their getSiblings.
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
    return(getDescendants(dG, getSiblings(bG, nodes)))
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
    VdG <- V(dG)
    VbG <- V(bG)
    m <- vcount(dG)
    if (is.null(VdG$names) || is.null(VbG$names) || any(VdG$names != 1:m) || any(VbG$names !=
        1:m)) {
        stop(paste("Input graphs to getMixedCompForNode must have vertices named",
            "1:m in order."))
    }

    ## This is a workaround of a graph.bfs bug in igraph
    ## versions 1.0.1 and before
    restricted = subNodes
    if (utils::packageVersion("igraph") <= "1.0.1") {
      restricted = restricted - 1
    }

    bidirectedComp <- as.integer(sort(graph.bfs(bG, root = node,
                                                restricted = restricted,
                                                mode = "total",
                                                unreachable = F)$order))

    bidirectedComp =
      as.integer(sort(graph.bfs(bG, root = node, restricted = restricted,
                                mode = "total", unreachable = F)$order))


    incomingNodes <- intersect(setdiff(getParents(dG, bidirectedComp), bidirectedComp),
        subNodes)
    return(list(biNodes = bidirectedComp, inNodes = incomingNodes))
}


#' Size of largest HT system Y satisfying the HTC for a node v except perhaps
#' having |getParents(v)| < |Y|.
#'
#' For an input mixed graph H, constructs the Gflow graph as described in Foygel
#' et al. (2012) for a subgraph G of H. A max flow algorithm is then run on
#' Gflow to determine the largest half-trek system in G to a particular node's
#' getParents given a set of allowed nodes. Here G should consist of a bidirected
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
    if (length(intersect(allowedNodes, c(node, which(O[node, ] == 1)))) != 0) {
        stop("Allowed nodes contained getSiblings of input node or the node itself.")
    }
    allowedNodes <- union(allowedNodes, inNodes)
    m <- length(biNodes) + length(inNodes)

    if (m == 1) {
        return(0)
    }

    oldNumToNewNum <- numeric(m)
    oldNumToNewNum[biNodes] <- 1:length(biNodes)
    oldNumToNewNum[inNodes] <- (length(biNodes) + 1):m

    nodesToUse <- c(biNodes, inNodes)
    allowedNodes <- intersect(allowedNodes, nodesToUse)
    L[c(inNodes, biNodes), inNodes] <- 0
    O[inNodes, inNodes] <- 0
    L <- L[nodesToUse, nodesToUse]
    O <- O[nodesToUse, nodesToUse]

    # 1 & 2 = source & target 2 + {1,...,m} = L(i) for i=1,...,m 2 + m + {1,...,m} =
    # R(i)-in for i=1,...,m 2 + 2*m + {1,...,m} = R(i)-out for i=1,...,m

    Cap.matrix <- matrix(0, 2 + 3 * m, 2 + 3 * m)

    for (i in 1:m) {
        # edge from L(i) to R(i)-in, and to R(j)-in for all getSiblings j of i
        Cap.matrix[2 + i, 2 + m + c(i, which(O[i, ] == 1))] <- 1
        # edge from R(i)-in to R(i)-out
        Cap.matrix[2 + m + i, 2 + 2 * m + i] <- 1
        # edge from R(i)-out to R(j)-in for all directed edges i->j
        Cap.matrix[2 + 2 * m + i, 2 + m + which(L[i, ] == 1)] <- 1
    }

    allowedNodes <- oldNumToNewNum[allowedNodes]
    node <- oldNumToNewNum[node]
    Cap.matrix[1, 2 + allowedNodes] <- 1
    Cap.matrix[2 + 2 * m + which(L[, node] == 1), 2] <- 1

    return(igraph::graph.maxflow(igraph::graph.adjacency(Cap.matrix), source = 1,
        target = 2)$value)
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

    dG <- igraph::graph.adjacency(L)
    newOrder <- as.integer(topological.sort(dG))
    L <- L[newOrder, newOrder]
    O <- O[newOrder, newOrder]

    dG <- igraph::graph.adjacency(L)
    bG <- igraph::graph.adjacency(O)
    V(dG)$names <- 1:m
    V(bG)$names <- 1:m

    # Generates a list where, for each node v, we have a vector corresponding to all
    # the nodes that could ever be in a half-trek system for v
    halfTrekSources <- vector("list", length = m)
    for (i in 1:m) {
        halfTrekSources[[i]] <- getSiblings(bG, getAncestors(dG, i))
    }

    # A matrix determining which nodes are half-trek reachable from each node
    Dependence.matrix <- O + diag(m)
    for (i in 1:m) {
        Dependence.matrix <- ((Dependence.matrix + Dependence.matrix %*% L) > 0)
    }

    Solved.nodes <- rep(0, m)
    Solved.nodes[which(colSums(L) == 0)] <- 1  # nodes with no getParents
    change <- 1
    count <- 1
    while (change == 1) {
        change <- 0

        Unsolved.nodes <- which(Solved.nodes == 0)
        for (i in Unsolved.nodes) {
            # A <- (Solved Nodes 'union' nodes not htr from i) \ ({i} 'union' sibs(i))
            A <- setdiff(c(which(Solved.nodes > 0), which(Dependence.matrix[i, ] ==
                0)), c(i, which(O[i, ] == 1)))
            # A <- A intersect (nodes that can ever be in a HT system for i)
            A <- intersect(A, halfTrekSources[[i]])

            mixedCompList <- getMixedCompForNode(dG, bG, getAncestors(dG, c(i, A)),
                i)
            flow <- getMaxFlow(L, O, A, mixedCompList$biNodes, mixedCompList$inNodes,
                i)
            if (flow == sum(L[, i])) {
                change <- 1
                count <- count + 1
                Solved.nodes[i] <- count
                next
            }

            mixedCompList <- getMixedCompForNode(dG, bG, getAncestors(dG, i), i)
            A <- intersect(A, unlist(mixedCompList))
            flow <- getMaxFlow(L, O, A, mixedCompList$biNodes, mixedCompList$inNodes,
                i)
            if (flow == sum(L[, i])) {
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
