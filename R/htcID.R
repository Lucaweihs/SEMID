#' Create an htc identification function.
#'
#' A helper function for \code{\link{htcIdentifyStep}}, creates an identifier
#' function based on its given parameters. This created identifier function will
#' identify the directed edges from 'targets' to 'node.'
#'
#' @param idFunc identification of edge coefficients often requires that other
#'        edge coefficients already be identified. This argument should be a
#'        function that produces all such identifications. The newly created
#'        identifier function will return these identifications along with its
#'        own.
#' @param sources the sources of the half-trek system.
#' @param targets the targets of the half-trek system (these should be the
#'        parents of \code{node}).
#' @param node the node for which all incoming edges are to be identified
#'        (the tails of which are targets).
#' @param htrSources the nodes in sources which are half-trek reachable from
#'        \code{node}. All incoming edges to these sources should be identified by
#'        \code{idFunc} for the newly created identification function to work.
#'
#' @return an identification function
#'
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
createHtcIdentifier <- function(idFunc, sources, targets, node, htrSources) {
    # Necessary redundant assignments
    idFunc <- idFunc
    sources <- sources
    targets <- targets
    node <- node
    htrSources <- htrSources
    return(function(Sigma) {
        m <- nrow(Sigma)
        identifiedParams <- idFunc(Sigma)
        Lambda <- identifiedParams$Lambda

        SigmaMinus <- Sigma
        for (source in htrSources) {
            SigmaMinus[source, ] <- Sigma[source, ] - t(Lambda[, source, drop = F]) %*%
                Sigma
        }

        if (abs(det(SigmaMinus[sources, targets, drop = F])) < 10^-10) {
            stop("In identification, found near-singular system. Is the input matrix generic?")
        }

        Lambda[targets, node] <- solve(SigmaMinus[sources, targets, drop = F], SigmaMinus[sources,
            node, drop = F])

        if (!any(is.na(Lambda))) {
            Omega <- t(diag(m) - Lambda) %*% Sigma %*% (diag(m) - Lambda)
            return(list(Lambda = Lambda, Omega = Omega))
        }
        return(list(Lambda = Lambda, Omega = identifiedParams$Omega))
    })
}

#' Perform one iteration of HTC identification.
#'
#' A function that does one step through all the nodes in a mixed graph
#' and tries to identify new edge coefficients using the existence of
#' half-trek systems as described in Foygel, Draisma, Drton (2012).
#'
#' @param mixedGraph a \code{\link{MixedGraph}} object representing
#'         the mixed graph.
#' @param unsolvedParents a list whose ith index is a vector of all the parents
#'        j of i in G which for which the edge j->i is not yet known to be
#'        generically identifiable.
#' @param solvedParents the complement of \code{unsolvedParents}, a list whose
#'        ith index is a vector of all parents j of i for which the edge i->j
#'        is known to be generically identifiable (perhaps by other algorithms).
#' @param identifier an identification function that must produce the
#'        identifications corresponding to those in solved parents. That is
#'        \code{identifier} should be a function taking a single argument Sigma
#'        (any generically generated covariance matrix corresponding
#'        to the mixed graph) and returns a list with two named arguments
#' \describe{
#'   \item{Lambda}{denote the number of nodes in \code{mixedGraph} as n. Then
#'                 Lambda is an nxn matrix whose i,jth entry
#' \enumerate{
#'   \item equals 0 if i is not a parent of j,
#'   \item equals NA if i is a parent of j but \code{identifier} cannot
#'         identify it generically,
#'   \item equals the (generically) unique value corresponding to the weight
#'         along the edge i->j that was used to produce Sigma.
#' }}
#'   \item{Omega}{just as Lambda but for the bidirected edges in the mixed
#'                graph}
#' }
#'        such that if j is in \code{solvedParents[[i]]} we must have that
#'        Lambda[j,i] is not NA.
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
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713
htcIdentifyStep <- function(mixedGraph, unsolvedParents, solvedParents, identifier) {
    identifiedEdges <- numeric(0)
    allNodes <- mixedGraph$nodes()
    solvedNodes <- which(sapply(unsolvedParents, FUN = function(x) {
        length(x) == 0
    }))
    for (i in setdiff(allNodes, solvedNodes)) {
        htrFromNode <- mixedGraph$htrFrom(i)
        allowedNodes <- setdiff(solvedNodes, mixedGraph$siblings(i))
        allowedNodes <- union(setdiff(allNodes, htrFromNode), allowedNodes)
        nodeParents <- mixedGraph$parents(i)
        if (length(allowedNodes) < length(nodeParents)) {
            next
        }
        halfTrekSystemResult <- mixedGraph$getHalfTrekSystem(allowedNodes, nodeParents)

        if (halfTrekSystemResult$systemExists) {
            identifiedEdges <- c(identifiedEdges, as.integer(rbind(nodeParents, i)))
            activeFrom <- halfTrekSystemResult$activeFrom
            identifier <- createHtcIdentifier(identifier, activeFrom, nodeParents,
                i, intersect(activeFrom, htrFromNode))
            solvedParents[[i]] <- nodeParents
            unsolvedParents[[i]] <- integer(0)
            solvedNodes <- c(i, solvedNodes)
        }
    }
    return(list(identifiedEdges = matrix(identifiedEdges, byrow = T, ncol = 2), unsolvedParents = unsolvedParents,
        solvedParents = solvedParents, identifier = identifier))
}

#' Determines which edges in a mixed graph are HTC-identifiable.
#'
#' Uses the half-trek criterion of Foygel, Draisma, and Drton (2012) determine
#' which edges in a mixed graph are generically identifiable. Depending on your
#' application it faster to use the \code{\link{graphID.htcID}} function
#' instead of this one, this function has the advantage of returning additional
#' information.
#'
#' @export
#'
#' @inheritParams generalGenericID
#' @inheritParams semID
#'
#' @return see the return value of \code{\link{generalGenericID}}.
#'
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
#'
#' Jin Tian. 2005. Identifying direct causal effects in linear models. In
#' \emph{Proceedings of the 20th national conference on Artificial intelligence
#' - Volume 1} (AAAI'05), Anthony Cohn (Ed.), Vol. 1. AAAI Press 346-352.
htcID <- function(mixedGraph, tianDecompose = T) {
    result <- generalGenericID(mixedGraph, list(htcIdentifyStep), tianDecompose = tianDecompose)
    result$call <- match.call()
    return(result)
}

#' Determines if a mixed graph is HTC-identifiable.
#'
#' Uses the half-trek criterion of Foygel, Draisma, and Drton (2013) to check
#' if an input mixed graph is generically identifiable.
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return The vector of HTC-identifiable nodes.
#'
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
graphID.htcID <- function(L, O) {
    m <- nrow(L)
    validateMatrices(L, O)
    O <- 1 * ((O + t(O)) != 0)

    # 1 & 2 = source & target 2 + {1,...,m} = L(i) for i=1,...,m 2+m + {1,...,m} =
    # R(i)-in for i=1,...,m 2+2*m + {1,...,m} = R(i)-out for i=1,...,m

    Cap.matrix.init <- matrix(0, 2 + 3 * m, 2 + 3 * m)
    for (i in 1:m) {
        # edge from L(i) to R(i)-in, and to R(j)-in for all siblings j of i
        Cap.matrix.init[2 + i, 2 + m + c(i, which(O[i, ] == 1))] <- 1
        # edge from R(i)-in to R(i)-out
        Cap.matrix.init[2 + m + i, 2 + 2 * m + i] <- 1
        # edge from R(i)-out to R(j)-in for all directed edges i->j
        Cap.matrix.init[2 + 2 * m + i, 2 + m + which(L[i, ] == 1)] <- 1
    }

    # when testing if a set A satisfies the HTC with respect to a node i, need to add
    # (1) edge from source to L(j) for all j in A and (2) edge from R(j)-out to
    # target for all parents j of i

    Dependence.matrix <- O + diag(m)
    for (i in 1:m) {
        Dependence.matrix <- (Dependence.matrix + Dependence.matrix %*% L > 0)
    }

    Solved.nodes <- rep(0, m)
    Solved.nodes[which(colSums(L) == 0)] <- 1  # nodes with no parents
    change <- 1
    count <- 1

    while (change == 1) {
        change <- 0

        for (i in which(Solved.nodes == 0)) {
            A <- setdiff(c(which(Solved.nodes > 0), which(Dependence.matrix[i, ] ==
                0)), c(i, which(O[i, ] == 1)))

            Cap.matrix <- Cap.matrix.init

            Cap.matrix[1, 2 + A] <- 1
            Cap.matrix[2 + 2 * m + which(L[, i] == 1), 2] <- 1

            flow <- graph.maxflow(graph.adjacency(Cap.matrix), source = 1, target = 2)$value

            if (flow == sum(L[, i])) {
                change <- 1
                count <- count + 1
                Solved.nodes[i] <- count
            }
        }
    }

    if (all(Solved.nodes == 0)) {
        Solved.nodes <- NULL
    } else {
        Solved.nodes <- order(Solved.nodes)[(1 + sum(Solved.nodes == 0)):m]
    }
    return(Solved.nodes)
}

#' Check for generic infinite-to-one via the half-trek criterion.
#'
#' Checks if a mixed graph is infinite-to-one using the half-trek criterion
#' presented by Foygel, Draisma, and Drton (2012).
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return TRUE if the graph could be determined to be generically
#' non-identifiable, FALSE if this test was inconclusive.
#'
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
graphID.nonHtcID <- function(L, O) {
    m <- nrow(L)
    validateMatrices(L, O)
    if (m == 1) {
        # Always identifiable in this case
        return(F)
    }
    O <- 1 * ((O + t(O)) != 0)

    # 1 & 2 = source & target 2 + {1,...,N} = L{i_n,j_n} for the n-th nonsibling
    # pair, n=1,...,N (2+N) + 2*m^2 = 2 copies of R_i(j) -- in copy & out copy where
    # (2+N) + (i-1)*m + j = R_i(j) in copy & (2+N+m^2) + (i-1)*m + j = R_i(j) out
    # copy

    nonsibs <- NULL
    N <- 0
    for (i in 1:(m - 1)) {
        for (j in (i + 1):m) {
            if (O[i, j] == 0) {
                N <- N + 1
                nonsibs <- rbind(nonsibs, c(i, j))
            }
        }
    }

    Cap.matrix <- matrix(0, 2 * m^2 + N + 2, 2 * m^2 + N + 2)
    if (N != 0) {
        Cap.matrix[1, 2 + (1:N)] <- 1  # edges from source to L{i,j} for each
        # nonsibling pair {i,j} = nonsibs[n,1:2] edge from L{i,j} to R_i(j), and to
        # R_i(k)-in for all siblings k of node j
        for (n in 1:N) {
            Cap.matrix[2 + n, 2 + N + (nonsibs[n, 1] - 1) * m + c(nonsibs[n, 2],
                which(O[nonsibs[n, 2], ] == 1))] <- 1
            # edge from L{i,j} to R_j(i), and to R_j(i)-in for all siblings k of node i
            Cap.matrix[2 + n, 2 + N + (nonsibs[n, 2] - 1) * m + c(nonsibs[n, 1],
                which(O[nonsibs[n, 1], ] == 1))] <- 1
        }
    }
    for (i in 1:m) {
        # edge from R_i(j)-out to target when j is a parent of i
        Cap.matrix[2 + N + m^2 + (i - 1) * m + which(L[, i] == 1), 2] <- 1
        for (j in 1:m) {
            # edge from R_i(j)-in to R_i(j)-out
            Cap.matrix[2 + N + (i - 1) * m + j, 2 + N + m^2 + (i - 1) * m + j] <- 1
            # edge from R_i(j)-out to R_i(k)-in where j->k is a directed edge
            Cap.matrix[2 + N + m^2 + (i - 1) * m + j, 2 + N + (i - 1) * m + which(L[j,
                ] == 1)] <- 1
        }
    }

    HTC.nonID <- graph.maxflow(igraph::graph.adjacency(Cap.matrix), source = 1, target = 2)$value
    return(HTC.nonID < sum(L))
}
