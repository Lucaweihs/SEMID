#' Create an htc identification function
#'
#' TODO: Add details
#'
#' @return an identification function
createHtcIdentifier <- function(idFunc, sources, targets, node,
                                htrSources) {
  # Necessary redundant assignments
  idFunc <- idFunc
  sources <- sources
  targets <- targets
  node <- node
  htrSources <- htrSources
  return(
    function(Sigma) {
      m <- nrow(Sigma)
      identifiedParams <- idFunc(Sigma)
      Lambda <- identifiedParams$Lambda

      SigmaMinus = Sigma
      for (source in htrSources) {
        SigmaMinus[source,] = Sigma[source,] - t(Lambda[,source,drop = F]) %*% Sigma
      }

      if (abs(det(SigmaMinus[sources, targets, drop = F])) < 10^-10) {
        stop("In identification, found near-singular system. Is the input matrix generic?")
      }

      Lambda[targets, node] <-
        solve(SigmaMinus[sources, targets, drop = F],
              SigmaMinus[sources, node, drop = F])

      if (!any(is.na(Lambda))) {
        Omega = t(diag(m) - Lambda) %*% Sigma %*% (diag(m) - Lambda)
        return(list(Lambda = Lambda, Omega = Omega))
      }
      return(list(Lambda = Lambda, Omega = identifiedParams$Omega))
    }
  )
}

#' Perform one iteration of htc identification.
#'
#' A function that does one step through all the nodes in a mixed graph
#' and tries to identify new edge coefficients using half-treks.
#'
#' @export
#'
#' @return a list
htcIdentifyStep = function(mixedGraph, unsolvedParents, solvedParents,
                           identifier) {
  identifiedEdges = c()
  m = mixedGraph$numNodes()
  solvedNodes = which(sapply(unsolvedParents, FUN = function(x) { length(x) == 0 }))
  for (i in setdiff(1:m, solvedNodes)) {
    htrFromNode = mixedGraph$htrFrom(i)
    allowedNodes = setdiff(solvedNodes, mixedGraph$allSiblings(i))
    allowedNodes = union(setdiff(1:m, htrFromNode),
                         intersect(allowedNodes, htrFromNode))
    nodeParents = mixedGraph$allParents(i)
    if (length(allowedNodes) < length(nodeParents)) {
      next
    }
    halfTrekSystemResult = mixedGraph$getHalfTrekSystem(allowedNodes,
                                                        nodeParents)

    if (halfTrekSystemResult$systemExists) {
      identifiedEdges = c(identifiedEdges, as.numeric(rbind(nodeParents, i)))
      activeFrom = halfTrekSystemResult$activeFrom
      identifier = createHtcIdentifier(identifier, activeFrom, nodeParents, i,
                                       intersect(activeFrom, htrFromNode))
      solvedParents[[i]] = nodeParents
      unsolvedParents[[i]] = integer()
      solvedNodes = c(i, solvedNodes)
    }
  }
  return(list(identifiedEdges = identifiedEdges, unsolvedParents = unsolvedParents,
              solvedParents = solvedParents, identifier = identifier))
}

#' Determines which edges in a mixed graph are HTC-identifiable.
#'
#' Uses the half-trek criterion of Foygel, Draisma, and Drton (2013) determine
#' which edges in a mixed graph are generically identifiable.
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return a list
htcID <- function(L, O) {
  return(generalGenericID(L, O, list(htcIdentifyStep)))
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
  #return(generalGenericID(L, O, list(htcIdentifyStep)))
  m <- nrow(L)
  validateMatrices(L, O)
  O <- 1 * ((O + t(O)) != 0)

  # 1 & 2 = source & target
  # 2 + {1,...,m} = L(i) for i=1,...,m
  # 2+m + {1,...,m} = R(i)-in for i=1,...,m
  # 2+2*m + {1,...,m} = R(i)-out for i=1,...,m

  Cap.matrix.init <- matrix(0, 2 + 3*m, 2 + 3*m)
  for (i in 1:m) {
    # edge from L(i) to R(i)-in, and to R(j)-in for all siblings j of i
    Cap.matrix.init[2 + i, 2 + m + c(i, which(O[i,] == 1))] <- 1
    # edge from R(i)-in to R(i)-out
    Cap.matrix.init[2 + m + i, 2 + 2*m + i] <- 1
    # edge from R(i)-out to R(j)-in for all directed edges i->j
    Cap.matrix.init[2 + 2*m + i, 2 + m + which(L[i,] == 1)] <- 1
  }

  # when testing if a set A satisfies the HTC with respect to a node i,
  #    need to add (1) edge from source to L(j) for all j in A
  #            and (2) edge from R(j)-out to target for all parents j of i

  Dependence.matrix <- O + diag(m)
  for (i in 1:m) {
    Dependence.matrix <- (Dependence.matrix + Dependence.matrix %*% L > 0)
  }

  Solved.nodes <- rep(0,m)
  Solved.nodes[which(colSums(L) == 0)] <- 1 # nodes with no parents
  change <- 1
  count <- 1

  while (change == 1) {
    change <- 0

    for (i in which(Solved.nodes == 0)) {
      A <- setdiff(c(which(Solved.nodes > 0),
                     which(Dependence.matrix[i, ] == 0)),
                   c(i, which(O[i,] == 1)))

      Cap.matrix <- Cap.matrix.init

      Cap.matrix[1, 2 + A] <- 1
      Cap.matrix[2 + 2*m + which(L[,i] == 1), 2] <- 1

      flow <-
        graph.maxflow(graph.adjacency(Cap.matrix), source = 1, target = 2)$value

      if (flow == sum(L[,i])) {
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

  # 1 & 2 = source & target
  # 2 + {1,...,N} = L{i_n,j_n} for the n-th nonsibling pair, n=1,...,N
  # (2+N) + 2*m^2 = 2 copies of R_i(j) -- in copy & out copy
  #         where (2+N) + (i-1)*m + j    = R_i(j) in copy
  #         & (2+N+m^2) + (i-1)*m + j    = R_i(j) out copy

  nonsibs <- NULL
  N <- 0
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      if (O[i,j] == 0) {
        N <- N + 1
        nonsibs <- rbind(nonsibs, c(i, j))
      }
    }
  }

  Cap.matrix <- matrix(0, 2*m^2 + N + 2, 2*m^2 + N + 2)
  if (N != 0) {
    Cap.matrix[1, 2 + (1:N)] <- 1 # edges from source to L{i,j} for each
    # nonsibling pair
    for (n in 1:N) {  #{i,j} = nonsibs[n,1:2]
      # edge from L{i,j} to R_i(j), and to R_i(k)-in for all siblings k of
      # node j
      Cap.matrix[2 + n, 2 + N +
                   (nonsibs[n,1] - 1)*m +
                   c(nonsibs[n,2], which(O[nonsibs[n,2], ] == 1))] <- 1
      # edge from L{i,j} to R_j(i), and to R_j(i)-in for all siblings k of
      # node i
      Cap.matrix[2 + n, 2 + N +
                   (nonsibs[n,2] - 1)*m +
                   c(nonsibs[n,1], which(O[nonsibs[n,1], ] == 1))] <- 1
    }
  }
  for (i in 1:m) {
    # edge from R_i(j)-out to target when j is a parent of i
    Cap.matrix[2 + N + m^2 + (i - 1)*m + which(L[, i] == 1), 2] <- 1
    for (j in 1:m) {
      # edge from R_i(j)-in to R_i(j)-out
      Cap.matrix[2 + N + (i - 1)*m + j,
                 2 + N + m^2 + (i - 1)*m + j] <- 1
      # edge from R_i(j)-out to R_i(k)-in where j->k is a directed edge
      Cap.matrix[2 + N + m^2 + (i - 1)*m + j,
                 2 + N + (i - 1)*m + which(L[j, ] == 1)] <- 1
    }
  }

  HTC.nonID <-
    graph.maxflow(igraph::graph.adjacency(Cap.matrix), source = 1, target = 2)$value
  return(HTC.nonID < sum(L))
}
