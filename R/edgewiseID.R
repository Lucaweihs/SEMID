createHtrGraph = function(L, O) {
  m = nrow(L)
  htrGraphAdjMatrix = matrix(0, 2 * m, 2 * m)
  for (i in 1:m) {
    # Left i points to right i
    htrGraphAdjMatrix[i, i + m] = 1
    for (j in 1:m) {
      if (O[i,j] == 1) {
        # If bidirected edge from i to j then
        # left i should point to right j
        htrGraphAdjMatrix[i, j + m] = 1
      }

      if (L[i,j] == 1) {
        # If directed edge from i to j then
        # right i should point to right j
        htrGraphAdjMatrix[i + m, j + m] = 1
      }
    }
  }
  return(graph.adjacency(htrGraphAdjMatrix, mode = "directed"))
}

flowGraphWithVertexCaps = function(L, O, vertexCap = 1) {
  m = nrow(L)
  adjMat = matrix(0, 3 * m + 2, 3 * m + 2)
  for (i in 1:m) {
    # Left-i points to right-i-in, right-i-in points to right-i-out
    adjMat[i, i + m] = 1 # Left-i -> right-i-in
    adjMat[i + m, i + 2 * m] = 1 # Right-i-in -> right-i-out

    for (j in 1:m) {
      if (i != j) {
        if (O[i,j] == 1) {
          # If bidirected edge from i to j then
          # left-i should point to right-j-in
          adjMat[i, j + m] = 1
        }

        if (L[i,j] == 1) {
          # If directed edge from i to j then
          # right-i-out should point to right-j-in
          adjMat[i + 2 * m, j + m] = 1
        }
      }
    }
  }
  adjMat[3 * m + 1, 1:m] = 1 # Source points to all lefts
  adjMat[(2 * m + 1):(3 * m), 3 * m + 2] = 1 # All right-outs point to target
  g = graph.adjacency(adjMat, mode = "directed")
  E(g)$capacity = vertexCap
  return(list(g = g, source = 3 * m + 1, target = 3 * m + 2))
}

createIdentifierBaseCase <- function(L) {
  return(function(Sigma) {
    Lambda <- matrix(NA, nrow(L), nrow(L))
    Lambda[L == 0] = 0
    return(Lambda)
  })
}

createLinearIdentifierFunc <- function(idFunc, sources, targets, node, sourceHtrParents) {
  # These assignments may seem redundent but they are necessary as they
  # assign these variables to the local environment of the function call.
  # This allows them to persist and still be usable by the returned function.
  idFunc <- idFunc
  sources <- sources
  targets <- targets
  node <- node
  sourceHtrParents <- sourceHtrParents
  return(
    function(Sigma) {
      m <- nrow(Sigma)
      Lambda <- idFunc(Sigma)

      SigmaMinus = Sigma
      for (sourceInd in 1:length(sources)) {
        source = sources[sourceInd]
        parents = sourceHtrParents[[sourceInd]]
        if (length(parents) != 0) {
          SigmaMinus[source,] <- Sigma[source, , drop = F] -
            t(Lambda[parents, source, drop = F]) %*% Sigma[parents,, drop = F]
        }
      }

      if (abs(det(SigmaMinus[sources, targets, drop = F])) < 10^-10) {
        stop("In identification, found near-singular system. Is the input matrix generic?")
      }

      Lambda[targets, node] <-
        solve(SigmaMinus[sources, targets, drop = F],
              SigmaMinus[sources, node, drop = F])

      return(Lambda)
    }
  )
}

edgewiseID = function(L, O) {
  validateMatrices(L, O)
  O = 1 * ((O + t(O)) != 0)
  diag(O) = 0
  m = nrow(L)

  htrGraph = createHtrGraph(L, O)
  dirGraph = graph.adjacency(L, mode = "directed")

  htr = neighborhood(htrGraph, order = 2 * m, nodes = 1:m, mode = "out", mindist = 1)
  unsolvedParents = neighborhood(dirGraph, order = 1, nodes = 1:m, mode = "in", mindist = 1)
  parents = unsolvedParents
  for (i in 1:m) {
    htr[[i]] = as.numeric(htr[[i]]) - m
    unsolvedParents[[i]] = as.numeric(unsolvedParents[[i]])
    parents[[i]] = unsolvedParents[[i]]
  }

  flowGraphList = flowGraphWithVertexCaps(L, O, 1)
  flowGraph = flowGraphList$g
  s = flowGraphList$source
  t = flowGraphList$target
  sOutEdgeIds = get.edge.ids(flowGraph, as.numeric(rbind(rep(s, m), 1:m)))
  tInEdgeIds = get.edge.ids(flowGraph, as.numeric(rbind((2 * m + 1):(3 * m), rep(t, m))))
  identifier = createIdentifierBaseCase(L)
  changeFlag = T
  solvedParents = rep(list(numeric(0)), m)
  while (changeFlag) {
    changeFlag = F
    for (i in 1:m) {
      unsolved = unsolvedParents[[i]]
      if (length(unsolved) != 0) {
        allowedNodes = logical(m)
        for (j in 1:m) {
          if (i != j &&
              O[i,j] != 1 &&
              length(intersect(htr[[j]], unsolved)) != 0 &&
              length(intersect(htr[[i]], unsolvedParents[[j]])) == 0) {
            allowedNodes[j] = TRUE
          }
        }
        if (all(!allowedNodes)) {
          next
        }

        subsetFound = F
        for (k in length(unsolved):1) {
          if (length(unsolved) == 1) {
            # Silly conditional required because the combn function
            # interprets inputs of size 1 as a length rather than an as
            # a vector from which to get subsets
            subsets = list(unsolved)
          } else {
            subsets = combn(unsolved, k, simplify = F)
          }
          for (subset in subsets) {
            E(flowGraph)$capacity[sOutEdgeIds] = 0
            for (a in which(allowedNodes)) {
              if (all(intersect(htr[[a]], unsolved) %in% subset)) {
                E(flowGraph)$capacity[sOutEdgeIds[a]] = 1
              }
            }
            E(flowGraph)$capacity[tInEdgeIds] = 0
            E(flowGraph)$capacity[tInEdgeIds[subset]] = 1
            flowResult = graph.maxflow(flowGraph, s, t)
            if (flowResult$value == k) {
              changeFlag = T
              subsetFound = T
              sources = which(flowResult$flow[sOutEdgeIds] == 1)
              identifier = createLinearIdentifierFunc(identifier,
                                                      sources,
                                                      subset,
                                                      i,
                                                      lapply(sources, function(x) {
                                                        intersect(htr[[i]], parents[[x]])
                                                      }))
              solvedParents[[i]] = sort(c(subset, solvedParents[[i]]))
              unsolvedParents[[i]] = setdiff(unsolved, subset)
              break
            }
          }
          if (subsetFound) {
            break
          }
        }
      }
    }
  }

  return(list(solvedParents = solvedParents,
              unsolvedParents = unsolvedParents,
              identifier = identifier))
}



