library(SEMID)
context("LatentDigraph object.")

source("helperFunctions.R")

test_that("Empty graph works properly", {
  g = LatentDigraph(matrix(0, 0, 0))

  expect_equal(g$numObserved(), 0)
  expect_equal(g$numLatents(), 0)
  expect_equal(g$numNodes(), 0)
  expect_equal(g$L(), matrix(0, 0, 0))
  expect_equal(g$latentNodes(), integer(0))
  expect_equal(g$observedNodes(), integer(0))

  expect_equal(g$parents(c()), integer(0))
  expect_error(g$parents(1))

  expect_equal(g$ancestors(c()), integer(0))
  expect_error(g$ancestors(1))

  expect_equal(g$descendants(c()), integer(0))
  expect_error(g$descendants(1))

  expect_error(g$trFrom(1, 1))
  expect_error(g$getTrekSystem(1, 1))
  expect_error(g$inducedSubgraph(1))

  expect_equal(g$observedParents(c()), integer(0))
  expect_equal(g$getMixedGraph()$L(), matrix(0, 0, 0))
  expect_equal(g$getMixedGraph()$O(), matrix(0, 0, 0))
})

test_that("Single node graph works properly", {
  L = matrix(0, 1, 1)

  expect_error(LatentDigraph(L + 1))
  expect_error(LatentDigraph(L, latentNodes = c(1)))
  expect_error(LatentDigraph(L, latentNodes = c(2)))
  expect_error(LatentDigraph(L, observedNodes = c(1,2)))

  g = LatentDigraph(L)

  expect_equal(g$numObserved(), 1)
  expect_equal(g$numLatents(), 0)
  expect_equal(g$L(), L)
  expect_equal(g$latentNodes(), integer(0))
  expect_equal(g$observedNodes(), c(1))

  expect_equal(g$parents(c()), integer(0))
  expect_equal(g$parents(1), integer(0))
  expect_error(g$parents(0))

  expect_equal(g$ancestors(c()), integer(0))
  expect_equal(g$ancestors(1), c(1))
  expect_error(g$ancestors(2))

  expect_equal(g$descendants(c()), integer(0))
  expect_equal(g$descendants(1), c(1))
  expect_error(g$descendants(2))

  expect_equal(g$trFrom(1), 1)
  expect_equal(g$trFrom(1, avoidLeftNodes = 1), integer(0))
  expect_equal(g$trFrom(1, avoidRightNodes = 1), 1)
  expect_error(g$trFrom(0))

  expect_equal(g$getTrekSystem(1, 1), list(systemExists = TRUE, activeFrom = 1))

  expect_equal(g$inducedSubgraph(1)$L(), L)

  expect_equal(g$observedParents(c()), integer(0))
  expect_equal(g$observedParents(1), integer(0))
  expect_equal(g$getMixedGraph()$L(), matrix(0, 1, 1))
  expect_equal(g$getMixedGraph()$O(), matrix(0, 1, 1))
})

test_that("Two nodes, no latents, graph works properly", {
  L = matrix(c(0,1,0,0), ncol = 2, byrow = T)

  expect_error(LatentDigraph(L + 1))
  expect_error(LatentDigraph(L, latentNodes = 1))
  expect_error(LatentDigraph(L, observedNodes = c(1,2,3)))

  g = LatentDigraph(L)

  expect_equal(g$numObserved(), 2)
  expect_equal(g$numLatents(), 0)
  expect_equal(g$L(), L)
  expect_equal(g$latentNodes(), integer(0))
  expect_equal(g$observedNodes(), c(1,2))

  expect_equal(g$parents(c()), integer(0))
  expect_equal(g$parents(1), integer(0))
  expect_equal(g$parents(2), c(1))
  expect_error(g$parents(0))

  expect_equal(g$ancestors(c()), integer(0))
  expect_equal(g$ancestors(1), c(1))
  expect_equal(sort(g$ancestors(2)), c(1,2))
  expect_error(g$ancestors(0))

  expect_equal(g$descendants(c()), integer(0))
  expect_equal(sort(g$descendants(1)), c(1,2))
  expect_equal(g$descendants(2), c(2))
  expect_error(g$descendants(0))

  expect_equal(sort(g$trFrom(1)), c(1,2))
  expect_equal(sort(g$trFrom(2)), c(1,2))
  expect_equal(sort(g$trFrom(2, avoidLeftNodes = 1)), c(2))
  expect_equal(sort(g$trFrom(1, avoidRightNodes = 2)), c(1))
  expect_error(g$trFrom(1, avoidNodesOnRight = 2))
  expect_error(g$trFrom(0))

  expect_equal(g$getTrekSystem(1, 1), list(systemExists = TRUE, activeFrom = 1))
  expect_equal(g$getTrekSystem(2, 1), list(systemExists = TRUE, activeFrom = 2))
  expect_equal(g$getTrekSystem(c(1,2), c(2,1)), list(systemExists = TRUE, activeFrom = c(1,2)))

  expect_equal(g$inducedSubgraph(c(1,2))$L(), L)
  expect_equal(g$inducedSubgraph(c(1))$L(), matrix(0,1,1))

  expect_equal(g$observedParents(c()), integer(0))
  expect_equal(g$observedParents(1), integer(0))
  expect_equal(g$observedParents(2), c(1))
  expect_equal(g$getMixedGraph()$L(), L)
  expect_equal(g$getMixedGraph()$O(), matrix(c(0, 0, 0, 0), nrow=2))
})

test_that("Two node (different labels), one latent, graph works properly", {
  L = matrix(c(0,1,0,
               0,0,0,
               1,1,0), ncol = 3, byrow = T)
  observedNodes = c(44, 12)
  latentNodes = c(3)
  expect_error(LatentDigraph(L + 1, observedNodes = observedNodes, latentNodes = latentNodes))
  expect_error(LatentDigraph(L, observedNodes = observedNodes, latentNodes = c(latentNodes, 1)))
  expect_error(LatentDigraph(L, observedNodes = c(12, observedNodes), latentNodes = latentNodes))

  g = LatentDigraph(L, observedNodes = observedNodes, latentNodes = latentNodes)

  expect_equal(g$numObserved(), 2)
  expect_equal(g$numLatents(), 1)
  expect_equal(g$L(), L)
  expect_equal(g$latentNodes(), latentNodes)
  expect_equal(g$observedNodes(), observedNodes)

  expect_equal(g$parents(c()), integer(0))
  expect_equal(g$parents(44), c(3))
  expect_equal(sort(g$parents(12)), c(3, 44))
  expect_equal(g$parents(3), integer(0))
  expect_error(g$parents(0))

  expect_equal(g$ancestors(c()), integer(0))
  expect_equal(sort(g$ancestors(44)), c(3, 44))
  expect_equal(sort(g$ancestors(12)), c(3, 12,44))
  expect_equal(sort(g$ancestors(3)), c(3))
  expect_error(g$ancestors(0))

  expect_equal(g$descendants(c()), integer(0))
  expect_equal(sort(g$descendants(3)), c(3, 12,44))
  expect_equal(sort(g$descendants(44)), c(12,44))
  expect_equal(g$descendants(12), c(12))
  expect_error(g$descendants(0))

  expect_equal(sort(g$trFrom(44)), c(3, 12, 44))
  expect_equal(sort(g$trFrom(12)), c(3, 12, 44))
  expect_equal(sort(g$trFrom(3)), c(3, 12, 44))
  expect_equal(sort(g$trFrom(44, avoidLeftNodes = 3)), c(12, 44))
  expect_equal(sort(g$trFrom(44, avoidLeftNodes = 12, avoidRightNodes = 3)), c(3, 12, 44))
  expect_equal(sort(g$trFrom(44,  avoidLeftNodes = 3, avoidRightNodes = 12)), c(44))
  expect_error(g$trFrom(0))

  expect_equal(g$getTrekSystem(44, 44), list(systemExists = TRUE, activeFrom = 44))
  expect_equal(g$getTrekSystem(12, 44), list(systemExists = TRUE, activeFrom = 12))
  expect_equal(g$getTrekSystem(c(44,12), c(12,44)), list(systemExists = TRUE, activeFrom = c(44,12)))
  expect_equal(g$getTrekSystem(c(3,12), c(12,44)), list(systemExists = TRUE, activeFrom = c(12,3)))

  expect_equal(g$inducedSubgraph(c(44,12,3))$L(), L)
  expect_equal(g$inducedSubgraph(c(44))$L(), matrix(0,1,1))
  expect_equal(g$inducedSubgraph(c(44))$latentNodes(), integer(0))
  expect_equal(g$inducedSubgraph(c(12))$L(), matrix(0,1,1))
  expect_equal(g$inducedSubgraph(c(12))$latentNodes(), integer(0))
  expect_equal(g$inducedSubgraph(c(44,3))$L(), matrix(c(0,0,1,0), ncol = 2, byrow = T))


  expect_equal(g$observedParents(c()), integer(0))
  expect_equal(g$observedParents(44), integer(0))
  expect_equal(g$observedParents(12), c(44))
  expect_equal(g$getMixedGraph()$L(), matrix(c(0,1,0,0), ncol = 2, byrow = T))
  expect_equal(g$getMixedGraph()$O(), matrix(c(0,1,1,0), ncol = 2, byrow = T))
})


test_that("More complicated graph, 5 observed, 2 latents, works properly", {
  observedNodes = c(2, 4, 6, 8, 10)
  latentNodes = c(11, 12)
  L = matrix(
    c(0,1,0,0,0,0,0,
      0,0,1,0,0,0,0,
      0,1,0,0,0,0,0,
      0,0,0,0,1,0,1,
      0,0,0,0,0,0,0,
      0,1,1,0,1,0,0,
      1,0,0,0,1,0,0), ncol = 7, byrow = T)
  expect_error(LatentDigraph(L + 1, observedNodes = observedNodes, latentNodes = latentNodes))
  expect_error(LatentDigraph(L, observedNodes = observedNodes, latentNodes = c(1, latentNodes)))
  expect_error(LatentDigraph(L, observedNodes = c(1, observedNodes), latentNodes = latentNodes))

  g = LatentDigraph(L, observedNodes = observedNodes, latentNodes = latentNodes)

  expect_equal(g$numObserved(), 5)
  expect_equal(g$numLatents(), 2)
  expect_equal(g$L(), L)
  expect_equal(g$latentNodes(), latentNodes)
  expect_equal(g$observedNodes(), observedNodes)

  expect_equal(g$parents(c()), integer(0))
  expect_equal(g$parents(2), c(12))
  expect_equal(g$parents(4), c(2,6,11))
  expect_equal(g$parents(6), c(4,11))
  expect_equal(g$parents(8), integer(0))
  expect_equal(g$parents(10), c(8,11,12))
  expect_equal(g$parents(c(2,4,8)), c(2,6,11,12))
  expect_equal(g$parents(c(2,4,6,10)), c(2,4,6,8,11,12))
  expect_equal(g$parents(c(12)), c(8))
  expect_error(g$parents(0))

  expect_equal(g$ancestors(c()), integer(0))
  expect_equal(sort(g$ancestors(2)), c(2,8,12))
  expect_equal(sort(g$ancestors(4)), c(2,4,6,8,11,12))
  expect_equal(sort(g$ancestors(6)), c(2,4,6,8,11,12))
  expect_equal(sort(g$ancestors(10)), c(8,10,11,12))
  expect_equal(sort(g$ancestors(c(2,10))), c(2,8,10,11,12))
  expect_equal(sort(g$ancestors(c(12))), c(8,12))
  expect_error(g$ancestors(0))

  expect_equal(g$descendants(c()), integer(0))
  expect_equal(sort(g$descendants(2)), c(2,4,6))
  expect_equal(sort(g$descendants(4)), c(4,6))
  expect_equal(sort(g$descendants(c(2,4))), c(2,4,6))
  expect_equal(sort(g$descendants(c(10,6))), c(4,6,10))
  expect_equal(sort(g$descendants(c(8))), c(2,4,6,8,10,12))
  expect_error(g$descendants(0))

  expect_equal(sort(g$trFrom(2)), c(2,4,6,8,10,12))
  expect_equal(sort(g$trFrom(2, avoidLeftNodes = 12)), c(2,4,6))
  expect_equal(sort(g$trFrom(6)), c(2,4,6,8,10,11,12))
  expect_equal(sort(g$trFrom(6, avoidLeftNodes = 12, avoidRightNodes = 4)), c(2,4,6,10,11))
  expect_equal(sort(g$trFrom(6, avoidLeftNodes = 2)), c(4,6,10,11))
  expect_equal(sort(g$trFrom(c(6,8), avoidRightNodes = 12, avoidLeftNodes = c(4, 11))), c(4,6,8,10))

  expect_equal(g$getTrekSystem(c(6,10), c(10, 8)),
               list(systemExists = TRUE, activeFrom = c(6,10)))
  expect_equal(g$getTrekSystem(c(2,4), c(2, 10)),
               list(systemExists = TRUE, activeFrom = c(2,4)))
  expect_equal(g$getTrekSystem(c(2,4), c(2, 10), avoidLeftNodes = 11),
               list(systemExists = FALSE, activeFrom = c(2)))
  expect_equal(g$getTrekSystem(c(2), c(10), avoidLeftNodes = 12),
               list(systemExists = FALSE, activeFrom = integer(0)))
  expect_equal(g$getTrekSystem(c(2), 8, avoidLeftNodes = 12),
               list(systemExists = FALSE, activeFrom = integer(0)))
  expect_equal(g$getTrekSystem(c(2), 10),
               list(systemExists = TRUE, activeFrom = 2))

  a = g$inducedSubgraph(c(10,4,6,11))
  expect_equal(a$L(),
               matrix(c(0,0,0,0,
                        0,0,1,0,
                        0,1,0,0,
                        1,1,1,0), ncol = 4, byrow = T))
  expect_equal(a$observedNodes(), c(10,4,6))
  expect_equal(a$latentNodes(), c(11))

  a = g$inducedSubgraph(c(11,12,2,10,4,6))
  expect_equal(a$L(),
               matrix(c(0,0,1,0,0,0,
                        0,0,0,0,0,0,
                        0,0,0,1,0,0,
                        0,0,1,0,0,0,
                        0,1,1,1,0,0,
                        1,1,0,0,0,0), ncol = 6, byrow = T))
  expect_equal(a$latentNodes(), c(11,12))
  expect_equal(a$observedNodes(), c(2,10,4,6))

  expect_equal(g$observedParents(c()), integer(0))
  expect_equal(g$observedParents(8), integer(0))
  expect_equal(g$observedParents(11), integer(0))
  expect_equal(g$observedParents(4), c(2,6))
  expect_equal(g$observedParents(12), c(8))
  expect_equal(g$observedParents(c(4,6)), c(2,4,6))
  expect_error(g$getMixedGraph())
})


test_that("getMixedGraph is working properly on larger graph where all latent nodes are source nodes", {
  observedNodes = c(2, 4, 6, 8, 10)
  latentNodes = c(11, 12)
  L = matrix(
    c(0,1,0,0,0,0,0,
      0,0,1,0,0,0,0,
      0,1,0,0,0,0,0,
      0,0,0,0,1,0,0,
      0,0,0,0,0,0,0,
      0,1,1,0,1,0,0,
      1,0,0,0,1,0,0), ncol = 7, byrow = T)

  g = LatentDigraph(L, observedNodes = observedNodes, latentNodes = latentNodes)

  expectedL = matrix(
    c(0,1,0,0,0,
      0,0,1,0,0,
      0,1,0,0,0,
      0,0,0,0,1,
      0,0,0,0,0), ncol = 5, byrow = T)

  expectedO = matrix(
    c(0,0,0,0,1,
      0,0,1,0,1,
      0,1,0,0,1,
      0,0,0,0,0,
      1,1,1,0,0), ncol = 5, byrow = T)

  expect_equal(g$getMixedGraph()$L(), expectedL)
  expect_equal(g$getMixedGraph()$O(), expectedO)

})

