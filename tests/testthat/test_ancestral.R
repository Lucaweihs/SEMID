library(SEMID)
context("Testing components related to generic identifiability by ancestor decomposition.")

test_that("ancestors function works as expected.", {
  ## Output should be empty
  expect_equal(ancestors(graph.empty(), c()), numeric(0))
  expect_equal(ancestors(graph.empty(), c(1,2,3)), numeric(0))
  expect_equal(ancestors(graph.star(5, mode="in"), c()), numeric(0))

  # Graph with no edges
  expect_equal(ancestors(graph.empty(5), 1), 1)
  expect_equal(ancestors(graph.empty(5), 2), 2)
  expect_equal(ancestors(graph.empty(5), c(2,4,5)), c(2,4,5))

  # Graphs with edges
  expect_equal(ancestors(graph.star(5, mode="in"), 1), 1:5)
  expect_equal(ancestors(graph.star(5, mode="in"), 2), 2)
  expect_equal(ancestors(graph.star(5, mode="out"), 1), 1)
  expect_equal(ancestors(graph.full(5, directed=T), 1), 1:5)
  g = graph.edgelist(matrix(c(1,2,2,3,3,5,5,7,2,4,4,6,6,7,4,5), ncol=2, byrow=T))
  expect_equal(ancestors(g, 1), 1)
  expect_equal(ancestors(g, 2), c(1,2))
  expect_equal(ancestors(g, 4), c(1,2,4))
  expect_equal(ancestors(g, 5), c(1,2,3,4,5))
  expect_equal(ancestors(g, 6), c(1,2,4,6))
  expect_equal(ancestors(g, c(6, 3)), c(1,2,3,4,6))
  expect_equal(ancestors(g, 7), 1:7)
})

test_that("parents function works as expected.", {
  ## Output should be empty
  expect_equal(parents(graph.empty(), c()), numeric(0))
  expect_equal(parents(graph.empty(), c(1,2,3)), numeric(0))
  expect_equal(parents(graph.star(5, mode="in"), c()), numeric(0))

  # Graph with no edges
  expect_equal(parents(graph.empty(5), 1), 1)
  expect_equal(parents(graph.empty(5), 2), 2)
  expect_equal(parents(graph.empty(5), c(2,4,5)), c(2,4,5))

  # Graphs with edges
  expect_equal(parents(graph.star(5, mode="in"), 1), 1:5)
  expect_equal(parents(graph.star(5, mode="in"), 2), 2)
  expect_equal(parents(graph.star(5, mode="out"), 1), 1)
  expect_equal(parents(graph.full(5, directed=T), 1), 1:5)
  g = graph.edgelist(matrix(c(1,2,2,3,3,5,5,7,2,4,4,6,6,7,4,5,1,7), ncol=2, byrow=T))
  expect_equal(parents(g, 1), 1)
  expect_equal(parents(g, 2), c(1,2))
  expect_equal(parents(g, 4), c(2,4))
  expect_equal(parents(g, 5), c(3,4,5))
  expect_equal(parents(g, 6), c(4,6))
  expect_equal(parents(g, c(6, 3)), c(2,3,4,6))
  expect_equal(parents(g, 7), c(1,5,6,7))
})

test_that("siblings function works as expected.", {
  ## Output should be empty
  expect_equal(siblings(graph.empty(), c()), numeric(0))
  expect_equal(siblings(graph.empty(), c(1,2,3)), numeric(0))
  expect_equal(siblings(graph.star(5, mode="in"), c()), numeric(0))

  # Graph with no edges
  expect_equal(siblings(graph.empty(5), 1), 1)
  expect_equal(siblings(graph.empty(5), 2), 2)
  expect_equal(siblings(graph.empty(5), c(2,4,5)), c(2,4,5))

  # Graphs with edges
  expect_equal(siblings(graph.star(5, mode="in"), 1), 1:5)
  expect_equal(siblings(graph.star(5, mode="in"), 2), c(1,2))
  expect_equal(siblings(graph.star(5, mode="out"), 1), 1:5)
  expect_equal(siblings(graph.star(5, mode = "undirected"), 1), 1:5)
  expect_equal(siblings(graph.full(5, directed=F), 1), 1:5)
  g = graph.edgelist(matrix(c(1,2,2,3,3,5,5,7,2,4,4,6,6,7,4,5,1,7), ncol=2, byrow=T), directed = F)
  expect_equal(siblings(g, 1), c(1,2,7))
  expect_equal(siblings(g, 2), c(1,2,3,4))
  expect_equal(siblings(g, c(6, 3)), c(2,3,4,5,6,7))
})

test_that("getMixedCompForNode function works as expected.", {
  ## Graph with single node
  dG = graph.empty(1, directed=T)
  bG = graph.empty(1, directed=F)
  expect_error(getMixedCompForNode(dG, bG, 1, 1)) # Since vertices are unnamed
  V(dG)$names = 1
  V(bG)$names = 1

  compList = getMixedCompForNode(dG, bG, 1, 1)
  expect_equal(compList, list(biNodes=1, inNodes=numeric(0)))

  compList = getMixedCompForNode(dG, bG, 1, 2)
  expect_equal(compList, list(biNodes=numeric(0), inNodes=numeric(0)))

  compList = getMixedCompForNode(dG, bG, 2, 1)
  expect_equal(compList, list(biNodes=numeric(0), inNodes=numeric(0)))

  # Graph with 2 nodes
  dG = graph.empty(2, directed=T)
  bG = graph.empty(2, directed=F)
  V(dG)$names = 1:2
  V(bG)$names = 2:1
  expect_error(getMixedCompForNode(dG, bG, 1, 1)) # Since vertices are wrongly named

  V(dG)$names = 1:2
  V(bG)$names = 1:2

  compList = getMixedCompForNode(dG, bG, 1, 1)
  expect_equal(compList, list(biNodes=1, inNodes=numeric(0)))

  compList = getMixedCompForNode(dG, bG, 1, 2)
  expect_equal(compList, list(biNodes=numeric(0), inNodes=numeric(0)))

  compList = getMixedCompForNode(dG, bG, 2, 1)
  expect_equal(compList, list(biNodes=numeric(0), inNodes=numeric(0)))

  # The sink marginalization example
  dG = graph.edgelist(matrix(c(1,2, 1,3, 1,6, 2,3, 2,4, 2,5, 2,6, 3,4, 4,5),
                             ncol=2, byrow=T))
  bG = graph.edgelist(matrix(c(1,6, 1,4, 2,3, 2,5, 2,6),
                             ncol=2, byrow=T), directed=F)
  V(dG)$names = 1:6
  V(bG)$names = 1:6
  compList = getMixedCompForNode(dG, bG, 1:6, 1)
  expect_equal(compList, list(biNodes=1:6, inNodes=numeric(0)))
  compList = getMixedCompForNode(dG, bG, c(4,5), 5)
  expect_equal(compList, list(biNodes=5, inNodes=4))
  compList = getMixedCompForNode(dG, bG, 1:5, 5)
  expect_equal(compList, list(biNodes=c(2,3,5), inNodes=c(1,4)))
  compList = getMixedCompForNode(dG, bG, c(2,3,4), 4)
  expect_equal(compList, list(biNodes=c(4), inNodes=c(2,3)))
  compList = getMixedCompForNode(dG, bG, rev(c(2,4,5,6)), 2)
  expect_equal(compList, list(biNodes=c(2,5,6), inNodes=c(4)))

  bG = delete.edges(bG, get.edge.ids(bG, c(1,6)))
  compList = getMixedCompForNode(dG, bG, 1:6, 1)
  expect_equal(compList, list(biNodes=c(1,4), inNodes=c(2,3)))
  compList = getMixedCompForNode(dG, bG, 1:6, 3)
  expect_equal(compList, list(biNodes=c(2,3,5,6), inNodes=c(1,4)))
})

test_that("getMaxFlow function works as expected.", {
  ## Graph with single node
  dG = graph.empty(1, directed=T)
  bG = graph.empty(1, directed=F)

  getMaxFlow(L, O, allowedNodes, biNodes, inNodes, node)
  expect_error(getMixedCompForNode(dG, bG, 1, 1)) # Since vertices are unnamed
  V(dG)$names = 1
  V(bG)$names = 1

  getMaxFlow(L, O, allowedNodes, biNodes, inNodes, node)

  compList = getMixedCompForNode(dG, bG, 1, 1)
  expect_equal(compList, list(biNodes=1, inNodes=numeric(0)))

  compList = getMixedCompForNode(dG, bG, 1, 2)
  expect_equal(compList, list(biNodes=numeric(0), inNodes=numeric(0)))

  compList = getMixedCompForNode(dG, bG, 2, 1)
  expect_equal(compList, list(biNodes=numeric(0), inNodes=numeric(0)))

  # Graph with 2 nodes
  dG = graph.empty(2, directed=T)
  bG = graph.empty(2, directed=F)
  V(dG)$names = 1:2
  V(bG)$names = 2:1
  expect_error(getMixedCompForNode(dG, bG, 1, 1)) # Since vertices are wrongly named

  V(dG)$names = 1:2
  V(bG)$names = 1:2

  compList = getMixedCompForNode(dG, bG, 1, 1)
  expect_equal(compList, list(biNodes=1, inNodes=numeric(0)))

  compList = getMixedCompForNode(dG, bG, 1, 2)
  expect_equal(compList, list(biNodes=numeric(0), inNodes=numeric(0)))

  compList = getMixedCompForNode(dG, bG, 2, 1)
  expect_equal(compList, list(biNodes=numeric(0), inNodes=numeric(0)))

  # The sink marginalization example
  dG = graph.edgelist(matrix(c(1,2, 1,3, 1,6, 2,3, 2,4, 2,5, 2,6, 3,4, 4,5),
                             ncol=2, byrow=T))
  bG = graph.edgelist(matrix(c(1,6, 1,4, 2,3, 2,5, 2,6),
                             ncol=2, byrow=T), directed=F)
  V(dG)$names = 1:6
  V(bG)$names = 1:6
  compList = getMixedCompForNode(dG, bG, 1:6, 1)
  expect_equal(compList, list(biNodes=1:6, inNodes=numeric(0)))
  compList = getMixedCompForNode(dG, bG, c(4,5), 5)
  expect_equal(compList, list(biNodes=5, inNodes=4))
  compList = getMixedCompForNode(dG, bG, 1:5, 5)
  expect_equal(compList, list(biNodes=c(2,3,5), inNodes=c(1,4)))
  compList = getMixedCompForNode(dG, bG, c(2,3,4), 4)
  expect_equal(compList, list(biNodes=c(4), inNodes=c(2,3)))
  compList = getMixedCompForNode(dG, bG, rev(c(2,4,5,6)), 2)
  expect_equal(compList, list(biNodes=c(2,5,6), inNodes=c(4)))

  bG = delete.edges(bG, get.edge.ids(bG, c(1,6)))
  compList = getMixedCompForNode(dG, bG, 1:6, 1)
  expect_equal(compList, list(biNodes=c(1,4), inNodes=c(2,3)))
  compList = getMixedCompForNode(dG, bG, 1:6, 3)
  expect_equal(compList, list(biNodes=c(2,3,5,6), inNodes=c(1,4)))
})