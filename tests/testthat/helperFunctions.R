###
# A collection of helper functions for testing, help create random examples.
###

rConnectedAdjMatrix = function(n,p) {
  weights = runif(n*(n-1)/2)
  g = minimum.spanning.tree(graph.full(n), weights=weights)
  adjMatrix = as.matrix(get.adjacency(g))
  adjMatrix = (upper.tri(matrix(0, n, n)) &
                 matrix(sample(c(T, F), n^2,
                               replace = T, prob = c(p, 1 - p)), ncol = n)) | adjMatrix
  adjMatrix = 1*(adjMatrix | t(adjMatrix))
  return(adjMatrix)
}

rAcyclicDirectedAdjMatrix = function(n,p) {
  return(1*(upper.tri(matrix(0, n, n)) & matrix(sample(c(T, F), n^2, replace = T,
                                                       prob = c(p, 1 - p)), ncol = n)))
}

rDirectedAdjMatrix = function(n,p) {
  return((1 - diag(n)) * matrix(sample(c(T, F), n^2, replace = T,
                                       prob = c(p, 1 - p)), ncol = n))
}

rUndirectedAdjMat = function(n, p) {
  mat = matrix(runif(n^2), ncol=n) < p
  mat = mat * upper.tri(mat)
  return(mat + t(mat))
}

getAdjMat = function(g) { as.matrix(get.adjacency(g)) }