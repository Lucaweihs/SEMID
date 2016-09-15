library(SEMID)
context("Components related to edgewise identification.")

solvedParentsToIndexMat = function(solvedParents) {
  nSolved = sum(sapply(solvedParents, length))
  mat = matrix(0, ncol = 2, nrow = nSolved)
  k = 1
  for (i in 1:length(solvedParents)) {
    solved = solvedParents[[i]]
    if (length(solved) == 0) {
      next
    }
    mat[k:(k + length(solved) - 1),] = cbind(solved, rep(i, length(solved)))
    k = k + 1
  }
  return(mat)
}

test_that("edgewiseID function works as expected.", {
  # Random test
  set.seed(2323)
  ps = c(.3, .4)
  sims = 10
  ns = c(4, 6)
  for (p in ps) {
    for (n in ns) {
      for (i in 1:sims) {
        L = rDirectedAdjMatrix(n, p)
        O = rUndirectedAdjMat(n, p)

        gid = graphID.htcID(L, O)
        gin = graphID.nonHtcID(L, O)
        eid = edgewiseID(L, O)

        extraSolved = 0
        for (j in 1:nrow(L)) {
          if (!(j %in% gid)) {
            extraSolved = extraSolved + length(eid$solvedParents[[j]])
          }
        }
        print(extraSolved)

        edgeSolvedNodes = which(sapply(eid$unsolvedParents,
                                function(x) { length(x) == 0 }))

        expect_true(length(gid) <= length(edgeSolvedNodes))
        expect_true(all(gid %in% edgeSolvedNodes))
        expect_true((length(edgeSolvedNodes) != n) || !gin)

        L1 = L * matrix(runif(n^2, .1, 1), n)
        O1 = (diag(n) + O) * matrix(runif(n^2, .1, 1), n)
        O1 = O1 + t(O1)
        S1 = t(solve(diag(n) - L1)) %*% O1 %*% solve(diag(n) - L1)

        toCheck = solvedParentsToIndexMat(eid$solvedParents)
        expect_true(all(abs(eid$identifier(S1)[toCheck] - L1[toCheck]) < 10^-6))
      }
    }
  }
})
