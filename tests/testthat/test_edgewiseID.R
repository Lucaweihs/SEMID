library(SEMID)
context("Components related to edgewise identification.")

source("helperFunctions.R")

edgewiseID = function(L, O, maxTrekSepSubsetSize = 3) {
  return(generalGenericID(L, O, list(edgewiseIdentifyStep)))
}

test_that("Edgewise identification does not identify edges erroneously.", {
  # Random test
  set.seed(2323)
  ps = c(.3)
  sims = 30
  ns = c(5, 6)
  for (p in ps) {
    for (n in ns) {
      for (i in 1:sims) {
        L = rDirectedAdjMatrix(n, p)
        O = rUndirectedAdjMat(n, p)
        gid = htcID(L, O)
        eid = edgewiseID(L, O)

        expect_true(all(sapply(1:n, FUN = function(x) {
          all(gid$solvedParents[[x]] %in% eid$solvedParents[[x]])
        })))

        randomIdentificationTest(eid$identifier, L, O, eid$solvedParents)
      }
    }
  }
})
