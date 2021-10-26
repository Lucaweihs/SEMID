library(SEMID)
context("Testing that latent-factor half-trek criterion function for generic identifiability works properly.")

source("graphExamples.R")

test_that("lfhtcID returns correct value for known examples.", {
  for (i in 1:length(digraphExamples)) {
    res = lfhtcID(digraphExamples[[i]]$graph)
    lfhtcid = (sum(sapply(res$unsolvedParents, length)) == 0)
    # TRUE = LF-HTC-identifiable, FALSE = not LF-HTC-identifiable
    expect_equal(lfhtcid, digraphExamples[[i]]$lfhtcid)
  }
})