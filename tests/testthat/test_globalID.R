library(SEMID)
context("Testing that function for global identifiability works properly.")

source("graphExamples.R")

test_that("globalID returns correct value for known examples.", {
    for (i in 1:length(graphExamples)) {
        graphExample <- graphExamples[[i]]
        L <- graphExample$L
        O <- graphExample$O
        globalId <- graphExample$globalId
        m <- nrow(L)

        result <- globalID(MixedGraph(L,O))
        if (globalId == 1) {
            expect_true(result)
        } else {
            expect_false(result)
        }
    }
})
