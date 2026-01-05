#' Check Extended M-Identifiability
#'
#' @param lambda adjacency matrix with number of cols = number of latent nodes,
#'        number of rows = number of observed nodes
#' @param maxCard (optional) maximum size of set W in Corollary 4.4
#'
#' @return a list consisting of a Boolean, whether the graph is
#'         sign-identifiable and if yes, a list consisting of the sets
#' @export
#'
#' @references
#' Sturma, N., Kranzlmüller, M., Portakal, I., and Drton, M. (2025) Matching
#' Criterion for Identifiability in Sparse Factor Analysis.
#' arXiv:2502.02986l
extmID <- function(lambda, maxCard = length(observedNodes)){

  input <- transformLambda(lambda)
  adjMatrix <- input[[1]]
  latentNodes <- input[[2]]
  observedNodes <- input[[3]]
  result <- list()
  result$latentNodes <- latentNodes
  result$observedNodes <- observedNodes
  result$call <- match.call()
  class(result) <- "extmIDresult"

  S <- {}

  # all nodes without children are identifiable and added to S
  for(latent in latentNodes){
    hasChildren <- FALSE
    for(observed in observedNodes){
      if(adjMatrix[latent, observed]==1){
        hasChildren = TRUE
        break
      }
    }

    if(!hasChildren){
      S <- union(S, latent)
    }
  }

  adjMatrix[S,] = 0

  observedWithoutParents <- which(colSums(adjMatrix[])==0)
  observedNodes <- setdiff(observedNodes,observedWithoutParents)

  # first, try to find nodes which are identifiable using the local BB criterion,
  # second, try to find nodes which are identifiable using the matching criterion

  latentNodes <- setdiff(latentNodes,S)
  notIdentifiedNodes <- latentNodes
  flowGraphAdjMatrix <- flowGraphMatrix(adjMatrix, latentNodes, observedNodes)
  tupleList <- list()
  while (! length(latentNodes)==0) {
    foundIdentifiableNode <- FALSE

    tupleForSolvedNodes <- checkLocalBBCriterion(adjMatrix, latentNodes, observedNodes)
    if(tupleForSolvedNodes$found){
      foundIdentifiableNode <- TRUE

      tuple <- list(list("criterion"="localBB", "S"=S, "new nodes in S"=tupleForSolvedNodes$newNodesInS, "U"=tupleForSolvedNodes$U))
      tupleList <- c(tupleList,tuple)
      latentNodes <- setdiff(latentNodes,tupleForSolvedNodes$newNodesInS)
      S <- union(S, tupleForSolvedNodes$newNodesInS)

      adjMatrix[tupleForSolvedNodes$newNodesInS,] = 0

      observedWithoutParents <- which(colSums(adjMatrix[])==0)
      observedNodes <- setdiff(observedNodes,observedWithoutParents)
    }

    if(!foundIdentifiableNode){
      for(h in notIdentifiedNodes){
        tupleForNode <- checkMatchingCriterion(flowGraphAdjMatrix, adjMatrix, h, latentNodes, observedNodes, maxCard)
        if(tupleForNode$found){
          foundIdentifiableNode <- TRUE
          latentNodes <- setdiff(latentNodes,h)
          notIdentifiedNodes <- notIdentifiedNodes[! notIdentifiedNodes %in% c(h)]
          tuple <- list(list("criterion"="matching", "h"=tupleForNode$h, "S"=S, "v"=tupleForNode$v, "W"=tupleForNode$W, "U"=tupleForNode$U))
          tupleList <- c(tupleList,tuple)
          S <- union(S,h)

          adjMatrix[h,] = 0

          observedWithoutParents <- which(colSums(adjMatrix[])==0)
          observedNodes <- setdiff(observedNodes,observedWithoutParents)
        }
      }
    }
    if(!foundIdentifiableNode){
      result$identifiable <- FALSE
      result$tupleList <- tupleList
    }
  }
  result$identifiable <- TRUE
  result$tupleList <- tupleList
  return(result)
}

#' Prints a extmIDresult object
#'
#' Prints a extmIDresult object as returned by
#' \code{\link{extmID}}. Invisibly returns its argument via
#' \code{\link{invisible}(x)} as most print functions do.
#'
#' @export
#'
#' @param x the extmIDresult object
#' @param ... optional parameters, currently unused.
print.extmIDresult <- function(x, ...) {
  cat("Call: ")
  print(x$call)

  cat("\nFactor Analysis Graph Info:\n")
  cat("latent nodes: ", x$latentNodes, "\n")
  cat("observed nodes: ", x$observedNodes, "\n\n")

  cat("Generic Sign-Identifiability Summary:\n")
  cat(sprintf("extM-identifiable:    %s\n", x$identifiable))
  cat("Tuple list:\n")

  for (i in seq_along(x$tupleList)) {
    cat("  Tuple", i, ":\n")
    for (nm in names(x$tupleList[[i]])) {
      val <- x$tupleList[[i]][[nm]]

      if (is.list(val)) {
        cat("    ", nm, ": {\n", sep = "")
        for (subnm in names(val)) {
          subval <- val[[subnm]]
          # recursively check if subval is list
          if (is.list(subval)) {
            cat("      ", subnm, ": { ... }\n", sep = "")
          } else {
            cat("      ", subnm, ": ", paste(subval, collapse = ", "), "\n", sep = "")
          }
        }
        cat("    }\n")
      } else {
        cat("    ", nm, ": ", paste(val, collapse = ", "), "\n", sep = "")
      }
    }
  }

  invisible(x)
}