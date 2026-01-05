#' Check M-Identifiability.
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
#' Sturma, N., Kranzlmüller, M., Portakal, I., and Drton, M.  (2025) Matching
#' Criterion for Identifiability in Sparse Factor Analysis.
#' arXiv:2502.02986
mID <- function(lambda, maxCard = length(observedNodes)){

  input <- transformLambda(lambda)
  adjMatrix <- input[[1]]
  latentNodes <- input[[2]]
  observedNodes <- input[[3]]
  result <- list()
  result$latentNodes <- latentNodes
  result$observedNodes <- observedNodes
  result$call <- match.call()
  class(result) <- "mIDresult"

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

  latentNodes <- setdiff(latentNodes, S)

  notIdentifiedNodes <- setdiff(latentNodes,S)
  flowGraphAdjMatrix <- flowGraphMatrix(adjMatrix, latentNodes, observedNodes)
  tupleList <- list()
  while (!(length(latentNodes)==0)) {
    foundIdentifiableNode <- FALSE
    for(h in notIdentifiedNodes){
      tupleForNode <- checkMatchingCriterion(flowGraphAdjMatrix, adjMatrix, h, latentNodes, observedNodes, maxCard)
      if(tupleForNode$found){
        foundIdentifiableNode <- TRUE
        latentNodes <- setdiff(latentNodes,h)
        notIdentifiedNodes <- notIdentifiedNodes[! notIdentifiedNodes %in% c(h)]
        tuple <- list(list("h"=tupleForNode$h, "S"=S, "v"=tupleForNode$v, "W"=tupleForNode$W, "U"=tupleForNode$U))
        tupleList <- c(tupleList,tuple)
        S <- union(S,h)
      }
    }

    if(!foundIdentifiableNode){
      result$identifiable <- FALSE
      result$tupleList <- tupleList
      return(result)
    }
  }
  result$identifiable <- TRUE
  result$tupleList <- tupleList
  return(result)
}

#' Prints a mIDresult object
#'
#' Prints a mIDresult object as returned by
#' \code{\link{mID}}. Invisibly returns its argument via
#' \code{\link{invisible}(x)} as most print functions do.
#'
#' @export
#'
#' @param x the mIDresult object
#' @param ... optional parameters, currently unused.
print.mIDresult <- function(x, ...) {
  cat("Call: ")
  print(x$call)

  cat("\nFactor Analysis Graph Info:\n")
  cat("latent nodes: ", x$latentNodes, "\n")
  cat("observed nodes: ", x$observedNodes, "\n\n")

  cat("Generic Sign-Identifiability Summary:\n")
  cat(sprintf("M-identifiable:    %s\n", x$identifiable))
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