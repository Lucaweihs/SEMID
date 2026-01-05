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

      tuple <- list(list("S"=S, "newNodesInS"=tupleForSolvedNodes$newNodesInS, "U"=tupleForSolvedNodes$U))
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
          S <- union(S,h)
          notIdentifiedNodes <- notIdentifiedNodes[! notIdentifiedNodes %in% c(h)]
          tuple <- list(list("h"=tupleForNode$h, "S"=S, "v"=tupleForNode$v, "W"=tupleForNode$W, "U"=tupleForNode$U))
          tupleList <- c(tupleList,tuple)

          adjMatrix[h,] = 0

          observedWithoutParents <- which(colSums(adjMatrix[])==0)
          observedNodes <- setdiff(observedNodes,observedWithoutParents)
        }
      }
    }
    if(!foundIdentifiableNode){
      return(list("identifiable" = FALSE, "tupleList" = list()))
    }
  }
  return(list("identifiable" = TRUE, "tupleList" = tupleList))
}