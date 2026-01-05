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
      print(tupleForNode)
      if(tupleForNode$found){
        foundIdentifiableNode <- TRUE
        latentNodes <- setdiff(latentNodes,h)
        S <- union(S,h)
        notIdentifiedNodes <- notIdentifiedNodes[! notIdentifiedNodes %in% c(h)]
        tuple <- list(list("h"=tupleForNode$h, "S"=S, "v"=tupleForNode$v, "W"=tupleForNode$W, "U"=tupleForNode$U))
        tupleList <- c(tupleList,tuple)
      }
    }

    if(!foundIdentifiableNode){
      return(list("identifiable" = FALSE, "tupleList" = list()))
    }
  }
  return(list("identifiable" = TRUE, "tupleList" = tupleList))
}
