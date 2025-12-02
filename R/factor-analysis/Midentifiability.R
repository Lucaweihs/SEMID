# check M-identifiability
# check sign-identifiability of a graph by iterating over all latent nodes
checkMidentifiability <- function(lambda, maxCard = length(observedNodes)){

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
