# extended M-identifiability check
# check sign-identifiability of a graph
checkExtendedMidentifiability <- function(lambda, maxCard = length(observedNodes)){

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