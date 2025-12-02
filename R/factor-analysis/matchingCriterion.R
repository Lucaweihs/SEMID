# implementation of the matching criterion
# for a given latent node try to find a tuple satisfying the matching criterion
checkMatchingCriterion <- function(flowGraphAdjMatrix, adjMatrix, h, latentNodes, observedNodes, maxCard = length(observedNodes)){
  for(v in observedNodes){
    # find parents of v
    vParents <- c()
    for(node in latentNodes){
      if (adjMatrix[node,v]==1){
        vParents <- union(vParents,node)
      }
    }

    if(setequal(vParents,h)){
      observedNodesWithoutV <- setdiff(observedNodes, v)
      maxSizeW <- min(min(length(observedNodesWithoutV)/2,length(latentNodes)),maxCard)
      for(W in powerSet(observedNodesWithoutV, maxSizeW)){
        if((length(W)>0) && (sum(adjMatrix[h,W])>=1)){
          observedNodesWithoutW <- setdiff(observedNodesWithoutV,W)
          possibleU = children(adjMatrix, parents(adjMatrix, W, latentNodes), observedNodesWithoutW)
          if (length(possibleU) >= length(W)){
            for(U in combn(possibleU,length(W),simplify=FALSE)){
              if(sum(adjMatrix[h,U])>=1){
                fulfillsMatchingCriterion <- matchingCriterion(flowGraphAdjMatrix, adjMatrix, v, W, U, h, vParents, latentNodes, observedNodes)
                if(fulfillsMatchingCriterion){
                  return(list(found = TRUE, h = h, v = v, W = W, U = U))
                }
              }
            }
          }
        }
      }
    }
  }
  return(list(found = FALSE))
}

# check for a given tuple whether the matching criterion is fulfilled
matchingCriterion <- function(flowGraphAdjMatrix, adjMatrix, v, W, U, h, vParents, latentNodes, observedNodes){

  # (i)
  if (any(v==union(W,U))){
    return(FALSE)
  }

  # (ii)
  if (length(intersect(W,U))>0 | !(length(W)==length(U)) | length(W)==0  ){
    return(FALSE)
  }
  # (iii)
  maxFlowThird <- maxFlowSTGraph(flowGraphAdjMatrix, adjMatrix, latentNodes, W, U)$value
  if (!(maxFlowThird==length(W))){
    return(FALSE)
  }

  # (iv)
  maxFlowFourth <- maxFlowSTGraph(flowGraphAdjMatrix, adjMatrix, latentNodes, union(W,v), union(U,v))$value
  if (maxFlowFourth==(length(W)+1)){
    return(FALSE)
  }
  return(TRUE)
}

# calculate generic adjacency matrix that is adapted for each flow graph
flowGraphMatrix <- function(adjMatrixGraph, latentNodes, observedNodes){
  m <- nrow(adjMatrixGraph)
  flowAdjMatrix <- matrix(0, 2*m+2, 2*m+2)

  s <- 2*m+1
  t <- 2*m+2

  for (x in observedNodes){
    # connect observed Nodes to s
    flowAdjMatrix[s, x] <- 1
    # connect observed Nodes to t
    flowAdjMatrix[m+x, t] <- 1
  }

  # connect latent nodes to copy of themselves
  for (y in latentNodes){
    flowAdjMatrix[y, m+y] <- 1
  }
  return(flowAdjMatrix)
}

# calculate maximum flow
maxFlowSTGraph <- function(flowGraphMatrix, adjMatrixGraph, latentNodes, W, U){
  flowAdjMatrix <- flowGraphMatrix
  m <- nrow(adjMatrixGraph)
  s <- 2*m+1
  t <- 2*m+2

  for (y in latentNodes){
    # connect nodes in W to latent nodes not in S
    for (x in W){
      if (adjMatrixGraph[y, x] == 1 ){
        flowAdjMatrix[x,y] <- 1
      }
    }
    # connect nodes in U to latent nodes not in S
    for (x in U){
      if (adjMatrixGraph[y, x] == 1 ){
        flowAdjMatrix[m+y, m+x] <- 1
      }
    }
  }
  stFlowGraph <- graph_from_adjacency_matrix(flowAdjMatrix, mode = "directed")
  value <- max_flow(stFlowGraph, s, t)
  return(value)
}