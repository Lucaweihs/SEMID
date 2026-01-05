#' Check Local BB-Criterion
#'
#' for a set of latent nodes, find a set of observed nodes U so that the tuple
#' satisfies the local BB-criterion and the inequality concerning the
#' cardinalities of the sets
#'
#' @param adjMatrix of the graph
#' @param latentNodes set of latent nodes
#' @param observedNodes set of observed nodes
#'
#' @returns a list consisting of a Boolean, whether the graph is
#'         sign-identifiable and if yes, a list consisting of the sets
#' @export
#'
#' @references
#' Sturma, N., Kranzlmüller, M., Portakal, I., and Drton, M.  (2025) Matching
#' Criterion for Identifiability in Sparse Factor Analysis.
#' arXiv:2502.02986
checkLocalBBCriterion <- function(adjMatrix, latentNodes, observedNodes){

  for(h in latentNodes){
    childrenOfH <- childrenOfNodes(adjMatrix, h, observedNodes)
    for(U in rje::powerSet(childrenOfH)){
      if(length(U)>2){
        jointParentsU <- jointParents(adjMatrix, U, latentNodes)
        p <- length(U)
        m <- length(jointParentsU)
        if(p*(m+1)- choose(m,2) < choose(p+1, 2)){
          if(fullFactorCriterion(adjMatrix, U, jointParentsU, latentNodes, observedNodes)){
            return(list(found= TRUE, newNodesInS = jointParentsU, U = U))
          }
        }
      }
    }
  }
  return(list(found = FALSE))
}

# check for a given tuple whether the local BB-criterion is fulfilled
fullFactorCriterion <- function(adjMatrix, U, jointParentsU, latentNodes, observedNodes){

  # criterion (i)
  inducedSubgraphMatrix <- adjMatrixInducedSubgraph(adjMatrix, U, jointParentsU, latentNodes, observedNodes)

  checkFullFactorZUTA <- checkFullFactorZUTA(inducedSubgraphMatrix, jointParentsU, U)
  isZUTA <- checkFullFactorZUTA$zuta
  if(!isZUTA){
    return(FALSE)
  }

  # criterion (ii)
  orderingZUTA <- checkFullFactorZUTA$ordering
  for(h in jointParentsU){
    setOfV <- setdiff(childrenOfNodes(adjMatrix, h, observedNodes), U)
    positionOfH <- match(h, orderingZUTA)
    setOfL <- orderingZUTA[1:positionOfH]

    remainingVtoCheck <- setOfV
    while(length(remainingVtoCheck)>0){
      foundU <- FALSE
      setUWithCheckedV <- U
      for(v in remainingVtoCheck){
        for(u in setUWithCheckedV){
          jointParentsVandU <- jointParents(adjMatrix, c(u,v), latentNodes)
          if(all(jointParentsVandU %in% setOfL)){
            foundU <- TRUE
            break
          }
        }
        if(foundU){
          remainingVtoCheck <- setdiff(remainingVtoCheck, v)
          setUWithCheckedV <- union(setUWithCheckedV, v)
          break
        }
      }
      if(!foundU){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

# check if a factor analysis graph satisfies the Zero Upper Triangular Assumption but has all other edges
checkFullFactorZUTA <- function(adjMatrix, latentNodes, observedNodes){

  # order latent nodes depending on their number of children
  numberOfChildren <- c()

  for(parent in latentNodes){
    counter <- 0
    for(child in observedNodes){
      if(adjMatrix[parent, child] == 1){
        counter <- counter+1
      }
    }
    numberOfChildren <- append(numberOfChildren, counter)
  }
  orderedLatentNodes <- latentNodes[order(numberOfChildren, decreasing = TRUE)]

  # check if node with most children has all children
  if(length(observedNodes)>max(numberOfChildren)){
    # not ZUTA
    return(list(zuta = FALSE))
  }
  # check if one number of children is missing e.g. 5, 3, 2
  if(!(all(c((length(observedNodes)-length(latentNodes)+1):length(observedNodes)) %in% numberOfChildren))){
    # not ZUTA
    return(list(zuta = FALSE))
  }
  # check if two nodes have same number of children
  if(length(latentNodes)>length(unique(numberOfChildren))){
    # not ZUTA
    return(list(zuta = FALSE))
  }

  remainingOrderedLatentNodes <- orderedLatentNodes

  for(node in orderedLatentNodes){
    remainingOrderedLatentNodes <- setdiff(remainingOrderedLatentNodes, node)
    if(length(remainingOrderedLatentNodes) == 0){
      break
    }
    secondNode <- remainingOrderedLatentNodes[1]
    optionOnlyChild <- c()

    # find the observed node, that is child of first, but not of second latent node
    for(observedNode in observedNodes){
      if(adjMatrix[node, observedNode]==1 & adjMatrix[secondNode, observedNode]==0){
        optionOnlyChild <- append(optionOnlyChild,observedNode)
      }
    }

    # check if there is exactly one node, that is child of first, but not of second node
    if(length(optionOnlyChild) == 0 | length(optionOnlyChild) > 1){
      # not ZUTA
      return(list(zuta = FALSE))
    }

    # check if this node is also not a child of all other latent nodes
    for(otherNode in setdiff(remainingOrderedLatentNodes, secondNode)){
      if(adjMatrix[node, optionOnlyChild] == 1 & adjMatrix[otherNode, optionOnlyChild] == 1){
        # not ZUTA
        return(list(zuta = FALSE))
      }
    }
  }
  return(list(zuta = TRUE, ordering = orderedLatentNodes))
}

# compute the adjacency matrix for the induced subgraph of a set U of observed nodes
# and the joint parents of U without S of latent nodes
adjMatrixInducedSubgraph <- function(adjMatrix, U, jointParentsU, latentNodes, observedNodes){
  adjMatrixSubgraph <- adjMatrix
  for(latentNode in latentNodes){
    if(latentNode %in% jointParentsU){
      for(observedNode in observedNodes){
        if(!(observedNode %in% U)){
          adjMatrixSubgraph[latentNode, observedNode] <- 0
        }
      }
    } else {
      for(observedNode in observedNodes){
        adjMatrixSubgraph[latentNode, observedNode] <- 0
      }
    }
  }
  return(adjMatrixSubgraph)
}
