# check whether a graph fulfills the Zero Upper Triangular Assumption
checkZUTA <- function(lambda){

  input <- transformLambda(lambda)
  adjMatrix <- input[1]
  latentNodes <- input[2]
  observedNodes <- input[3]

  # generate matrix with only latent rows and observed columns and no zero-rows or zero-columns
  numberOfNodes <- nrow(adjMatrix)

  allRows <- c(1:numberOfNodes)
  notLatentRows <- (setdiff(allRows,latentNodes)) * (-1)
  notObservedColumns <- (setdiff(allRows,observedNodes)) * (-1)

  cleanMatrix <- adjMatrix[notLatentRows, notObservedColumns]

  if (sum(rowSums(cleanMatrix[])>0)==1){
    return(TRUE)
  }
  if (sum(colSums(cleanMatrix[])>0)==1){
    if (sum(cleanMatrix)==1){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  cleanMatrix <- cleanMatrix[rowSums(cleanMatrix[])>0,colSums(cleanMatrix[])>0]


  if(nrow(cleanMatrix)>1){
    if(findColumnsWith1(cleanMatrix)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(TRUE)
  }
}

# find a column with sum=1,  delete the row with the 1 in that column, check smaller matrix
findColumnsWith1 <- function(cleanMatrix){
  if(any(colSums(cleanMatrix[])==1)){
    if(nrow(cleanMatrix) > 2){
      for(column in which(colSums(cleanMatrix[])==1)){
        row <- which(cleanMatrix[,column]==1)
        newMatrix <- cleanMatrix[-row, ]
        if(findColumnsWith1(newMatrix)){
          isZUTA <- TRUE
          return(TRUE)
        }
      }
      return(FALSE)
    } else {
      return (TRUE)
    }
  } else {
    return(FALSE)
  }
}