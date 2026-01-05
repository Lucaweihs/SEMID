#' Check the Zero Upper Triangular Assumption
#'
#' @param lambda adjacency matrix with number of cols = number of latent nodes,
#'        number of rows = number of observed nodes
#'
#' @returns a Boolean, whether the graph fulfills ZUTA
#' @export
#'
#' @references
#' Sturma, N., Kranzlmüller, M., Portakal, I., and Drton, M.  (2025) Matching
#' Criterion for Identifiability in Sparse Factor Analysis.
#' arXiv:2502.02986
ZUTA <- function(lambda){

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
    if(findColumnsWithSumOne(cleanMatrix)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(TRUE)
  }
}

#' A Helper Function for Check ZUTA
#'
#' iterative function that tries to find a column with sum=1, delete the row
#' with the 1 in that column, check smaller matrix
#' @param cleanMatrix a matrix without zero columns and rows
#'
#' @returns a Boolean whether it finds a column with sum=1 and the new matrix,
#'          without the row in which the entry is, fulfills ZUTA
findColumnsWithSumOne <- function(cleanMatrix){
  if(any(colSums(cleanMatrix[])==1)){
    if(nrow(cleanMatrix) > 2){
      for(column in which(colSums(cleanMatrix[])==1)){
        row <- which(cleanMatrix[,column]==1)
        newMatrix <- cleanMatrix[-row, ]
        if(findColumnsWithSumOne(newMatrix)){
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