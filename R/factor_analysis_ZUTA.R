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
  adjMatrix <- input[[1]]
  latentNodes <- input[[2]]
  observedNodes <- input[[3]]

  result <- list()
  result$latentNodes <- latentNodes
  result$observedNodes <- observedNodes
  result$call <- match.call()
  class(result) <- "ZUTAresult"

  # generate matrix with only latent rows and observed columns and no zero-rows or zero-columns
  numberOfNodes <- length(latentNodes) + length(observedNodes)

  allRows <- c(1:numberOfNodes)
  notLatentRows <- (setdiff(allRows,latentNodes)) * (-1)
  notObservedColumns <- (setdiff(allRows,observedNodes)) * (-1)

  cleanMatrix <- adjMatrix[notLatentRows, notObservedColumns]

  if (sum(rowSums(cleanMatrix[])>0)==1){
    result$zuta <- TRUE
    return(result)
  }
  if (sum(colSums(cleanMatrix[])>0)==1){
    if (sum(cleanMatrix)==1){
      result$zuta <- TRUE
      return(result)
    } else {
      result$zuta <- FALSE
      return(result)
    }
  }

  cleanMatrix <- cleanMatrix[rowSums(cleanMatrix[])>0,colSums(cleanMatrix[])>0]


  if(nrow(cleanMatrix)>1){
    if(findColumnsWithSumOne(cleanMatrix)){
      result$zuta <- TRUE
      return(result)
    } else {
      result$zuta <- FALSE
      return(result)
    }
  } else {
    result$zuta <- TRUE
    return(result)
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

#' Prints a ZUTAresult object
#'
#' Prints a ZUTAresult object as returned by
#' \code{\link{ZUTA}}. Invisibly returns its argument via
#' \code{\link{invisible}(x)} as most print functions do.
#'
#' @export
#'
#' @param x the ZUTAresult object
#' @param ... optional parameters, currently unused.
print.ZUTAresult <- function(x, ...) {
  cat("Call: ")
  print(x$call)

  cat("\nFactor Analysis Graph Info:\n")
  cat("latent nodes: ", x$latentNodes, "\n")
  cat("observed nodes: ", x$observedNodes, "\n\n")

  cat(sprintf("ZUTA:    %s\n", x$zuta))

  invisible(x)
}