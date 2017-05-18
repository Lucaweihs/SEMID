#' Identifiability of linear structural equation models.
#'
#' @description
#' NOTE: \code{graphID} has been deprecated, use \code{\link{semID}} instead.
#'
#' This function checks global and generic identifiability of linear
#' structural equation models. For generic identifiability the function
#' checks a sufficient criterion as well as a necessary criterion but this
#' check may be inconclusive.
#'
#' @export
#'
#' @param L Adjacency matrix for the directed part of the path
#' diagram/mixed graph; an edge pointing from i to j is encoded as L[i,j]=1 and
#' the lack of an edge between i and j is encoded as L[i,j]=0. There should be
#' no directed self loops, i.e. no i such that L[i,i]=1.
#' @param O Adjacency matrix for the bidirected part of the path diagram/mixed
#' graph. Edges are encoded as for the L parameter. Again there should be no
#' self loops. Also this matrix will be coerced to be symmetric so it is only
#' necessary to specify an edge once, i.e. if O[i,j]=1 you may, but are not
#' required to, also have O[j,i]=1.
#' @param output.type A character string indicating whether output is
#' printed ('matrix'), saved to a file ('file'), or returned as a list
#' ('list') for further processing in R.
#' @param file.name A character string naming the output file.
#' @param decomp.if.acyclic A logical value indicating whether an input graph
#' that is acyclic is to be decomposed before applying identifiability criteria.
#' @param test.globalID A logical value indicating whether or not global
#' identifiability is checked.
#' @param test.genericID A logical value indicating whether or not a sufficient
#' condition for generic identifiability is checked.
#' @param test.nonID A logical value indicating whether or not a condition
#' implying generic non-identifiability is checked.
#'
#' @return
#'   A list or printed matrix indicating the identifiability status of the
#'   linear SEM given by the input graph.  Optionally the graph's
#'   components are listed.
#'
#'   With output.type = 'list', the function returns a list of components
#'   for the graph.  Each list entry is again a list that indicates first
#'   which nodes form the component and second whether the component forms
#'   a mixed graph that is acyclic.  The next entries in the list show
#'   HTC-identifiable nodes, meaning nodes v for which the coefficients for
#'   all the directed edges pointing to v can be identified using the
#'   methods from Foygel et al. (2012).  The HTC-identifiable nodes are
#'   listed in the order in which they are found by the recursive
#'   identification algorithm.  The last three list entries are
#'   logical values that indicate whether or not the graph component is
#'   generically identifiable, globally identifiable or not identifiable;
#'   compare Drton et al. (2011) and Foygel et al. (2012).  In the latter
#'   case the Jacobian of the parametrization does not have full rank.
#'
#'   With output.type = 'matrix', a summary of the above
#'   information is printed.
#'
#' @references Drton, M., Foygel, R., and Sullivant, S.  (2011) Global
#' identifiability of linear structural equation models. \emph{Ann. Statist.}
#' 39(2): 865-886.
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
#'
#' @examples
#' \dontrun{
#' L = t(matrix(
#'   c(0, 1, 0, 0, 0,
#'     0, 0, 1, 0, 0,
#'     0, 0, 0, 1, 0,
#'     0, 0, 0, 0, 1,
#'     0, 0, 0, 0, 0), 5, 5))
#' O = t(matrix(
#'   c(0, 0, 1, 1, 0,
#'     0, 0, 0, 1, 1,
#'     0, 0, 0, 0, 0,
#'     0, 0, 0, 0, 0,
#'     0, 0, 0, 0, 0), 5, 5))
#' O=O+t(O)
#' graphID(L,O)
#'
#'
#' ## Examples from Foygel, Draisma & Drton (2012)
#' demo(SEMID)
#' }
graphID <- function(L, O, output.type = "matrix", file.name = NULL, decomp.if.acyclic = TRUE, 
    test.globalID = TRUE, test.genericID = TRUE, test.nonID = TRUE) {
    .Deprecated("semID", package = "SEMID")
    if (!is.matrix(L) || !is.matrix(O)) {
        print("L and O must be two square matrices of the same size")
        return(NULL)
    } else {
        m <- unique(c(dim(L), dim(O)))
        if (length(m) > 1) {
            print("L and O must be two square matrices of the same size")
            return(NULL)
        } else {
            if (!is.character(file.name) && output.type == "file") {
                print("need to input a file name for output.type = 'file'")
                return(NULL)
            } else {
                L <- (L != 0)
                diag(L) <- 0
                O <- O + t(O)
                O <- (O != 0)
                diag(O) <- 0
                Graph.output <- graphID.decompose(L, O, decomp.if.acyclic, test.globalID, 
                  test.genericID, test.nonID)
                if (all(c(test.globalID, test.genericID, test.nonID))) {
                  graphID.output.alltests(Graph.output$Components, Graph.output$Decomp, 
                    output.type, file.name)
                } else {
                  graphID.output.notalltests(Graph.output$Components, Graph.output$Decomp, 
                    output.type, file.name, test.globalID, test.genericID, test.nonID)
                }
            }
        }
    }
}

### Function to format the output when all tests have been run.
graphID.output.alltests <- function(Graph.components, Decomp, output.type, file.name) {
    if (output.type == "list") {
        return(Graph.components)
    } else {
        Graph.output.matrix <- NULL
        for (comp in 1:length(Graph.components)) {
            if (Graph.components[[comp]]$acyclic) {
                if (Graph.components[[comp]]$GlobalID) {
                  Comp.status <- "Globally ID'able"
                } else {
                  if (Graph.components[[comp]]$GenericID) {
                    Comp.status <- "Generically ID'able"
                  } else {
                    if (Graph.components[[comp]]$NonID) {
                      Comp.status <- "Generically non-ID'able"
                    } else {
                      Comp.status <- "HTC-inconclusive"
                    }
                  }
                }
            } else {
                if (Graph.components[[comp]]$GenericID) {
                  Comp.status <- "Generically ID'able"
                } else {
                  if (Graph.components[[comp]]$NonID) {
                    Comp.status <- "Generically non-ID'able"
                  } else {
                    Comp.status <- "HTC-inconclusive"
                  }
                }
            }
            
            if (length(Graph.components[[comp]]$HTC.ID.nodes) == 0) {
                Graph.components[[comp]]$HTC.ID.nodes <- "none"
            }
            
            Comp.row <- c(comp, paste(Graph.components[[comp]]$Nodes, collapse = ","), 
                paste(Graph.components[[comp]]$HTC.ID.nodes, collapse = ","), Comp.status)
            Graph.output.matrix <- rbind(Graph.output.matrix, Comp.row)
        }
        
        Graph.output.matrix <- rbind(c("Component", "Nodes", "HTC-ID'able nodes", 
            "Identifiability"), Graph.output.matrix)
        colnames(Graph.output.matrix) <- rep("", 4)
        rownames(Graph.output.matrix) <- rep("", length(Graph.components) + 1)
        
        if (!Decomp) {
            Graph.output.matrix <- Graph.output.matrix[, -(1:2)]
        }
        
        Max.string <- max(nchar(Graph.output.matrix))
        for (i in 1:dim(Graph.output.matrix)[1]) {
            for (j in 1:dim(Graph.output.matrix)[2]) {
                Graph.output.matrix[i, j] <- paste(c(Graph.output.matrix[i, j], rep(" ", 
                  5 + Max.string - nchar(Graph.output.matrix[i, j]))), collapse = "")
            }
        }
        
        if (output.type == "file") {
            write.table(Graph.output.matrix, file = file.name, quote = FALSE)
        } else {
            print(Graph.output.matrix, quote = FALSE)
        }
    }
}

## Format output when not all tests have been run.
graphID.output.notalltests <- function(Graph.components, Decomp, output.type, file.name, 
    test.globalID, test.genericID, test.nonID) {
    if (output.type == "list") {
        return(Graph.components)
    } else {
        Graph.output.matrix <- NULL
        for (comp in 1:length(Graph.components)) {
            if (Graph.components[[comp]]$acyclic && test.globalID) {
                if (Graph.components[[comp]]$GlobalID) {
                  Comp.status.GlobalID <- "yes"
                } else {
                  Comp.status.GlobalID <- "no"
                }
            } else {
                Comp.status.GlobalID <- "not tested"
            }
            if (test.genericID) {
                if (Graph.components[[comp]]$GenericID) {
                  Comp.status.GenericID <- "yes"
                } else {
                  Comp.status.GenericID <- "no"
                }
            } else {
                Comp.status.GenericID <- "not tested"
            }
            if (test.nonID) {
                if (Graph.components[[comp]]$NonID) {
                  Comp.status.NonID <- "yes"
                } else {
                  Comp.status.NonID <- "no"
                }
            } else {
                Comp.status.NonID <- "not tested"
            }
            
            if (length(Graph.components[[comp]]$HTC.ID.nodes) == 0) {
                if (test.genericID) {
                  Graph.components[[comp]]$HTC.ID.nodes <- "none"
                } else {
                  Graph.components[[comp]]$HTC.ID.nodes <- "not tested"
                }
            }
            
            Comp.row <- c(comp, paste(Graph.components[[comp]]$Nodes, collapse = ","), 
                paste(Graph.components[[comp]]$HTC.ID.nodes, collapse = ","), Comp.status.GlobalID, 
                Comp.status.GenericID, Comp.status.NonID)
            Graph.output.matrix <- rbind(Graph.output.matrix, Comp.row)
        }
        
        Graph.output.matrix <- rbind(c("Component", "Nodes", "HTC-ID'able nodes", 
            "Globally ID'able?", "Generically ID'able?", "Generically non-ID'able?"), 
            Graph.output.matrix)
        colnames(Graph.output.matrix) <- rep("", 6)
        rownames(Graph.output.matrix) <- rep("", length(Graph.components) + 1)
        
        if (!Decomp) {
            Graph.output.matrix <- Graph.output.matrix[, -(1:2)]
        }
        
        Max.string <- max(nchar(Graph.output.matrix))
        for (i in 1:dim(Graph.output.matrix)[1]) {
            for (j in 1:dim(Graph.output.matrix)[2]) {
                Graph.output.matrix[i, j] <- paste(c(Graph.output.matrix[i, j], rep(" ", 
                  3 + Max.string - nchar(Graph.output.matrix[i, j]))), collapse = "")
            }
        }
        
        
        if (output.type == "file") {
            write.table(Graph.output.matrix, file = file.name, quote = FALSE)
        } else {
            print(Graph.output.matrix, quote = FALSE)
        }
    }
}


#' Determine generic identifiability by Tian Decomposition and HTC
#'
#' Split a graph into mixed Tian components and solve each separately
#' using the HTC.
#'
#' @inheritParams graphID
#'
#' @return A list with two named components:
#'
#'   1. Components - a list of lists. Each list represents one mixed Tian component
#'       of the graph. Each list contains named components corresponding to which
#'       nodes are in the component and results of various tests of
#'       identifiability on the component (see the parameter descriptions).
#'
#'  2. Decomp - true if a decomposition occured, false if not.
graphID.decompose <- function(L, O, decomp.if.acyclic = TRUE, test.globalID = TRUE, 
    test.genericID = TRUE, test.nonID = TRUE) {
    m <- nrow(L)
    L <- (L != 0)
    diag(L) <- 0
    O <- O + t(O)
    O <- (O != 0)
    diag(O) <- 0
    
    Decomp <- FALSE
    if (decomp.if.acyclic) {
        test.acyclic <- TRUE
        nodes.acyclic <- 1:m
        while (test.acyclic && !Decomp) {
            if (any(colSums(L[nodes.acyclic, nodes.acyclic, drop = FALSE]) == 0)) {
                nodes.acyclic <- nodes.acyclic[-which(colSums(L[nodes.acyclic, nodes.acyclic, 
                  drop = FALSE]) == 0)]
                if (length(nodes.acyclic) == 0) {
                  Decomp <- TRUE
                }
            } else {
                test.acyclic <- FALSE
            }
        }
    }
    
    if (Decomp) {
        Bidir.comp <- diag(m) + O
        for (i in 1:(m - 1)) {
            Bidir.comp <- (Bidir.comp %*% (diag(m) + O) > 0) + 0
        }
        Components <- list()
        V <- 1:m
        num.comp <- 0
        while (length(V) > 0) {
            num.comp <- num.comp + 1
            i <- min(V)
            Components[[num.comp]] <- which(Bidir.comp[i, ] == 1)
            V <- setdiff(V, Components[[num.comp]])
        }
    } else {
        Components <- list()
        Components[[1]] <- 1:m
        num.comp <- 1
    }
    
    for (comp in 1:num.comp) {
        Component <- Components[[comp]]
        Component.parents <- sort(setdiff(which(rowSums(L[, Component, drop = FALSE]) > 
            0), Component))
        if (length(Components[[comp]]) == 1) {
            Components[[comp]] <- list()
            Components[[comp]]$Nodes <- Component
            Components[[comp]]$acyclic <- TRUE
            if (test.genericID) {
                Components[[comp]]$HTC.ID.nodes <- Component
                Components[[comp]]$GenericID <- TRUE
            }
            if (test.globalID) {
                Components[[comp]]$GlobalID <- TRUE
            }
            if (test.nonID) {
                Components[[comp]]$NonID <- FALSE
            }
        } else {
            m1 <- length(Component)
            m2 <- length(Component.parents)
            L.Component <- cbind(L[c(Component, Component.parents), Component], matrix(0, 
                m1 + m2, m2))
            O.Component <- cbind(rbind(O[Component, Component], matrix(0, m2, m1)), 
                matrix(0, m1 + m2, m2))
            Components[[comp]] <- c(list(Nodes = Component), graphID.main(L.Component, 
                O.Component, test.globalID, test.genericID, test.nonID))
            Components[[comp]]$HTC.ID.nodes <- Component[intersect(Components[[comp]]$HTC.ID.nodes, 
                1:length(Component))]
        }
    }
    Graph.output <- list()
    Graph.output$Components <- Components
    Graph.output$Decomp <- Decomp
    return(Graph.output)
}

#' Helper function to handle a graph component.
#'
#' Calls the other functions that determine identifiability status.
#'
#' @inheritParams graphID
#'
#' @return A list containing named components of the results of various tests
#' desired based on the input parameters.
graphID.main <- function(L, O, test.globalID = TRUE, test.genericID = TRUE, test.nonID = TRUE) {
    m <- nrow(L)
    L <- (L != 0)
    diag(L) <- 0
    O <- O + t(O)
    O <- (O != 0)
    diag(O) <- 0
    
    ILinv <- diag(m)
    for (i in 1:m) {
        ILinv <- 0 + (diag(m) + ILinv %*% L > 0)
    }
    
    Output <- list()
    
    Output$acyclic <- (max(ILinv + t(ILinv) - diag(m)) == 1)
    
    if (test.genericID) {
        Output$HTC.ID.nodes <- graphID.htcID(L, O)
        
        if (length(Output$HTC.ID.nodes) < m) {
            Output$GenericID <- FALSE
            if (Output$acyclic && test.globalID) {
                Output$GlobalID <- FALSE
            }
            if (test.nonID) {
                Output$NonID <- graphID.nonHtcID(L, O)
            }
        } else {
            Output$GenericID <- TRUE
            if (Output$acyclic && test.globalID) {
                Output$GlobalID <- graphID.globalID(L, O)
            }
            if (test.nonID) {
                Output$NonID <- FALSE
            }
        }
    } else {
        if (test.nonID) {
            Output$NonID <- graphID.nonHtcID(L, O)
            if (!Output$NonID) {
                if (Output$acyclic && test.globalID) {
                  Output$GlobalID <- graphID.globalID(L, O)
                }
            }
        } else {
            if (Output$acyclic && test.globalID) {
                Output$GlobalID <- graphID.globalID(L, O)
            }
        }
    }
    return(Output)
}

#' Check for global identifiability of a mixed graph.
#'
#' Checks for the global identifiability of a mixed graph using techniques
#' presented in Drton, Foygel, Sullivant (2011).
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return TRUE if the graph was globally identifiable, FALSE otherwise.
#'
#' @references
#' Drton, Mathias; Foygel, Rina; Sullivant, Seth. Global identifiability of
#' linear structural equation models. \emph{Ann. Statist.}  39 (2011), no. 2,
#' 865--886.
graphID.globalID <- function(L, O) {
    m <- nrow(L)
    validateMatrices(L, O)
    O <- 1 * ((O + t(O)) != 0)
    
    if (any(L + O == 2) || !igraph::is.dag(igraph::graph.adjacency(L))) {
        return(F)
    }
    
    for (i in 1:m) {
        dirGraph <- igraph::graph.adjacency(L)
        biGraph <- igraph::graph.adjacency(O, mode = "undirected")
        ancestors <- as.integer(igraph::neighborhood(dirGraph, m, nodes = i, mode = "in")[[1]])
        biComponent <- as.integer(igraph::neighborhood(biGraph, m, nodes = i)[[1]])
        validVertices <- intersect(ancestors, biComponent)
        
        newL <- 0 * L
        newL[validVertices, validVertices] <- L[validVertices, validVertices]
        newDirGraph <- igraph::graph.adjacency(newL)
        
        for (j in setdiff(validVertices, i)) {
            if (igraph::graph.maxflow(newDirGraph, source = j, target = i)$value != 
                0) {
                return(FALSE)
            }
        }
    }
    
    return(TRUE)
}

#' Determine generic identifiability of a mixed graph.
#'
#' If directed part of input graph is cyclic then will check for generic
#' identifiability using the half-trek criterion. Otherwise will use the a
#' slightly stronger version of the half-trek criterion using ancestor
#' decompositions.
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return The vector of nodes that could be determined to be generically
#' identifiable.
#'
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
#'
#' @references
#' {Drton}, M. and {Weihs}, L. (2015) Generic Identifiability of Linear
#' Structural Equation Models by Ancestor Decomposition. arXiv 1504.02992
graphID.genericID <- function(L, O) {
    .Deprecated("semID", package = "SEMID")
    if (is.dag(igraph::graph.adjacency(L, mode = "directed"))) {
        return(graphID.ancestralID(L, O))
    } else {
        return(graphID.htcID(L, O))
    }
}
