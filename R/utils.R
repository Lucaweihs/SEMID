#' Checks a MixedGraph has appropriate node numbering
#'
#' Checks that the input mixed graph has vertices are numbered from 1
#' to mixedGraph$numNodes(). Throws an error if they are not.
#'
#' @param mixedGraph the mixed graph object
mixedGraphHasSimpleNumbering <- function(mixedGraph) {
    nodes <- mixedGraph$nodes()
    if (any(nodes != 1:length(nodes))) {
        stop(paste("Currently only mixed graphs whose vertices are numbered from 1",
            "to mixedGraph$numNodes() in order are supported. The inputted mixed",
            "graph as vertices numbered as mixedGraph$nodes() == c(", paste(nodes,
                collapse = ", "), ")."))

    }
}

#' Identifiability of linear structural equation models.
#'
#' This function can be used to check global and generic identifiability of
#' linear structural equation models (L-SEMs). In particular, this function
#' takes a \code{\link{MixedGraph}} object corresponding to the L-SEM and
#' checks different conditions known for global and generic identifiability.
#'
#' @export
#'
#' @param mixedGraph a \code{\link{MixedGraph}} object representing the L-SEM.
#' @param testGlobalID TRUE or FALSE if the graph should be tested for global
#'        identifiability. This uses the \code{\link{graphID.globalID}}
#'        function.
#' @param testGenericNonID TRUE of FALSE if the graph should be tested for
#'        generic non-identifiability, that is, if for every generic choice
#'        of parameters for the L-SEM there are infinitely many
#'        other choices that lead to the same covariance matrix. This currently
#'        uses the \code{\link{graphID.nonHtcID}} function.
#' @param genericIdStepFunctions a list of the generic identifier step functions
#'        that should be used for testing generic identifiability. See
#'        \code{\link{generalGenericID}} for a discussion of such functions. If
#'        this list is empty then generic identifiability is not tested. By
#'        default this will (only) run the half-trek criterion (see
#'        \code{\link{htcIdentifyStep}}) for generic identifiability.
#' @param tianDecompose TRUE or FALSE if the mixed graph should be Tian
#'                      decomposed before running the identification algorithms
#'                      (when appropriate). In general letting this be TRUE will
#'                      make the algorithm faster and more powerful. Note that
#'                      this is a version of the Tian decomposition that works
#'                      also with cyclic graphs.
#' @inheritParams generalGenericID
#'
#' @return returns an object of \link{class} '\code{SEMIDResult},' this
#'         object is just a list with 6 components:
#' \describe{
#'   \item{\code{isGlobalID}}{If testGlobalID == TRUE, then TRUE or FALSE if
#'   the graph is globally identifiable. If testGlobalID == FALSE then NA.}
#'   \item{\code{isGenericNonID}}{If testGenericNonID == TRUE, then TRUE if the
#'   graph is generically non-identifiable or FALSE the test is inconclusive.
#'   If testGenericNonID == FALSE then NA.}
#'   \item{\code{genericIDResult}}{If length(genericIdStepFunctions) != 0 then
#'   a \code{GenericIDResult} object as returned by
#'   \code{\link{generalGenericID}}. Otherwise a list of length 0.}
#'   \item{\code{mixedGraph}}{the inputted mixed graph object.}
#'   \item{\code{tianDecompose}}{the argument tianDecompose.}
#'   \item{\code{call}}{the call made to this function.}
#' }
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
#' O = O + t(O)
#' graph = MixedGraph(L,O)
#' semID(graph)
#'
#' ## Examples from Foygel, Draisma & Drton (2012)
#' demo(SEMID)
#' }
semID <- function(mixedGraph, testGlobalID = TRUE, testGenericNonID = TRUE, genericIdStepFunctions = list(htcIdentifyStep),
    tianDecompose = TRUE) {
    isGlobalID <- NA
    if (testGlobalID) {
        isGlobalID <- graphID.globalID(mixedGraph$L(), mixedGraph$O())
    }

    isGenericNonID <- NA
    if (testGenericNonID) {
        if (tianDecompose) {
            isGenericNonID <- FALSE
            for (cComponent in mixedGraph$tianDecompose()) {
                isGenericNonID <- isGenericNonID || graphID.nonHtcID(cComponent$L,
                  cComponent$O)
                if (isGenericNonID) {
                  break
                }
            }
        } else {
            isGenericNonID <- graphID.nonHtcID(mixedGraph$L(), mixedGraph$O())
        }
    }

    genericIDResult <- list()
    if (length(genericIdStepFunctions) != 0) {
        genericIDResult <- generalGenericID(mixedGraph, idStepFunctions = genericIdStepFunctions,
            tianDecompose = tianDecompose)
    }

    result <- list(isGlobalID = isGlobalID, isGenericNonID = isGenericNonID, genericIDResult = genericIDResult,
        mixedGraph = mixedGraph, tianDecompose = tianDecompose, call = match.call())
    class(result) <- "SEMIDResult"
    return(result)
}

#' Prints a SEMIDResult object
#'
#' Prints a SEMIDResult object as returned by
#' \code{\link{semID}}. Invisibly returns its argument via
#' \code{\link{invisible}(x)} as most print functions do.
#'
#' @export
#'
#' @param x the SEMIDResult object
#' @param ... optional parameters, currently unused.
print.SEMIDResult <- function(x, ...) {
    cat("Call: ")
    print(x$call)

    cat("\nAttempted Tian decomposition?\n")
    cat(x$tianDecompose, "\n")

    if (!is.na(x$isGlobalID)) {
        cat(paste("\nIs globally identifiable?:\n", x$isGlobalID, "\n", sep = ""))
    }

    if (!is.na(x$isGenericNonID)) {
        cat(paste("\nHas a generically infinite-to-one parameterization?:\n", if (x$isGenericNonID) {
            "TRUE"
        } else if (length(x$genericIDResult) != 0 && length(unlist(x$genericIDResult$unsolvedParents)) ==
            0) {
            "FALSE"
        } else {
            "INCONCLUSIVE"
        }, "\n", sep = ""))
    }

    if (length(x$genericIDResult) != 0) {
        cat(paste("\nNumber of parameters shown generically identifiable:\n"))
        cat(paste("Directed edges:", length(unlist(x$genericIDResult$solvedParents)),
            "out of", sum(x$genericIDResult$mixedGraph$L()), "\n"))
        cat(paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2,
            "out of", sum(x$genericIDResult$mixedGraph$O())/2, "\n"))
    }

    invisible(x)
}

#' A helper function to validate input matrices.
#'
#' This helper function validates that the two input matrices, L and O, are of
#' the appropriate form to be interpreted by the other functions. In particular
#' they should be square matrices of 1's and 0's with all 0's along their
#' diagonals. We do not require O to be symmetric here.
#'
#' @param L See above description.
#' @param O See above description.
#'
#' @return This function has no return value.
validateMatrices <- function(L, O) {
    if (!is.matrix(L) || !is.matrix(O)) {
        stop("L and O must be matrices.")
    } else if (length(unique(c(dim(L), dim(O)))) != 1) {
        stop("L and O must both be square matrices of the same dimensions.")
    }
    t1 <- all(L %in% c(0, 1))
    t2 <- all(O %in% c(0, 1))
    if (!t1 || !t2) {
        stop("L and O must contain only 1's and 0's.")
    } else if (any(diag(L) != 0) || any(diag(O) != 0)) {
        stop("L and O must have 0's along their diagonals.")
    }
}

#' Identify bidirected edges if all directed edges are identified
#'
#' Creates an identifier function that assumes that all directed edges have
#' already been identified and then is able to identify all bidirected edges
#' simultaneously.
#'
#' @param idFunc an identifier function that identifies all directed edges
#'
#' @return a new identifier function that identifies everything.
createSimpleBiDirIdentifier <- function(idFunc) {
    idFunc <- idFunc
    return(function(Sigma) {
        m <- nrow(Sigma)
        identifiedParams <- idFunc(Sigma)
        Lambda <- identifiedParams$Lambda
        if (!any(is.na(Lambda))) {
            Omega <- t(diag(m) - Lambda) %*% Sigma %*% (diag(m) - Lambda)
            return(list(Lambda = Lambda, Omega = Omega))
        }
        return(list(Lambda = Lambda, Omega = identifiedParams$Omega))
    })
}

#' Globally identify the covariance matrix of a C-component
#'
#' The Tian decomposition of a mixed graph G allows one to globally identify
#' the covariance matrices Sigma' of special subgraphs of G called c-components.
#' This function takes the covariance matrix Sigma corresponding to G and
#' a collection of node sets which specify the c-component, and returns the
#' Sigma' corresponding to the c-component.
#'
#' @param Sigma the covariance matrix for the mixed graph G
#' @param internal an integer vector corresponding to the vertices of the
#'        C-component that are in the bidirected equivalence classes (if the
#'        graph is not-acyclic then these equivalence classes must be enlarged
#'        by combining two bidirected components if there are two vertices, one
#'        in each component, that are simultaneously on the same directed cycle).
#' @param incoming the parents of vertices in internal that are not in the set
#'        internal themselves
#' @param topOrder a topological ordering of c(internal, incoming) with respect
#'        to the graph G. For vertices in a strongly connected component the
#'        ordering is allowed to be arbitrary.
#'
#' @export
#'
#' @return the new Sigma corresponding to the c-component
tianSigmaForComponent <- function(Sigma, internal, incoming, topOrder) {
    if (length(incoming) == 0) {
        return(Sigma[topOrder, topOrder, drop = F])
    }
    newSigmaInv <- matrix(0, length(topOrder), length(topOrder))
    for (j in 1:length(topOrder)) {
        node <- topOrder[j]
        if (node %in% internal) {
            if (j == 1) {
                newSigmaInv[j, j] <- newSigmaInv[j, j] + 1/Sigma[node, node, drop = F]
            } else {
                inds <- topOrder[1:(j - 1)]

                SigmaIndsInv <- solve(Sigma[inds, inds, drop = F])
                schurInv <- solve(Sigma[node, node, drop = F] - Sigma[node, inds,
                  drop = F] %*% SigmaIndsInv %*% Sigma[inds, node, drop = F])
                newSigmaInv[j, j] <- newSigmaInv[j, j] + schurInv
                meanMat <- as.numeric(solve(Sigma[inds, inds]) %*% Sigma[inds, node] %*%
                  schurInv)
                newSigmaInv[j, 1:(j - 1)] <- newSigmaInv[j, 1:(j - 1)] - meanMat
                newSigmaInv[1:(j - 1), j] <- newSigmaInv[j, 1:(j - 1)]
                newSigmaInv[1:(j - 1), 1:(j - 1)] <- newSigmaInv[1:(j - 1), 1:(j -
                  1)] + SigmaIndsInv %*% Sigma[inds, node] %*% meanMat
            }
        }
    }
    newSigmaInv[topOrder %in% incoming, topOrder %in% incoming] <- newSigmaInv[topOrder %in%
        incoming, topOrder %in% incoming] + diag(length(incoming))
    newSigma <- solve(newSigmaInv)
    return(newSigma)
}

#' Identifies components in a tian decomposition
#'
#' Creates an identification function which combines the identification
#' functions created on a collection of c-components into a identification
#' for the full mixed graph.
#'
#' @param idFuncs a list of identifier functions for the c-components
#' @param cComponents the c-components of the mixed graph as returned by
#'                    \code{\link{tianDecompose}}.
#'
#' @return a new identifier function
tianIdentifier <- function(idFuncs, cComponents) {
    idFuncs <- idFuncs
    cComponents <- cComponents
    return(function(Sigma) {
        Lambda <- matrix(NA, ncol(Sigma), ncol(Sigma))
        Omega <- matrix(NA, ncol(Sigma), ncol(Sigma))

        for (i in 1:length(cComponents)) {
            internal <- cComponents[[i]]$internal
            incoming <- cComponents[[i]]$incoming
            topOrder <- cComponents[[i]]$topOrder

            newSigma <- tianSigmaForComponent(Sigma, internal, incoming, topOrder)

            result <- idFuncs[[i]](newSigma)
            internalInds <- which(topOrder %in% internal)
            Lambda[topOrder, internal] <- result$Lambda[, internalInds]
            Lambda[-topOrder, internal] <- 0
            Omega[internal, internal] <- result$Omega[internalInds, internalInds]
            Omega[-internal, internal] <- 0
            Omega[internal, -internal] <- 0
        }
        return(list(Lambda = Lambda, Omega = Omega))
    })
}

#' A general generic identification algorithm template.
#'
#' A function that encapsulates the general structure of our algorithms for
#' testing generic identifiability. Allows for various identification algorithms
#' to be used in concert, in particular it will use the identifier functions
#' in the list \code{idStepFunctions} sequentially until it can find no more
#' identifications. The step functions that are currently available for use are
#' in \code{idStepFunctions}
#' \enumerate{
#'   \item htcIdentifyStep
#'   \item ancestralIdentifyStep
#'   \item edgewiseIdentifyStep
#'   \item trekSeparationIdentifyStep
#' }
#'
#' @export
#'
#' @inheritParams semID
#' @param tianDecompose TRUE or FALSE determining whether or not the Tian
#'                      decomposition should be used before running the
#'                      current generic identification algorithm. In general
#'                      letting this be TRUE will make the algorithm faster and
#'                      more powerful.
#' @param idStepFunctions a list of identification step functions
#'
#' @return returns an object of \link{class} '\code{GenericIDResult},' this
#'         object is just a list with 9 components:
#' \describe{
#'   \item{\code{solvedParents}}{a list whose ith element contains a vector
#'   containing the subsets of parents of node i for which the edge j->i could
#'   be shown to be generically identifiable.}
#'   \item{\code{unsolvedParents}}{as for \code{solvedParents} but for the
#'   unsolved parents.}
#'   \item{\code{solvedSiblings}}{as for \code{solvedParents} but for the
#'   siblings of node i (i.e. the bidirected neighbors of i).}
#'   \item{\code{unsolvedSiblings}}{as for \code{solvedSilbings} but for the
#'   unsolved siblings of node i (i.e. the bidirected neighbors of i).}
#'   \item{\code{identifier}}{a function that takes a (generic) covariance
#'   matrix corresponding to the graph and identifies the edges parameters
#'   from solvedParents and solvedSiblings. See \code{\link{htcIdentifyStep}}
#'   for a more in-depth discussion of identifier functions.}
#'   \item{\code{mixedGraph}}{a mixed graph object of the graph.}
#'   \item{\code{idStepFunctions}}{a list of functions used to generically
#'   identify parameters. For instance, htcID uses the function
#'   \code{\link{htcIdentifyStep}} to identify edges.}
#'   \item{\code{tianDecompose}}{the argument tianDecompose.}
#'   \item{\code{call}}{the call made to this function.}
#' }
generalGenericID <- function(mixedGraph, idStepFunctions, tianDecompose = T) {
    mixedGraphHasSimpleNumbering(mixedGraph)
    m <- mixedGraph$numNodes()
    unsolvedParents <- lapply(1:m, function(node) {
        mixedGraph$parents(node)
    })
    solvedParents <- rep(list(numeric(0)), m)

    if (!tianDecompose) {
        identifier <- createIdentifierBaseCase(mixedGraph$L(), mixedGraph$O())

        changeFlag <- T
        while (changeFlag) {
            for (idStepFunction in idStepFunctions) {
                idResult <- idStepFunction(mixedGraph, unsolvedParents, solvedParents,
                  identifier)
                changeFlag <- length(idResult$identifiedEdges) != 0
                unsolvedParents <- idResult$unsolvedParents
                solvedParents <- idResult$solvedParents
                identifier <- idResult$identifier
                if (changeFlag) {
                  break
                }
            }
        }
        if (length(unlist(unsolvedParents)) == 0) {
            identifier <- createSimpleBiDirIdentifier(identifier)
            solvedSiblings <- lapply(1:m, FUN = function(x) {
                mixedGraph$siblings(x)
            })
            unsolvedSiblings <- rep(list(integer(0)), m)
        } else {
            solvedSiblings <- rep(list(integer(0)), m)
            unsolvedSiblings <- lapply(1:m, FUN = function(x) {
                mixedGraph$siblings(x)
            })
        }
    } else {
        solvedSiblings <- rep(list(integer(0)), m)
        cComps <- mixedGraph$tianDecompose()

        compResults <- vector("list", length(cComps))
        identifiers <- vector("list", length(cComps))

        for (i in 1:length(cComps)) {
            result <- generalGenericID(MixedGraph(cComps[[i]]$L, cComps[[i]]$O),
                idStepFunctions, tianDecompose = F)
            topOrder <- cComps[[i]]$topOrder
            compResults[[i]] <- result
            for (j in 1:length(topOrder)) {
                solvedParents[[topOrder[j]]] <- c(solvedParents[[topOrder[j]]], topOrder[result$solvedParents[[j]]])
                solvedSiblings[[topOrder[j]]] <- c(solvedSiblings[[topOrder[j]]],
                  topOrder[result$solvedSiblings[[j]]])
            }

            identifiers[[i]] <- result$identifier
        }

        unsolvedParents <- lapply(1:m, FUN = function(x) {
            setdiff(mixedGraph$parents(x), solvedParents[[x]])
        })
        unsolvedSiblings <- lapply(1:m, FUN = function(x) {
            setdiff(mixedGraph$siblings(x), solvedSiblings[[x]])
        })
        identifier <- tianIdentifier(identifiers, cComps)
    }
    result <- list()
    class(result) <- "GenericIDResult"
    result$solvedParents <- solvedParents
    result$unsolvedParents <- unsolvedParents
    result$solvedSiblings <- solvedSiblings
    result$unsolvedSiblings <- unsolvedSiblings
    result$identifier <- identifier
    result$mixedGraph <- mixedGraph
    result$idStepFunctions <- idStepFunctions
    result$tianDecompose <- tianDecompose
    result$call <- match.call()
    return(result)
}

#' Prints a GenericIDResult object
#'
#' Prints a GenericIDResult object as returned by
#' \code{\link{generalGenericID}}. Invisibly returns its argument via
#' \code{\link{invisible}(x)} as most print functions do.
#'
#' @export
#'
#' @param x the GenericIDResult object
#' @param ... optional parameters, currently unused.
print.GenericIDResult <- function(x, ...) {
    cat("Call: ")
    print(x$call)

    solvedParents <- x$solvedParents
    solvedSiblings <- x$solvedSiblings
    n <- length(x$solvedParents)

    cat(paste("\nMixed Graph Info.\n"))
    cat(paste("# nodes:", n, "\n"))
    cat(paste("# dir. edges:", sum(x$mixedGraph$L()), "\n"))
    cat(paste("# bi. edges:", sum(x$mixedGraph$O())/2, "\n"))

    cat(paste("\nGeneric Identifiability Summary\n"))
    cat(paste("# dir. edges shown gen. identifiable:", length(unlist(solvedParents)),
        "\n"))
    cat(paste("# bi. edges shown gen. identifiable:", length(unlist(solvedSiblings))/2,
        "\n"))

    cat("\nGenerically identifiable dir. edges:\n")
    edges <- character(min(length(unlist(solvedParents))/2, 11))
    k <- 0
    for (i in 1:n) {
        if (length(solvedParents[[i]]) != 0) {
            for (j in solvedParents[[i]]) {
                k <- k + 1
                edges[k] <- paste(j, "->", i, sep = "")

                if (k == 10) {
                  edges[11] <- "..."
                  break
                }
            }
            if (k == 10) {
                break
            }
        }
    }
    if (length(edges) == 0) {
        cat("None\n")
    } else {
        cat(paste(paste(edges, collapse = ", "), "\n"))
    }

    cat("\nGenerically identifiable bi. edges:\n")
    edges <- character(min(length(unlist(solvedSiblings))/2, 11))
    k <- 0
    for (i in 1:n) {
        solvedSibs <- setdiff(solvedSiblings[[i]], 1:i)
        if (length(solvedSibs) != 0) {
            for (j in solvedSibs) {
                k <- k + 1
                edges[k] <- paste(i, "<->", j, sep = "")

                if (k == 10) {
                  edges[11] <- "..."
                  break
                }
            }
            if (k == 10) {
                break
            }
        }
    }
    if (length(edges) == 0) {
        cat("None\n")
    } else {
        cat(paste(paste(edges, collapse = ", "), "\n"))
    }
    invisible(x)
}

#' Returns all subsets of a certain size
#'
#' For an input vector x, returns in a list, the collection of all subsets
#' of x of size k.
#'
#' @param x a vector from which to get subsets
#' @param k the size of the subsets returned
#'
#' @return a list of all subsets of x of a given size k
subsetsOfSize <- function(x, k) {
    if (k > length(x)) {
        return(list())
    }
    if (k == 0) {
        return(list(numeric(0)))
    }
    if (length(x) == k) {
        return(list(x))
    }
    return(combn(x, k, simplify = F))
}

#' Create an identifier base case
#'
#' Identifiers are functions that take as input a covariance matrix Sigma
#' corresponding to some mixed graph G and, from that covariance matrix,
#' identify some subset of the coefficients in the mixed graph G. This function
#' takes as input the matrices, L and O, defining G and creates an identifier
#' that does not identify any of the coefficients of G. This is useful as a
#' base case when building more complex identification functions.
#'
#' @inheritParams graphID
#'
#' @return a function that takes as input a covariance matrix compatible with
#'         the mixed graph defined by L/O and returns a list with two
#'         named components:
#'         Lambda - a matrix equal to L but with NA values instead of 1s,
#'         Omega - a matrix equal to O but with NA values instead of 1s.
#'         When building more complex identifiers these NAs will be replaced
#'         by the value that can be identified from Sigma.
createIdentifierBaseCase <- function(L, O) {
    # Redundant assignment puts L into the environment of the, below returned,
    # function
    validateMatrices(L, O)
    L <- L
    O <- O
    return(function(Sigma) {
        Lambda <- matrix(NA, nrow(L), nrow(L))
        Lambda[L == 0] <- 0
        Omega <- matrix(NA, nrow(O), nrow(O))
        Omega[(O == 0) & !(diag(nrow(O)))] <- 0
        return(list(Lambda = Lambda, Omega = Omega))
    })
}
