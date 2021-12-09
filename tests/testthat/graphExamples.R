source("helperFunctions.R")
graphExamples <- list()

## Empty graph
L <- t(matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 1, genId = 1,
    htcId = 1)))

## Verma graph
L <- t(matrix(c(0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 4, 4))
O <- t(matrix(c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 4, 4))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 1, genId = 1,
    htcId = 1)))


## Ex. 3a HTC id
L <- t(matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = 1)))


## Ex 3b HTC id
L <- t(matrix(c(0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = 1)))


## Ex 3c HTC inc
L <- t(matrix(c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = -1)))


## Ex 3d HTC id
L <- t(matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = 1)))

## Ex 3e HTC inc
L <- t(matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = -1)))

## Ex 4a HTC nonid
L <- t(matrix(c(0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 0,
    htcId = 0)))


## Ex 4b HTC inc
L <- t(matrix(c(0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 0,
    htcId = -1)))


## Ex 4c HTC nonid
L <- t(matrix(c(0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 0,
    htcId = 0)))


## Ex 4d HTC inc
L <- t(matrix(c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 0,
    htcId = -1)))


## Ex 5a HTC inc
L <- t(matrix(c(0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 0,
    htcId = -1)))


## Ex 5b HTC inc
L <- t(matrix(c(0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 0,
    htcId = -1)))


## Ex 5c HTC inc
L <- t(matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 0,
    htcId = -1)))

## Ex 5d HTC inc

L <- t(matrix(c(0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 0,
    htcId = -1)))


## Ex 7a (Fig 9a)
L <- t(matrix(c(0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0), 5, 5))
O <- t(matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0), 5, 5))
O <- O + t(O)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = -1, tianId = 1)))

## Sink node marginalization example from ancestral decomposition paper
dG <- igraph::graph.edgelist(matrix(c(1, 2, 1, 3, 1, 6, 2, 3, 2, 4, 2, 5, 2, 6, 3,
    4, 4, 5), ncol = 2, byrow = T))
bG <- igraph::graph.edgelist(matrix(c(1, 6, 1, 4, 2, 3, 2, 5, 2, 6), ncol = 2, byrow = T),
    directed = F)
L <- as.matrix(igraph::get.adjacency(dG))
O <- as.matrix(igraph::get.adjacency(bG))
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = -1, tianId = -1, ancId = 1)))

## Simple 2 node unidentifiable
L <- t(matrix(c(0, 1, 0, 0), 2, 2))
O <- t(matrix(c(0, 1, 1, 0), 2, 2))
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 0,
    htcId = 0)))

## Instrumental variable model
L <- t(matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3))
O <- t(matrix(c(0, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3))
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = 1)))

## Several examples where ancestor decomposition was seen to be helpful
n <- 6
p <- 0.1
p1 <- 0.2
set.seed(176796)
O <- rConnectedAdjMatrix(n, p)
L <- rAcyclicDirectedAdjMatrix(n, p1)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = -1, ancId = 1)))
set.seed(335911)
O <- rConnectedAdjMatrix(n, p)
L <- rAcyclicDirectedAdjMatrix(n, p1)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = -1, ancId = 1)))

set.seed(762097)
O <- rConnectedAdjMatrix(n, p)
L <- rAcyclicDirectedAdjMatrix(n, p1)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = -1, ancId = 1)))

n <- 8
p <- 0.3
p1 <- 0.4
set.seed(501)
O <- rConnectedAdjMatrix(n, p)
L <- rAcyclicDirectedAdjMatrix(n, p1)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = -1, ancId = 1)))

n <- 10
p <- 0.2
p1 <- 0.5
set.seed(3178)
O <- rConnectedAdjMatrix(n, p)
L <- rAcyclicDirectedAdjMatrix(n, p1)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = 1,
    htcId = -1, ancId = 1)))

set.seed(5536)
O <- rConnectedAdjMatrix(n, p)
L <- rAcyclicDirectedAdjMatrix(n, p1)
graphExamples <- c(graphExamples, list(list(L = L, O = O, globalId = 0, genId = -1,
    htcId = -1, ancId = -1)))


#############################################################
## latent digraphs where all latent nodes are source nodes ##
#############################################################
digraphExamples <- list()

## Insturmental variables
L = matrix(c(0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 0,
             0, 1, 1, 0), 4, 4, byrow=TRUE)
observedNodes = seq(1,3)
latentNodes = c(4)
g = LatentDigraph(L, observedNodes, latentNodes)
digraphExamples <- c(digraphExamples, list(list(graph = g, lfhtcid = TRUE)))

## One global latent variable
L = matrix(c(0, 1, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 0,
             0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0,
             1, 1, 1, 1, 1, 0), 6, 6, byrow=TRUE)
observedNodes = seq(1,5)
latentNodes = c(6)
g = LatentDigraph(L, observedNodes, latentNodes)
digraphExamples <- c(digraphExamples, list(list(graph = g, lfhtcid = TRUE)))

## Figure 1
L = matrix(c(0, 0, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 0,
             0, 0, 0, 1, 1, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0,
             1, 1, 1, 1, 1, 0), 6, 6, byrow=TRUE)
observedNodes = seq(1,5)
latentNodes = c(6)
g = LatentDigraph(L, observedNodes, latentNodes)
digraphExamples <- c(digraphExamples, list(list(graph = g, lfhtcid = TRUE)))


## Figure 2
L = matrix(c(0, 1, 1, 0, 0, 0,
             0, 0, 1, 0, 0, 0,
             0, 0, 0, 1, 1, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0,
             1, 1, 1, 1, 1, 0), 6, 6, byrow=TRUE)
observedNodes = seq(1,5)
latentNodes = c(6)
g = LatentDigraph(L, observedNodes, latentNodes)
digraphExamples <- c(digraphExamples, list(list(graph = g, lfhtcid = FALSE)))


## Figure 3 a)
L = matrix(c(0, 0, 1, 0, 0, 0,
             0, 0, 1, 0, 0, 0,
             0, 0, 0, 1, 0, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0,
             1, 1, 1, 1, 1, 0), 6, 6, byrow=TRUE)
observedNodes = seq(1,5)
latentNodes = c(6)
g = LatentDigraph(L, observedNodes, latentNodes)
digraphExamples <- c(digraphExamples, list(list(graph = g, lfhtcid = FALSE)))


## Figure 3 b)
L = matrix(c(0, 1, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 0,
             0, 0, 0, 1, 0, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0,
             1, 1, 1, 1, 1, 0), 6, 6, byrow=TRUE)
observedNodes = seq(1,5)
latentNodes = c(6)
g = LatentDigraph(L, observedNodes, latentNodes)
digraphExamples <- c(digraphExamples, list(list(graph = g, lfhtcid = FALSE)))


## Figure 5
L = matrix(c(0, 0, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 0,
             0, 0, 0, 1, 1, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0,
             1, 0, 1, 1, 0, 0), 6, 6, byrow=TRUE)
observedNodes = seq(1,5)
latentNodes = c(6)
g = LatentDigraph(L, observedNodes, latentNodes)
digraphExamples <- c(digraphExamples, list(list(graph = g, lfhtcid = TRUE)))


## Figure 6
L = matrix(c(0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 1, 0, 0, 0,
             0, 0, 0, 0, 1, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,
             1, 1, 1, 1, 0, 0, 0, 0,
             0, 0, 0, 1, 1, 0, 0, 0,
             0, 0, 1, 0, 1, 0, 0, 0), 8, 8, byrow=TRUE)
observedNodes = seq(1,5)
latentNodes = seq(6,8)
g = LatentDigraph(L, observedNodes, latentNodes)
digraphExamples <- c(digraphExamples, list(list(graph = g, lfhtcid = FALSE)))


## Figure 10
L = matrix(c(0, 1, 0, 0, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 1, 1, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,
             1, 1, 1, 1, 0, 0, 0, 0,
             0, 0, 0, 1, 1, 1, 0, 0), 8, 8, byrow=TRUE)
observedNodes = seq(1,6)
latentNodes = seq(7,8)
g = LatentDigraph(L, observedNodes, latentNodes)
digraphExamples <- c(digraphExamples, list(list(graph = g, lfhtcid = TRUE)))


## Not enough to iterate over L subset of pa_V(v)
L = matrix(c(0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0,
             1, 1, 1, 1, 1, 1, 0, 0,
             1, 0, 0, 1, 1, 1, 0, 0), 8, 8, byrow=TRUE)
observedNodes = seq(1,6)
latentNodes = seq(7,8)
g = LatentDigraph(L, observedNodes, latentNodes)
digraphExamples <- c(digraphExamples, list(list(graph = g, lfhtcid = TRUE)))