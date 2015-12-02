###
# Wrapper function --- this is the one that gets called
###
graphID <- function(L, O, output.type = 'matrix', file.name = NULL, decomp.if.acyclic = TRUE, test.globalID = TRUE, test.genericID = TRUE, test.nonID = TRUE){
  ## L = directed edge matrix
  ## O = bidirected edge matrix

  if(!is.matrix(L) || !is.matrix(O)){
    print('L and O must be two square matrices of the same size'); return(NULL)
  }else{
    m=unique(c(dim(L),dim(O)))
    if(length(m) > 1){
      print('L and O must be two square matrices of the same size'); return(NULL)
    }else{
      if(!is.character(file.name) && output.type == 'file'){
        print('need to input a file name for output.type = \'file\''); return(NULL)
      }else{
        L <- (L != 0) ; diag(L) <- 0
        O <- O + t(O) ; O <- (O != 0) ; diag(O) <- 0
        Graph.output <- graphID.decompose(L, O, decomp.if.acyclic, test.globalID, test.genericID, test.nonID)
        if (all(c(test.globalID, test.genericID, test.nonID))){
          graphID.output.alltests(Graph.output$Components, Graph.output$Decomp, output.type, file.name)
        }else{
          graphID.output.notalltests(Graph.output$Components, Graph.output$Decomp, output.type, file.name, test.globalID, test.genericID, test.nonID)
        }

      }
    }
  }
}

###
# Function to format the output when all tests have been run.
###
graphID.output.alltests <- function(Graph.components, Decomp, output.type, file.name){
  if(output.type == 'list'){
    return(Graph.components)
  }else{
    Graph.output.matrix <- NULL
    for(comp in 1:length(Graph.components)){
      if(Graph.components[[comp]]$acyclic){
        if(Graph.components[[comp]]$GlobalID) {Comp.status <- 'Globally ID\'able'}else{
          if(Graph.components[[comp]]$GenericID) {Comp.status <- 'Generically ID\'able'}else{
            if(Graph.components[[comp]]$NonID) {Comp.status <- 'Generically non-ID\'able'}else{
              Comp.status <- 'HTC-inconclusive'
            }
          }
        }
      }else{
        if(Graph.components[[comp]]$GenericID) {Comp.status <- 'Generically ID\'able'}else{
          if(Graph.components[[comp]]$NonID) {Comp.status <- 'Generically non-ID\'able'}else{
            Comp.status <- 'HTC-inconclusive'
          }
        }
      }

      if (length(Graph.components[[comp]]$HTC.ID.nodes) == 0){
        Graph.components[[comp]]$HTC.ID.nodes <- 'none'
      }

      Comp.row <- c(comp,
                    paste(Graph.components[[comp]]$Nodes,collapse=','),
                    paste(Graph.components[[comp]]$HTC.ID.nodes,collapse=','),							Comp.status)
      Graph.output.matrix <- rbind(Graph.output.matrix, Comp.row)
    }

    Graph.output.matrix <- rbind(c('Component', 'Nodes', 'HTC-ID\'able nodes', 'Identifiability'), Graph.output.matrix)
    colnames(Graph.output.matrix) <- rep('',4)
    rownames(Graph.output.matrix) <- rep('',length(Graph.components)+1)

    if(!Decomp){
      Graph.output.matrix <- Graph.output.matrix[, -(1:2)]
    }

    Max.string <- max(nchar(Graph.output.matrix))
    for(i in 1:dim(Graph.output.matrix)[1]){
      for(j in 1:dim(Graph.output.matrix)[2]){
        Graph.output.matrix[i,j] <- paste(c(Graph.output.matrix[i,j], rep(' ', 5 + Max.string - nchar(Graph.output.matrix[i,j]))), collapse='')
      }
    }

    if(output.type=='file'){
      write.table(Graph.output.matrix, file=file.name, quote=FALSE)
    }else{
      print(Graph.output.matrix, quote=FALSE)
    }
  }
}

###
# Function to format the output when not all tests have been run
###
graphID.output.notalltests <- function(Graph.components, Decomp, output.type, file.name, test.globalID, test.genericID, test.nonID){
  if(output.type == 'list'){
    return(Graph.components)
  }else{
    Graph.output.matrix <- NULL
    for(comp in 1:length(Graph.components)){
      if(Graph.components[[comp]]$acyclic && test.globalID){
        if(Graph.components[[comp]]$GlobalID){
          Comp.status.GlobalID <- 'yes'
        }else{
          Comp.status.GlobalID <- 'no'
        }
      }else{
        Comp.status.GlobalID <- 'not tested'
      }
      if(test.genericID){
        if(Graph.components[[comp]]$GenericID){
          Comp.status.GenericID <- 'yes'
        }else{
          Comp.status.GenericID <- 'no'
        }
      }else{
        Comp.status.GenericID <- 'not tested'
      }
      if(test.nonID){
        if(Graph.components[[comp]]$NonID){
          Comp.status.NonID <- 'yes'
        }else{
          Comp.status.NonID <- 'no'
        }
      }else{
        Comp.status.NonID <- 'not tested'
      }


      if (length(Graph.components[[comp]]$HTC.ID.nodes) == 0){
        if(test.genericID){
          Graph.components[[comp]]$HTC.ID.nodes <- 'none'
        }else{
          Graph.components[[comp]]$HTC.ID.nodes <- 'not tested'
        }
      }


      Comp.row <- c(comp,
                    paste(Graph.components[[comp]]$Nodes,collapse=','),
                    paste(Graph.components[[comp]]$HTC.ID.nodes,collapse=','), Comp.status.GlobalID, Comp.status.GenericID, Comp.status.NonID)
      Graph.output.matrix <- rbind(Graph.output.matrix, Comp.row)
    }

    Graph.output.matrix <- rbind(c('Component', 'Nodes', 'HTC-ID\'able nodes', 'Globally ID\'able?', 'Generically ID\'able?', 'Generically non-ID\'able?'), Graph.output.matrix)
    colnames(Graph.output.matrix) <- rep('',6)
    rownames(Graph.output.matrix) <- rep('',length(Graph.components)+1)

    if(!Decomp){
      Graph.output.matrix <- Graph.output.matrix[, -(1:2)]
    }

    Max.string <- max(nchar(Graph.output.matrix))
    for(i in 1:dim(Graph.output.matrix)[1]){
      for(j in 1:dim(Graph.output.matrix)[2]){
        Graph.output.matrix[i,j] <- paste(c(Graph.output.matrix[i,j], rep(' ', 3 + Max.string - nchar(Graph.output.matrix[i,j]))), collapse='')
      }
    }


    if(output.type=='file'){
      write.table(Graph.output.matrix, file=file.name, quote=FALSE)
    }else{
      print(Graph.output.matrix, quote=FALSE)
    }
  }
}

###
# Tian decomposition: Split a graph into bidirected connected components & then
# solve each separately
###
graphID.decompose <- function(L, O, decomp.if.acyclic = TRUE, test.globalID = TRUE, test.genericID = TRUE, test.nonID = TRUE){
  m <- dim(L)[1]
  L <- (L != 0) ; diag(L) <- 0
  O <- O + t(O) ; O <- (O != 0) ; diag(O) <- 0

  Decomp <- FALSE
  if (decomp.if.acyclic) {
    test.acyclic <- TRUE
    nodes.acyclic <- 1:m
    while (test.acyclic && !Decomp) {
      if (any(colSums(L[nodes.acyclic, nodes.acyclic, drop=FALSE]) == 0)) {
        nodes.acyclic <- nodes.acyclic[-which(colSums(L[nodes.acyclic, nodes.acyclic, drop=FALSE]) == 0)]
        if (length(nodes.acyclic) == 0) {
          Decomp <- TRUE
        }
      }else{
        test.acyclic <- FALSE
      }
    }
  }

  if (Decomp) {
    Bidir.comp <- diag(m) + O
    for(i in 1:(m-1)){
      Bidir.comp <- (Bidir.comp %*% (diag(m) + O)>0) + 0
    }
    Components <- list()
    V <- 1:m
    num.comp=0
    while (length(V) > 0) {
      num.comp <- num.comp + 1
      i <- min(V)
      Components[[num.comp]] <- which(Bidir.comp[i,]==1)
      V <- setdiff(V,Components[[num.comp]])
    }
  }else{
    Components <- list()
    Components[[1]] <- 1:m
    num.comp <- 1
  }

  for (comp in 1:num.comp) {
    Component <- Components[[comp]]
    Component.parents <- sort(setdiff(which(rowSums(L[,Component,drop=FALSE])>0),Component))
    if(length(Components[[comp]]) == 1){
      Components[[comp]] <- list()
      Components[[comp]]$Nodes <- Component
      Components[[comp]]$acyclic <- TRUE
      if (test.genericID){
        Components[[comp]]$HTC.ID.nodes <- Component
        Components[[comp]]$GenericID <- TRUE
      }
      if (test.globalID){
        Components[[comp]]$GlobalID <- TRUE
      }
      if(test.nonID){
        Components[[comp]]$NonID <- FALSE
      }
    }else{
      m1 <- length(Component) ; m2 <- length(Component.parents)
      L.Component <- cbind( L[c(Component,Component.parents),Component], matrix(0, m1 + m2, m2) )
      O.Component <- cbind( rbind( O[Component,Component], matrix(0, m2, m1) ), matrix(0, m1 + m2, m2) )
      Components[[comp]] <- c(list(Nodes=Component),graphID.main(L.Component, O.Component, test.globalID, test.genericID, test.nonID))
      Components[[comp]]$HTC.ID.nodes <- Component[intersect(Components[[comp]]$HTC.ID.nodes,1:length(Component))]
    }
  }
  Graph.output=list()
  Graph.output$Components <- Components
  Graph.output$Decomp <- Decomp
  return(Graph.output)
}

###
# Main function to handle a graph component: calls the other functions that
# determine identifiability status
###
graphID.main <- function(L, O, test.globalID = TRUE, test.genericID = TRUE, test.nonID = TRUE){
  m <- dim(L)[1]
  L <- (L != 0) ; diag(L) <- 0
  O <- O + t(O) ; O <- (O != 0) ; diag(O) <- 0

  ILinv <- diag(m)
  for (i in 1:m){
    ILinv <- 0 + (diag(m) + ILinv %*% L>0)
  }

  Output <- list()

  Output$acyclic <- (max(ILinv + t(ILinv) - diag(m)) == 1)


  if (test.genericID){

    Output$HTC.ID.nodes <- graphID.genericID(L,O)

    if (length(Output$HTC.ID.nodes) < m){
      Output$GenericID <- FALSE
      if (Output$acyclic && test.globalID){
        Output$GlobalID <- FALSE
      }
      if (test.nonID){
        Output$NonID <- graphID.nonID(L,O)
      }
    }else{
      Output$GenericID <- TRUE
      if (Output$acyclic && test.globalID){
        Output$GlobalID <- graphID.globalID(L,O)
      }
      if (test.nonID){
        Output$NonID <- FALSE
      }
    }
  }else{
    if (test.nonID){
      Output$NonID <- graphID.nonID(L,O)
      if (!Output$NonID){
        if (Output$acyclic && test.globalID){
          Output$GlobalID <- graphID.globalID(L,O)
        }
      }
    }else{
      if (Output$acyclic && test.globalID){
        Output$GlobalID <- graphID.globalID(L,O)
      }
    }
  }
  return(Output)
}

###
# Check for global identifiability (see arXiv:1003:1146).
###
graphID.globalID <- function(L, O){  # for acyclic graphs only
  m <- dim(L)[1]
  L <- (L != 0) ; diag(L) <- 0
  O <- O + t(O) ; O <- (O != 0) ; diag(O) <- 0

  Global.ID=TRUE
  i <- 0
  while (Global.ID == 1 && i < m){
    i <- i+1
    S <- 1:m
    change <- 1
    while (change == 1){
      change <- 0
      S.old <- S
      for (s in setdiff(S.old, i)){
        if (graph.maxflow(graph.adjacency(L), source=s, target=i)$value == 0){
          S <- setdiff(S, s)
          change <- 1
        }
      }
      S.old <- S
      for (s in setdiff(S.old, i)){
        if (graph.maxflow(graph.adjacency(O, mode='undirected'), source=s, target=i)$value == 0){
          S <- setdiff(S, s)
          change <- 1
        }
      }
    }
    if (length(S) > 1){
      Global.ID <- FALSE
    }
  }

  return(Global.ID)
}

###
# Check for generic infinite-to-one via the half-trek criterion
# (see arXiv:1107:5552).
###
graphID.nonID <- function(L, O){
  m <- dim(L)[1]
  L <- (L != 0) ; diag(L) <- 0
  O <- O + t(O) ; O <- (O != 0) ; diag(O) <- 0

                                        # 1 & 2 = source & target
                                        # 2 + {1,...,N} = L{i_n,j_n} for the n-th nonsibling pair, n=1,...,N
                                        # (2+N) + 2*m^2 = 2 copies of R_i(j) -- in copy & out copy
                                        #         where (2+N) + (i-1)*m + j    = R_i(j) in copy
                                        #         & (2+N+m^2) + (i-1)*m + j    = R_i(j) out copy

  nonsibs <- NULL; N <- 0
  for (i in 1:(m-1)) { for (j in (i+1):m) {
    if (O[i,j] == 0) {
      N <- N+1
      nonsibs <- rbind(nonsibs, c(i, j))
    }
  }}

  Cap.matrix <- matrix(0,2*m^2+N+2,2*m^2+N+2)
  Cap.matrix[1, 2+(1:N)] <- 1 # edges from source to L{i,j} for each nonsibling pair
  for (n in 1:N) {  #{i,j} = nonsibs[n,1:2]
                                        # edge from L{i,j} to R_i(j), and to R_i(k)-in for all siblings k of node j
    Cap.matrix[2+n, 2+N + (nonsibs[n,1]-1)*m + c(nonsibs[n,2], which(O[nonsibs[n,2], ] == 1))] <- 1
                                        # edge from L{i,j} to R_j(i), and to R_j(i)-in for all siblings k of node i
    Cap.matrix[2+n, 2+N + (nonsibs[n,2]-1)*m + c(nonsibs[n,1], which(O[nonsibs[n,1], ] == 1))] <- 1
  }
  for (i in 1:m) {
                                        # edge from R_i(j)-out to target when j is a parent of i
    Cap.matrix[2+N+m^2 + (i-1)*m + which(L[, i] == 1), 2] <- 1
    for (j in 1:m) {
                                        # edge from R_i(j)-in to R_i(j)-out
      Cap.matrix[2+N + (i-1)*m + j,2+N+m^2 + (i-1)*m + j] <- 1
                                        # edge from R_i(j)-out to R_i(k)-in where j->k is a directed edge
      Cap.matrix[2+N+m^2 + (i-1)*m + j,2+N + (i-1)*m + which(L[j, ] == 1)] <- 1
    }}



  HTC.nonID <- (graph.maxflow( graph.adjacency(Cap.matrix), source=1, target=2)$value < sum(L) )

  return( HTC.nonID )

}




###
# If directed part of input graph is cyclic then will check for generic
# identifiability using the half-trek criterion (see arXiv:1107:5552). Otherwise
# will use the a slightly stronger version of the half-trek criterion using
# ancestor decompositions (see arXiv:1504.02992).
###
graphID.genericID <- function(L, O, ILinv=0){
  m <- nrow(L)
  L <- (L != 0)
  diag(L) <- 0
  if (is.dag(graph.adjacency(L, mode="directed"))) {
    return(graphID.ancestral(L, O))
  }
  O <- O + t(O)
  O <- (O != 0)
  diag(O) <- 0

                                        # 1 & 2 = source & target
                                        # 2 + {1,...,m} = L(i) for i=1,...,m
                                        # 2+m + {1,...,m} = R(i)-in for i=1,...,m
                                        # 2+2*m + {1,...,m} = R(i)-out for i=1,...,m


  Cap.matrix.init <- matrix(0, 2+3*m, 2+3*m)

  for (i in 1:m){
                                        # edge from L(i) to R(i)-in, and to R(j)-in for all siblings j of i
    Cap.matrix.init[2+i, 2+m+c(i, which(O[i,] == 1))] <- 1
                                        # edge from R(i)-in to R(i)-out
    Cap.matrix.init[2+m+i, 2+2*m+i] <- 1
                                        # edge from R(i)-out to R(j)-in for all directed edges i->j
    Cap.matrix.init[2+2*m+i, 2+m + which(L[i,] == 1)] <- 1
  }

                                        # when testing if a set A satisfies the HTC with respect to a node i,
                                        #    need to add (1) edge from source to L(j) for all j in A
                                        #            and (2) edge from R(j)-out to target for all parents j of i

  Dependence.matrix <- O + diag(m)
  if(is.matrix(ILinv)){
    Dependence.matrix <- Dependence.matrix %*% ILinv
  }else{
    for (i in 1:m){
      Dependence.matrix <- (Dependence.matrix + Dependence.matrix %*% L > 0)
    }
  }

  Solved.nodes <- rep(0,m)
  Solved.nodes[which(colSums(L) == 0)] <- 1 # nodes with no parents
  change <- 1
  count <- 1

  while (change == 1){
    change <- 0

    for (i in which(Solved.nodes == 0)){
      A <- setdiff(c(which(Solved.nodes > 0), which(Dependence.matrix[i, ] == 0)), c(i, which(O[i,] == 1)))

      Cap.matrix <- Cap.matrix.init

      Cap.matrix[1, 2+A] <- 1
      Cap.matrix[2+2*m + which(L[,i] == 1), 2] <- 1

      flow <- graph.maxflow (graph.adjacency(Cap.matrix), source=1, target=2)$value

      if (flow == sum(L[,i])){
        change <- 1
        count <- count+1
        Solved.nodes[i] <- count
      }

    }

  }

  if(all(Solved.nodes == 0)){
    Solved.nodes <- NULL
  }else{
    Solved.nodes <- order(Solved.nodes)[(1+sum(Solved.nodes == 0)):m]
  }

  return(Solved.nodes)
}

###
# Get the ancestors of a collection of nodes in a graph g, the ancestors DO
# include the the nodes themselves.
###
ancestors <- function(g, nodes) {
  if (length(V(g)) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  sort(unique(unlist(neighborhood(g, length(V(g)), nodes=nodes, mode="in"))))
}

###
# Get the parents of a collection of nodes in a graph g, the parents DO include
# the input nodes themselves.
###
parents <- function(g, nodes) {
  if (length(V(g)) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  sort(unique(unlist(neighborhood(g, 1, nodes=nodes, mode="in"))))
}

###
# Get the siblings of a collection of nodes in a graph g, the siblings DO
# include the input nodes themselves.
###
siblings = function(g, nodes) {
  if (length(V(g)) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  sort(unique(unlist(neighborhood(g, 1, nodes=nodes, mode="all"))))
}

###
# For an input mixed graph H and set of nodes A, let GA be the subgraph of
# H on the nodes A. This function returns the mixed component of GA containing
# a specified node.
#
# @param dG a directed graph representing the directed part of the mixed graph.
#        Here dG is expected to have vertices with names 1,...,m where m is the
#        number of vertices in the graph.
# @param bG an undirected graph representing the undirected part of the mixed
#        graph. The names of vertices in bG should correspond to those in dG
# @param subNodes an ancestral set of nodes in the mixed graph, this set should
#        include the special node nodeName
# @param nodeName the node for which the mixed component is found.
# @returns a list with two named elements:
#          biNodes - the nodes of the mixed graph in the biDirected component
#                    containing nodeName w.r.t the ancestral set of nodes
#          inNodes - the nodes in the graph which are not part of biNodes
#                    but which are a parent of some node in biNodes.
###
getMixedCompForNode <- function(dG, bG, subNodes, nodeName) {
  m = length(V(dG))
  if (is.null(V(dG)$names) || is.null(V(bG)) ||
      any(V(dG)$names != 1:m) || any(V(bG)$names != 1:m)) {
    stop("Input graphs to getMixedCompForNode must have vertices named 1:m in order.")
  }

  nodesToDelete = setdiff(1:m, subNodes)
  bG = delete.vertices(bG, nodesToDelete)

  nodeNum = which(V(bG)$names == nodeName)
  clustMember = clusters(bG)$membership
  nodeMember = clustMember[nodeNum]

  bidirectedComp = which(clustMember == nodeMember)
  bidirectedComp = (V(bG)$names)[bidirectedComp]

  incomingNodes = intersect(setdiff(parents(dG, bidirectedComp),
                                    bidirectedComp),
                            subNodes)

  return(list(biNodes=bidirectedComp, inNodes=incomingNodes))
}

###
# For an input mixed graph H, constructs the Gflow graph as described in
# arXiv:1107:5552 for a subgraph G of H. A max flow algorithm is then run on
# Gflow to determine the largest half-trek system in G to a particular node's
# parents given a set of allowed nodes. Here G should consist of a bidirected
# part and nodes which are not in the bidirected part but are a parent of some
# node in the bidirected part. G should contain the node for which to compute
# the max flow.
#
# @param L an adjacency matrix representing the directed part of the mixed graph.
# @param O an adjacency matrix representing the bidirected part of the mixed graph.
# @param allowedNodes the set of allowed nodes.
# @param biNodes the set of nodes in the subgraph G which are part of the
#        bidirected part.
# @param inNodes the nodes of the subgraph G which are not in the bidirected
#        part but are a parent of some node in the bidirected component.
# @param node the node (as an integer) for which the maxflow should be computed.
###
getMaxFlow = function(L, O, allowedNodes, biNodes, inNodes, node) {
  if (!(node %in% biNodes) || any(O[biNodes, inNodes] != 0)) {
    stop("When getting max flow either some in-nodes were connected to some bi-nodes or node was not in biNodes.")
  }
  m = length(biNodes) + length(inNodes)

  oldNumToNewNum = numeric(m)
  oldNumToNewNum[biNodes] = 1:length(biNodes)
  oldNumToNewNum[inNodes] = (length(biNodes)+1):m

  nodesToUse = c(biNodes, inNodes)
  allowedNodes = intersect(allowedNodes, nodesToUse)
  L[c(inNodes, biNodes), inNodes] = 0
  O[inNodes, inNodes] = 0
  L = L[nodesToUse, nodesToUse]
  O = O[nodesToUse, nodesToUse]

  # 1 & 2 = source & target
  # 2 + {1,...,m} = L(i) for i=1,...,m
  # 2 + m + {1,...,m} = R(i)-in for i=1,...,m
  # 2 + 2*m + {1,...,m} = R(i)-out for i=1,...,m

  Cap.matrix <- matrix(0, 2 + 3*m, 2 + 3*m)

  for (i in 1:m){
    # edge from L(i) to R(i)-in, and to R(j)-in for all siblings j of i
    Cap.matrix[2 + i, 2 + m + c(i, which(O[i,] == 1))] <- 1
    # edge from R(i)-in to R(i)-out
    Cap.matrix[2 + m + i, 2 + 2*m + i] <- 1
    # edge from R(i)-out to R(j)-in for all directed edges i->j
    Cap.matrix[2 + 2*m + i, 2 + m + which(L[i,] == 1)] <- 1
  }

  allowedNodes = oldNumToNewNum[allowedNodes]
  node = oldNumToNewNum[node]
  Cap.matrix[1, 2 + allowedNodes] = 1
  Cap.matrix[2 + 2*m + which(L[,node] == 1), 2] = 1

  return(graph.maxflow(graph.adjacency(Cap.matrix), source=1, target=2)$value)
}

###
# For an input, acyclic, mixed graph attempts to determine if the graph is
# generically identifiable using ancestral subsets (see arXiv:1504.02992).
###
graphID.ancestral <- function(L, O) {
  m <- nrow(L)
  L <- (L != 0)
  diag(L) <- 0
  O <- (O + t(O)) != 0
  diag(O) <- 0

  dG = graph.adjacency(L)
  newOrder = topological.sort(dG)
  L = L[newOrder,newOrder]
  O = O[newOrder,newOrder]

  dG = graph.adjacency(L)
  bG = graph.adjacency(O)
  V(dG)$names = 1:m
  V(bG)$names = 1:m

  # Generates a list where, for each node v, we have a vector
  # corresponding to all the nodes that could ever be in a
  # half-trek system for v
  halfTrekSources = vector("list", length=m)
  for(i in 1:m) {
    halfTrekSources[[i]] = siblings(bG, ancestors(dG, i))
  }

  # A matrix determining which nodes are htr from each node
  Dependence.matrix <- O + diag(m)
  for (i in 1:m) {
    Dependence.matrix <- ((Dependence.matrix + Dependence.matrix %*% L) > 0)
  }

  Solved.nodes <- rep(0,m)
  Solved.nodes[which(colSums(L) == 0)] <- 1 # nodes with no parents
  change <- 1
  count <- 1
  while (change == 1){
    change <- 0

    Unsolved.nodes <- which(Solved.nodes == 0)
    for (i in Unsolved.nodes) {
      # A <- (Solved Nodes 'union' nodes not htr from i) \ ({i} 'union' siblings(i))
      A = setdiff(c(which(Solved.nodes > 0), which(Dependence.matrix[i, ] == 0)),
                  c(i, which(O[i,] == 1)))
      # A <- A intersect (nodes that can ever be in a HT system for i)
      A = intersect(A, halfTrekSources[[i]])

      mixedCompList = getMixedCompForNode(dG, bG, ancestors(dG, c(i,A)), i)
      flow = getMaxFlow(L, O, A, mixedCompList$biNodes, mixedCompList$inNodes, i)
      if (flow == sum(L[,i])){
        change <- 1
        count <- count + 1
        Solved.nodes[i] <- count
        next
      }

      mixedCompList = getMixedCompForNode(dG, bG, ancestors(dG, i), i)
      A = union(intersect(A, unlist(mixedCompList)), mixedCompList$inNodes)
      flow = getMaxFlow(L, O, A, mixedCompList$biNodes, mixedCompList$inNodes, i)
      if (flow == sum(L[,i])){
        change <- 1
        count <- count + 1
        Solved.nodes[i] <- count
        next
      }
    }
  }

  if(all(Solved.nodes == 0)){
    Solved.nodes <- NULL
  }else{
    Solved.nodes <- order(Solved.nodes)[(1+sum(Solved.nodes == 0)):m]
  }

  return(newOrder[Solved.nodes])
}
