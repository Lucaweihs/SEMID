## Empty graph
L = t(matrix(
  c(0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex. 3a HTC id
L = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 1, 1, 0,
    0, 0, 0, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 3b HTC id

L = t(matrix(
  c(0, 1, 0, 0, 1,
    0, 0, 1, 1, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 1, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 3c HTC inc
L = t(matrix(
  c(0, 1, 1, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 1, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 3d HTC id
L = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 1, 0, 1, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 0, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")

## Ex 3e HTC inc
L = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 1,
    0, 0, 0, 0, 1,
    1, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 1, 0, 0,
    0, 0, 1, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 4a HTC nonid
L = t(matrix(
  c(0, 1, 1, 1, 1,
    0, 0, 1, 1, 0,
    0, 0, 0, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 0, 0, 0,
    0, 0, 1, 0, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 4b HTC inc
L = t(matrix(
  c(0, 1, 1, 1, 1,
    0, 0, 1, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 1, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 4c HTC nonid
L = t(matrix(
  c(0, 1, 1, 0, 0,
    1, 0, 0, 0, 0,
    1, 0, 0, 0, 0,
    1, 0, 0, 0, 0,
    1, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 4d HTC inc
L = t(matrix(
  c(0, 0, 0, 1, 1,
    0, 0, 1, 0, 1,
    1, 1, 0, 1, 0,
    0, 0, 0, 0, 1,
    0, 1, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 5a HTC inc
L = t(matrix(
  c(0, 1, 1, 1, 1,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 1, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 5b HTC inc
L = t(matrix(
  c(0, 1, 0, 1, 1,
    0, 0, 1, 0, 1,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 1, 0, 1,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 5c HTC inc
L = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1,
    1, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 5d HTC inc
L = t(matrix(
  c(0, 1, 0, 1, 1,
    0, 0, 1, 0, 1,
    1, 0, 0, 1, 0,
    0, 1, 0, 0, 0,
    1, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
readline("Press any key to continue to next example")


## Ex 7a (Fig 9a)
L = t(matrix(
  c(0, 1, 1, 0, 0,
    0, 0, 1, 1, 1,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 0, 1, 0,
    0, 0, 1, 0, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5)); O = O + t(O)
g = MixedGraph(L,O) # Create graph
g$plot() # Plot grapth
semID(g) # Test for global/generic identifiability
