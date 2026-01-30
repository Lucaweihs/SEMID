---
editor_options: 
  markdown: 
    wrap: 72
---

# SEMID Package

## Purpose

This package offers a number of functions for determining parameter
identifiability in different classes of linear structural equation
models (SEMs) with latent variables. Each model is defined by a directed
graph or by a mixed graph, depending on the modeling assumptions. The
following sections highlight the primary ways in which the package can
be used.

## Linear SEMs given by Mixed Graphs

In the `SEMID` package, we represent mixed graphs via the MixedGraph
class.

```         
> # Mixed graphs are specified by their directed adjacency matrix L and
> # bidirected adjacency matrix O.
> library(SEMID)
> L = t(matrix(
+ c(0, 1, 0, 0, 0,
+   0, 0, 0, 1, 1,
+   0, 0, 0, 1, 0,
+   0, 1, 0, 0, 1,
+   0, 0, 0, 1, 0), 5, 5))
>
> O = t(matrix(
+ c(0, 0, 0, 0, 0,
+   0, 0, 1, 0, 1,
+   0, 0, 0, 1, 0,
+   0, 0, 0, 0, 0,
+   0, 0, 0, 0, 0), 5, 5)); O=O+t(O)
>
> # Create the mixed graph object corresponding to L and O
> g = MixedGraph(L, O)
>
> # Plot the mixed graph
> g$plot()
```

See the documentation for the MixedGraph class `?MixedGraph` for more
information.

### Global Identifiability

For deciding global identifiability in mixed graphs, there exists an ‘if
and only if’ graphical criterion developed by

Drton, M., Foygel, R., and Sullivant, S. (2011) Global identifiability
of linear structural equation models. *Ann. Statist.* 39(2): 865-886.
<https://doi.org/10.1214/10-AOS859>.

This criterion can be accessed through the function `globalID`.

```         
> # Check global identifiability
> globalID(g)
[1] FALSE
```

### Generic Identifiability

There still do not exist any ‘if and only if’ graphical conditions for
testing whether or not a mixed graph is generically identifiable.
However, there do exist sufficient and necessary conditions. The `SEMID`
package contains implementations of various sufficient conditions.

-   The half-trek criterion:

Rina Foygel, Jan Draisma, Mathias Drton (2012). Half-trek criterion for
generic identifiability of linear structural equation models. *Ann.
Statist.* 40(3):1682--1713. <https://doi.org/10.1214/12-AOS1012>.

-   Ancestor decomposition techniques:

Mathias Drton, Luca Weihs (2016). Generic Identifiability of Linear
Structural Equation Models by Ancestor Decomposition. *Scand. J.
Statist.* 43:1035--1045. <https://doi.org/10.1111/sjos.12227>.

-   Edgewise and determinantal criteria:

Luca Weih, Bill Robinson, Emilie Dufresne, Jennifer Kenkel, Kaie Kubjas,
Reginald McGee II,Nhan Nguyen, Elina Robeva, Mathias Drton (2017).
Determinantal Generalizations of Instrumental Variables. *J. Causal
Inference* 6(1). <https://doi.org/10.1515/jci-2017-0009>.

```         
> # Check generic identifiability using different criteria
> # Start with the half-trek criterion
> htcID(g)
Call: SEMID::htcID(mixedGraph = g)

Mixed Graph Info.
# nodes: 5 
# dir. edges: 7 
# bi. edges: 3 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 2 
# bi. edges shown gen. identifiable: 0 

Generically identifiable dir. edges:
2->5, 4->5 

Generically identifiable bi. edges:
None

> # Ancestor decomposition techniques:
> ancestralID(g)
Call: ancestralID(mixedGraph = g)

Mixed Graph Info.
# nodes: 5 
# dir. edges: 7 
# bi. edges: 3 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 2 
# bi. edges shown gen. identifiable: 0 

Generically identifiable dir. edges:
2->5, 4->5 

Generically identifiable bi. edges:
None

> # Edgewise identification algorithm:
> edgewiseID(g)
Call: edgewiseID(mixedGraph = g)

Mixed Graph Info.
# nodes: 5 
# dir. edges: 7 
# bi. edges: 3 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 4 
# bi. edges shown gen. identifiable: 0 

Generically identifiable dir. edges:
2->4, 5->4, 2->5, 4->5 

Generically identifiable bi. edges:
None

> # Edgewise identification algorithm leveraging trek-separation relations:
> edgewiseTSID(g)
Call: edgewiseTSID(mixedGraph = g)

Mixed Graph Info.
# nodes: 5 
# dir. edges: 7 
# bi. edges: 3 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 4 
# bi. edges shown gen. identifiable: 0 

Generically identifiable dir. edges:
2->4, 5->4, 2->5, 4->5 

Generically identifiable bi. edges:
None
```

Note that, by default, all strategies first apply a Tian decomposition
and then check identifiability on each of the components. This yields
faster computations as described in Section 8 of Foygel, Draisma, and
Drton (2012). It is also possible to apply different identification
strategies repeatedly until no further edges can be identified. This is
possible via the function `generalGenericID`.

```         
> # Check generic identifiability by repeatedly applying different criteria
> generalGenericID(mixedGraph = g, 
+                   idStepFunctions = list(htcIdentifyStep,
+                                          ancestralIdentifyStep, 
+                                          edgewiseIdentifyStep, 
+                                          trekSeparationIdentifyStep), 
+                   tianDecompose = TRUE)
Call: generalGenericID(mixedGraph = g, idStepFunctions = list(htcIdentifyStep, 
    ancestralIdentifyStep, edgewiseIdentifyStep, trekSeparationIdentifyStep), 
    tianDecompose = TRUE)

Mixed Graph Info.
# nodes: 5 
# dir. edges: 7 
# bi. edges: 3 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 4 
# bi. edges shown gen. identifiable: 0 

Generically identifiable dir. edges:
2->4, 5->4, 2->5, 4->5 

Generically identifiable bi. edges:
None
```

In this example, we do not get additional edges certified to be
generically identifiability. Therefore, we check the necessary condition
from Foygel, Draisma, and Drton (2012) for generic identifiability of
the whole graph, which is also implemented in `SEMID`.

```         
> graphID.nonHtcID(g$L(), g$O())
[1] TRUE
```

This means that the given graph is infinite-to-one and, in particular,
not generically identifiable.

## Linear SEMs given by Latent-Factor Graphs

The latent-factor half-trek criterion (LF-HTC) by Barber, Drton, Sturma
and Weihs (2022) is a sufficient criterion to check generic
identifiability in directed graphical models with explicitly modeled
latent variables. These models correspond to latent-factor graphs, which
we represent via the LatentDigraph class.

```         
> # Latent digraphs are specified by their directed adjacency matrix L
> library(SEMID)
> L = matrix(c(0, 1, 0, 0, 0, 0,
+              0, 0, 1, 0, 0, 0,
+              0, 0, 0, 0, 0, 0,
+              0, 0, 0, 0, 1, 0,
+              0, 0, 0, 0, 0, 0,
+              1, 1, 1, 1, 1, 0), 6, 6, byrow=TRUE)
> observedNodes = seq(1,5)
> latentNodes = c(6)
>
> # Create the latent digraph object corresponding to L
> g = LatentDigraph(L, observedNodes, latentNodes)
>
> # Plot latent digraph
> plot(g)
```

The function `lfhtcID` implements the algorithm to check
LF-HTC-identifiability as presented in

Rina Foygel Barber, Mathias Drton, Nils Sturma, Luca Weihs (2022).
Half-Trek Criterion for Identifiability of Latent Variable Models. *Ann.
Statist.* 50(6):3174--3196. <https://doi.org/doi:10.1214/22-AOS2221>.

The LF-HTC is applicable to all graphs where the latent nodes are source
nodes.

```         
> lfhtcID(g)
Call: lfhtcID(graph = g)

Latent Digraph Info
# observed nodes: 5 
# latent nodes: 1 
# total nr. of edges between observed nodes: 3 

Generic Identifiability Summary
# nr. of edges between observed nodes shown gen. identifiable: 3 
# gen. identifiable edges: 1->2, 2->3, 4->5
```

Note that the corresponding mixed graph obtained from a latent
projection is not identifiable; see Section 4 in Barber et al. (2022).

```         
> # Get a mixed graph via latent projection
> gMixed <- g$getMixedGraph()
> gMixed$plot()

> # Check the original half-trek criterion on the mixed graph
> htcID(gMixed)
Call: htcID(mixedGraph = gMixed)

Mixed Graph Info.
# nodes: 5 
# dir. edges: 3 
# bi. edges: 10 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 0 
# bi. edges shown gen. identifiable: 0 

Generically identifiable dir. edges:
None

Generically identifiable bi. edges:
None
```

### Estimating Direct Causal Effects

If a graph is generically identifiable, we can use the identification
formulas to obtain estimators of the direct causal effects. For an
example, see the more detailed description
<https://st-mardi.quarto.pub/gmci/chapters/notebook_gallery/notebooks/GMCI-notebook-SEMID/notebook.html>.

## Identifiability in Sparse Factor Analysis

The matching criterion is a sufficient condition for generic
identification of the factor loading matrix (up to column sign) in
factor analysis. It is developed in the following paper:

Nils Sturma, Miriam Kranzlmüller, Irem Portakal, Mathias Drton (2025).
Matching Criterion for Identifiability in Sparse Factor Analysis. *arXiv
preprint* arXiv:2502.02986

We represent sparse factor analysis graphs via the adjacency matrix
`lambda`, where the columns represent latent nodes and the rows
represent the observed nodes.

```         
> # The factor analysis graph is specified by the matrix lambda
> library(SEMID)
> lambda = matrix(c(1, 0, 0,
+                   1, 1, 0,
+                   0, 1, 1,
+                   1, 0, 1,
+                   0, 1, 0,
+                   0, 0, 1), 6, 3, byrow=TRUE)
> # The latent nodes are nodes 1, 2, and 3, while the observed nodes are the 
> # nodes 4, 5, 6, 7, 8, and 9.
```

The function `mID` implements an algorithm to check M-identifiability:

```         
> mID(lambda)
Call: mID(lambda = lambda)

Factor Analysis Graph Info:
latent nodes:  1 2 3 
observed nodes:  4 5 6 7 8 9 

Generic Sign-Identifiability Summary:
M-identifiable:    TRUE
Tuple list:
  Tuple 1 
    h: 1
    S: 
    v: 4
    W: 5
    U: 7
  Tuple 2 
    h: 2
    S: 1
    v: 5
    W: 6
    U: 8
  Tuple 3 
    h: 3
    S: 1, 2
    v: 6
    W: 7
    U: 9
```

M-identifiability can only establish identifiability of graphs that
satisfy the Zero Upper Triangular Assumption (ZUTA). Via the function
`ZUTA` we can check this assumption.

```         
> ZUTA(lambda)
Call: ZUTA(lambda = lambda)

Factor Analysis Graph Info:
latent nodes:  1 2 3 4 5 
observed nodes:  6 7 8 9 10 11 12 13 14 15 

ZUTA:    TRUE
```

Sturma et al. (2025) also provide an extended, more powerful sufficient
condition. We can check 'extended M-identifiability' as follows.

```         
> # The factor analysis graph is specified by the matrix lambda
> library(SEMID)
> lambda = matrix(c(1, 0, 0, 0, 0,
+                   1, 1, 0, 0, 0,
+                   1, 1, 1, 0, 0,
+                   1, 1, 1, 1, 0,
+                   1, 1, 1, 1, 0,
+                   1, 1, 1, 1, 0,
+                   1, 1, 1, 1, 1,
+                   1, 1, 1, 1, 0,
+                   0, 0, 0, 0, 1,
+                   0, 0, 0, 0, 1), 10, 5, byrow=TRUE)
> # The latent nodes are nodes 1, 2, 3, 4, and 5, while the observed nodes are the 
> # nodes 6, 7, 8, 9, 10, 11, 12, 13, 14, and 15.
         
> extmID(lambda)
Call: extmID(lambda = lambda)

Factor Analysis Graph Info:
latent nodes:  1 2 3 4 5 
observed nodes:  6 7 8 9 10 11 12 13 14 15 

Generic Sign-Identifiability Summary:
extM-identifiable:    TRUE
Tuple list:
  Tuple 1 
    criterion: localBB
    S: 
    new nodes in S: 1, 2, 3, 4
    U: 6, 7, 8, 9, 10, 11, 12, 13
  Tuple 2 
    criterion: matching
    h: 5
    S: 1, 2, 3, 4
    v: 12
    W: 14
    U: 15
```
