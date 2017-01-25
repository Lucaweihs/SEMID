# SEMID Package

## Purpose

This package offers a number of different functions for determining
global and generic identifiability of path diagrams / mixed graphs. The
following sections highlight the primary ways in which the package can be used.

## Global Identifiability

Drton, Foygel, and Sullivant (2011) showed that there exist if and only if
graphical conditions for testing whether or not the parameters in a mixed graph
are globally identifiable. This criterion can be accessed through the function
`graphID.globalID`.

 Drton, M., Foygel, R., and Sullivant, S.  (2011) Global
identifiability of linear structural equation models. _Ann. Statist._
39(2): 865-886.

## Generic Identiability of Parameters

There still do not exist any 'if and only if' graphical conditios for testing
whether or not certain parameters in a mixed graph are generically identifiable. There do, howeover, exist some necessary and some sufficient conditions which
work for a large collection of graphs.

### Sufficient Conditions

Until recently, criterions for generic identifiability, like the half-trek
criterion of Foygel, Draisma, and Drton (2012), had to show that all edges incoming to a node where generically identifiable simultaenously and thus, if
any single such edge incoming to a node was generically nonidentifiable, the
criterion would fail. The recent work of Weihs, Robeva, Robinson, et al. (2017)
develops new criteria that are able to identify subsets of edges coming into a
node substantially improving upon prior methods at the cost of computational
efficiency. We list both the older algorithms (available in prior versions
of this package) and the newer algorithms below.

#### In prior version of the package

* `graphID.htcID` implements the half-trek criterion of

Foygel, Rina; Draisma, Jan; Drton, Mathias. Half-trek criterion for generic identifiability of linear structural equation models. Ann. Statist. 40 (2012), no. 3, 1682--1713. doi:10.1214/12-AOS1012.

* `graphID.ancestralID` implements the ancestor decomposition techniques of 

Drton, M., and Weihs, L. (2016) Generic Identifiability of Linear Structural Equation Models by Ancestor Decomposition. Scand J Statist, 43: 1035–1045. doi: 10.1111/sjos.12227.

* `graphID` - a single function that gives access to a number of the other
identification strategies at once, see the documentation for more information.

#### New additions

* `htcID` - the half-trek critierion from `graphID.htcID` but whose returned
values are more consistent with a focus on identifying individual edges.

* `ancestralID` - an updated version of the `graphID.ancestralID` which works on
cyclic mixed graphs.

* `edgewiseID` - an edgewise identification algorithm strictly improving upon
the half-trek criterion.

* `edgewiseTSID` - an edgewise identification algorithm leveraging 
trek-separation relations to identify even more edges than `edgewiseID`, this
has the downside being computationally expensive.

* `generalGenericID` - an generic identification algorithm template that allows
you to mix and match different identification techniques to find the right 
balance between computational efficiency and exhaustiveness of the procedure. 
See the below examples for more details.

### Necessary Conditions

* `graphID.nonHtcID` - the best known test for whether any edges in a mixed 
graph are generically non-identifiabile. This comes Foygel, Draisma, and Drton (2012).

### Examples

Mixed graphs are specified by their adjacency matrix \code{L} and bidirected
adjacency matrix \code{O}. Lets use the `graphID` function to check global

```
> library(SEMID)
> L = t(matrix(
+   c(0, 1, 1, 0, 0,
+     0, 0, 1, 1, 1,
+     0, 0, 0, 1, 0,
+     0, 0, 0, 0, 1,
+     0, 0, 0, 0, 0), 5, 5))

> O = t(matrix(
+   c(0, 0, 0, 1, 0,
+     0, 0, 1, 0, 1,
+     0, 0, 0, 0, 0,
+     0, 0, 0, 0, 0,
+     0, 0, 0, 0, 0), 5, 5)); O=O+t(O)

> # Without using decomposition techniques we can't identify all nodes
> # just using the half-trek criterion
> graphID(L, O, decomp.if.acyclic = F)
                                              
 HTC-ID'able nodes      Identifiability       
 1,2                    HTC-inconclusive
>
> # If we use decomposition techniques then we can
> graphID(L, O, decomp.if.acyclic = T)
                                                                           
 Component                Nodes                    HTC-ID'able nodes       
 1                        1,4                      1,4                     
 2                        2,3,5                    2,5,3                   
                         
 Identifiability         
 Globally ID'able        
 Generically ID'able 
>
> # The edgewiseTSID function can show that all edges are generically
> # identifiable without proprocessing with decomposition techniques
> result = edgewiseTSID(L, O, tianDecompose = F)
> result$unsolvedParents
[[1]]
integer(0)

[[2]]
integer(0)

[[3]]
integer(0)

[[4]]
integer(0)

[[5]]
integer(0)
> # The above shows that each node in the graph has no unsolved parents
> # after running the edgewiseTSID algorithm. See the help of edgewiseTSID
> # to find out more information about what else is returned by edgewiseTSID.
```

Using the `generalGenericId` method we can also mix and match different
identification strategies. Lets say we wanted to first try to identify
everything using the half-trek criterion but then, if there are still things
that cant be shown generically identifiable, we want to use the edgewise
criterion by limiting the edgesets it looks at to be a small size. We can
do this as follows:

```
> library(SEMID)
> # Lets first define some matrices for a mixed graph
> L = t(matrix(
+   c(0, 1, 0, 0, 0,
+     0, 0, 0, 1, 1,
+     0, 0, 0, 1, 0,
+     0, 1, 0, 0, 1,
+     0, 0, 0, 1, 0), 5, 5))

> O = t(matrix(
+   c(0, 0, 0, 0, 0,
+     0, 0, 1, 0, 1,
+     0, 0, 0, 1, 0,
+     0, 0, 0, 0, 0,
+     0, 0, 0, 0, 0), 5, 5)); O=O+t(O)

> # Now lets define an "identification step" function corresponding to using
> # the edgewise identification algorithm but with subsets controlled by 1.
> restrictedEdgewiseIdentifyStep <- function(mixedGraph, unsolvedParents,
+                                              solvedParents, identifier) {
+   return(edgewiseIdentifyStep(mixedGraph, unsolvedParents, solvedParents,
+                               identifier, subsetSizeControl = 1))
+ }

> # Now we run an identification algorithm that iterates between the htc
> # and the "restricted" edgewise identification algorithm
> generalGenericID(L, O, list(htcIdentifyStep, restrictedEdgewiseIdentifyStep),
+                   tianDecompose = F)$unsolvedParents
[[1]]
integer(0)

[[2]]
[1] 1 4

[[3]]
integer(0)

[[4]]
[1] 2 3 5

[[5]]
integer(0)

> # We can do better (fewer unsolvd parents) if we don't restrict the edgewise 
> # identifier algorithm as much
> generalGenericID(L, O, list(htcIdentifyStep, edgewiseIdentifyStep),
+                  tianDecompose = F)$unsolvedParents
[[1]]
integer(0)

[[2]]
[1] 1 4

[[3]]
integer(0)

[[4]]
[1] 3

[[5]]
integer(0)
```

## Experiments with mixed graphs

To be able to implement the above algorithms we also developed a MixedGraph
class using the `R.oo` package of

 Bengtsson, H. (2003)The R.oo package - Object-Oriented Programming with References Using Standard R Code, Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003), ISSN 1609-395X, Hornik, K.; Leisch, F. & Zeileis, A. (eds.) URL https://www.r-project.org/conferences/DSC-2003/Proceedings/Bengtsson.pdf

This class can make it much easier to represent a mixed graph and run
experiments with them. For instance we can create a mixed graph, plot it, and
test if there is a half-trek system between two sets of vertices very easily:

```
> library(SEMID)
> L = t(matrix(
+       c(0, 1, 0, 0, 0,
+         0, 0, 0, 1, 1,
+         0, 0, 0, 1, 0,
+         0, 1, 0, 0, 1,
+         0, 0, 0, 1, 0), 5, 5))
> 
> O = t(matrix(
+       c(0, 0, 0, 0, 0,
+         0, 0, 1, 0, 1,
+         0, 0, 0, 1, 0,
+         0, 0, 0, 0, 0,
+         0, 0, 0, 0, 0), 5, 5)); O=O+t(O)
>      
> # Create the mixed graph object corresponding to L and O
> g = MixedGraph(L, O)
> 
> # Plot the mixed graph
> g$plot()
> 
> # Test whether or not there is a half-trek system from the nodes
> # 1,2 to 3,4
> g$getHalfTrekSystem(c(1,2), c(3,4))
$systemExists
[1] TRUE

$activeFrom
[1] 1 2
```

See the documentation for the MixedGraph class `?MixedGraph` for more 
information.
