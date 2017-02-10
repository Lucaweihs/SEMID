# SEMID Package

## Purpose

This package offers a number of different functions for determining
global and generic identifiability of path diagrams / mixed
graphs. The following sections highlight the primary ways in which the
package can be used. Much of the package's functionality can be accessed
through the wrapper function `semID`.

## The MixedGraph class

To be able to implement the different algorithms described below we
created a MixedGraph class using the `R.oo` package of

 Bengtsson, H. (2003)The R.oo package - Object-Oriented Programming
 with References Using Standard R Code, Proceedings of the 3rd
 International Workshop on Distributed Statistical Computing (DSC
 2003), ISSN 1609-395X, Hornik, K.; Leisch, F. & Zeileis, A. (eds.)
 URL
 https://www.r-project.org/conferences/DSC-2003/Proceedings/Bengtsson.pdf

This class can make it much easier to represent a mixed graph and run
experiments with them. For instance we can create a mixed graph, plot
it, and test if there is a half-trek system between two sets of
vertices very easily:

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

## Global Identifiability

Drton, Foygel, and Sullivant (2011) showed that there exist if and
only if graphical conditions for testing whether or not the parameters
in a mixed graph are globally identifiable. This criterion can be
accessed through the function `graphID.globalID`.

 Drton, M., Foygel, R., and Sullivant, S.  (2011) Global
identifiability of linear structural equation models. _Ann. Statist._
39(2): 865-886.

## Generic Identiability of Parameters

There still do not exist any 'if and only if' graphical conditions for
testing whether or not certain parameters in a mixed graph are
generically identifiable. There do, howeover, exist some necessary and
some sufficient conditions which work for a large collection of
graphs.

### Sufficient Conditions

Until recently, criterions for generic identifiability, like the
half-trek criterion of Foygel, Draisma, and Drton (2012), had to show
that all edges incoming to a node where generically identifiable
simultaenously and thus, if any single such edge incoming to a node
was generically nonidentifiable, the criterion would fail. The recent
work of Weihs, Robeva, Robinson, et al. (2017) develops new criteria
that are able to identify subsets of edges coming into a node
substantially improving upon prior methods at the cost of
computational efficiency. We list both the older algorithms (available
in prior versions of this package) and the newer algorithms below.

* `htcID` - implements the half-trek critierion of

Foygel, Rina; Draisma, Jan; Drton, Mathias. Half-trek criterion for
generic identifiability of linear structural equation
models. Ann. Statist. 40 (2012), no. 3,
1682--1713. doi:10.1214/12-AOS1012.

and is an updated version of the (deprecated) function
`graphID.htcID`.

* `ancestralID` - implements the ancestor decomposition techniques of

Drton, M., and Weihs, L. (2016) Generic Identifiability of Linear
Structural Equation Models by Ancestor Decomposition. Scand J Statist,
43: 1035–1045. doi: 10.1111/sjos.12227.

and is an updated version of the (deprecated) function
`graphID.ancestralID`. This new version of the function works also on
cyclic graphs by using an updated version of the Tian decomposition.

* `edgewiseID` - an edgewise identification algorithm strictly
improving upon the half-trek criterion.

* `edgewiseTSID` - an edgewise identification algorithm leveraging
trek-separation relations to identify even more edges than
`edgewiseID`, this has the downside being computationally expensive.

* `generalGenericID` - an generic identification algorithm template
that allows you to mix and match different identification techniques
to find the right balance between computational efficiency and
exhaustiveness of the procedure.  See the below examples for more
details.

### Necessary Conditions

* `graphID.nonHtcID` - the best known test for whether any edges in a
mixed graph are generically non-identifiabile. This comes Foygel,
Draisma, and Drton (2012).

### Examples

Lets use a few of the above functions to check the generic identifiability of 
parameters in a mixed graph.

```
> library(SEMID)
> # Mixed graphs are specified by their directed adjacency matrix L and
> # bidirected adjacency matrix O.
> L = t(matrix(
+ c(0, 1, 1, 0, 0,
+   0, 0, 1, 1, 1,
+   0, 0, 0, 1, 0,
+   0, 0, 0, 0, 1,
+   0, 0, 0, 0, 0), 5, 5))
>
> O = t(matrix(
+ c(0, 0, 0, 1, 0,
+   0, 0, 1, 0, 1,
+   0, 0, 0, 0, 0,
+   0, 0, 0, 0, 0,
+   0, 0, 0, 0, 0), 5, 5)); O=O+t(O)
>
> # Create a mixed graph object
> graph = MixedGraph(L, O)
>
> # Without using decomposition techniques we can't identify all nodes
> # just using the half-trek criterion
> htcID(graph, tianDecompose = F)
Call: htcID(mixedGraph = graph, tianDecompose = F)

Mixed Graph Info.
# nodes: 5 
# dir. edges: 7 
# bi. edges: 3 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 1 
# bi. edges shown gen. identifiable: 0 

Generically identifiable dir. edges:
1->2 

Generically identifiable bi. edges:
None
>
> # The edgewiseTSID function can show that all edges are generically
> # identifiable without proprocessing with decomposition techniques
> edgewiseTSID(graph, tianDecompose = F)
Call: edgewiseTSID(mixedGraph = graph, tianDecompose = F)

Mixed Graph Info.
# nodes: 5 
# dir. edges: 7 
# bi. edges: 3 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 7 
# bi. edges shown gen. identifiable: 3 

Generically identifiable dir. edges:
1->2, 1->3, 2->3, 2->4, 3->4, 2->5, 4->5 

Generically identifiable bi. edges:
1<->4, 2<->3, 2<->5 
>
> # The above shows that all edges in the graph are generically identifiable.
> # See the help of edgewiseTSID to find out more information about what
> # else is returned by edgewiseTSID.
```

Using the `generalGenericId` method we can also mix and match
different identification strategies. Lets say we wanted to first try
to identify everything using the half-trek criterion but then, if
there are still things that cant be shown generically identifiable, we
want to use the edgewise criterion by limiting the edgesets it looks
at to be a small size. We can do this as follows:

```
> library(SEMID)
> # Lets first define some matrices for a mixed graph
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
> # Create a mixed graph object
> graph = MixedGraph(L, O)
>
> # Now lets define an "identification step" function corresponding to
> # using the edgewise identification algorithm but with subsets
> # controlled by 1.
> restrictedEdgewiseIdentifyStep <- function(mixedGraph,
+                                            unsolvedParents,
+                                            solvedParents, 
+                                            identifier) {
+     return(edgewiseIdentifyStep(mixedGraph, unsolvedParents,
+                                 solvedParents, identifier, 
+                                 subsetSizeControl = 1))
+ }
>
> # Now we run an identification algorithm that iterates between the
> # htc and the "restricted" edgewise identification algorithm
> generalGenericID(graph, list(htcIdentifyStep,
+                               restrictedEdgewiseIdentifyStep),
+	                 tianDecompose = F)
Call: generalGenericID(mixedGraph = graph, idStepFunctions = list(htcIdentifyStep, 
    restrictedEdgewiseIdentifyStep), tianDecompose = F)

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
>
> # We can do better (fewer unsolvd parents) if we don't restrict the edgewise 
> # identifier algorithm as much
> generalGenericID(graph, list(htcIdentifyStep, edgewiseIdentifyStep),
+                  tianDecompose = F)
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
