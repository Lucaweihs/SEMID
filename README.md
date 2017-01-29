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
> library(SEMID)
> L = t(matrix(
+ c(0, 1, 0, 0, 0,
+   0, 0, 0, 1, 1,
+   0, 0, 0, 1, 0,
+   0, 1, 0, 0, 1,
+   0, 0, 0, 1, 0), 5, 5))
>
> O =t(matrix(
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

Mixed graphs are specified by their adjacency matrix \code{L} and
bidirected adjacency matrix \code{O}. Lets use the `graphID` function
to check global

``` > library(SEMID)
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
                               restrictedEdgewiseIdentifyStep),
	           tianDecompose = F)

R version 3.3.2 (2016-10-31) -- "Sincere Pumpkin Patch"
Copyright (C) 2016 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin13.4.0 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[Workspace loaded from ~/Dropbox/Research/Half Trek/SEMID/.RData]

> semID(MixedGraph(L, O))
Error: could not find function "semID"

Restarting R session...

> library(SEMID)
> semID(MixedGraph(L, O))
Call: semID(mixedGraph = MixedGraph(L, O))

Is globally identifiable?:
  

Has a generically infinite-to-one parameterization?:
 INCONCLUSIVE 
Number of parameters shown generically identifiable:
Directed edges: 0 out of 6 
Called from: eval(expr, envir, enclos)
Browse[1]> n
debug at /Users/lucaweihs/Dropbox/Research/Half Trek/SEMID/R/utils.R#130: cat(paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2, 
    "out of", sum(x$genericIDResult$mixedGraph$O()))/2, "\n")
Browse[2]> n
Error in paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2,  : 
  non-numeric argument to binary operator
In addition: Warning message:
In if (!is.na(x$genericIDResult)) { :
  the condition has length > 1 and only the first element will be used
Called from: cat(paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2, 
    "out of", sum(x$genericIDResult$mixedGraph$O()))/2, "\n")
Browse[2]> n
> semID(MixedGraph(L, O))
Call: semID(mixedGraph = MixedGraph(L, O))

Is globally identifiable?:
  

Has a generically infinite-to-one parameterization?:
 INCONCLUSIVE 
Number of parameters shown generically identifiable:
Directed edges: 0 out of 6 
Called from: eval(expr, envir, enclos)
Browse[1]> n
debug at /Users/lucaweihs/Dropbox/Research/Half Trek/SEMID/R/utils.R#130: cat(paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2, 
    "out of", sum(x$genericIDResult$mixedGraph$O()))/2, "\n")
Browse[2]> n
Error in paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2,  : 
  non-numeric argument to binary operator
In addition: Warning message:
In if (!is.na(x$genericIDResult)) { :
  the condition has length > 1 and only the first element will be used
Called from: cat(paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2, 
    "out of", sum(x$genericIDResult$mixedGraph$O()))/2, "\n")
Browse[2]> x
 [1]  0.5557319 -0.7436013  0.9768292  0.9552242  0.7984791
 [6] -0.6042898  0.2636222  0.7721992 -0.3662376 -0.8263269
Browse[2]> x
 [1]  0.5557319 -0.7436013  0.9768292  0.9552242  0.7984791
 [6] -0.6042898  0.2636222  0.7721992 -0.3662376 -0.8263269
Browse[2]> sum(x$genericIDResult$mixedGraph$L())
Error during wrapup: $ operator is invalid for atomic vectors
Browse[2]> Q
> rm(x)
> semID(MixedGraph(L, O))
Call: semID(mixedGraph = MixedGraph(L, O))

Is globally identifiable?:
  

Has a generically infinite-to-one parameterization?:
 INCONCLUSIVE 
Number of parameters shown generically identifiable:
Directed edges: 0 out of 6 
Called from: eval(expr, envir, enclos)
Browse[1]> n
debug at /Users/lucaweihs/Dropbox/Research/Half Trek/SEMID/R/utils.R#130: cat(paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2, 
    "out of", sum(x$genericIDResult$mixedGraph$O()))/2, "\n")
Browse[2]> n
Error in paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2,  : 
  non-numeric argument to binary operator
In addition: Warning message:
In if (!is.na(x$genericIDResult)) { :
  the condition has length > 1 and only the first element will be used
Called from: cat(paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2, 
    "out of", sum(x$genericIDResult$mixedGraph$O()))/2, "\n")
Browse[2]> x$genericIDResult
Error during wrapup: object 'x' not found
Browse[2]> Q

Restarting R session...

> library(SEMID)
> semID(MixedGraph(L, O))
Call: semID(mixedGraph = MixedGraph(L, O))

Is globally identifiable?:
  

Has a generically infinite-to-one parameterization?:
 INCONCLUSIVE 
Called from: eval(expr, envir, enclos)
Browse[1]> n
debug at /Users/lucaweihs/Dropbox/Research/Half Trek/SEMID/R/utils.R#126: cat(paste("Number of parameters shown generically identifiable:\n"))
Browse[2]> n
Number of parameters shown generically identifiable:
debug at /Users/lucaweihs/Dropbox/Research/Half Trek/SEMID/R/utils.R#127: cat(paste("Directed edges:", length(unlist(x$genericIDResult$solvedParents)), 
    "out of", sum(x$genericIDResult$mixedGraph$L())), "\n")
Browse[2]> x$genericIDResult
Call: generalGenericID(mixedGraph = mixedGraph, idStepFunctions = genericIdStepFunctions, 
    tianDecompose = tianDecompose)

Mixed Graph Info.
# nodes: 6 
# dir. edges: 6 
# bi. edges: 5 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 0 
# bi. edges shown gen. identifiable: 0 

Generically identifiable dir. edges:
None

Generically identifiable bi. edges:
None
Warning message:
In if (!is.na(x$genericIDResult)) { :
  the condition has length > 1 and only the first element will be used
Browse[2]> x$genericIDResult
Call: generalGenericID(mixedGraph = mixedGraph, idStepFunctions = genericIdStepFunctions, 
    tianDecompose = tianDecompose)

Mixed Graph Info.
# nodes: 6 
# dir. edges: 6 
# bi. edges: 5 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 0 
# bi. edges shown gen. identifiable: 0 

Generically identifiable dir. edges:
None

Generically identifiable bi. edges:
None
Browse[2]> Q

Restarting R session...

> library(SEMID)

Restarting R session...

> library(SEMID)
> semID(MixedGraph(L, O))
Call: semID(mixedGraph = MixedGraph(L, O))

Is globally identifiable?:
  

Has a generically infinite-to-one parameterization?:
 INCONCLUSIVE 
Called from: eval(expr, envir, enclos)
Browse[1]> n
debug at /Users/lucaweihs/Dropbox/Research/Half Trek/SEMID/R/utils.R#126: cat(paste("Number of parameters shown generically identifiable:\n"))
Browse[2]> n
Number of parameters shown generically identifiable:
debug at /Users/lucaweihs/Dropbox/Research/Half Trek/SEMID/R/utils.R#127: cat(paste("Directed edges:", length(unlist(x$genericIDResult$solvedParents)), 
    "out of", sum(x$genericIDResult$mixedGraph$L())), "\n")
Browse[2]> 1+1
[1] 2
Browse[2]> x$genericIDResult$solvedParents
[[1]]
numeric(0)

[[2]]
numeric(0)

[[3]]
numeric(0)

[[4]]
numeric(0)

[[5]]
numeric(0)

[[6]]
numeric(0)

Browse[2]> length(unlist(x$genericIDResult$solvedParents))
[1] 0
Browse[2]> x$genericIDResult$mixedGraph$L())
Error: unexpected ')' in "x$genericIDResult$mixedGraph$L())"
Browse[4]> x$genericIDResult$mixedGraph$L()
     [,1] [,2] [,3] [,4] [,5] [,6]
[1,]    0    1    1    1    0    1
[2,]    0    0    0    0    0    0
[3,]    0    0    0    0    0    0
[4,]    0    0    0    0    1    0
[5,]    0    0    0    0    0    0
[6,]    0    0    0    0    1    0
Browse[4]> sum(x$genericIDResult$mixedGraph$L())
[1] 6
Browse[4]> n
Browse[2]> n
Directed edges: 0 out of 6 
debug at /Users/lucaweihs/Dropbox/Research/Half Trek/SEMID/R/utils.R#130: cat(paste("Bidirected edges:", length(unlist(x$genericIDResult$solvedSiblings))/2, 
    "out of", sum(x$genericIDResult$mixedGraph$O()))/2, "\n")
Browse[2]> x$genericIDResult$solvedSiblings
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

[[6]]
integer(0)

Browse[2]> length(unlist(x$genericIDResult$solvedSiblings)) / 2
[1] 0
Browse[2]> sum(x$genericIDResult$mixedGraph$O())) / 2
Error: unexpected ')' in "sum(x$genericIDResult$mixedGraph$O()))"
Browse[4]> sum(x$genericIDResult$mixedGraph$O())) / 
Error during wrapup: unexpected ')' in "sum(x$genericIDResult$mixedGraph$O()))"
Browse[4]> sum(x$genericIDResult$mixedGraph$O())
[1] 10
Browse[4]> sum(x$genericIDResult$mixedGraph$O()) / 2
[1] 5
Browse[4]> Q

Restarting R session...

> library(SEMID)
> semID(MixedGraph(L, O))
Call: semID(mixedGraph = MixedGraph(L, O))

Is globally identifiable?:
  

Has a generically infinite-to-one parameterization?:
 INCONCLUSIVE 
Number of parameters shown generically identifiable:
Directed edges: 0 out of 6 
Bidirected edges: 0 out of 5 

Restarting R session...

> library(SEMID)
> semID(MixedGraph(L, O))
Call: semID(mixedGraph = MixedGraph(L, O))

Is globally identifiable?:
 FALSE 

Has a generically infinite-to-one parameterization?:
 INCONCLUSIVE 
Number of parameters shown generically identifiable:
Directed edges: 0 out of 6 
Bidirected edges: 0 out of 5 

Restarting R session...

> library(SEMID)

Restarting R session...

> library(SEMID)
> semID(MixedGraph(L, O))
Call: semID(mixedGraph = MixedGraph(L, O))

Is globally identifiable?:
FALSE

Has a generically infinite-to-one parameterization?:
 INCONCLUSIVE

Number of parameters shown generically identifiable:
Directed edges: 0 out of 6 
Bidirected edges: 0 out of 5 

Restarting R session...

> library(SEMID)
> semID(MixedGraph(L, O))
Call: semID(mixedGraph = MixedGraph(L, O))

Is globally identifiable?:
FALSE

Has a generically infinite-to-one parameterization?:
INCONCLUSIVE

Number of parameters shown generically identifiable:
Directed edges: 0 out of 6 
Bidirected edges: 0 out of 5 
> semID(MixedGraph(L, O), genericIdStepFunctions = list(edgewiseIdentifyStep, trekSeparationIdentifyStep))
Call: semID(mixedGraph = MixedGraph(L, O), genericIdStepFunctions = list(edgewiseIdentifyStep, 
    trekSeparationIdentifyStep))

Is globally identifiable?:
FALSE

Has a generically infinite-to-one parameterization?:
INCONCLUSIVE

Number of parameters shown generically identifiable:
Directed edges: 6 out of 6 
Bidirected edges: 5 out of 5 

Restarting R session...

> library(SEMID)
> semID(MixedGraph(L, O), genericIdStepFunctions = list(edgewiseIdentifyStep, trekSeparationIdentifyStep))
Call: semID(mixedGraph = MixedGraph(L, O), genericIdStepFunctions = list(edgewiseIdentifyStep, 
    trekSeparationIdentifyStep))

Is globally identifiable?:
FALSE

Has a generically infinite-to-one parameterization?:
FALSE

Number of parameters shown generically identifiable:
Directed edges: 6 out of 6 
Bidirected edges: 5 out of 5 

Restarting R session...

> library(SEMID)
> semID(MixedGraph(L, O), genericIdStepFunctions = list(edgewiseIdentifyStep, trekSeparationIdentifyStep))
Call: semID(mixedGraph = MixedGraph(L, O), genericIdStepFunctions = list(edgewiseIdentifyStep, 
    trekSeparationIdentifyStep))

Did Tian decomposition?
TRUE 

Is globally identifiable?:
FALSE

Has a generically infinite-to-one parameterization?:
FALSE

Number of parameters shown generically identifiable:
Directed edges: 6 out of 6 
Bidirected edges: 5 out of 5 
> plotMixedGraph(L, O)
> htcID(L, O, tianDecompose = F)
Error in htcID(L, O, tianDecompose = F) : unused argument (O)
> L
     [,1] [,2] [,3] [,4] [,5] [,6]
[1,]    0    1    1    1    0    1
[2,]    0    0    0    0    0    0
[3,]    0    0    0    0    0    0
[4,]    0    0    0    0    1    0
[5,]    0    0    0    0    0    0
[6,]    0    0    0    0    1    0
> htcID(MixedGraph(L, O), tianDecompose = F)
Call: htcID(mixedGraph = MixedGraph(L, O), tianDecompose = F)

Mixed Graph Info.
# nodes: 6 
# dir. edges: 6 
# bi. edges: 5 

Generic Identifiability Summary
# dir. edges shown gen. identifiable: 0 
# bi. edges shown gen. identifiable: 0 

Generically identifiable dir. edges:
None

Generically identifiable bi. edges:
None
> L = t(matrix(
+     c(0, 1, 1, 0, 0,
+       0, 0, 1, 1, 1,
+       0, 0, 0, 1, 0,
+       0, 0, 0, 0, 1,
+       0, 0, 0, 0, 0), 5, 5))
> 
> O = t(matrix(
+     c(0, 0, 0, 1, 0,
+       0, 0, 1, 0, 1,
+       0, 0, 0, 0, 0,
+       0, 0, 0, 0, 0,
+       0, 0, 0, 0, 0), 5, 5)); O=O+t(O)
> graph = MixedGraph(L, O)
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
> L
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    1    1    0    0
[2,]    0    0    1    1    1
[3,]    0    0    0    1    0
[4,]    0    0    0    0    1
[5,]    0    0    0    0    0
> O
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    0    0    1    0
[2,]    0    0    1    0    1
[3,]    0    1    0    0    0
[4,]    1    0    0    0    0
[5,]    0    1    0    0    0
> g = MixedGraph(L, O)
> g$plot()
> g$tianDecompose()
[[1]]
[[1]]$internal
[1] 1 4

[[1]]$incoming
[1] 2 3

[[1]]$topOrder
[1] 1 2 3 4

[[1]]$L
     [,1] [,2] [,3] [,4]
[1,]    0    0    0    0
[2,]    0    0    0    1
[3,]    0    0    0    1
[4,]    0    0    0    0

[[1]]$O
     [,1] [,2] [,3] [,4]
[1,]    0    0    0    1
[2,]    0    0    0    0
[3,]    0    0    0    0
[4,]    1    0    0    0


[[2]]
[[2]]$internal
[1] 2 3 5

[[2]]$incoming
[1] 1 4

[[2]]$topOrder
[1] 1 2 3 4 5

[[2]]$L
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    1    1    0    0
[2,]    0    0    1    0    1
[3,]    0    0    0    0    0
[4,]    0    0    0    0    1
[5,]    0    0    0    0    0

[[2]]$O
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    0    0    0    0
[2,]    0    0    1    0    1
[3,]    0    1    0    0    0
[4,]    0    0    0    0    0
[5,]    0    1    0    0    0


> semID(L, O, genericIdStepFunctions = list(edgewiseIdentifyStep, trekSeparationIdentifyStep))
Error in mixedGraph$tianDecompose : 
  $ operator is invalid for atomic vectors
In addition: Warning message:
In if (testGlobalID) { :
  the condition has length > 1 and only the first element will be used
> semID(L, O, genericIdStepFunctions = list(edgewiseIdentifyStep, trekSeparationIdentifyStep))
Error in mixedGraph$tianDecompose : 
  $ operator is invalid for atomic vectors
In addition: Warning message:
In if (testGlobalID) { :
  the condition has length > 1 and only the first element will be used
> semID(g, genericIdStepFunctions = list(edgewiseIdentifyStep, trekSeparationIdentifyStep))
Call: semID(mixedGraph = g, genericIdStepFunctions = list(edgewiseIdentifyStep, 
    trekSeparationIdentifyStep))

Did Tian decomposition?
TRUE 

Is globally identifiable?:
FALSE

Has a generically infinite-to-one parameterization?:
FALSE

Number of parameters shown generically identifiable:
Directed edges: 7 out of 7 
Bidirected edges: 3 out of 3 
> L
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    1    1    0    0
[2,]    0    0    1    1    1
[3,]    0    0    0    1    0
[4,]    0    0    0    0    1
[5,]    0    0    0    0    0
> O
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    0    0    1    0
[2,]    0    0    1    0    1
[3,]    0    1    0    0    0
[4,]    1    0    0    0    0
[5,]    0    1    0    0    0
> L = t(matrix(
+     c(0, 1, 0, 0, 0,
+       0, 0, 0, 1, 1,
+       0, 0, 0, 1, 0,
+       0, 1, 0, 0, 1,
+       0, 0, 0, 1, 0), 5, 5))
> O = t(matrix(
+     c(0, 0, 0, 0, 0,
+       0, 0, 1, 0, 1,
+       0, 0, 0, 1, 0,
+       0, 0, 0, 0, 0,
+       0, 0, 0, 0, 0), 5, 5)); O=O+t(O)
> restrictedEdgewiseIdentifyStep <- function(mixedGraph,
+                                            unsolvedParents,
+                                            solvedParents, 
+                                            identifier) {
+     return(edgewiseIdentifyStep(mixedGraph, unsolvedParents,
+                                 solvedParents, identifier, 
+                                 subsetSizeControl = 1))
+ }
> graph = MixedGraph(L, O)
> generalGenericID(graph, list(htcIdentifyStep,
+                              restrictedEdgewiseIdentifyStep),
+                  tianDecompose = F)
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
> # We can do better (fewer unsolvd parents) if we don't restrict the
> # edgewise identifier algorithm as much
> generalGenericID(graph, list(htcIdentifyStep,
+                              edgewiseIdentifyStep),
+                  tianDecompose = F)
Call: generalGenericID(mixedGraph = graph, idStepFunctions = list(htcIdentifyStep, 
    edgewiseIdentifyStep), tianDecompose = F)

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
