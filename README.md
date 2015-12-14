# SEMID Package

## Purpose

This package offers a number of different functions for determining
global and generic identifiability of path diagrams / mixed graphs. The primary
functionality of the package can be accessed through the function `graphID`
which allows one to perform multiple tests simultaneously and formats the output
nicely. The different tests of identifiability can, however, also be accessed
directly. In particular we have the following functions:

* `graphID.htcID` and `graphID.nonHtcID` which implement the two algorithms from

 Drton, M., Foygel, R., and Sullivant, S.  (2011) Global
identifiability of linear structural equation models. _Ann. Statist._
39(2): 865-886.

* `graphID.globalID` implementing the criterion of

 Drton, Mathias; Foygel, Rina; Sullivant, Seth. Global identifiability of
linear structural equation models. _Ann. Statist._  39 (2011), no. 2,
865--886.

* `graphID.ancestralID` implementing the algorithm from

 Drton, M. and Weihs, L. (2015) Generic Identifiability of Linear
Structural Equation Models by Ancestor Decomposition. arXiv 1504.02992

## Example

A number of examples can be found by loading the package and running
`demo(SEMID)`, one particular such example follows:

```
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

> graphID(L,O)
                                                                           
 Component                Nodes                    HTC-ID'able nodes       
 1                        1,4                      1,4                     
 2                        2,3,5                    2,5,3                   
                         
 Identifiability         
 Globally ID'able        
 Generically ID'able 
```

Replace `graphID` in the example with the other above functions to see their
output.
