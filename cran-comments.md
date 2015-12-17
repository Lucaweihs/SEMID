## Test environments
* local OS X install, R 3.2.3
* Using devtools::build_win()

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking R code for possible problems ... NOTE
graphID.output.alltests: no visible global function definition for
  'write.table'
graphID.output.notalltests: no visible global function definition for
  'write.table'
Undefined global functions or variables:
  write.table
Consider adding
  importFrom("utils", "write.table")
  
  This note only occurs when checking the package using the unstable development
  version of R (2015-12-13 r69768) when running the devtools::build_win()
  command. As the write.table function is global and I could not find any 
  documentation noting that this was to be deprecated I can only imagine this
  is not a problem of this package.

## Downstream dependencies
There are currently no downstream dependencies for this package.
