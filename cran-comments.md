Fixes segfaults on R-devel builds caused by careless use of uninitialized 
Rcpp::NumericVectors.

## Test environments
* ubuntu 20.04 (devel, release, oldrel) 
* macOs catalina (release) 
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

## Downstream dependencies
No downstream dependencies yet.
