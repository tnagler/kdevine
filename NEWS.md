kdevine 0.4.3
---------

BUG FIXES

   * removed all calls to `closeAllConnections()`.


kdevine 0.4.2
---------

The kdevine package is no longer actively developed. Consider using   

   - the 'kde1d' package for marginal estimation,   
   - the functions `vine()` and `vinecop()` from the 'rvinecopulib' 
     package as replacements for `kdevine()` and `kdevinecop()`.
     
BUG FIXES

   * fixed `rkde1d()` for discrete data.


kdevine 0.4.1
---------

DEPENDENCIES

   * Removed dependency on `ks` package (and thereby its dependencies like rgl)
     by calling `KernSmooth::dpik()` directly instead of `ks::hpi()`.

BUG FIXES

   * fixed `plot.kde1d` type for continuous data (was shown as histogram).


kdevine 0.4.0
---------

NEW FEATURES

   * support for discrete variables based on continuous convolution (using R 
     package cctools).

   * new bandwidth default `mult_1d = log(1 + d)` in `kdevine()`.


kdevine 0.3.0
---------

First release on CRAN.
