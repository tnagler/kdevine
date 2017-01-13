kdevine
=========

> Multivariate kernel density estimation with vine copulas

[![Build status Linux](https://travis-ci.org/tnagler/kdevine.svg?branch=master)](https://travis-ci.org/tnagler/kdevine) 
[![Build status Windows](https://ci.appveyor.com/api/projects/status/epfs987wspjqkwlk?svg=true)](https://ci.appveyor.com/project/tnagler/kdevine) [![CRAN version](http://www.r-pkg.org/badges/version/kdevine)](https://cran.r-project.org/package=kdevine) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/kdevine)](https://cran.r-project.org/package=kdevine)

This package implements a vine copula based kernel density estimator. The 
estimator does not suffer from the curse of dimensionality and is therefore well 
suited for high-dimensional applications (see, Nagler and Czado, 2016). 

You can install the latest development version as follows:

``` r
devtools::install_github("tnagler/kdevine")
```


References
----------


Nagler, T., Czado, C. (2016)  
Evading the curse of dimensionality in nonparametric density estimation with simplified vine copulas  
*Journal of Multivariate Analysis 151, 69-89* ([doi:10.1016/j.jmva.2016.07.003](http://dx.doi.org/10.1016/j.jmva.2016.07.003))  
[preprint](http://arxiv.org/abs/1503.03305)

