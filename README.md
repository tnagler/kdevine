kdevine
=========

> Multivariate kernel density estimation with vine copulas

[![Build status Linux](https://travis-ci.org/tnagler/kdevine.svg?branch=master)](https://travis-ci.org/tnagler/kdevine) 

This package implements a vine copula based kernel density estimator. The 
estimator does not suffer from the curse of dimensionality and is therefore well 
suited for high-dimensional applications (see, Nagler and Czado, 2015). 

You can install the latest development version as follows:

``` r
devtools::install_github("tnagler/kdevine")
```


References
----------

Nagler, T., Czado, C. (2015) 
Evading the curse of dimensionality in nonparametric density estimation. 
[*arXiv:1503.03305 [stat.ME]*](http://arxiv.org/abs/1503.03305)

