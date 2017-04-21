blopmatch: Matching Estimator based on a Bilevel Optimization Problem
================

[![Travis-CI Build Status](https://travis-ci.org/gvegayon/blopmatch.svg?branch=master)](https://travis-ci.org/gvegayon/blopmatch) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/gvegayon/blopmatch?branch=master&svg=true)](https://ci.appveyor.com/project/gvegayon/blopmatch) [![Coverage Status](https://img.shields.io/codecov/c/github/gvegayon/blopmatch/master.svg)](https://codecov.io/github/gvegayon/blopmatch?branch=master)

An implementation of Díaz, Rau and Rivera (2015) matching estimator for causal inference. The authors propose a matching estimator based on a Bilevel Optimization Problem. In raw terms, the two problems are (1) finding a convex combination that (2) using the closets neighbors possible. The solution to this problem allows computing Treatment Effect estimators that significantly improve balance in case-control studies, and furthermore, can be used for data imputation.

<!-- README.md is generated from README.Rmd. Please edit that file -->
Installation
------------

``` r
devtools::install_github("gvegayon/blopmatch")
```

Examples
--------

``` r
# Loading the package
library(blopmatch)

# Simulating data
set.seed(1331)
X <- matrix(rnorm(20*2), ncol=2)

# Matching individual 5 to the rest
ans <- blopi_glpk(X[5,,drop=FALSE], X[-5,,drop=FALSE])

# Resulting weights (matches)
ans$lambda
#> 1 x 19 sparse Matrix of class "dgCMatrix"
#>                                                                   
#> [1,] . . . . . 0.1793899 . 0.2995506 . . . 0.5210595 . . . . . . .

# Target vs Projected
X[5,]
#> [1] -0.1544987  0.7678420
ans$lambda %*% X[-5,]
#> 1 x 2 Matrix of class "dgeMatrix"
#>            [,1]     [,2]
#> [1,] -0.1544987 0.767842
```
