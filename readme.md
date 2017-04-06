
blopmatch: Matching Estimator based on a Bilevel Optimization Problem
=====================================================================

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
ans <- blopi_glpk(5, X)

# Resulting weights (matches)
ans$lambda
#>  [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.1793899 0.0000000
#>  [8] 0.2995506 0.0000000 0.0000000 0.0000000 0.5210595 0.0000000 0.0000000
#> [15] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000

# Target vs Projected
X[5,]
#> [1] -0.1544987  0.7678420
colSums(ans$lambda*X[-5,])
#> [1] -0.1544987  0.7678420
```
