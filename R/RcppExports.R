# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

matching_group_cpp <- function(Treat, exact, zeroindex = 1L) {
    .Call('_blopmatch_matching_group_cpp', PACKAGE = 'blopmatch', Treat, exact, zeroindex)
}

#' Weighted Norm
#' 
#' @template mat
#' @templateVar X 1
#' @templateVar W 1
#' @templateVar p 1
#' @export
#' 
#' @details
#' Computes
#' \deqn{
#' \left[(\mathbf{x}_i - \mathbf{x}_j)  \times \mathbf{W} \times (\mathbf{x}_i - \mathbf{x}_j)^\mathbf{t}\right]^\frac{p}{2},
#' \quad\forall i,j
#' }{
#' [(x_i - x_j) \%*\% W \%*\% (x_i - x_j)]^(p/2), for all i, j
#' }
#' 
#' @return A square matrix of size \eqn{n\times n}{n * n}.
#' @examples
#' 
#' set.seed(1231)
#' X <- matrix(rnorm(20*2), ncol=2)
#' W <- diag(2)
#' 
#' weighted_norm(X, W, 1.0)
#' 
weighted_norm <- function(X, W, p = 1.0) {
    .Call('_blopmatch_weighted_norm', PACKAGE = 'blopmatch', X, W, p)
}

