% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{weighted_norm}
\alias{weighted_norm}
\title{Weighted Norm}
\usage{
weighted_norm(X, W, p = 1)
}
\arguments{
\item{X}{Numeric matrix of size \eqn{n\times K}{n*K}.}

\item{W}{Numeric matrix of size \eqn{K\times K}{K*K}. When \code{W} is equal to \code{solve(var(X))}, it is the Mahalanobis norm.}

\item{p}{Numeric scalar. See \code{\link{weighted_norm}}.}
}
\value{
A square matrix of size \eqn{n\times n}{n * n}.
}
\description{
Weighted Norm

\code{}
}
\details{
Computes
\deqn{
\left[(\mathbf{x}_i - \mathbf{x}_j)  \times \mathbf{W} \times (\mathbf{x}_i - \mathbf{x}_j)^\mathbf{t}\right]^\frac{p}{2},
\quad\forall i,j
}{
[(x_i - x_j) \%*\% W \%*\% (x_i - x_j)]^(p/2), for all i, j
}
}
\examples{

set.seed(1231)
X <- matrix(rnorm(20*2), ncol=2)
W <- diag(2)

weighted_norm(X, W, 1.0)

}
