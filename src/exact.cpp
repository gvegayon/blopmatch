#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List matching_group_cpp(
  const arma::ivec & Treat,
  const arma::mat & exact,
  int zeroindex = 1
) {
  
  unsigned int n = exact.n_rows, i, j;
  unsigned int K = exact.n_cols, k;

  std::vector< std::vector< arma::uword > > ans0(n);

  for (i = 0u; i < n; i++) 
    for (j = i + 1u; j < n; j++) {
      
      // Initializing current search
      bool match;    
  
      // Comparing i to j  
      match = true;
      for (k = 0u; k < K; k++) {
        
        // Opposite groups?
        if ((Treat.at(i) >= 0) && (Treat.at(i) == Treat.at(j))) {
          match = false;
          break;
        }
        
        // If any attribute in E is different, then 
        // no exact and return false
        if (exact.at(i, k) != exact.at(j, k)) {
          match = false;
          break;
        }
        
      }
      
      // Checking if exact
      if (match) {
        ans0.at(i).push_back(j + zeroindex);
        ans0.at(j).push_back(i + zeroindex);
      }
    }
    
  // Coercing into a list
  List ans(n);
  for (i = 0u; i < n; i++)
    if (ans0.at(i).size() > 0)
      ans.at(i) = arma::conv_to< arma::uvec >::from(ans0.at(i));
    
  return ans;
}

//' Weighted Norm
//' 
//' @template mat
//' @templateVar X 1
//' @templateVar W 1
//' @templateVar p 1
//' @export
//' 
//' @details
//' Computes
//' \deqn{
//' \left[(\mathbf{x}_i - \mathbf{x}_j)  \times \mathbf{W} \times (\mathbf{x}_i - \mathbf{x}_j)^\mathbf{t}\right]^\frac{p}{2},
//' \quad\forall i,j
//' }{
//' [(x_i - x_j) \%*\% W \%*\% (x_i - x_j)]^(p/2), for all i, j
//' }
//' 
//' @return A square matrix of size \eqn{n\times n}{n * n}.
//' @examples
//' 
//' set.seed(1231)
//' X <- matrix(rnorm(20*2), ncol=2)
//' W <- diag(2)
//' 
//' weighted_norm(X, W, 1.0)
//' 
// [[Rcpp::export]]
arma::mat weighted_norm(
    const arma::mat & X, 
    const arma::mat & W,
    double p = 1.0
  ) {
  
  unsigned int n = X.n_rows, i, j, K = X.n_cols;
  arma::mat ans(n, n);
  
  // Checking size
  if ( (K != W.n_rows) || (K != W.n_cols) )
    stop("W must be a matrix of size %ix%i", K, K);
  
  for (i = 0u; i < n; i++)
    for (j = i; j < n; j++) {
      ans.at(i,j) = std::pow(arma::dot((X.row(i) - X.row(j))* W, (X.row(i) - X.row(j))), .5*p);
      if (i != j) ans.at(j,i) = ans.at(i,j);
    }
    
    return ans;
}
