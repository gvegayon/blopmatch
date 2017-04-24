#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Finds matching group
//' @template mat
//' @templateVar X 1
//' @templateVar Treat 1
//' @param zeroindex Integer scalar. Added to the resulting index.
//' @return A list of length \eqn{n} specifying for each i to which individuals
//' each j != i will be matched.
//' @export 
//' 
//' @examples
//' # Asigning matching groups for lalonde
//' 
//' data(lalonde, package = "MatchIt")
//' dat <- lalonde
//' 
//' # Match (no exact)
//' m_noexact <- matching_group(cbind(rep(1, nrow(dat))), dat$treat)
//' 
//' # How many matches?
//' table(sapply(m_noexact, length))
//' table(dat$treat)
//' 
//' # What if we ask exact match on black
//' 
//' m_exact <- matching_group(
//'   as.matrix(dat[,"black",drop=FALSE]),
//'   dat$treat)
//' 
//' table(sapply(m_exact, length))
//' with(dat, table(treat, black))
// [[Rcpp::export]]
List matching_group(
  const arma::mat & X,
  const arma::ivec & Treat,
  int zeroindex = 1
) {
  
  unsigned int n = X.n_rows, i, j;
  unsigned int K = X.n_cols, k;

  std::vector< std::vector< arma::uword > > ans0(n);

  for (i = 0u; i < n; i++) 
    for (j = i + 1u; j < n; j++) {
      
      // Initializing current search
      bool match;    
  
      // Comparing i to j  
      match = true;
      for (k = 0u; k < K; k++) {
        
        // Opposit groups?
        if (Treat.at(i) >= 0 && Treat.at(i) == Treat.at(j)) {
          match = false;
          break;
        }
        
        // If any attribute in E is different, then 
        // no exact and return false
        if (X.at(i, k) != X.at(j, k)) {
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

//' Computes Quadratic Form Distance
//' @template mat
//' @templateVar X 1
//' @templateVar W 1
//' @param p Numeric scalar.
//' @export
//' 
//' 
// [[Rcpp::export]]
arma::mat generalized_norm(
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
