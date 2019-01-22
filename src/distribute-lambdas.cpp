#include <Rcpp.h>
using namespace Rcpp;

//' Preprocess data identifying duplicates
//' This function identifies duplicates, groups them up and returns a vector
//' of indices of the observations to be used in the model, i.e. a new dataset
//' without duplicates.
//' @param x Numeric matrix. Covariates over which duplicates should be tagged.
//' @param zeroindex Logical. Internal use.
//' @return
//' A list with the following elements:
//' - `groups` A list with the ids of duplicates grouped up.
//' - `selected_ids` An integer vector with the ids to be used.
//' - `n` `length(selected_ids)`.
//' @export
// [[Rcpp::export]]
List tag_duplicates(const NumericMatrix & x, bool zeroindex = true) {
  
  // Creating the output
  int n = x.nrow();
  std::vector< std::vector< int >> ans(x.nrow());
  std::vector< int > marked;
  marked.reserve(n);
  
  // Creating the index
  std::vector<int> idx(n);
  for (int i=0; i< n; i++)
    idx.at(i) = i;
  
  int last = n;
  for (int i = 0; i < n; ++i) {
    
    // Assigning
    marked.push_back(idx[i] + (zeroindex? 1 : 0));
    
    // Are we in the last?
    if (last == i)
      break;
    
    for (int j = 0; j < last; ++j) {
      
      if (idx[i] == idx[j])
        continue;
      
      // Are these equal?
      LogicalVector test = all(x.row(idx[i]) == x.row(idx[j]));
      if (test.at(0)) 
        ans.at(idx[i]).push_back(idx[j] + (zeroindex? 1 : 0));
      
    }
    
    // Moving them to the tail
    for (int j = 0; j < ans.at(idx[i]).size(); ++j) {
      
      std::swap(
        idx.at(ans.at(idx[i]).at(j) - (zeroindex? 1 : 0)),
        idx.at((last-- - 1))
        );
      
    }
    
  }
  
  // Reallocating memory, i.e. removing unused capacity.
  marked.shrink_to_fit();
  
  return List::create(
    _["groups"] = wrap(ans),
    _["selected_ids"]  = wrap(marked),
    _["n"]       = n
  );
  
}


