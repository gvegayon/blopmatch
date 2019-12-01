#include <Rcpp.h>
using namespace Rcpp;

//' Preprocess data identifying duplicates
//' 
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
  
  std::vector< std::vector< int >> group(n);
  std::vector< int > marked;
  
  marked.reserve(n);
  
  std::vector< bool > checked(n);
  std::fill(checked.begin(), checked.end(), false);
    
  
  
  // Creating the index
  std::vector< int > idx(n);
  for (int i=0; i < n; ++i)
    idx.at(i) = i;
  
  int last = n;
  for (int i = 0; i < n; ++i) {
    
    // Has been checked?
    if (checked[i])
      continue;
    group[i].push_back(idx[i] + (zeroindex? 1 : 0));
    checked[idx[i]] = true;
    
    // Assigning
    marked.push_back(idx[i] + (zeroindex? 1 : 0));
    
    // Are we in the last?
    if (last == 0)
      break;
    
    // Comparing all the members available
    for (int j = 0; j < last; ++j) {
      
      // If seen before, then cont
      if (checked[idx[j]])
        continue;
      
      // // If the same, then cont
      // if (idx[i] == idx[j])
      //   continue;
      
      // Are these equal?
      LogicalVector test = all(x.row(idx[i]) == x.row(idx[j]));
      if (test.at(0)) {
        
        // Adding it to the group, and marking it as checked
        group[i].push_back(idx[j] + (zeroindex? 1 : 0));
        checked[idx[j]] = true;
        
        
      }
      
    }

    // Moving them to the tail
    for (int j = 0; j < group[i].size(); ++j) {
      
      // Rprintf("i: %i, j: %i, last: %i\n", i, j, last);
      // print(wrap(ans.at(idx[i])));
      
      // Are we in the last?
      if (last == 0)
        break;
      
      std::swap(
        idx.at(group[i][j] - (zeroindex? 1 : 0)),
        idx.at((last-- - 1))
        );
      
    }
    
  }
  
  // Reallocating memory, i.e. removing unused capacity.
  marked.shrink_to_fit();
  group.resize(marked.size());
  
  return List::create(
    _["groups"]       = wrap(group),
    _["selected_ids"] = wrap(marked),
    _["n"]            = n
  );
  
}
