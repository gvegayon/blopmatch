#include <Rcpp.h>
#include <glpk.h>

using namespace Rcpp;

// Function to arrange a numeric matrix as the desired form of array for GLPK.
void as_glpk_array(
    const NumericMatrix & x,
    int *ai, int *aj,
    double *av
) {
  
  int n = x.nrow();
  int k = x.ncol();
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < k; j++) {
      ai[i + 1 + j*k] = i + 1;
      aj[i + 1 + j*k] = j + 1;
      av[i + 1 + j*k] = x.at(i,j);
    }
    
    
    return;
}

class glpkObj {
  
public:
  
  // Declaring the constructor
  glpkObj(
    const NumericVector & obj,
    const NumericMatrix & subj_lhs,
    const NumericVector & subj_rhs
  );
  
  // Cleaning up
  ~glpkObj() {
    delete [] ia;
    delete [] ja;
    delete [] ar;
    glp_delete_prob(lp);
  };
  
  // Solver
  void simplex() {glp_simplex(lp, &param);}
  
  // Getter
  List getSol();
  
  
  
private:
  glp_prob *lp;
  glp_smcp param;
  int * ia;
  int * ja;
  double * ar;
  int k;
  
};

// Defining the constructor
glpkObj::glpkObj(
  const NumericVector & obj,
  const NumericMatrix & subj_lhs,
  const NumericVector & subj_rhs
) {
  
  // Creating arrays of the given size
  k = subj_lhs.ncol();
  ia = new int[subj_lhs.size() + 1];
  ja = new int[subj_lhs.size() + 1];
  ar = new double[subj_lhs.size() + 1];
  
  // We want to maximize
  glp_set_obj_dir(lp, GLP_MAX);
  
  glp_add_rows(lp, subj_lhs.nrow());
  
  int i = 0;
  for (NumericVector::const_iterator it = subj_rhs.begin(); it != subj_rhs.end(); ++it)
    glp_set_row_bnds(lp, ++i, GLP_UP, 0.0, *it);
  
  glp_add_cols(lp, obj.length());
  
  i = 0;
  for (NumericVector::const_iterator it = obj.begin(); it != obj.end(); ++it) {
    // glp_set_col_name(lp, 1, "x1");
    glp_set_col_bnds(lp, ++i, GLP_LO, 0.0, 0.0);
    glp_set_obj_coef(lp, i, *it);
  }
  
  // Creating arrays from the LHS.
  as_glpk_array(subj_lhs, ia, ja, ar);
  glp_load_matrix(lp, subj_lhs.ncol()*subj_lhs.nrow(), ia, ja, ar);
  
}

// Returns the solutions
List glpkObj::getSol() {
  
  double val = glp_get_obj_val(lp); 
  
  NumericVector par(k);
  for (int i =0;i<k; i++)
    par.at(i) = glp_get_col_prim(lp, i+1);
  
  return List::create(
    _["val"] = val,
    _["par"] = par
  );
  
}

