#include <Rcpp.h>
#include <glpk.h>

using namespace Rcpp;

typedef NumericVector::const_iterator nveciter;

// "Row indices of each element are stored in the array ia, column indices are stored in
// the array ja, and numerical values of corresponding elements are stored in the array ar."
// So... data must be presented as:
// (Row, Column, Data) arrays of vectors starting at 1 (not 0)
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

//' A mininmal example
//' @param obj Numeric vector of size \eqn{K}. Coeficients in the objective function.
//' @param subj_lhs Numeric matrix of size \eqn{K\times m}. Constraints.
//' @param pname Character scalar. Name of the LP.
//' @export
//' @examples
//' obj      <- c(10, 6, 4)
//' subj_lhs <- matrix(c(1, 10, 2, 1, 4, 2, 1, 5, 6), ncol = 3)
//' subj_rhs <- (100, 600, 300)
//' glpk_example(obj, subj_lhs, subj_rhs)
// [[Rcpp::export]]
List glpk_example(
    const NumericVector & obj,
    const NumericMatrix & subj_lhs,
    const NumericVector & subj_rhs,
    const StringVector pname = "sample")
{
  // Checking dimensions
  if (obj.size() != subj_lhs.ncol())
    stop("The objective function doesn't have as many columns in there are in constraints.");
  
  if (subj_lhs.nrow() != subj_rhs.size())
    stop("The constraints don't have as many rows as the RHS.");
  
  // Initializing objects and parameters
  glp_prob *lp;
  glp_smcp param;
  glp_init_smcp(&param);
  param.msg_lev = GLP_MSG_ERR;
  
  // Creating arrays of the given size
  int    * ia = new int[subj_lhs.size() + 1];
  int    * ja = new int[subj_lhs.size() + 1];
  double * ar = new double[subj_lhs.size() + 1];
  
  double z, x1, x2, x3;
  
  // Setting up the problem
  lp = glp_create_prob();
  glp_set_prob_name(lp, pname[0]);
  
  // We want to maximize
  glp_set_obj_dir(lp, GLP_MAX);
  
  // We have 3 constriants: First specify the slack variables
  /* type of auxiliary/structural variable: */
// #define GLP_FR             1  /* free (unbounded) variable */
// #define GLP_LO             2  /* variable with lower bound */
// #define GLP_UP             3  /* variable with upper bound */
// #define GLP_DB             4  /* double-bounded variable */
// #define GLP_FX             5  /* fixed variable */
  
  glp_add_rows(lp, subj_lhs.nrow());
  
  int i = 0;
  for (nveciter it = subj_rhs.begin(); it != subj_rhs.end(); ++it)
    glp_set_row_bnds(lp, ++i, GLP_UP, 0.0, *it);
  
  glp_add_cols(lp, obj.length());
  
  i = 0;
  for (nveciter it = obj.begin(); it != obj.end(); ++it) {
    // glp_set_col_name(lp, 1, "x1");
    glp_set_col_bnds(lp, ++i, GLP_LO, 0.0, 0.0);
    glp_set_obj_coef(lp, i, *it);
  }
  
  // Creating arrays from the LHS.
  as_glpk_array(subj_lhs, ia, ja, ar);
  glp_load_matrix(lp, subj_lhs.ncol()*subj_lhs.nrow(), ia, ja, ar);
  
  // Solving the problem and retrieving the obj
  glp_simplex(lp, &param);
  
  NumericVector ans(obj.length());
  
  z = glp_get_obj_val(lp);
  for (int i =0;i<ans.length(); i++)
    ans.at(i) = glp_get_col_prim(lp, i+1);

  
  // Returning and deleting the problem
  // Rprintf("\nz = %g; x1 = %g; x2 = %g\n", z, x1, x2, x3);
  glp_delete_prob(lp);
  delete [] ia;
  delete [] ja;
  delete [] ar;
  
  
  return List::create(
    _["z"] = z,
    _["x"] = ans
  );
}