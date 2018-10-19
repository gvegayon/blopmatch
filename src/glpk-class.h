#include <Rcpp.h>
#include <glpk.h>

using namespace Rcpp;

// "Row indices of each element are stored in the array ia, column indices are stored in
// the array ja, and numerical values of corresponding elements are stored in the array ar."
// So... data must be presented as:
// (Row, Column, Data) arrays of vectors starting at 1 (not 0)

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
    const NumericVector & coef,
    const NumericMatrix & subj_lhs,
    const NumericVector & subj_rhs,
    int DIR = GLP_MAX
  );
  
  // Cleaning up
  ~glpkObj() {
    delete [] ia;
    delete [] ja;
    delete [] ar;
    glp_delete_prob(lp);
  };
  
  // void set_col_bnds() {glp_set_col_bnds(lp, ++i, GLP_LO, 0.0, 0.0);}
  template<typename T> void set_obj_coef(T & coef);
  
  // Column bounds
  template<typename Ttype, typename Tlb, typename Tub>
  void set_col_bnds(Ttype & type, Tlb & lb, Tub & ub);
  
  void set_col_bnds(int type, int lb, int ub);
  
  // Row bounds
  template<typename Ttype, typename Tlb, typename Tub>
  void set_row_bnds(Ttype & type, Tlb & lb, Tub & ub);
  
  template<typename Tub> void set_row_bnds(int type, int lb, Tub & ub);
  template<typename Tlb> void set_row_bnds(int type, Tlb & lb, int ub);
  
  void set_row_bnds(int type, int lb, int ub);
  
  // Simple methods ------------------------------------------------------------
  
  int add_rows(int nrs) {return glp_add_rows(lp, nrs);}
  int add_cols(int ncs) {return glp_add_cols(lp, ncs);}
  void simplex() {glp_simplex(lp, &param);}
  void set_obj_dir(int DIR) {glp_set_obj_dir(lp, DIR);}
  
  // Retrieving data -----------------------------------------------------------
  List getSol();
  
  
private:
  glp_prob *lp;
  glp_smcp param;
  int * ia;
  int * ja;
  double * ar;
  int nvars;
  int nconstraints;
  
};

// Setting the coefficients of the objective function --------------------------
template<typename T> void glpkObj::set_obj_coef(T & coef) {
  
  int i = 0;
  for (typename T::const_iterator it = coef.begin(); it != coef.end(); ++it)
    glp_set_obj_coef(lp, ++i, *it);
  
}

// Setting the colum bounds ----------------------------------------------------
template<typename Ttype, typename Tlb, typename Tub>
void glpkObj::set_col_bnds(Ttype & type, Tlb & lb, Tub & ub) {
  
  for (int i = 0; i < nvars; ++i) 
    glp_set_col_bnds(lp, i + 1, type[i], lb[i], ub[i]);
  
}

void glpkObj::set_col_bnds(int type, int lb, int ub) {
  
  for (int i = 0; i < nvars; ++i) 
    glp_set_col_bnds(lp, i + 1, type, lb, ub);
  
}

// Setting the row bounds ------------------------------------------------------
template<typename Ttype, typename Tlb, typename Tub>
void glpkObj::set_row_bnds(Ttype & type, Tlb & lb, Tub & ub) {
  
  for (int i = 0; i < nconstraints; ++i) 
    glp_set_row_bnds(lp, i + 1, type[i], lb[i], ub[i]);
  
}

template<typename Tub>
void glpkObj::set_row_bnds(int type, int lb, Tub & ub) {
  
  for (int i = 0; i < nconstraints; ++i) 
    glp_set_row_bnds(lp, i + 1, type, lb, ub[i]);
  
}

template<typename Tlb>
void glpkObj::set_row_bnds(int type, Tlb & lb, int ub) {
  
  for (int i = 0; i < nconstraints; ++i) 
    glp_set_row_bnds(lp, i + 1, type, lb[i], ub);
  
}

void glpkObj::set_row_bnds(int type, int lb, int ub) {
  
  for (int i = 0; i < nconstraints; ++i) 
    glp_set_row_bnds(lp, i + 1, type, lb, ub);
  
}

// Defining the constructor ---------------------------------------------------
glpkObj::glpkObj(
  const NumericVector & coef,
  const NumericMatrix & subj_lhs,
  const NumericVector & subj_rhs,
  int DIR
) {
  
  // Creating arrays of the given size
  nvars        = subj_lhs.ncol();
  nconstraints = subj_lhs.nrow();
  
  ia = new int[subj_lhs.size() + 1];
  ja = new int[subj_lhs.size() + 1];
  ar = new double[subj_lhs.size() + 1];
  
  // Initializing the problem
  lp = glp_create_prob();
  glp_set_obj_dir(lp, DIR);
  glp_init_smcp(&param);
  param.msg_lev = GLP_MSG_ERR;
  
  // We want to maximize
  this->add_rows(nconstraints);
  this->set_row_bnds(GLP_UP, 0.0, subj_rhs); 
  this->add_cols(nvars);
  
  // Setting bounds and objective coefficient
  this->set_col_bnds(GLP_LO, 0.0, 0.0);
  this->set_obj_coef(coef);
  
  
  // Creating arrays from the LHS.
  as_glpk_array(subj_lhs, ia, ja, ar);
  glp_load_matrix(lp, subj_lhs.ncol()*subj_lhs.nrow(), ia, ja, ar);
  
}

// Returns the solutions
List glpkObj::getSol() {
  
  double val = glp_get_obj_val(lp); 
  
  NumericVector par(nvars);
  for (int i =0;i<nvars; i++)
    par.at(i) = glp_get_col_prim(lp, i+1);
  
  return List::create(
    _["val"] = val,
    _["par"] = par
  );
  
}

