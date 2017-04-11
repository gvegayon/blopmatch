# Solvers

#' @importFrom Rglpk Rglpk_solve_LP
#' @importFrom lpSolveAPI make.lp set.column set.objfn set.constr.type set.rhs
#' set.bounds solve.lpExtPtr get.objective get.variables get.constraints
NULL

#' @useDynLib blopmatch
#' @importFrom Rcpp sourceCpp
NULL

# Importing from the Matrix pkg ------------------------------------------------

#' @import methods
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix Matrix rowSums colSums
#' @importMethodsFrom Matrix t
NULL

# Imports from R CODE

#' @importFrom stats dist
#' @importFrom graphics plot text
NULL