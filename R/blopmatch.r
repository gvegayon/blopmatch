#' @title Bilevel Optimization Problem Matching Solvers
#' 
#' @description 
#' \code{blop} matches all rows \eqn{i} in \code{X} to \eqn{j \neq i}{j!=i}.
#' This function is a wrapper of both \code{blopi_lpsolve} and \code{blopi_glpk}.
#' 
#' @template mat
#' @templateVar X 1
#' @templateVar Treat 1
#' @templateVar W 1
#' @templateVar p 1
#' @templateVar exact 1
#' @param solver A character scalar. Either \code{"glpk"} or \code{"lpsolve"}.
#' @param xi Numeric vector of length \eqn{k}. Ith individual covariates.
#' @param D Numeric vector. Distances vector (internal use only).
#' @author George G. Vega Yon
#' 
#' @details Both \code{W} and \code{p} are passed to \code{\link{weighted_norm}}.
#' 
#' @return In the case of \code{blop}, a list of class \code{blopmatch_match}:
#' 
#' \item{matches}{A list of size \code{N} with the BLOP solutions}
#' \item{status}{An integer vector of length \eqn{N}. A value equal to one indicates
#' that the problem was binded.}
#' \item{X}{Original matrix \code{X} passed to the function.}
#' 
#' @export
#' @examples 
#' set.seed(123)
#' 
#' X <- matrix(rnorm(10*2), ncol=2)
#' ans <- blop(X)
#' 
#' plot(ans)
#' 
#' data(lalonde, package="MatchIt")
#' 
blop <- function(
  X, 
  Treat = rep(-1, nrow(X)), 
  exact = NULL, 
  solver="glpk",
  W = diag(ncol(X)),
  p = 1
  ) {
  # Solving the BLOP for each one
  iseq <- 1L:nrow(X)
  
  # Selecting the matching group
  groups <- matching_group(Treat, exact)
  
  # Computing Distances
  D <- weighted_norm(X = X, W = W, p = p)
  
  # Solving the problem
  ans <- if (solver == "glpk")
    lapply(iseq, function(i, covars) {
      
      # Checking if it can be matched or not
      if (!length(groups[[i]])) {
        warning("Row ", i, " couldn't be matched.")
        return(NULL)
      }
      
      # Solving the blop
      c(
        blopi_glpk(
          xi = covars[i,,drop=TRUE],
          X  = covars[groups[[i]],,drop=FALSE],
          D = D[i, groups[[i]]]
          ),
        list(against = groups[[i]])
      )
      
      
      
      }, covars=X)
  else if (solver == "lpsolve")
    lapply(iseq, function(i, covars) {
      
      # Checking if it can be matched or not
      if (!length(groups[[i]])) {
        warning("Row ", i, " couldn't be matched.")
        return(NULL)
      }
      
      # Solving the blop
      c(
        blopi_lpsolve(
          xi = covars[i,,drop=TRUE],
          X  = covars[groups[[i]],,drop=FALSE],
          D  = D[i, groups[[i]]]
          ),
        list(against = groups[[i]])
      )
      
      }, covars=X)
  else stop("-solver- should be either 'glpk' or 'lpsolve'.")
  
  # Checking out which ones were't solved perfectly
  status <- sapply(ans, function(x) {
    if (!length(x)) -1
    else if (any(x[["slack"]] != 0)) 1
    else 0
  })
  
  structure(
    list(
      matches = ans,
      status  = status,
      X       = X,
      Treat   = Treat,
      exact   = exact,
      solver  = solver,
      W       = W,
      p       = p
    ),
    class = "blopmatch_match"
  )
}


#' \code{blopi_lpsolve} solves the problem using the \pkg{lpSolverAPI} R package
#' which implements the lp_solver library.
#' @rdname blop
#' @export
blopi_lpsolve <- function(xi, X, D = NULL) {
  # Constraints:
  #  sum(lambda) = 1 : 1
  #  lambda*xk' = xk : k 
  #  TOTAL: 1 + k
  #
  # Variables:
  #  lambda: n
  #  slack variables: in the case of infeasibility: k + 1
  
  N <- nrow(X)
  K <- ncol(X)
  
  my.lp <- lpSolveAPI::make.lp(1 + K, N + (K + 1))
  
  # Setting columns ------------------------------------------------------------
  # lambda columns
  for (j in 1:N)
    lpSolveAPI::set.column(my.lp, j, c(
      1, X[j, , drop=TRUE]
    ))
  
  # Slack variables columns 
  for (j in (N + 1):(N + (K + 1)))
    lpSolveAPI::set.column(my.lp, j, c(
      rep(-1, K + 1)
    ))
  
  # Objective function ---------------------------------------------------------
  #  lambda*distance
  if (!length(D))
    D <- apply(X, 1, function(x) dist(rbind(x, xi)))
  
  lpSolveAPI::set.objfn(
    my.lp, c(D, rep(1, K + 1))
  )
  
  # Constraints ----------------------------------------------------------------
  # Recall: We have 
  # sum(lambda) == 1: 1
  # lambda*x_k' = x_k == k
  # {1 = "<=", 2 = ">=", 3 = "="}
  lpSolveAPI::set.constr.type(my.lp, c(3, rep(3, K)))
  lpSolveAPI::set.rhs(my.lp, c(1, xi))
  lpSolveAPI::set.bounds(my.lp, rep(0, N + (K + 1)))
  
  # Solving the problem
  ans <- lpSolveAPI::solve.lpExtPtr(my.lp)
  
  structure(
    list(
      obj    = lpSolveAPI::get.objective(my.lp),
      lambda = methods::as(
        matrix(lpSolveAPI::get.variables(my.lp)[1:N], nrow=1),
        "dgCMatrix"),
      slack  = lpSolveAPI::get.variables(my.lp)[(N + 1):(N + K + 1)],
      constr = lpSolveAPI::get.constraints(my.lp),
      status = ans,
      xi     = xi
    ), class = "blopmatch_matchi"
  )
  
}

#' \code{blopi_glpk} solves the problem using the \pkg{Rglpk} R package
#' which implements the GNU Linear Programming Kit.
#' @rdname blop
#' @export
blopi_glpk <- function(xi, X, D = NULL) {
  # Constraints:
  #  sum(lambda) = 1 : 1
  #  lambda*xk' = xk : k 
  #  TOTAL: 1 + k
  #
  # Variables:
  #  lambda: n - 1 
  #  slack variables: in the case of infeasibility: k + 1
  
  N <- nrow(X)
  K <- ncol(X)
  
  # Objective function ---------------------------------------------------------
  #  lambda*distance
  if (!length(D))
    D <- apply(X, 1, function(x) stats::dist(rbind(xi, x)))
  
  # Solving the LP
  ans <- Rglpk::Rglpk_solve_LP(
    obj   = c(D, rep(1, K + 1)),               # |lambda| + |slack|
    mat   = rbind(
      c(rep(1, N), rep(-1, K + 1)),        # sum(lambda + slack) = 1
      cbind(
        t(X),                    # Proj(X) = X
        matrix(-1, nrow = K, ncol = (K + 1))     # Slack vars
      )
    ),
    dir    = rep("==", K + 1),
    rhs    = c(1, xi),
    bounds = list(
      lower = list(ind = 1L:N, val = rep(0, N))
    ),
    control = list(presolve = FALSE, tm_limit = 500)
  )
  
  structure(
    list(
      obj    = ans$objval,
      lambda = methods::as(
        matrix(ans$solution[1L:N], nrow=1),
        "dgCMatrix"
        ),
      slack  = ans$solution[(N + 1L):(N + K + 1L)],
      constr = NA,
      status = ans$status,
      xi     = xi
    ), class = "blopmatch_matchi"
  )
  
}


#' @export
print.blopmatch_match <- function(x, ...) {
  
  N      <- nrow(x$X)
  binded <- sum(x$status) 
  K      <- ncol(x$X)
  
  
  cat(sprintf("BILEVEL OPTIMIZATION MATCHING PROBLEM\n"))
  cat(
    sep="",
    sprintf("%% of perfect matches: %.2f%%\n", (1 - binded/N) *100),
    sprintf("N: %i, K: %i\n", N, K)
  )
}

#' @export
plot.blopmatch_match <- function(x, y=1:min(2, ncol(x$X)), ...) {
  
  plot(x$X[,y,drop=FALSE], pch=20, col="lightgray")
  binded <- which(x$status == 1)
  if (length(binded))
    text(x$X[binded,y,drop=FALSE], labels = binded, col="red")
  
}
