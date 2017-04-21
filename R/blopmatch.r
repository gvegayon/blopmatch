#' @title Bilevel Optimization Problem Matching Solvers
#' 
#' @description 
#' \code{blop} matches all rows \eqn{i} in \code{X} to \eqn{j \neq i}{j!=i}.
#' This function is a wrapper of both \code{blopi_lpsolve} and \code{blopi_glpk}.
#' 
#' @param X A matrix of size \eqn{N\times K}{N * K}. Matching parameters.
#' @param solver A character scalar. Either \code{"glpk"} or \code{"lpsolve"}.
#' @param xi Numeric vector of length \eqn{k}. Ith individual covariates.
#' @author George G. Vega Yon
#' 
#' @return In the case of \code{blop}, a list of class \code{blopmatch_match}:
#' 
#' \item{matches}{A list of size \code{N} with the BLOP solutions}
#' \item{relaxed}{An integer vector indicating in which of the \eqn{n} the
#' constraints were not binding (infeasibility of the problem).}
#' \item{X}{Original matrix \code{X} passed to the function.}
#' \item{X_pred}{A matrix of size \eqn{N\times K}{N * K} with the predicted values.}
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
blop <- function(X, solver="glpk") {
  # Solving the BLOP for each one
  iseq <- 1L:nrow(X)
  
  # Solving the problem
  ans <- if (solver == "glpk")
    lapply(iseq, function(i, covars) blopi_glpk(covars[i,,drop=TRUE], covars[-i,,drop=FALSE]), covars=X)
  else if (solver == "lpsolve")
    lapply(iseq, function(i, covars) blopi_lpsolve(covars[i,,drop=TRUE], covars[-i,,drop=FALSE]), covars=X)
  else stop("-solver- should be either 'glpk' or 'lpsolve'.")
  
  # Checking out which ones were't solved perfectly
  relaxed <- which(colSums(sapply(ans, "[[", "slack")) != 0)
  
  structure(
    list(
      matches = ans,
      relaxed = relaxed,
      X       = X,
      X_pred  = do.call("rbind", lapply(iseq, function(i) {
        l <- ans[[i]]$lambda
        matrix((l/sum(l)) %*% X[-i,,drop=FALSE], nrow=1)
      }))
    ),
    class = "blopmatch_match"
  )
}


#' \code{blopi_lpsolve} solves the problem using the \pkg{lpSolverAPI} R package
#' which implements the lp_solver library.
#' @rdname blop
#' @export
blopi_lpsolve <- function(xi, X) {
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
  D <- apply(X, 1, function(x) dist(rbind(x, xi)))
  lpSolveAPI::set.objfn(
    my.lp, c(D, rep(1e3, K + 1))
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
blopi_glpk <- function(xi, X) {
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
  D <- apply(X, 1, function(x) stats::dist(rbind(xi, x)))
  
  # Solving the LP
  ans <- Rglpk::Rglpk_solve_LP(
    obj   = c(D, rep(1e3, K + 1)),               # |lambda| + |slack|
    mat   = rbind(
      c(rep(1, N), rep(-1, K + 1)),        # sum(lambda + slack) = 1
      cbind(
        t(X),                    # Proj(X) = X
        matrix(-1, nrow = K, ncol = (K + 1))     # Slack vars
      )
    ),
    dir    = rep("==", 3),
    rhs    = c(1, xi),
    bounds = list(
      lower = list(ind = 1L:N, val = rep(0, N))
    )
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
  
  N        <- nrow(x$X)
  Nrelaxed <- length(x$relaxed) 
  K        <- ncol(x$X)
  
  
  cat(sprintf("BILEVEL OPTIMIZATION MATCHING PROBLEM\n"))
  cat(
    sep="",
    sprintf("%% of perfect matches: %.2f%%\n", (1 - Nrelaxed/N) *100),
    sprintf("N: %i, K: %i\n", N, K)
  )
}

#' @export
plot.blopmatch_match <- function(x, y=1:min(2, ncol(x$X)), ...) {
  
  plot(x$X[,y,drop=FALSE], pch=20, col="lightgray")
  if (length(x$relaxed))
    text(x$X[x$relaxed,y,drop=FALSE], labels = x$relaxed, col="red")
  
}
