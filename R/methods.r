#' @export
#' @rdname blop
as.matrix.blopmatch_match <- function(x, ...) {
  # Making space
  n <- as.integer(nrow(x$X))
  m <- methods::new("dgCMatrix", Dim=c(n, n), p=rep(0L, n + 1L))
  
  g <- lapply(x$matches, "[[", "against")
  l <- lapply(x$matches, "[[", "lambda")
  
  for (i in 1:n) {
    m[i, as.vector(g[[i]]), drop=FALSE] <- l[[i]]
  }
  
  m
}

#' Treatment Effect Estimators using BLOP
#' 
#' @param x An object of class \code{\link[blop]{blopmatch_match}}
#' @template mat
#' @templateVar y 1
#' @templateVar X 1
#' @templateVar Treat 1
#' @author George G. Vega Yon
#' @export
teffect <- function(X, y, effect = "ate", ...) {
  UseMethod("teffect")
}

#' @export
#' @rdname teffect
teffect.default <- function(X, y, Treat, effect = "ate", ...) {
  
}

#' @export
#' @rdname teffect
teffect.blopmatch_match <- function(x, y, effect = "ate", tol = 1e-20, X = NULL, ...) {
  
  # Subseting
  ids       <- which(x$Treat != -1)
  x$X       <- x$X[ids,,drop=FALSE]
  x$matches <- x$matches[ids]
  x$Treat   <- x$Treat[ids]
  
  # Checking teffect
  if (effect == "att") {
    ids <- which(x$Treat == 1)
  } else if (effect == "atc") {
    ids <- which(x$Treat == 0)
  } else if (effect == "ate") {
    ids <- 1:nrow(x$X)
  } else 
    stop("Invalid effect.")
  
  # Retrieving matrix
  M <- as.matrix(x)
  M@x[M@x <= tol] <- 0
  M   <- M[ids , , drop=FALSE]
  M <- M / Matrix::rowSums(M, na.rm = TRUE)
  
  yhat <- as.vector(M %*% cbind(y))
  
  estimator <- mean((y[ids] - yhat)*ifelse(x$Treat[ids] == 1, 1, -1), na.rm=TRUE)
  
  estimator
}
