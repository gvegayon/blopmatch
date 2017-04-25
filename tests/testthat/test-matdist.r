context("Mat dist")

# Generating fake data
X <- cbind(
  c(0:5, 0),
  c(1:6, 1)
)

X[3,] <- X[1,]
X[4, 1] <- X[1,1]

# ------------------------------------------------------------------------------
test_that("Exact match", {
  ans <- matching_group(rep(-1, 7), X)
  expect_equivalent(ans[[1]][,1], c(3, 7))
  expect_equivalent(ans[[3]][,1], c(1, 7))
  expect_equivalent(ans[[7]][,1], c(1, 3))
})

# ------------------------------------------------------------------------------
test_that("Quadform distance", {
  # Diagonal
  ans0 <- weighted_norm(X, diag(2), 1)
  ans1 <- as.matrix(dist(X))
  
  expect_equivalent(ans0, ans1)
  
  # Mahalanobis
  generalized_normR <- function(X, W) {
    d <- matrix(NA, nrow=nrow(X), ncol=nrow(X))
    for (i in 1:nrow(X)) {
      for (j in i:nrow(X)) {
        val <- X[i,,drop=FALSE] - X[j,,drop=FALSE]
        d[i,j] <- sqrt(val %*% W %*% t(val))
        d[j,i] <- d[i,j]
      }
    }
    
    d
  }
  
  W <- solve(var(X))
  
  ans0 <- weighted_norm(X, W)
  ans1 <- weighted_norm(X, W)
  
  expect_equivalent(ans0, ans1)
  
})


