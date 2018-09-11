library(magrittr)
library(blopmatch)

set.seed(1222)

k <- 2
n <- 10
xi <- runif(k)
X  <- matrix(runif(k*n), ncol=k) 
D  <- Map(
  function(xj) weighted_norm(rbind(xi, xj), diag(k))[1,2],
  xj = lapply(1:n, function(i) X[i,])
  ) %>% unlist

          
ans <- blopi_lpsolve(xi, X, D)

ans$lambda %*% X %>% as.matrix() %>% rbind(., xi) %>% dist

plot(X)
points(xi[1], xi[2], col="red")

# Objective vs matched
ans$lambda %*% X
xi
