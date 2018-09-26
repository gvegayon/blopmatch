library(magrittr)
library(blopmatch)

set.seed(1144)

k <- 2
n <- 100
xi <- runif(k) # cbind(1,1)
X  <- matrix(runif(k*n), ncol=k) # matrix(c(.5,.8,.5,0), ncol=2, byrow = TRUE)
D  <- Map(
  function(xj) weighted_norm(rbind(xi, xj), diag(k))[1,2],
  xj = lapply(1:n, function(i) X[i,])
  ) %>% unlist

          
ans <- blopi_lpsolve(xi, X, D)

f <- Rglpk::Rglpk_read_file("misc/example1.lp", "MPS_free")
ans2 <- Rglpk::Rglpk_solve_LP(
  f$objective, f$constraints[[1]], f$constraints[[2]],
  f$constraints[[3]], f$bounds, f$types, f$maximum)

ans$lambda %*% X %>% as.matrix() %>% rbind(., xi) %>% dist
ans2$solution %*% X %>% as.matrix() %>% rbind(., xi) %>% dist

# Objective vs matched
xihat <- ans$lambda %*% X

plot(x=X[,1], y=X[,2], xlim = c(0,1), ylim = c(0,1))
points(x=xi[1], y=xi[2], col="red", pch=2)
points(x=xihat[1], y=xihat[2], col="blue", pch=3)
points(x=ans$xi_feasible[1], y=ans$xi_feasible[2], col="green", pch=4)


weighted_norm(
  as.matrix(rbind(rbind(xi, xihat), X)),
  diag(k))[1,]
