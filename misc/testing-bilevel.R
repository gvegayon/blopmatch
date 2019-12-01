library(magrittr)
library(blopmatch)

set.seed(1144)

k <- 2
n <- 20
xi <- runif(k) # cbind(1,1)
X  <- matrix(runif(k*n), ncol=k) # matrix(c(.5,.8,.5,0), ncol=2, byrow = TRUE)
D  <- Map(
  function(xj) weighted_norm(rbind(xi, xj), diag(k))[1,2],
  xj = lapply(1:n, function(i) X[i,])
  ) %>% unlist

          
ans0 <- blopi_lpsolve(xi, X, D)
ans1 <- blopi_glpk(xi, X, D)

# microbenchmark::microbenchmark(
#   lpsolve = blopi_lpsolve(xi, X, D),
#   glpk = blopi_glpk(xi, X, D)
# )

ans0$lambda %*% X %>% as.matrix() %>% rbind(., xi) %>% dist
ans1$lambda %*% X %>% as.matrix() %>% rbind(., xi) %>% dist

# Objective vs matched
xihat <- ans0$lambda %*% X

plot(x=X[,1], y=X[,2], xlim = c(0,1), ylim = c(0,1))
points(x=xi[1], y=xi[2], col="red", pch=2)
points(x=xihat[1], y=xihat[2], col="blue", pch=3)
points(x=ans1$xi_feasible[1], y=ans1$xi_feasible[2], col="green", pch=4)

ans0$constr
ans1$constr

weighted_norm(
  as.matrix(rbind(rbind(xi, xihat), X)),
  diag(k))[1,]


data(nsw)
# nsw <- lalonde
nsw$re75 <- scale(nsw$re75)
X <- as.matrix(subset(
  nsw, select = c("age", "education", "black", "hispanic", "married", "nodegree","re75"))
)
bsol <- blop(
  X     = X,
  Treat = nsw$treat,
  p     = 2,
  # W     = solve(var(X)),
  # solver = "lpsolve"
  )


teffect(bsol, nsw$re78, "ate")
teffect(bsol, nsw$re78, "att")

# > teffect(bsol, nsw$re78, "ate")
# [1] 938.4912
# > teffect(bsol, nsw$re78, "att")
# [1] 614.415

