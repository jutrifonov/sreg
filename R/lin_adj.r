#-------------------------------------------------------------------
# %#     Function that implements the calculation of \hat{\mu} --
# %#     i.e., calculates linear adjustments
#-------------------------------------------------------------------
lin.adj.sreg <- function(a, S, X, model)
#-------------------------------------------------------------------
{
  data <- data.frame(S, X)
  theta.mtrx <- model[[a + 1]]
  theta.vec.matched <- theta.mtrx[data$S, ]
  mu.hat <- diag(as.matrix(X) %*% t(theta.vec.matched))
  return(mu.hat)
}
#-------------------------------------------------------------------
lin.adj.creg <- function(a, data, model)
#-------------------------------------------------------------------
{
  X.data <- data[, 6:ncol(data)]
  theta.mtrx <- model$theta.list[[a + 1]]
  theta.vec.matched <- theta.mtrx[data$S, ]
  mu.hat <- diag(as.matrix(X.data) %*% t(theta.vec.matched))

  return(mu.hat)
}
