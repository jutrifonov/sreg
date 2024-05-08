#-------------------------------------------------------------------
# %#     Auxiliary function providing the appropriate data.frame
# %#     for the subsequent iterative OLS estimation. Takes into
# %#     account the number of observations and creates indicators.
#-------------------------------------------------------------------
subsample.ols.sreg <- function(Y, S, D, X, s, d)
#-------------------------------------------------------------------
{
  X <- as.matrix(X)
  data <- data.frame(Y, S, D, X)
  keep.s <- s
  keep.d <- d
  filtered.data <- data[data$D %in% keep.d & data$S %in% keep.s, ]
  data.ols <- filtered.data
  return(data.ols)
}
#-------------------------------------------------------------------
subsample.ols.creg <- function(data, s, d)
#-------------------------------------------------------------------
{
  keep.s <- s
  keep.d <- d
  filtered.data <- data[data$D == keep.d & data$S == keep.s, ]
  data.ols <- filtered.data
  return(data.ols)
}
