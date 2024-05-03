#-------------------------------------------------------------------
# %#     Function that implements the OLS estimation of the
# %#     fully-saturated regression via lm() on the corresponding
# %#     subsamples generated via filter.ols.sreg()
#-------------------------------------------------------------------
lm.iter.sreg <- function(Y, S, D, X)
#-------------------------------------------------------------------
{
  theta.list <- rep(list(matrix(NA, ncol = ncol(X), nrow = max(S))), (max(D) + 1))

  for (d in 0:max(D))
  {
    for (s in 1:max(S))
    {
      data.filtered <- subsample.ols.sreg(Y, S, D, X, s, d)
      data.X <- data.filtered[, 4:(4 + ncol(X) - 1)]
      data.filtered.adj <- data.frame(Y = data.filtered$Y, data.X)
      result <- tryCatch(
        {
          lm(Y ~ ., data = data.filtered.adj)
        },
        error = function(e) {
          NULL
        }
      )
      if (is.null(result)) {
        theta.list[[d + 1]][s, ] <- rep(NA, ncol(X))
      } else {
        theta.list[[d + 1]][s, ] <- coef(result)[2:(2 + ncol(X) - 1)]
      }
    }
  }
  list.rtrn <- theta.list
  return(list.rtrn)
}
#-------------------------------------------------------------------
lm.iter.creg <- function(Y, S, D, G.id, Ng, X)
#-------------------------------------------------------------------
{
  if (!is.null(Ng)) {
    working.df <- data.frame(Y, S, D, G.id, Ng, X)
  } else {
    working.df <- data.frame(Y, S, D, G.id, X)
    working.df <- working.df %>%
      group_by(G.id) %>%
      mutate(Ng = n()) %>%
      ungroup() %>%
      select(Y, S, D, G.id, Ng, all_of(names(X)))
    working.df <- as.data.frame(working.df)
  }
  Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
  cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", names(working.df)[6:ncol(working.df)])])
  cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
  data <- cl.lvl.data
  theta.list <- rep(list(matrix(NA, ncol = ncol(X), nrow = max(S))), (max(D) + 1))
  for (d in 0:max(D))
  {
    for (s in 1:max(S))
    {
      data.filtered <- subsample.ols.creg(data, s, d)
      data.X <- data.filtered[, 6:(6 + ncol(X) - 1)]
      data.filtered.adj <- data.frame(Y.bar.Ng = data.filtered$Y.bar * data.filtered$Ng, data.X)
      result <- tryCatch(
        {
          lm(Y.bar.Ng ~ ., data = data.filtered.adj)
        },
        error = function(e) {
          NULL
        }
      )
      if (is.null(result)) {
        theta.list[[d + 1]][s, ] <- rep(NA, ncol(X))
      } else {
        theta.list[[d + 1]][s, ] <- coef(result)[2:(1 + ncol(X))]
      }
    }
  }
  list.rtrn <- list(
    "theta.list"  = theta.list,
    "cl.lvl.data" = data
  )
  return(list.rtrn)
}
