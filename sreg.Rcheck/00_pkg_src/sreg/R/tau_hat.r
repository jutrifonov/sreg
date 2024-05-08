#-------------------------------------------------------------------
# %#     Function that implements \hat{\tau} --
# %#     i.e. the ATE estimator
#-------------------------------------------------------------------
tau.hat.sreg <- function(Y, S, D, X=NULL, model=NULL)
#-------------------------------------------------------------------
{
  tau.hat <- numeric(max(D))
  for (d in 1:max(D))
  {
    if (!is.null(X)) {
      data <- data.frame(Y, S, D, X)
      data$pi <- pi.hat.sreg(S, D)[, d]
      data$pi.0 <- pi.hat.sreg(S, D, inverse = T)[, 1]
      data$A <- ifelse(D == d, 1, ifelse(D == 0, 0, -999999))
      data$I <- as.numeric(data$A != -999999)

      mu.hat.d <- lin.adj.sreg(d, data$S, data[4:(4 + ncol(X) - 1)], model)
      mu.hat.0 <- lin.adj.sreg(0, data$S, data[4:(4 + ncol(X) - 1)], model)

      Ksi.vec <- data$I * (((data$A * (data$Y - mu.hat.d)) / data$pi) -
        (((1 - data$A) * (data$Y - mu.hat.0)) / (data$pi.0))) +
        (mu.hat.d - mu.hat.0)

      tau.hat[d] <- mean(Ksi.vec)
    } else {
      data <- data.frame(Y, S, D)
      data$pi <- pi.hat.sreg(S, D)[, d]
      data$pi.0 <- pi.hat.sreg(S, D, inverse = T)[, 1]
      data$A <- ifelse(D == d, 1, ifelse(D == 0, 0, -999999))
      data$I <- as.numeric(data$A != -999999)
      mu.hat.d <- 0
      mu.hat.0 <- 0

      Ksi.vec <- data$I * (((data$A * (data$Y - mu.hat.d)) / data$pi) -
        (((1 - data$A) * (data$Y - mu.hat.0)) / (data$pi.0))) +
        (mu.hat.d - mu.hat.0)

      tau.hat[d] <- mean(Ksi.vec)
    }
  }
  return(tau.hat)
}
#-------------------------------------------------------------------
tau.hat.creg <- function(Y, S, D, G.id, Ng, X=NULL, model=NULL)
#-------------------------------------------------------------------
{
  tau.hat.vec <- numeric(max(D))
  Y.bar.g.list <- rep(list(NA), max(D))
  mu.hat.list <- rep(list(NA), max(D))
  pi.hat.list <- rep(list(NA), max(D))
  data.list <- rep(list(NA), max(D))
  if (!is.null(X)) {
    cl.lvl.data <- model$cl.lvl.data
    data <- cl.lvl.data
    Ng.full <- data$Ng
    Y.bar.full <- data$Y.bar
    for (d in 1:max(D))
    {
      data$pi <- pi.hat.creg(data$S, data$D)[, d]
      data$pi.0 <- pi.hat.creg(data$S, data$D, inverse = T)[, 1]
      data$A <- ifelse(data$D == d, 1, ifelse(data$D == 0, 0, -999999))
      data$I <- as.numeric(data$A != -999999)
      data.list[[d]] <- data
      pi.hat.list[[d]] <- data$pi

      mu.hat.d <- lin.adj.creg(d, data = cl.lvl.data, model)
      mu.hat.0 <- lin.adj.creg(0, data = cl.lvl.data, model)

      Xi.g <- data$I * (((data$A * (Y.bar.full * data$Ng - mu.hat.d)) / data$pi) -
        (((1 - data$A) * (Y.bar.full * data$Ng - mu.hat.0)) / data$pi.0)) +
        (mu.hat.d - mu.hat.0)

      mu.hat.list[[d]] <- as.matrix(cbind(mu.hat.0, mu.hat.d), ncol = 2)

      tau.hat <- mean(Xi.g) / mean(Ng.full)
      tau.hat.vec[d] <- tau.hat
    }
    rtrn.list <- list(
      "tau.hat"   = tau.hat.vec,
      "mu.hat"    = mu.hat.list,
      "pi.hat"    = pi.hat.list,
      "pi.hat.0"  = data$pi.0,
      "data.list" = data.list,
      "Y.bar.g"   = Y.bar.full,
      "Ng"        = Ng.full
    )
  } else {
    if (!is.null(Ng)) {
      working.df <- data.frame(Y, S, D, G.id, Ng)
    } else {
      working.df <- data.frame(Y, S, D, G.id)
      working.df <- working.df %>%
        group_by(G.id) %>%
        mutate(Ng = n()) %>%
        ungroup() %>%
        select(Y, S, D, G.id, Ng)
      working.df <- as.data.frame(working.df)
    }
    Y.bar.full <- aggregate(Y ~ G.id, working.df, mean)$Y
    cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng")])
    Ng.full <- cl.lvl.data$Ng

    for (d in 1:max(D))
    {
      data <- cl.lvl.data
      data$pi <- pi.hat.creg(data$S, data$D)[, d]
      data$pi.0 <- pi.hat.creg(data$S, data$D, inverse = T)[, 1]
      data$A <- ifelse(data$D == d, 1, ifelse(data$D == 0, 0, -999999))
      data$I <- as.numeric(data$A != -999999)
      data.list[[d]] <- data
      pi.hat.list[[d]] <- data$pi

      mu.hat.d <- 0
      mu.hat.0 <- 0

      Xi.g <- data$I * (((data$A * (Y.bar.full * data$Ng - mu.hat.d)) / data$pi) -
        (((1 - data$A) * (Y.bar.full * data$Ng - mu.hat.0)) / data$pi.0)) +
        (mu.hat.d - mu.hat.0)

      mu.hat.list[[d]] <- as.matrix(cbind(mu.hat.0, mu.hat.d), ncol = 2)

      tau.hat <- mean(Xi.g) / mean(Ng.full)
      tau.hat.vec[d] <- tau.hat
    }
    rtrn.list <- list(
      "tau.hat"   = tau.hat.vec,
      "mu.hat"    = mu.hat.list,
      "pi.hat"    = pi.hat.list,
      "pi.hat.0"  = data$pi.0,
      "data.list" = data.list,
      "Y.bar.g"   = Y.bar.full,
      "Ng"        = Ng.full
    )
  }
  return(rtrn.list)
}
