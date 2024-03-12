# %##%##%##%###%##%##%##%###%##%##%##%###%##
### This R file provides the collection ####
### of functions to estimate the ATE    ####
### under CAR with multiple treatments  ####
###        & cluster-level treatment    ####
# %##%##%##%###%##%##%##%###%##%##%##%###%##
####      The code is developed by      ####
####      @Juri Trifonov, UChicago      ####
####            Supervisors:            ####
####      @Azeem Shaikh, UChicago       ####
####    @Max Tabord-Meehan, UChicago    ####
# %##%##%##%###%##%##%##%###%##%##%##%###%##
# %##%##%##%##
# %# v.1.2.5
# %##%##%##%##
#-------------------------------------------------------------------
# %##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##
# %##%##%##%##%##%#      I. ATE  estimator     #%##%##%##%##%##%##%#
# %##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##
#-------------------------------------------------------------------
#------------------------------------------------------------------
# %# (1) Auxiliary function providing the appropriate data.frame
# %#     for the subsequent iterative OLS estimation. Accounts for
# %#     the number of observations and creates indicators.
# %source function for theta.est.str()
#-------------------------------------------------------------------
filter.ols.creg <- function(data, s, d)
#-------------------------------------------------------------------
{
  keep.s <- s
  keep.d <- d
  filtered.data <- data[data$D == keep.d & data$S == keep.s, ]
  data.ols <- filtered.data
  return(data.ols)
}
#-------------------------------------------------------------------
lm.iter.creg <- function(Y, S, D, G.id, Ng, X=NULL, Ng.cov = FALSE)
#-------------------------------------------------------------------
{
  if (Ng.cov == TRUE) {
    working.df <- data.frame(Y, S, D, G.id, Ng)
    Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
    cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng")])
    cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
    data <- cl.lvl.data
    data$X <- NULL
    theta.list <- rep(list(matrix(NA, ncol = 1, nrow = max(S))), (max(D) + 1))
    for (d in 0:max(D))
    {
      for (s in 1:max(S))
      {
        data.filtered <- filter.ols.creg(data, s, d)
        data.filtered.adj <- data.frame(Y.bar.Ng = data.filtered$Y.bar * data.filtered$Ng, Ng = data.filtered$Ng)
        result <- lm(Y.bar.Ng ~ ., data = data.filtered.adj)
        theta.list[[d + 1]][s, ] <- coef(result)[2]
      }
    }
  } else {
    working.df <- data.frame(Y, S, D, G.id, Ng, X)
    Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
    cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", names(working.df)[6:ncol(working.df)])])
    cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
    data <- cl.lvl.data
    theta.list <- rep(list(matrix(NA, ncol = ncol(X) + 1, nrow = max(S))), (max(D) + 1))
    for (d in 0:max(D))
    {
      for (s in 1:max(S))
      {
        data.filtered <- filter.ols.creg(data, s, d)
        data.X <- data.filtered[, 6:(6 + ncol(X) - 1)]
        data.filtered.adj <- data.frame(Y.bar.Ng = data.filtered$Y.bar * data.filtered$Ng, Ng = data.filtered$Ng, data.X)
        result <- lm(Y.bar.Ng ~ ., data = data.filtered.adj)
        theta.list[[d + 1]][s, ] <- coef(result)[2:(2 + ncol(X))]
      }
    }
  }
  list.rtrn <- list(
    "theta.list"  = theta.list,
    "cl.lvl.data" = data
  )
  return(list.rtrn)
}
#-------------------------------------------------------------------
lin.adj.creg <- function(a, data, model, Ng.cov = FALSE)
#-------------------------------------------------------------------
{
  if (Ng.cov == TRUE) {
    theta.mtrx <- model$theta.list[[a + 1]]
    theta.vec.matched <- theta.mtrx[data$S, ]
    Ng.hat <- theta.vec.matched * data$Ng
    mu.hat <- Ng.hat
  } else {
    X.data <- data[, 6:ncol(data)]
    theta.mtrx <- model$theta.list[[a + 1]]
    theta.vec.matched <- theta.mtrx[data$S, ]
    Ng.hat <- theta.vec.matched[, 1] * data$Ng
    X.hat <- diag(as.matrix(X.data) %*% t(theta.vec.matched[, -1]))
    mu.hat <- Ng.hat + X.hat
  }

  return(mu.hat)
}
#-------------------------------------------------------------------
pi.hat.creg <- function(S, D, inverse = FALSE)
#-------------------------------------------------------------------
{
  n <- length(D)
  data <- data.frame(S, D)
  counts <- data %>%
    group_by(S, D) %>%
    summarise(n = n())
  scount <- data %>%
    group_by(S) %>%
    summarise(ns = n())

  j <- left_join(counts, scount, by = join_by(S == S))
  j$pi_hat <- j$n / j$ns
  pi_hat_all <- j %>%
    select(c(S, D, pi_hat)) %>%
    spread(key = D, value = pi_hat)
  if (inverse) {
    n_repeat <- max(counts$D)
    ret_df <- matrix(replicate(n_repeat, pi_hat_all$"0"), nrow = nrow(pi_hat_all))
  } else {
    pi.hat.df <- select(data.frame(pi_hat_all), -c(1, 2))
    ret_df <- as.matrix(pi.hat.df)
  }
  return(as.matrix(ret_df[S, ]))
}
#-------------------------------------------------------------------
tau.hat.creg <- function(Y, S, D, G.id, Ng, X=NULL, model=NULL, Ng.cov = FALSE)
#-------------------------------------------------------------------
{
  tau.hat.vec <- rep(NA, max(D))
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
    if (Ng.cov == TRUE) {
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

        mu.hat.d <- lin.adj.creg(d, data = cl.lvl.data, model, Ng.cov = TRUE)
        mu.hat.0 <- lin.adj.creg(0, data = cl.lvl.data, model, Ng.cov = TRUE)

        Xi.g <- data$I * (((data$A * (Y.bar.full * data$Ng - mu.hat.d)) / data$pi) -
          (((1 - data$A) * (Y.bar.full * data$Ng - mu.hat.0)) / data$pi.0)) +
          (mu.hat.d - mu.hat.0)

        mu.hat.list[[d]] <- as.matrix(cbind(mu.hat.0, mu.hat.d), ncol = 2)

        tau.hat <- mean(Xi.g) / mean(Ng)

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
      working.df <- data.frame(Y, S, D, G.id, Ng)
      Y.bar.full <- aggregate(Y ~ G.id, working.df, mean)$Y
      cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng")]) # created data on a cluster level for estimating pi.hat(s)
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

        tau.hat <- mean(Xi.g) / mean(Ng)

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
  }
  return(rtrn.list)
}

#-------------------------------------------------------------------
# %##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##
# %##%##%##%##%##%#    II. Variance  estimator     #%##%##%##%##%##%
# %##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# Variance Estimator
#-------------------------------------------------------------------
as.var.creg <- function(model = NULL, fit, HC1) {
  var.vec <- rep(NA, length(fit$tau.hat))
  n.vec <- rep(NA, length(fit$tau.hat))

  if (!is.null(model)) {
    for (d in 1:length(fit$tau.hat))
    {
      Y.bar.g <- fit$Y.bar.g
      Ng <- fit$Ng
      mu.hat.0 <- fit$mu.hat[[d]][, 1]
      mu.hat.d <- fit$mu.hat[[d]][, 2]
      tau.est <- fit$tau.hat
      pi.hat <- fit$pi.hat[[d]]
      pi.hat.0 <- fit$pi.hat.0
      data <- fit$data.list[[d]]
      n <- length(Y.bar.g)

      Xi.tilde.1 <- (mu.hat.d - mu.hat.0) +
        (Ng * Y.bar.g - mu.hat.d) / pi.hat

      Xi.tilde.0 <- (mu.hat.d - mu.hat.0) -
        (Ng * Y.bar.g - mu.hat.0) / pi.hat.0

      data <- data.frame(data, Xi.tilde.1, Xi.tilde.0, Y.Ng = Y.bar.g * Ng)

      count.Xi.1 <- data %>%
        group_by(S, A) %>%
        summarise(Xi.mean.1 = mean(Xi.tilde.1)) %>%
        filter(A != -999999)
      count.Xi.0 <- data %>%
        group_by(S, A) %>%
        summarise(Xi.mean.0 = mean(Xi.tilde.0)) %>%
        filter(A != -999999)
      count.Y <- data %>%
        group_by(S, A) %>%
        summarise(Y.bar = mean(Y.Ng)) %>%
        filter(A != -999999)
      count.Ng <- data %>%
        group_by(S) %>%
        summarise(Ng.bar = mean(Ng))

      j <- left_join(count.Xi.1, count.Xi.0, by = join_by(S == S, A == A)) %>%
        left_join(count.Y, by = join_by(S == S, A == A)) %>%
        left_join(count.Ng, by = join_by(S == S))

      Xi.tilde.1.all <- j %>%
        select(c(S, A, Xi.mean.1)) %>%
        spread(key = A, value = Xi.mean.1)
      Xi.tilde.0.all <- j %>%
        select(c(S, A, Xi.mean.0)) %>%
        spread(key = A, value = Xi.mean.0)
      Y.Ng.all <- j %>%
        select(c(S, A, Y.bar)) %>%
        spread(key = A, value = Y.bar)
      Ng.bar.all <- j %>%
        select(c(S, A, Ng.bar)) %>%
        spread(key = A, value = Ng.bar)

      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.Ng.mean <- as.matrix(select(data.frame(Y.Ng.all), -1))
      Ng.bar.mean <- as.matrix(select(data.frame(Ng.bar.all), -1))

      Xi.1.mean <- Xi.tilde.1.mean[data$S, 2]
      Xi.0.mean <- Xi.tilde.0.mean[data$S, 1]
      Y.g.bar.cl.1 <- Y.Ng.mean[data$S, 2]
      Y.g.bar.cl.0 <- Y.Ng.mean[data$S, 1]
      N.g.bar.cl <- Ng.bar.mean[data$S, 1]

      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.2 <- Y.g.bar.cl.1 - Y.g.bar.cl.0 - tau.est[d] * N.g.bar.cl

      sigma.hat.sq <- mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2) + Xi.hat.2^2) / (mean(Ng))^2

      if (HC1 == TRUE) {
        var.vec[d] <- ((mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2))) * (n / (n - (max(S) + max(D) * max(S)))) +
          mean(Xi.hat.2^2)) / (mean(Ng))^2
      } else {
        var.vec[d] <- sigma.hat.sq
      }
      n.vec[d] <- n
    }
  } else {
    for (d in 1:length(fit$tau.hat))
    {
      Y.bar.g <- fit$Y.bar.g
      Ng <- fit$Ng
      tau.est <- fit$tau.hat
      pi.hat <- fit$pi.hat[[d]]
      pi.hat.0 <- fit$pi.hat.0
      data <- fit$data.list[[d]]
      n <- length(Y.bar.g)

      mu.hat.0 <- 0
      mu.hat.d <- 0

      Xi.tilde.1 <- (mu.hat.d - mu.hat.0) +
        (Ng * Y.bar.g - mu.hat.d) / pi.hat

      Xi.tilde.0 <- (mu.hat.d - mu.hat.0) -
        (Ng * Y.bar.g - mu.hat.0) / pi.hat.0

      data <- data.frame(data, Xi.tilde.1, Xi.tilde.0, Y.Ng = Y.bar.g * Ng)

      count.Xi.1 <- data %>%
        group_by(S, A) %>%
        summarise(Xi.mean.1 = mean(Xi.tilde.1)) %>%
        filter(A != -999999)
      count.Xi.0 <- data %>%
        group_by(S, A) %>%
        summarise(Xi.mean.0 = mean(Xi.tilde.0)) %>%
        filter(A != -999999)
      count.Y <- data %>%
        group_by(S, A) %>%
        summarise(Y.bar = mean(Y.Ng)) %>%
        filter(A != -999999)
      count.Ng <- data %>%
        group_by(S) %>%
        summarise(Ng.bar = mean(Ng))

      j <- left_join(count.Xi.1, count.Xi.0, by = join_by(S == S, A == A)) %>%
        left_join(count.Y, by = join_by(S == S, A == A)) %>%
        left_join(count.Ng, by = join_by(S == S))

      Xi.tilde.1.all <- j %>%
        select(c(S, A, Xi.mean.1)) %>%
        spread(key = A, value = Xi.mean.1)
      Xi.tilde.0.all <- j %>%
        select(c(S, A, Xi.mean.0)) %>%
        spread(key = A, value = Xi.mean.0)
      Y.Ng.all <- j %>%
        select(c(S, A, Y.bar)) %>%
        spread(key = A, value = Y.bar)
      Ng.bar.all <- j %>%
        select(c(S, A, Ng.bar)) %>%
        spread(key = A, value = Ng.bar)

      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.Ng.mean <- as.matrix(select(data.frame(Y.Ng.all), -1))
      Ng.bar.mean <- as.matrix(select(data.frame(Ng.bar.all), -1))

      Xi.1.mean <- Xi.tilde.1.mean[data$S, 2]
      Xi.0.mean <- Xi.tilde.0.mean[data$S, 1]
      Y.g.bar.cl.1 <- Y.Ng.mean[data$S, 2]
      Y.g.bar.cl.0 <- Y.Ng.mean[data$S, 1]
      N.g.bar.cl <- Ng.bar.mean[data$S, 1]

      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.2 <- Y.g.bar.cl.1 - Y.g.bar.cl.0 - tau.est[d] * N.g.bar.cl

      sigma.hat.sq <- mean(data$I * (data$A * (Xi.hat.1)^2 + (1 - data$A) * (Xi.hat.0)^2) + Xi.hat.2^2) / (mean(Ng))^2

      if (HC1 == TRUE) {
        var.vec[d] <- ((mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2))) * (n / (n - (max(S) + max(D) * max(S)))) +
          mean(Xi.hat.2^2)) / mean(Ng)^2
      } else {
        var.vec[d] <- sigma.hat.sq
      }
      n.vec[d] <- n
    }
  }
  se.vec <- sqrt(var.vec / n.vec)
  return(se.vec)
}

#-------------------------------------------------------------------
# %# (10) The core function. It provides estimates of ATE, their s.e.,
# %#     calculates t-stats and corresponding p-values
#-------------------------------------------------------------------
res.creg <- function(Y, S, D, G.id, Ng, X, Ng.cov = FALSE, HC1)
#-------------------------------------------------------------------
{
  n <- length(Y)

  if(is.null(S)){
    S <- rep(1, n)
  }

  if (!is.null(X)) {
    model <- lm.iter.creg(Y, S, D, G.id, Ng, X)
    fit <- tau.hat.creg(Y, S, D, G.id, Ng, X, model)
    tau.est <- fit$tau.hat
    se.rob <- as.var.creg(model, fit, HC1)
    t.stat <- tau.est / se.rob
    p.value <- 2 * pmin(pnorm(t.stat), 1 - pnorm(t.stat))
    CI.left <- tau.est - qnorm(0.975) * se.rob
    CI.right <- tau.est + qnorm(0.975) * se.rob

    res.list <- list(
      "tau.hat"  = tau.est,
      "se.rob"   = se.rob,
      "t.stat"   = t.stat,
      "p.value"  = p.value,
      "as.CI"    = c(CI.left, CI.right),
      "CI.left"  = CI.left,
      "CI.right" = CI.right,
      "data"     = data.frame(Y, S, D, G.id, Ng, X),
      "lin.adj"  = data.frame(Ng, X)
    )
  } else {
    if (Ng.cov == T) {
      model <- lm.iter.creg(Y, S, D, G.id, Ng, X = NULL, Ng.cov = T)
      fit <- tau.hat.creg(Y, S, D, G.id, Ng, X = NULL, model, Ng.cov = T)
      tau.est <- fit$tau.hat
      se.rob <- as.var.creg(model, fit, HC1)
      t.stat <- tau.est / se.rob
      p.value <- 2 * pmin(pnorm(t.stat), 1 - pnorm(t.stat))
      CI.left <- tau.est - qnorm(0.975) * se.rob
      CI.right <- tau.est + qnorm(0.975) * se.rob
      res.list <- list(
        "tau.hat"  = tau.est,
        "se.rob"   = se.rob,
        "t.stat"   = t.stat,
        "p.value"  = p.value,
        "as.CI"    = c(CI.left, CI.right),
        "CI.left"  = CI.left,
        "CI.right" = CI.right,
        "data"     = data.frame(Y, S, D, G.id, Ng),
        "lin.adj"  = data.frame(Ng)
      )
      
    } else {
      fit <- tau.hat.creg(Y, S, D, G.id, Ng, X = NULL, model = NULL)
      tau.est <- fit$tau.hat
      se.rob <- as.var.creg(model = NULL, fit, HC1)
      t.stat <- tau.est / se.rob
      p.value <- 2 * pmin(pnorm(t.stat), 1 - pnorm(t.stat))
      CI.left <- tau.est - qnorm(0.975) * se.rob
      CI.right <- tau.est + qnorm(0.975) * se.rob
      res.list <- list(
        "tau.hat"  = tau.est,
        "se.rob"   = se.rob,
        "t.stat"   = t.stat,
        "p.value"  = p.value,
        "as.CI"    = c(CI.left, CI.right),
        "CI.left"  = CI.left,
        "CI.right" = CI.right,
        "data"     = data.frame(Y, S, D, G.id, Ng),
        "lin.adj"  = NULL
      )
    }
  }
  return(res.list)
}

#-------------------------------------------------------------------
# %# (11) Summary method for sreg(). Provide the output table.
#-------------------------------------------------------------------
summary.creg <- function(model)
#-------------------------------------------------------------------
{
  n <- length(model$data$Y)
  G <- length(unique(model$data$G.id))
  tau.hat <- as.vector(model$tau.hat)
  se.rob <- as.vector(model$se.rob)
  t.stat <- as.vector(model$t.stat)
  p.value <- as.vector(model$p.value)
  CI.left <- as.vector(model$CI.left)
  CI.right <- as.vector(model$CI.right)

  if (!is.null(model$lin.adj)) {
    cat("Saturated Model Estimation Results under CAR with clusters and linear adjustments\n")
  } else {
    cat("Saturated Model Estimation Results under CAR with clusters\n")
  }
  cat(paste("Observations:", n, "\n"))
  cat(paste("Clusters:", G, "\n"))
  cat(paste(
    "Number of treatments:",
    max(model$data$D), "\n"
  ))
  cat(paste(
    "Number of strata:",
    max(model$data$S), "\n"
  ))
  cat(paste0(
    "Covariates used in linear adjustments: ",
    paste(names(model$lin.adj), collapse = ", "), "\n"
  ))

  cat("---\n")

  cat("Coefficients:\n")

  m <- length(tau.hat)

  stars <- rep("", m)
  stars[p.value <= 0.001] <- "***"
  stars[(p.value > 0.001) & (p.value < 0.01)] <- "**"
  stars[(p.value > 0.01) & (p.value <= 0.05)] <- "*"
  stars[(p.value > 0.05) & (p.value <= 0.1)] <- "."

  df <- data.frame(
    "Tau"           = tau.hat,
    "As.se"         = se.rob,
    "T-stat"        = t.stat,
    "P-value"       = p.value,
    "CI.left(95%)"  = CI.left,
    "CI.right(95%)" = CI.right,
    "Significance"  = stars, 
    check.names     = FALSE
  )
  is.df.num.col <- sapply(df, is.numeric)
  df[, is.df.num.col] <- round(df[, is.df.num.col], 5)
  print(df)
  cat("---\n")
  cat(paste(
    "Signif. codes:  0 ‘***’ 0.001 ‘**’",
    "0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n"
  ))
}

#-------------------------------------------------------------------------------
# %##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%#%##%##%##%#
# %##%##%##%##%##%## III. DGP functions for simulations #%##%##%##%##%##%##%##%#%
# %##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%#%##%##%##%#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------
# %# (1) Cluster sizes generation
#-------------------------------------------------------------------
gen.cluster.sizes <- function(G, max.support)
#-------------------------------------------------------------------
{
  sample <- 10 * (rbbinom(G, max.support, alpha = 1, beta = 1) + 1)
  return(sample)
}

#-------------------------------------------------------------------
# %# (2) Potential outcomes generation
#-------------------------------------------------------------------
#------------------------------------------------------------------
dgp.po.creg <- function(Ng, G, tau.vec, sigma1 = sqrt(2),
                        gamma.vec = c(0.4, 0.2, 1), n.treat)
#------------------------------------------------------------------
{
  for (a in seq_along(tau.vec))
  {
    assign(paste("mu.", a, sep = ""), tau.vec[a]) # create mu.a, where a \in \mathbb{A}
  }

  beta.rv <- rbeta(G, 2, 2)
  Z.g.2 <- (beta.rv - 0.5) * sqrt(20)
  x_1 <- (rnorm(G, mean = 5, sd = 2) - 5) / 2
  x_2 <- (rnorm(G, mean = 2, sd = 1) - 2) / 1
  X <- data.frame(x_1, x_2)

  cluster.indicator <- rep(c(1:G), Ng)
  cl.id <- cluster.indicator
  total.sample <- length(cluster.indicator)

  for (a in 1:n.treat)
  {
    assign(paste("epsilon.ig.", a, sep = ""), rnorm(total.sample, 0, sigma1))
  }

  epsilon.ig.0 <- rnorm(total.sample, 0, 1)

  m.0 <- gamma.vec[1] * Z.g.2 + gamma.vec[2] * x_1 + gamma.vec[3] * x_2

  for (a in 1:n.treat) # create m() functions m.a for every treatment a \in \mathbb{A}
  {
    assign(paste("m.", a, sep = ""), m.0)
  }

  Yig.0 <- rep(m.0, Ng) + epsilon.ig.0

  for (a in 1:n.treat)
  {
    formula <- paste0("mu.", a, "+rep(m.", a, ",Ng)", "+epsilon.ig.", a)
    result <- eval(parse(text = formula))
    assign(paste("Yig.", a, sep = ""), result)
  }

  ret.names <- c(
    paste("Yig.", 0:n.treat, sep = ""),
    "Z.g.2", "X", "G", "Ng", "cl.id", paste("m.", 0:n.treat, sep = ""),
    paste("mu.", 1:n.treat, sep = "")
  )

  ret.list <- mget(ret.names)
  return(ret.list)
}

#-------------------------------------------------------------------
# %# (3) Random Treatment Assignment
# %source function for dgp.obs()
#-------------------------------------------------------------------
gen.treat.creg <- function(pi.matr.w, ns, k)
#-------------------------------------------------------------------
{
  rows <- nrow(pi.matr.w)
  code.elements <- character(rows + 1)

  for (i in 1:rows)
  {
    code.elements[i] <- paste0(
      "rep(", i,
      ", floor(pi.matr.w[", i, ",", k, "]*ns))"
    )
  }

  code.elements[rows + 1] <- paste0(
    "rep(0, (ns - ",
    paste0("floor(pi.matr.w[", 1:rows,
      ",", k, "]*ns)",
      collapse = " - "
    ), "))"
  )

  code <- paste(code.elements, collapse = ", ")

  result <- eval(parse(text = paste("sample(c(", code, "))")))

  return(result)
}

#-------------------------------------------------------------------
# %# (4) Generate the formula for Y.obs (Rubin model)
# %source function for dgp.obs()
#-------------------------------------------------------------------
gen.rubin.formula.creg <- function(n.treat)
#-------------------------------------------------------------------
{
  # Create a sequence of A values from 0 to max.A
  A.values <- 0:n.treat

  # Initialize an empty formula string
  formula <- "Y.obs = "

  # Generate the formula dynamically with indicators
  for (a in A.values)
  {
    if (a == 0) {
      formula <- paste(formula, paste0("Y.", a, " * (A == 0)"))
    } else {
      formula <- paste(formula, paste0("Y.", a, " * (A == ", a, ")"))
    }

    if (a < n.treat) {
      formula <- paste(formula, " + ")
    }
  }
  return(formula)
}

#-------------------------------------------------------------------
# %# (5) Generate observed outcomes,
# %#     by taking as input the potential outcomes,
# %#     matrix of strata assignments, pi.vec, and
# %#     number of treatments
#-------------------------------------------------------------------
dgp.obs.creg <- function(baseline, I.S, pi.vec, n.treat)
#-------------------------------------------------------------------
{
  if (n.treat != length(pi.vec)) {
    return(print("The number of treatments doesn't
                 match the length of vector pi.vec"))
  }
  num.strata <- ncol(I.S)
  n <- baseline$G
  A <- cbind(rep(0, n)) # treatment Assignment
  l.seq <- num.strata / 2

  pi.matr <- matrix(1, ncol = num.strata, nrow = n.treat)
  pi.matr.w <- pi.matr * pi.vec

  for (k in 1:num.strata)
  {
    index <- which(I.S[, k] == 1)
    ns <- length(index)

    # pick a random permutation of elements in \mathbb{A} and 0
    A[index] <- gen.treat.creg(pi.matr.w, ns, k)
  }
  strata.set <- data.frame(I.S)
  strata.set$S <- max.col(strata.set)
  cluster.indicator <- baseline$cl.id
  G.seq <- seq(c(1:baseline$G))
  data.short <- data.frame(
    "cl.id" = G.seq, A, S = strata.set$S, Ng = baseline$Ng,
    baseline$X
  )
  data.long <- data.frame("cl.id" = cluster.indicator)
  merged.data <- merge(data.long, data.short, by = "cl.id")
  length(merged.data$A)
  A <- merged.data$A
  S <- merged.data$S
  X <- merged.data[5:ncol(merged.data)]
  Ng <- merged.data$Ng
  # now we need to generate observed outcomes via the Rubin model
  for (a in 0:n.treat)
  {
    assign(paste("Y.", a, sep = ""), baseline[[paste("Yig.", a, sep = "")]])
  }
  formula <- gen.rubin.formula.creg(n.treat)
  Y.obs <- eval(parse(text = formula))

  ret.list <- list(
    "Y"           = Y.obs,
    "D"           = A,
    "S"           = S,
    "Z.2"         = baseline$Z.g.2,
    "X"           = X,
    "Ng"          = Ng,
    "G.id"        = cluster.indicator,
    "cl.lvl.data" = data.short
  )
  return(ret.list)
}

#---------------------------------------------------------------------
# %# (6) Function taken from Ivan's website
# %# to generate the strata from the observed covariates
# %# NB! Works only if we form strata from one W.
#---------------------------------------------------------------------
form.strata.creg <- function(baseline, num.strata)
#---------------------------------------------------------------------
{
  #-------------------------------------------------------------------
  # X-plain: Generates strata indicators from covariates (W).
  #         - W: covariates (must be scalar)
  #         - num.strata: the number of strata to be formed
  #         - model: model in baseline.
  #-------------------------------------------------------------------
  n <- baseline$G
  W <- baseline$Z.g.2
  bounds <- seq(min(W), max(W), length.out = num.strata + 1) # careful with bounds
  I.S <- matrix(0, n, num.strata)
  for (s in 1:num.strata)
  {
    I.S[, s] <- (W > bounds[s]) * (W <= bounds[s + 1])
  }

  return(I.S)
}
