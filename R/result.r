#-------------------------------------------------------------------
# %#     The core function. It provides estimates of ATE, their s.e.,
# %#     calculates t-stats and corresponding p-values
#-------------------------------------------------------------------
res.sreg <- function(Y, S=NULL, D, X=NULL, HC1)
#-------------------------------------------------------------------
{
  n <- length(Y)
  if (is.null(S)) {
    S <- rep(1, n)
  }
  if (!is.null(X)) {
    model <- lm.iter.sreg(Y, S, D, X)
    tau.est <- tau.hat.sreg(Y, S, D, X, model)
    se.rob <- as.var.sreg(Y, S, D, X, model, tau.est, HC1)
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
      "ols.iter" = model,
      "CI.left"  = CI.left,
      "CI.right" = CI.right,
      "data"     = data.frame(Y, S, D, X),
      "lin.adj"  = data.frame(X)
    )
  } else {
    tau.est <- tau.hat.sreg(Y, S, D, X = NULL, model = NULL)
    se.rob <- as.var.sreg(Y, S, D, X = NULL, model = NULL, tau.est, HC1)
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
      "ols.iter" = NULL,
      "CI.left"  = CI.left,
      "CI.right" = CI.right,
      "data"     = data.frame(Y, S, D),
      "lin.adj"  = NULL
    )
  }
  class(res.list) <- "sreg"
  return(res.list)
}
#-------------------------------------------------------------------
res.creg <- function(Y, S, D, G.id, Ng, X, HC1)
#-------------------------------------------------------------------
{
  n <- length(Y)

  if (is.null(S)) {
    S <- rep(1, n)
  }
  if (!is.null(X)) {
    df <- data.frame(G.id, X)
    if (!check.cluster(df)) {
      X.0 <- X
      df.mod <- as.data.frame(df %>%
        group_by(G.id) %>%
        mutate(across(everything(), ~ if (is.numeric(.)) mean(.x, na.rm = TRUE) else .x)) %>%
        ungroup())
      X <- df.mod[, 2:ncol(df.mod)]
    } else {
      X.0 <- X
    }
    model <- lm.iter.creg(Y, S, D, G.id, Ng, X)
    fit <- tau.hat.creg(Y, S, D, G.id, Ng, X, model)
    tau.est <- fit$tau.hat
    se.rob <- as.var.creg(model, fit, HC1)
    t.stat <- tau.est / se.rob
    p.value <- 2 * pmin(pnorm(t.stat), 1 - pnorm(t.stat))
    CI.left <- tau.est - qnorm(0.975) * se.rob
    CI.right <- tau.est + qnorm(0.975) * se.rob
    if (!is.null(Ng)) {
      res.list <- list(
        "tau.hat"  = tau.est,
        "se.rob"   = se.rob,
        "t.stat"   = t.stat,
        "p.value"  = p.value,
        "as.CI"    = c(CI.left, CI.right),
        "ols.iter" = model$theta.list,
        "CI.left"  = CI.left,
        "CI.right" = CI.right,
        "data"     = data.frame(Y, S, D, G.id, Ng, X.0),
        "lin.adj"  = data.frame(X.0)
      )
    } else {
      res.list <- list(
        "tau.hat"  = tau.est,
        "se.rob"   = se.rob,
        "t.stat"   = t.stat,
        "p.value"  = p.value,
        "as.CI"    = c(CI.left, CI.right),
        "ols.iter" = model$theta.list,
        "CI.left"  = CI.left,
        "CI.right" = CI.right,
        "data"     = data.frame(Y, S, D, G.id, X.0),
        "lin.adj"  = data.frame(X.0)
      )
    }
  } else {
    fit <- tau.hat.creg(Y, S, D, G.id, Ng, X = NULL, model = NULL)
    tau.est <- fit$tau.hat
    se.rob <- as.var.creg(model = NULL, fit, HC1)
    t.stat <- tau.est / se.rob
    p.value <- 2 * pmin(pnorm(t.stat), 1 - pnorm(t.stat))
    CI.left <- tau.est - qnorm(0.975) * se.rob
    CI.right <- tau.est + qnorm(0.975) * se.rob
    if (!is.null(Ng)) {
      res.list <- list(
        "tau.hat"  = tau.est,
        "se.rob"   = se.rob,
        "t.stat"   = t.stat,
        "p.value"  = p.value,
        "as.CI"    = c(CI.left, CI.right),
        "ols.iter" = NULL,
        "CI.left"  = CI.left,
        "CI.right" = CI.right,
        "data"     = data.frame(Y, S, D, G.id, Ng),
        "lin.adj"  = NULL
      )
    } else {
      res.list <- list(
        "tau.hat"  = tau.est,
        "se.rob"   = se.rob,
        "t.stat"   = t.stat,
        "p.value"  = p.value,
        "as.CI"    = c(CI.left, CI.right),
        "ols.iter" = NULL,
        "CI.left"  = CI.left,
        "CI.right" = CI.right,
        "data"     = data.frame(Y, S, D, G.id),
        "lin.adj"  = NULL
      )
    }
  }
  class(res.list) <- "sreg"
  return(res.list)
}
