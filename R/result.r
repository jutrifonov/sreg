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
#-------------------------------------------------------------------
res.sreg.ss <- function(Y, S, D, X, HC1 = TRUE) 
#-------------------------------------------------------------------
{
  N <- length(Y)
  if (!is.null(X)) {
    model <- tau.hat.sreg.ss(Y, D, X, S)
    tau.est <- model$tau.hat
    var.est <- as.var.sreg.ss(Y, D, X, S, fit = model, HC1) / N
    se.rob <- sqrt(var.est)
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
      "beta.hat" = model$beta.hat,
      "CI.left"  = CI.left,
      "CI.right" = CI.right,
      "data"     = data.frame(Y, S, D, X),
      "lin.adj"  = data.frame(X)
    )
  } else {
    model <- tau.hat.sreg.ss(Y, D, X = NULL, S)
    tau.est <- model$tau.hat
    var.est <- as.var.sreg.ss(Y, D, X = NULL, S, fit = NULL, HC1) / N
    se.rob <- sqrt(var.est)
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
      "beta.hat" = NULL,
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
res.creg.ss <- function(Y, S, D, G.id, Ng, X = NULL, HC1 = TRUE) 
#-------------------------------------------------------------------
{
  N <- length(unique(G.id))
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
    X <- data.frame(X)
    model <- tau.hat.creg.ss(Y, D, X, S, G.id, Ng)
    tau.est <- model$tau.hat
    var.est <- as.var.creg.ss(Y, D, X, S, G.id, Ng, model, HC1) / N
    se.rob <- sqrt(var.est)
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
        "beta.hat" = model$beta.hat,
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
        "beta.hat" = model$beta.hat,
        "CI.left"  = CI.left,
        "CI.right" = CI.right,
        "data"     = data.frame(Y, S, D, G.id, X.0),
        "lin.adj"  = data.frame(X.0)
      )
    }
  } else {
    model <- tau.hat.creg.ss(Y, D, X = NULL, S, G.id, Ng)
    tau.est <- model$tau.hat
    var.est <- as.var.creg.ss(Y, D, X = NULL, S, G.id, Ng, fit = NULL, HC1) / N
    se.rob <- sqrt(var.est)
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
        "beta.hat" = NULL,
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
        "beta.hat" = NULL,
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

res.sreg.mixed <- function(Y, S, D, X = NULL, HC1 = TRUE, small.strata = TRUE) {
  # Step 1: Prepare data and classify strata
  data <- data.frame(Y = Y, S = S, D = D)
  if (!is.null(X)) {
    data <- cbind(data, X)
  }

  data_all <- design.classifier(data, S = S, small.strata = small.strata)

  # Step 2: Split data
  data_small <- dplyr::filter(data_all, stratum_type == "small")
  data_big   <- dplyr::filter(data_all, stratum_type == "big")

  # Step 3: Extract covariates
  X_names <- if (!is.null(X)) colnames(X) else character(0)
  X_small <- if (length(X_names) > 0) data_small[, X_names, drop = FALSE] else NULL
  X_big   <- if (length(X_names) > 0) data_big[, X_names, drop = FALSE] else NULL

  # Step 4: Run estimators
  res_small <- res.sreg.ss(
    Y = data_small$Y,
    D = data_small$D,
    S = data_small$S,
    X = X_small,
    HC1 = HC1
  )

  res_big <- res.sreg(
    Y = data_big$Y,
    D = data_big$D,
    S = data_big$S,
    X = X_big,
    HC1 = HC1
  )

  # Step 5: Combine estimates
  N_small <- nrow(data_small)
  N_big   <- nrow(data_big)
  N_total <- N_small + N_big

  tau_hat <- (N_small / N_total) * res_small$tau.hat + (N_big / N_total) * res_big$tau.hat
  se_combined <- sqrt((N_small / N_total)^2 * res_small$se.rob^2 +
                      (N_big   / N_total)^2 * res_big$se.rob^2)

  t.stat <- tau_hat / se_combined
  p.value <- 2 * pmin(pnorm(t.stat), 1 - pnorm(t.stat))
  CI.left <- tau_hat - qnorm(0.975) * se_combined
  CI.right <- tau_hat + qnorm(0.975) * se_combined
  res.list <- list(
    "tau.hat"  = tau_hat,
    "se.rob"   = se_combined,
    "t.stat"   = t.stat,
    "p.value"  = p.value,
    "as.CI"    = c(CI.left, CI.right),
    "beta.hat" = res_small$beta.hat,
    "CI.left"  = CI.left,
    "CI.right" = CI.right,
    "ols.iter" = res_big$ols.iter,
    "lin.adj" = res_small$lin.adj,
    "data" = data_all,
    "res.small" = res_small,
    "res.big" = res_big,
    "mixed.design" = TRUE
  )

  class(res.list) <- "sreg"
  return(res.list)
}

res.creg.mixed <- function(Y, S, D, G.id, Ng, X = NULL, HC1 = TRUE, small.strata = TRUE) {
  data <- data.frame(Y = Y, D = D, S = S, G.id = G.id, Ng = Ng)
  if (!is.null(X)) {
    data <- cbind(data, X)
  }

  data_all <- design.classifier(data, S = S, G.id = G.id, small.strata = small.strata)
  data_small <- dplyr::filter(data_all, stratum_type == "small")
  data_big   <- dplyr::filter(data_all, stratum_type == "big")

  X_names <- if (!is.null(X)) colnames(X) else character(0)
  X_small <- if (length(X_names) > 0) data_small[, X_names, drop = FALSE] else NULL
  X_big   <- if (length(X_names) > 0) data_big[, X_names, drop = FALSE] else NULL

  res_small <- res.creg.ss(Y = data_small$Y, D = data_small$D, S = data_small$S,
                           G.id = data_small$G.id, Ng = data_small$Ng,
                           X = X_small, HC1 = HC1)
  res_big <- res.creg(Y = data_big$Y, D = data_big$D, S = data_big$S,
                      G.id = data_big$G.id, Ng = data_big$Ng,
                      X = NULL, HC1 = FALSE)

  N_small <- nrow(data_small)
  N_big   <- nrow(data_big)
  N_total <- N_small + N_big

  tau_hat <- (N_small / N_total) * res_small$tau.hat + (N_big / N_total) * res_big$tau.hat
  se_combined <- sqrt((N_small / N_total)^2 * res_small$se.rob^2 +
                      (N_big   / N_total)^2 * res_big$se.rob^2)

  t.stat <- tau_hat / se_combined
  p.value <- 2 * pmin(pnorm(t.stat), 1 - pnorm(t.stat))
  CI.left <- tau_hat - qnorm(0.975) * se_combined
  CI.right <- tau_hat + qnorm(0.975) * se_combined

  ret.list <- list(
    "tau.hat"      = tau_hat,
    "se.rob"       = se_combined,
    "t.stat"       = t.stat,
    "p.value"      = p.value,
    "as.CI"        = c(CI.left, CI.right),
    "beta.hat"     = res_small$beta.hat,
    "CI.left"      = CI.left,
    "CI.right"     = CI.right,
    "ols.iter"     = res_big$ols.iter,
    "lin.adj"      = res_small$lin.adj,
    "data"         = data_all,
    "res.small"    = res_small,
    "res.big"      = res_big,
    "mixed.design" = TRUE
  )
  class(ret.list) = "sreg"
  return(ret.list)
}
