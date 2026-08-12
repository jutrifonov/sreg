#-------------------------------------------------------------------
# %#     Function that implements \hat{\sigma^2} --
# %#     i.e. the variance estimator
#-------------------------------------------------------------------
as.var.sreg <- function(Y, S, D, X = NULL, model = NULL, tau, HC1)
#-------------------------------------------------------------------
{
  var.vec <- numeric(max(D))
  n.vec <- numeric(max(D))

  if (!is.null(X)) {
    for (d in 1:max(D))
    {
      data <- data.frame(Y, S, D, X)
      data$pi <- pi.hat.sreg(S, D)[, d]
      data$pi.0 <- pi.hat.sreg(S, D, inverse = TRUE)[, 1]
      n <- length(Y)

      data$A <- ifelse(D == d, 1, ifelse(D == 0, 0, -999999))
      data$I <- as.numeric(data$A != -999999)

      # Indicator for observations assigned to treatment arms
      # other than d and control
      data$I.other <- as.numeric(data$A == -999999)

      mu.hat.d <- lin.adj.sreg(d, data$S, data[4:(4 + ncol(X) - 1)], model)
      mu.hat.0 <- lin.adj.sreg(0, data$S, data[4:(4 + ncol(X) - 1)], model)

      Xi.tilde.1 <- (mu.hat.d - mu.hat.0) +
        (data$Y - mu.hat.d) / data$pi

      Xi.tilde.0 <- (mu.hat.d - mu.hat.0) -
        (data$Y - mu.hat.0) / data$pi.0

      # Contribution for observations in treatment arms other than d and control
      Xi.tilde.other <- mu.hat.d - mu.hat.0

      data <- data.frame(
        data,
        Xi.tilde.1,
        Xi.tilde.0,
        Xi.tilde.other,
        Y.tau.D = data$Y - tau[d] * data$A * data$I
      )

      count.Xi.1 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.1 = mean(.data$Xi.tilde.1), .groups = "drop") %>%
        filter(.data$A != -999999)

      count.Xi.0 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.0 = mean(.data$Xi.tilde.0), .groups = "drop") %>%
        filter(.data$A != -999999)

      count.Y <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Y.tau = mean(.data$Y.tau.D), .groups = "drop") %>%
        filter(.data$A != -999999)

      j <- left_join(
        count.Xi.1,
        count.Xi.0,
        by = join_by("S" == "S", "A" == "A")
      ) %>%
        left_join(count.Y, by = join_by("S" == "S", "A" == "A"))

      Xi.tilde.1.all <- j %>%
        select(c("S", "A", "Xi.mean.1")) %>%
        spread(key = "A", value = "Xi.mean.1")

      Xi.tilde.0.all <- j %>%
        select(c("S", "A", "Xi.mean.0")) %>%
        spread(key = "A", value = "Xi.mean.0")

      Y.tau.D.all <- j %>%
        select(c("S", "A", "Y.tau")) %>%
        spread(key = "A", value = "Y.tau")

      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.tau.D.mean <- as.matrix(select(data.frame(Y.tau.D.all), -1))

      S_reset <- as.integer(factor(S, levels = Xi.tilde.1.all$S))
      Xi.1.mean <- Xi.tilde.1.mean[S_reset, 2]

      S_reset <- as.integer(factor(S, levels = Xi.tilde.0.all$S))
      Xi.0.mean <- Xi.tilde.0.mean[S_reset, 1]

      S_reset <- as.integer(factor(S, levels = Y.tau.D.all$S))
      Y.tau.D.1.mean <- Y.tau.D.mean[S_reset, 2]
      Y.tau.D.0.mean <- Y.tau.D.mean[S_reset, 1]

      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean
      Xi.hat.2 <- Y.tau.D.1.mean - Y.tau.D.0.mean

      # Center the other-treatment add-back term within each stratum
      # using all observations in the stratum
      Xi.hat.other <- data$Xi.tilde.other -
        stats::ave(data$Xi.tilde.other, data$S, FUN = mean)

      # This term contributes only for observations assigned to other treatment arms
      Xi.hat.other[data$I.other == 0] <- 0

      # Within-stratum variance includes treatment, control, and other-arm terms
      within.var <- mean(
        as.numeric(data$A == 1) * Xi.hat.1^2 +
          as.numeric(data$A == 0) * Xi.hat.0^2 +
          data$I.other * Xi.hat.other^2
      )

      sigma.hat.sq <- within.var + mean(Xi.hat.2^2)

      if (HC1 == TRUE) {
        S_reset <- as.integer(factor(S, levels = Y.tau.D.all$S))
        adj_factor_denom <- n - (max(S_reset) + max(D) * max(S_reset))

        if (adj_factor_denom <= 0 || is.nan(adj_factor_denom)) {
          warning("HC1 adjustment unstable or undefined due to degenerate strata-treatment structure; reverting to unadjusted estimator.")
          var.vec[d] <- sigma.hat.sq
        } else {
          adj_factor <- n / adj_factor_denom
          var.vec[d] <- within.var * adj_factor + mean(Xi.hat.2^2)
        }
      } else {
        var.vec[d] <- sigma.hat.sq
      }

      n.vec[d] <- n
    }
  } else {
    for (d in 1:max(D))
    {
      data <- data.frame(Y, S, D)
      data$pi <- pi.hat.sreg(S, D)[, d]
      data$pi.0 <- pi.hat.sreg(S, D, inverse = TRUE)[, 1]
      n <- length(Y)

      data$A <- ifelse(D == d, 1, ifelse(D == 0, 0, -999999))
      data$I <- as.numeric(data$A != -999999)

      data$I.other <- as.numeric(data$A == -999999)

      mu.hat.d <- 0
      mu.hat.0 <- 0

      Xi.tilde.1 <- (mu.hat.d - mu.hat.0) +
        (data$Y - mu.hat.d) / data$pi

      Xi.tilde.0 <- (mu.hat.d - mu.hat.0) -
        (data$Y - mu.hat.0) / data$pi.0

      Xi.tilde.other <- mu.hat.d - mu.hat.0

      data <- data.frame(
        data,
        Xi.tilde.1,
        Xi.tilde.0,
        Xi.tilde.other,
        Y.tau.D = data$Y - tau[d] * data$A * data$I
      )

      count.Xi.1 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.1 = mean(.data$Xi.tilde.1), .groups = "drop") %>%
        filter(.data$A != -999999)

      count.Xi.0 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.0 = mean(.data$Xi.tilde.0), .groups = "drop") %>%
        filter(.data$A != -999999)

      count.Y <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Y.tau = mean(.data$Y.tau.D), .groups = "drop") %>%
        filter(.data$A != -999999)

      j <- left_join(
        count.Xi.1,
        count.Xi.0,
        by = join_by("S" == "S", "A" == "A")
      ) %>%
        left_join(count.Y, by = join_by("S" == "S", "A" == "A"))

      Xi.tilde.1.all <- j %>%
        select(c("S", "A", "Xi.mean.1")) %>%
        spread(key = "A", value = "Xi.mean.1")

      Xi.tilde.0.all <- j %>%
        select(c("S", "A", "Xi.mean.0")) %>%
        spread(key = "A", value = "Xi.mean.0")

      Y.tau.D.all <- j %>%
        select(c("S", "A", "Y.tau")) %>%
        spread(key = "A", value = "Y.tau")

      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.tau.D.mean <- as.matrix(select(data.frame(Y.tau.D.all), -1))

      S_reset <- as.integer(factor(S, levels = Xi.tilde.1.all$S))
      Xi.1.mean <- Xi.tilde.1.mean[S_reset, 2]

      S_reset <- as.integer(factor(S, levels = Xi.tilde.0.all$S))
      Xi.0.mean <- Xi.tilde.0.mean[S_reset, 1]

      S_reset <- as.integer(factor(S, levels = Y.tau.D.all$S))
      Y.tau.D.1.mean <- Y.tau.D.mean[S_reset, 2]
      Y.tau.D.0.mean <- Y.tau.D.mean[S_reset, 1]

      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean
      Xi.hat.2 <- Y.tau.D.1.mean - Y.tau.D.0.mean

      Xi.hat.other <- data$Xi.tilde.other -
        stats::ave(data$Xi.tilde.other, data$S, FUN = mean)

      Xi.hat.other[data$I.other == 0] <- 0

      within.var <- mean(
        as.numeric(data$A == 1) * Xi.hat.1^2 +
          as.numeric(data$A == 0) * Xi.hat.0^2 +
          data$I.other * Xi.hat.other^2
      )

      sigma.hat.sq <- within.var + mean(Xi.hat.2^2)

      if (HC1 == TRUE) {
        S_reset <- as.integer(factor(S, levels = Y.tau.D.all$S))
        adj_factor_denom <- n - (max(S_reset) + max(D) * max(S_reset))

        if (adj_factor_denom <= 0 || is.nan(adj_factor_denom)) {
          warning("HC1 adjustment unstable or undefined due to degenerate strata-treatment structure; reverting to unadjusted estimator.")
          var.vec[d] <- sigma.hat.sq
        } else {
          adj_factor <- n / adj_factor_denom
          var.vec[d] <- within.var * adj_factor + mean(Xi.hat.2^2)
        }
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
as.var.creg <- function(model = NULL, fit, HC1)
#-------------------------------------------------------------------
{
  var.vec <- numeric(length(fit$tau.hat))
  n.vec <- numeric(length(fit$tau.hat))

  if (!is.null(model)) {
    for (d in seq_along(fit$tau.hat))
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

      Xi.tilde.other <- mu.hat.d - mu.hat.0
      data$I.other <- as.numeric(data$A == -999999)

      data <- data.frame(data, Xi.tilde.1, Xi.tilde.0, Xi.tilde.other,
        Y.Ng = Y.bar.g * Ng)

      count.Xi.1 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.1 = mean(.data$Xi.tilde.1)) %>%
        filter(.data$A != -999999)
      count.Xi.0 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.0 = mean(.data$Xi.tilde.0)) %>%
        filter(.data$A != -999999)
      count.Y <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Y.bar = mean(.data$Y.Ng)) %>%
        filter(.data$A != -999999)
      count.Ng <- data %>%
        group_by(.data$S) %>%
        summarise(Ng.bar = mean(.data$Ng))

      j <- left_join(count.Xi.1, count.Xi.0, by = join_by("S" == "S", "A" == "A")) %>%
        left_join(count.Y, by = join_by("S" == "S", "A" == "A")) %>%
        left_join(count.Ng, by = join_by("S" == "S"))

      Xi.tilde.1.all <- j %>%
        select(c("S", "A", "Xi.mean.1")) %>%
        spread(key = "A", value = "Xi.mean.1")
      Xi.tilde.0.all <- j %>%
        select(c("S", "A", "Xi.mean.0")) %>%
        spread(key = "A", value = "Xi.mean.0")
      Y.Ng.all <- j %>%
        select(c("S", "A", "Y.bar")) %>%
        spread(key = "A", value = "Y.bar")
      Ng.bar.all <- j %>%
        select(c("S", "A", "Ng.bar")) %>%
        spread(key = "A", value = "Ng.bar")

      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.Ng.mean <- as.matrix(select(data.frame(Y.Ng.all), -1))
      Ng.bar.mean <- as.matrix(select(data.frame(Ng.bar.all), -1))

      S_reset <- as.integer(factor(data$S, levels = Xi.tilde.1.all$S))
      Xi.1.mean <- Xi.tilde.1.mean[S_reset, 2]
      Xi.0.mean <- Xi.tilde.0.mean[S_reset, 1]
      Y.g.bar.cl.1 <- Y.Ng.mean[S_reset, 2]
      Y.g.bar.cl.0 <- Y.Ng.mean[S_reset, 1]
      N.g.bar.cl <- Ng.bar.mean[S_reset, 1]

      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.other <- Xi.tilde.other -
        stats::ave(Xi.tilde.other, data$S, FUN = mean) -
        tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.2 <- Y.g.bar.cl.1 - Y.g.bar.cl.0 - tau.est[d] * N.g.bar.cl

      within.var <- mean(
        as.numeric(data$A == 1) * Xi.hat.1^2 +
          as.numeric(data$A == 0) * Xi.hat.0^2 +
          data$I.other * Xi.hat.other^2
      )
      sigma.hat.sq <- (within.var + mean(Xi.hat.2^2)) / (mean(Ng))^2

      if (HC1 == TRUE) {
        adj_denom <- n - (max(data$S) + max(data$D) * max(data$S))
        if (adj_denom <= 0 || is.nan(adj_denom)) {
          warning("HC1 adjustment unstable or undefined due to degenerate strata-treatment structure; reverting to unadjusted estimator.")
          var.vec[d] <- sigma.hat.sq
        } else {
          var.vec[d] <- (within.var * (n / adj_denom) +
            mean(Xi.hat.2^2)) / (mean(Ng))^2
        }
      } else {
        var.vec[d] <- sigma.hat.sq
      }
      n.vec[d] <- n
    }
  } else {
    for (d in seq_along(fit$tau.hat))
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

      Xi.tilde.other <- rep(0, n)
      data$I.other <- as.numeric(data$A == -999999)

      data <- data.frame(data, Xi.tilde.1, Xi.tilde.0, Xi.tilde.other,
        Y.Ng = Y.bar.g * Ng)

      count.Xi.1 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.1 = mean(.data$Xi.tilde.1)) %>%
        filter(.data$A != -999999)
      count.Xi.0 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.0 = mean(.data$Xi.tilde.0)) %>%
        filter(.data$A != -999999)
      count.Y <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Y.bar = mean(.data$Y.Ng)) %>%
        filter(.data$A != -999999)
      count.Ng <- data %>%
        group_by(.data$S) %>%
        summarise(Ng.bar = mean(.data$Ng))

      j <- left_join(count.Xi.1, count.Xi.0, by = join_by("S" == "S", "A" == "A")) %>%
        left_join(count.Y, by = join_by("S" == "S", "A" == "A")) %>%
        left_join(count.Ng, by = join_by("S" == "S"))

      Xi.tilde.1.all <- j %>%
        select(c("S", "A", "Xi.mean.1")) %>%
        spread(key = "A", value = "Xi.mean.1")
      Xi.tilde.0.all <- j %>%
        select(c("S", "A", "Xi.mean.0")) %>%
        spread(key = "A", value = "Xi.mean.0")
      Y.Ng.all <- j %>%
        select(c("S", "A", "Y.bar")) %>%
        spread(key = "A", value = "Y.bar")
      Ng.bar.all <- j %>%
        select(c("S", "A", "Ng.bar")) %>%
        spread(key = "A", value = "Ng.bar")

      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.Ng.mean <- as.matrix(select(data.frame(Y.Ng.all), -1))
      Ng.bar.mean <- as.matrix(select(data.frame(Ng.bar.all), -1))

      S_reset <- as.integer(factor(data$S, levels = Xi.tilde.1.all$S))
      Xi.1.mean <- Xi.tilde.1.mean[S_reset, 2]
      Xi.0.mean <- Xi.tilde.0.mean[S_reset, 1]
      Y.g.bar.cl.1 <- Y.Ng.mean[S_reset, 2]
      Y.g.bar.cl.0 <- Y.Ng.mean[S_reset, 1]
      N.g.bar.cl <- Ng.bar.mean[S_reset, 1]

      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.other <- Xi.tilde.other -
        stats::ave(Xi.tilde.other, data$S, FUN = mean) -
        tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.2 <- Y.g.bar.cl.1 - Y.g.bar.cl.0 - tau.est[d] * N.g.bar.cl

      within.var <- mean(
        as.numeric(data$A == 1) * Xi.hat.1^2 +
          as.numeric(data$A == 0) * Xi.hat.0^2 +
          data$I.other * Xi.hat.other^2
      )
      sigma.hat.sq <- (within.var + mean(Xi.hat.2^2)) / (mean(Ng))^2

      if (HC1 == TRUE) {
        adj_denom <- n - (max(data$S) + max(data$D) * max(data$S))
        if (adj_denom <= 0 || is.nan(adj_denom)) {
          warning("HC1 adjustment unstable or undefined due to degenerate strata-treatment structure; reverting to unadjusted estimator.")
          var.vec[d] <- sigma.hat.sq
        } else {
          adj_factor <- n / adj_denom
          var.vec[d] <- (within.var * adj_factor +
            mean(Xi.hat.2^2)) / (mean(Ng))^2
        }
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
as.var.sreg.ss <- function(Y, D, X = NULL, S, fit = NULL, HC1 = TRUE)
#-------------------------------------------------------------------
{
  # n = number of blocks
  n <- max(S)
  pi_hat_vec <- pi.hat.sreg(S, D, vector = TRUE)
  pi_hat_0 <- pi.hat.sreg(S, D, vector = TRUE, inverse = TRUE)[1]
  V <- numeric(max(D))

  if (!is.null(X)) {
    # Center X and compute the augmented outcome Y_a
    X_bar <- colMeans(X)
    X_dem <- sweep(as.matrix(X), 2, X_bar)
  } else {
    X_dem <- 0
  }
  for (d in 1:max(D))
  {
    if (!is.null(X)) {
      beta_hat <- fit$beta.hat[d, ]
      Y_a <- Y - X_dem %*% beta_hat # check carefully here and in the cluster function what is wrong with the transpose sign?
    } else {
      beta_hat <- 0
      Y_a <- Y
    }

    l <- sum(D == d) / n
    q <- sum(D == 0) / n

    pi_hat <- pi_hat_vec[d]
    # Compute Gamma_hat_1 and Gamma_hat_0
    Gamma_hat_1 <- sum(Y_a[D == d]) * (1 / sum(D == d))
    Gamma_hat_0 <- sum(Y_a[D == 0]) * (1 / sum(D == 0))

    # Precompute sums of Y_a for treated & untreated in each block
    sums_treated <- tapply(Y_a * (D == d), S, sum)
    sums_untreated <- tapply(Y_a * (D == 0), S, sum)


    #----------------------------------------
    # Compute rho_hat_00 and rho_hat_11
    # We consider pairs of adjacent blocks: (1,2), (3,4), ...
    #----------------------------------------
    # Indices of pairs
    idx1 <- seq(1, n, 2)
    idx2 <- seq(2, n, 2)

    # zeta_0 = sum of products of untreated across pairs of blocks
    zeta_0 <- sum(sums_untreated[idx1] * sums_untreated[idx2]) / (q^2)

    # zeta_1 = sum of products of treated across pairs of blocks
    zeta_1 <- sum(sums_treated[idx1] * sums_treated[idx2]) / (l^2)
    # Multiply each by (2/n) to get rho_00 and rho_11
    rho_hat_00 <- zeta_0 * (2 / n)
    rho_hat_11 <- zeta_1 * (2 / n)

    #----------------------------------------
    # Compute rho_hat_10
    # sum_rho_10 = sum over j of ( (sum of treated)*(sum of untreated) / (l*(k-l)) )
    # Then divide by n
    #----------------------------------------
    sum_rho_10 <- sum((sums_treated * sums_untreated) / (l * q))
    rho_hat_10 <- sum_rho_10 / n

    #----------------------------------------
    # Compute sigma_hat_1 and sigma_hat_0
    #----------------------------------------
    sigma_hat_1 <- sum((Y_a - Gamma_hat_1)^2 * (D == d)) * (1 / (n * l))
    sigma_hat_0 <- sum((Y_a - Gamma_hat_0)^2 * (D == 0)) * (1 / (n * q))
    #----------------------------------------
    # Compute the final variance components
    #----------------------------------------
    # v_hat_1_1 and v_hat_1_0
    v_hat_1_1 <- sigma_hat_1 - (rho_hat_11 - Gamma_hat_1^2)
    v_hat_1_0 <- sigma_hat_0 - (rho_hat_00 - Gamma_hat_0^2)

    # v_hat_2_11, v_hat_2_00, v_hat_2_10
    v_hat_2_11 <- rho_hat_11 - Gamma_hat_1 * Gamma_hat_1

    v_hat_2_00 <- rho_hat_00 - Gamma_hat_0 * Gamma_hat_0

    v_hat_2_10 <- rho_hat_10 - Gamma_hat_1 * Gamma_hat_0

    # Final V
    if (!is.null(X)) {
      if (HC1 == TRUE) {
        beta_hat <- fit$beta.hat[d, ]
        K <- length(beta_hat) + 1
        V_d <- (1 / pi_hat) * (n / (n - K)) * v_hat_1_1 +
          (1 / pi_hat_0) * (n / (n - K)) * v_hat_1_0 +
          v_hat_2_11 + v_hat_2_00 -
          2 * v_hat_2_10
      } else {
        V_d <- (1 / pi_hat) * v_hat_1_1 +
          (1 / pi_hat_0) * v_hat_1_0 +
          v_hat_2_11 + v_hat_2_00 -
          2 * v_hat_2_10
      }
    } else {
      V_d <- (1 / pi_hat) * v_hat_1_1 +
        (1 / pi_hat_0) * v_hat_1_0 +
        v_hat_2_11 + v_hat_2_00 -
        2 * v_hat_2_10
    }
    V[d] <- V_d
  }

  return(V)
}
#-------------------------------------------------------------------
# Small-strata cluster variance with the common-denominator correction.
as.var.creg.ss <- function(Y, D, X = NULL, S, G.id, Ng, fit = NULL, HC1 = TRUE) {
  data <- .creg_ss_cluster_data(Y, S, D, G.id, Ng, X)
  strata <- sort(unique(data$S))
  n <- length(strata)
  if (n %% 2L != 0L) stop("The paired-strata variance estimator requires an even number of strata.")
  nbar <- mean(data$Ng)
  arms <- sort(unique(data$D))
  treatments <- setdiff(arms, 0)
  pi <- vapply(arms, function(r) mean(data$D == r), numeric(1))
  names(pi) <- as.character(arms)
  idx1 <- seq(1, n, 2)
  idx2 <- seq(2, n, 2)
  V <- numeric(max(treatments))

  if (is.null(fit)) fit <- tau.hat.creg.ss(Y, D, X, S, G.id, Ng)

  for (d in treatments) {
    tau <- fit$tau.hat[d]
    if (!is.null(X)) {
      x_names <- paste0(".creg_x_", seq_len(ncol(X)))
      x_mat <- as.matrix(data[, x_names, drop = FALSE])
      x_dem <- sweep(x_mat, 2, colMeans(x_mat))
      beta <- as.numeric(fit$beta.hat[d, ])
      residual_total <- data$T - as.numeric(x_dem %*% beta)
      K <- length(beta) + 1L
    } else {
      residual_total <- data$T
      K <- 1L
    }

    # Signed arm-specific terms in the linearization of Qhat / Nbar.
    W <- -tau * pi[as.character(data$D)] * data$Ng / nbar
    W[data$D == d] <- W[data$D == d] + residual_total[data$D == d] / nbar
    W[data$D == 0] <- W[data$D == 0] - residual_total[data$D == 0] / nbar

    gamma <- sigma <- numeric(length(arms))
    names(gamma) <- names(sigma) <- as.character(arms)
    arm_sums <- matrix(0, nrow = n, ncol = length(arms),
      dimnames = list(as.character(strata), as.character(arms)))
    k_arm <- numeric(length(arms))
    names(k_arm) <- as.character(arms)

    for (r in arms) {
      key <- as.character(r)
      Wr <- W[data$D == r]
      gamma[key] <- mean(Wr)
      sigma[key] <- mean((Wr - gamma[key])^2)
      sums <- tapply(W * (data$D == r), factor(data$S, levels = strata), sum)
      arm_sums[, key] <- as.numeric(sums)
      k_arm[key] <- sum(data$D == r) / n
    }

    rho <- matrix(NA_real_, length(arms), length(arms),
      dimnames = list(as.character(arms), as.character(arms)))
    for (r in arms) {
      rkey <- as.character(r)
      rho[rkey, rkey] <- (2 / n) *
        sum(arm_sums[idx1, rkey] * arm_sums[idx2, rkey]) / k_arm[rkey]^2
    }
    if (length(arms) > 1L) {
      for (i in seq_along(arms)[-length(arms)]) {
        for (j in (i + 1L):length(arms)) {
          rkey <- as.character(arms[i])
          tkey <- as.character(arms[j])
          rho[rkey, tkey] <- rho[tkey, rkey] <-
            mean(arm_sums[, rkey] * arm_sums[, tkey]) /
            (k_arm[rkey] * k_arm[tkey])
        }
      }
    }

    v2 <- rho - outer(gamma, gamma)
    v1 <- sigma - diag(v2)
    hc <- if (isTRUE(HC1)) n / (n - K) else 1
    V[d] <- sum(hc * v1 / pi) + sum(v2)
  }
  V
}
