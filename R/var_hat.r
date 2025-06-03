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
      data$pi.0 <- pi.hat.sreg(S, D, inverse = T)[, 1]
      n <- length(Y)
      data$A <- ifelse(D == d, 1, ifelse(D == 0, 0, -999999))
      data$I <- as.numeric(data$A != -999999)

      mu.hat.d <- lin.adj.sreg(d, data$S, data[4:(4 + ncol(X) - 1)], model)
      mu.hat.0 <- lin.adj.sreg(0, data$S, data[4:(4 + ncol(X) - 1)], model)

      Xi.tilde.1 <- (mu.hat.d - mu.hat.0) +
        (data$Y - mu.hat.d) / data$pi

      Xi.tilde.0 <- (mu.hat.d - mu.hat.0) -
        (data$Y - mu.hat.0) / data$pi.0

      data <- data.frame(data, Xi.tilde.1, Xi.tilde.0, Y.tau.D = data$Y - tau[d] * data$A * data$I)

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
        summarise(Y.tau = mean(.data$Y.tau.D)) %>%
        filter(.data$A != -999999)

      j <- left_join(count.Xi.1, count.Xi.0, by = join_by("S" == "S", "A" == "A")) %>% left_join(count.Y, by = join_by("S" == "S", "A" == "A"))

      Xi.tilde.1.all <- j %>%
        select(c("S", "A", "Xi.mean.1")) %>%
        spread(key = "A", value = "Xi.mean.1")
      Xi.tilde.0.all <- j %>%
        select(c("S", "A", "Xi.mean.0")) %>%
        spread(key = "A", value = "Xi.mean.0")
      Y.tau.D.all <- j %>%
        select(c("S", "A", "Y.tau")) %>%
        spread(key = "A", value = "Y.tau")
      # print(data.frame(Xi.tilde.1.all))
      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.tau.D.mean <- as.matrix(select(data.frame(Y.tau.D.all), -1))
      # print(data.frame(Y.tau.D.all))

      # Xi.1.mean <- Xi.tilde.1.mean[S, 2]
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

      sigma.hat.sq <- mean(data$I * (data$A * (Xi.hat.1)^2 + (1 - data$A) * (Xi.hat.0)^2) + Xi.hat.2^2)

      if (HC1 == TRUE) {
        S_reset <- as.integer(factor(S, levels = Y.tau.D.all$S))
        var.vec[d] <- (mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2))) * (n / (n - (max(S_reset) + max(D) * max(S_reset)))) +
          mean(Xi.hat.2^2)
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
      data$pi.0 <- pi.hat.sreg(S, D, inverse = T)[, 1]
      n <- length(Y)
      data$A <- ifelse(D == d, 1, ifelse(D == 0, 0, -999999))
      data$I <- as.numeric(data$A != -999999)

      mu.hat.d <- 0
      mu.hat.0 <- 0

      Xi.tilde.1 <- (mu.hat.d - mu.hat.0) +
        (data$Y - mu.hat.d) / data$pi

      Xi.tilde.0 <- (mu.hat.d - mu.hat.0) -
        (data$Y - mu.hat.0) / data$pi.0

      data <- data.frame(data, Xi.tilde.1, Xi.tilde.0, Y.tau.D = data$Y - tau[d] * data$A * data$I)

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
        summarise(Y.tau = mean(.data$Y.tau.D)) %>%
        filter(.data$A != -999999)

      j <- left_join(count.Xi.1, count.Xi.0, by = join_by("S" == "S", "A" == "A")) %>% left_join(count.Y, by = join_by("S" == "S", "A" == "A"))

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

      # Use S_reset as in the X block
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

      sigma.hat.sq <- mean(data$I * (data$A * (Xi.hat.1)^2 + (1 - data$A) * (Xi.hat.0)^2) + Xi.hat.2^2)
      if (HC1 == TRUE) {
        S_reset <- as.integer(factor(S, levels = Y.tau.D.all$S))
        adj_factor_denom <- n - (max(S_reset) + max(D) * max(S_reset))
        if (adj_factor_denom <= 0 || is.nan(adj_factor_denom)) {
          warning("HC1 adjustment unstable or undefined due to degenerate strata-treatment structure; reverting to unadjusted estimator.")
          var.vec[d] <- sigma.hat.sq
        } else {
          adj_factor <- n / adj_factor_denom
          var.vec[d] <- mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2)) * adj_factor +
            mean(Xi.hat.2^2)
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
      Xi.hat.2 <- Y.g.bar.cl.1 - Y.g.bar.cl.0 - tau.est[d] * N.g.bar.cl

      sigma.hat.sq <- mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2) + Xi.hat.2^2) / (mean(Ng))^2

      if (HC1 == TRUE) {
        var.vec[d] <- ((mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2))) * (n / (n - (max(data$S) + max(data$D) * max(data$S)))) +
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
      Xi.hat.2 <- Y.g.bar.cl.1 - Y.g.bar.cl.0 - tau.est[d] * N.g.bar.cl

      if (HC1 == TRUE) {
        var.vec[d] <- ((mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2))) * (n / (n - (max(data$S) + max(data$D) * max(data$S)))) +
          mean(Xi.hat.2^2)) / mean(Ng)^2
      } else {
        sigma.hat.sq <- mean(data$I * (data$A * (Xi.hat.1)^2 + (1 - data$A) * (Xi.hat.0)^2) + Xi.hat.2^2) / (mean(Ng))^2
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
  # N = number of observations
  N <- length(Y)

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
    # print(beta_hat)
    # print(X_dem)
    # Y_a     <- Y - beta_hat * X_dem
    # Y_a <- Y - X_dem %*% beta_hat # check carefully here and in the cluster function what is wrong with the transpose sign?

    l <- sum(D == d) / n
    q <- sum(D == 0) / n

    pi_hat <- pi_hat_vec[d]
    # Compute Gamma_hat_1 and Gamma_hat_0
    Gamma_hat_1 <- sum(Y_a[D == d]) * (1 / sum(D == d))
    Gamma_hat_0 <- sum(Y_a[D == 0]) * (1 / sum(D == 0))

    # Precompute sums of Y_a for treated & untreated in each block
    sums_treated <- tapply(Y_a * (D == d), S, sum)
    sums_untreated <- tapply(Y_a * (D == 0), S, sum)

    # print(head(sums_treated))
    # print(head(sums_untreated))

    #----------------------------------------
    # Compute rho_hat_00 and rho_hat_11
    # We consider pairs of adjacent blocks: (1,2), (3,4), ...
    #----------------------------------------
    # Indices of pairs
    idx1 <- seq(1, n, 2)
    idx2 <- seq(2, n, 2)

    # zeta_0 = sum of products of untreated across pairs of blocks
    zeta_0 <- sum(sums_untreated[idx1] * sums_untreated[idx2]) / (q^2)
    # print(head(zeta_0))

    # zeta_1 = sum of products of treated across pairs of blocks
    zeta_1 <- sum(sums_treated[idx1] * sums_treated[idx2]) / (l^2)
    # print(head(zeta_1))
    # Multiply each by (2/n) to get rho_00 and rho_11
    rho_hat_00 <- zeta_0 * (2 / n)
    rho_hat_11 <- zeta_1 * (2 / n)
    # print(rho_hat_11)

    #----------------------------------------
    # Compute rho_hat_10
    # sum_rho_10 = sum over j of ( (sum of treated)*(sum of untreated) / (l*(k-l)) )
    # Then divide by n
    #----------------------------------------
    sum_rho_10 <- sum((sums_treated * sums_untreated) / (l * q))
    rho_hat_10 <- sum_rho_10 / n
    # print(sum_rho_10)
    # print(rho_hat_10)

    #----------------------------------------
    # Compute sigma_hat_1 and sigma_hat_0
    #----------------------------------------
    sigma_hat_1 <- sum((Y_a - Gamma_hat_1)^2 * (D == d)) * (1 / (n * l))
    sigma_hat_0 <- sum((Y_a - Gamma_hat_0)^2 * (D == 0)) * (1 / (n * q))
    # print(s)
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
as.var.creg.ss <- function(Y, D, X = NULL, S, G.id, Ng, fit = NULL, HC1 = TRUE)
#-------------------------------------------------------------------
{
  n <- max(S)
  if (!is.null(X)) {
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
    cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", setdiff(names(working.df), c("Y", "S", "D", "G.id", "Ng")))])
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
    Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
    cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng")])
  }
  cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
  data <- cl.lvl.data
  N.bar.G <- mean(data$Ng) # ??? Why is this so weird?
  # n = number of blocks
  if (!is.null(X)) {
    # X_mat <- as.matrix(data[, grepl("^x_", names(data))])
    covariate_cols <- names(X)
    X_mat <- as.matrix(data[, covariate_cols, drop = FALSE])
    X_bar <- colMeans(X_mat)
    X_dem <- sweep(X_mat, 2, X_bar)
  }
  pi_hat_vec <- pi.hat.creg(data$S, data$D, vector = TRUE)
  pi_hat_0 <- pi.hat.creg(data$S, data$D, vector = TRUE, inverse = TRUE)[1]
  V <- numeric(max(data$D))

  for (d in 1:max(data$D))
  {
    if (!is.null(X)) {
      beta_hat <- fit$beta.hat[d, ]
      Y_a <- (data$Ng / mean(data$Ng)) * data$Y.bar - X_dem %*% beta_hat * (1 / N.bar.G)
    } else {
      beta_hat <- 0
      Y_a <- (data$Ng / mean(data$Ng)) * data$Y.bar
    }
    l <- sum(data$D == d) / n
    q <- sum(data$D == 0) / n
    # print(l)
    # print(q)
    pi_hat <- pi_hat_vec[d]
    # Compute Gamma_hat_1 and Gamma_hat_0
    Gamma_hat_1 <- sum(Y_a[data$D == d]) * (1 / sum(data$D == d))
    Gamma_hat_0 <- sum(Y_a[data$D == 0]) * (1 / sum(data$D == 0))

    # Precompute sums of Y_a for treated & untreated in each block
    sums_treated <- as.numeric(tapply(Y_a * (data$D == d), data$S, sum))
    sums_untreated <- as.numeric(tapply(Y_a * (data$D == 0), data$S, sum))



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
    sigma_hat_1 <- sum((Y_a - Gamma_hat_1)^2 * (data$D == d)) * (1 / (n * l))
    sigma_hat_0 <- sum((Y_a - Gamma_hat_0)^2 * (data$D == 0)) * (1 / (n * q))

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
