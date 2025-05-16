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
#-------------------------------------------------------------------
tau.hat.sreg.ss <- function(Y, D, X = NULL, S) 
#-------------------------------------------------------------------
{
  if (!is.null(X)) {
    tau.hat <- numeric(max(D))
    beta.hat <- numeric(max(D))
    beta.hat <- matrix(ncol = ncol(X), nrow = max(D))
    data_full <- data.frame(Y, D, X, S)
    for (d in 1:max(D))
    {
      agg_data <- data_full %>%
        group_split(S) %>%
        map_dfr(~ {
          df <- .
          covariate_cols <- setdiff(names(df), c("Y", "D", "S"))
          list(
            S = df$S[1],
            Y_treated = mean(df$Y[df$D == d], na.rm = TRUE),
            Y_control = mean(df$Y[df$D == 0], na.rm = TRUE),
            X_treated = list(colMeans(df[df$D == d, covariate_cols, drop = FALSE], na.rm = TRUE)),
            X_control = list(colMeans(df[df$D == 0, covariate_cols, drop = FALSE], na.rm = TRUE)),
            k = nrow(df),
            l = sum(df$D == d),
            q = sum(df$D == 0)
          )
        })

      # Create Y_diff vector
      Y_diff <- with(agg_data, Y_treated - Y_control)
      # Compute differences correctly: row-wise for each pair of treated/control vectors
      X_diff_mat <- map2(
        agg_data$X_treated, agg_data$X_control,
        ~ .x - .y
      ) %>%
        do.call(rbind, .)
      # print(X_diff_mat)
      # print(Y_diff)



      data_decomp <- as.data.frame(agg_data)
      X_treated_mat <- do.call(rbind, data_decomp$X_treated)
      X_control_mat <- do.call(rbind, data_decomp$X_control)


      lm_model <- lm(Y_diff ~ ., data = as.data.frame(cbind(Y_diff, X_diff_mat)))

      # print(summary(lm_model))
      beta_hat <- unname(lm_model$coefficients[-1])

      covariate_cols <- setdiff(names(data_full), c("Y", "D", "S"))
      X_mat <- as.matrix(data_full[, covariate_cols, drop = FALSE])
      X_bar <- colMeans(X_mat)
      X_dem <- sweep(X_mat, 2, X_bar)

      # adjusted estimator:
      theta_hat_adj <- sum((data_full$Y * (data_full$D == d))) / sum((data_full$D == d)) -
        sum((data_full$Y) * (data_full$D == 0)) / sum(data_full$D == 0) -
        as.numeric(t(colMeans(X_dem[data_full$D == d, , drop = FALSE]) -
          colMeans(X_dem[data_full$D == 0, , drop = FALSE])) %*% beta_hat)
      tau.hat[d] <- theta_hat_adj
      beta.hat[d, ] <- beta_hat
    }
  } else {
    tau.hat <- numeric(max(D))
    # beta.hat <- numeric(max(D))
    # beta.hat <- matrix(ncol = ncol(X), nrow = max(D))
    data_full <- data.frame(Y, D, S)
    for (d in 1:max(D))
    {
      agg_data <- data_full %>%
        group_split(S) %>%
        map_dfr(~ {
          df <- .
          list(
            S = df$S[1],
            Y_treated = mean(df$Y[df$D == d], na.rm = TRUE),
            Y_control = mean(df$Y[df$D == 0], na.rm = TRUE),
            k = nrow(df),
            l = sum(df$D == d),
            q = sum(df$D == 0)
          )
        })

      # Create Y_diff vector
      Y_diff <- with(agg_data, Y_treated - Y_control)

      data_decomp <- as.data.frame(agg_data)

      # adjusted estimator:
      theta_hat <- sum((data_full$Y * (data_full$D == d))) / sum((data_full$D == d)) -
        sum((data_full$Y) * (data_full$D == 0)) / sum(data_full$D == 0)

      tau.hat[d] <- theta_hat
      beta.hat <- NULL
    }
  }

  ret_list <- list(
    tau.hat = tau.hat,
    beta.hat = beta.hat
  )
  return(ret_list)
}
#-------------------------------------------------------------------
tau.hat.creg.ss <- function(Y, D, X = NULL, S, G.id, Ng) 
#-------------------------------------------------------------------
{
  if (!is.null(X)) {
    tau.hat <- numeric(max(D))

    beta.hat <- matrix(ncol = ncol(X), nrow = max(D))

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

    cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", setdiff(names(working.df), c("Y", "S", "D", "G.id", "Ng"))) ])
    cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
    # print(cl.lvl.data)
    data <- cl.lvl.data
    N.bar.G <- mean(data$Ng) # ??? Why is this so weird?

    covariate_cols <- names(X)

    for (d in 1:max(D))
    {
      # Compute averages for treated and control within each stratum
      agg_data <- data %>%
        group_split(S) %>%
        map_dfr(~ {
          df <- .
          covariate_cols <- names(X)
          list(
            S = df$S[1],
            Y_treated = mean(df$Y.bar[df$D == d] * N.bar.G, na.rm = TRUE),
            Y_control = mean(df$Y.bar[df$D == 0] * N.bar.G, na.rm = TRUE),
            X_treated = list(colMeans(df[df$D == d, covariate_cols, drop = FALSE], na.rm = TRUE)),
            X_control = list(colMeans(df[df$D == 0, covariate_cols, drop = FALSE], na.rm = TRUE)),
            k = nrow(df), # Total units in stratum (should be 2)
            l = sum(df$D == d),
            q = sum(df$D == 0)
          )
        })

      # Create Y_diff vector
      Y_diff <- with(agg_data, Y_treated - Y_control)

      # Compute differences correctly: row-wise for each pair of treated/control vectors
      X_diff_mat <- map2(
        agg_data$X_treated, agg_data$X_control,
        ~ .x - .y
      ) %>%
        do.call(rbind, .)

      # print(Y_diff)
      # print(X_diff_mat)
      data_decomp <- as.data.frame(agg_data)
      X_treated_mat <- do.call(rbind, data_decomp$X_treated)
      X_control_mat <- do.call(rbind, data_decomp$X_control)

      # run the linear model for covariate adjustments
      lm_model <- lm(Y_diff ~ ., data = as.data.frame(cbind(Y_diff, X_diff_mat)))
      
      beta_hat <- unname(lm_model$coefficients[-1])

      covariate_cols <- names(X)
      X_mat <- as.matrix(data[, covariate_cols, drop = FALSE])

      X_bar <- colMeans(X_mat)
      X_dem <- sweep(X_mat, 2, X_bar)

      # X_bar <- mean(data$X)
      # X_dem <- data$X - X_bar
      # adjusted estimator:
      G_1 <- sum((data$D == d) * data$Ng)
      G_0 <- sum((data$D == 0) * data$Ng)
      # print(X_dem)
      # print(beta_hat)
      # print(as.numeric(t(colSums(X_dem[data$D == d, , drop = FALSE]) / G_1 -
      #                             colSums(X_dem[data$D == 0, , drop = FALSE]) / G_0)))
      theta_hat_adj <- sum((data$Y.bar * data$Ng * (data$D == d))) / sum((data$D == d) * data$Ng) -
        sum(data$Y.bar * data$Ng * (data$D == 0)) / sum((data$D == 0) * data$Ng) -
        as.numeric(t(colSums(X_dem[data$D == d, , drop = FALSE]) / G_1 -
          colSums(X_dem[data$D == 0, , drop = FALSE]) / G_0) %*% beta_hat)


      # (sum(X_dem * (data$D == d)) / sum((data$D == d) * data$Ng) -
      # sum(X_dem * (data$D == 0)) / sum((data$D == 0) * data$Ng)) * beta_hat
      tau.hat[d] <- theta_hat_adj
      beta.hat[d, ] <- beta_hat
    }
  } else {
    tau.hat <- numeric(max(D))

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
    cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
    # print(cl.lvl.data)
    data <- cl.lvl.data
    N.bar.G <- mean(data$Ng) # ??? Why is this so weird?

    for (d in 1:max(D))
    {
      # Compute averages for treated and control within each stratum
      agg_data <- data %>%
        group_split(S) %>%
        map_dfr(~ {
          df <- .
          list(
            S = df$S[1],
            Y_treated = mean(df$Y.bar[df$D == d] * N.bar.G, na.rm = TRUE),
            Y_control = mean(df$Y.bar[df$D == 0] * N.bar.G, na.rm = TRUE),
            k = nrow(df), # Total units in stratum (should be 2)
            l = sum(df$D == d),
            q = sum(df$D == 0)
          )
        })

      # Create Y_diff vector
      Y_diff <- with(agg_data, Y_treated - Y_control)

      # print(Y_diff)
      # print(X_diff_mat)
      data_decomp <- as.data.frame(agg_data)


      # run the linear model for covariate adjustments
      # print(summary(lm_model))
      # print(summary(lm_model))
      beta_hat <- 0

      # X_bar <- mean(data$X)
      # X_dem <- data$X - X_bar
      # adjusted estimator:
      G_1 <- sum((data$D == d) * data$Ng)
      G_0 <- sum((data$D == 0) * data$Ng)
      # print(X_dem)
      # print(beta_hat)
      # print(as.numeric(t(colSums(X_dem[data$D == d, , drop = FALSE]) / G_1 -
      #                             colSums(X_dem[data$D == 0, , drop = FALSE]) / G_0)))
      theta_hat <- sum((data$Y.bar * data$Ng * (data$D == d))) / sum((data$D == d) * data$Ng) -
        sum(data$Y.bar * data$Ng * (data$D == 0)) / sum((data$D == 0) * data$Ng)

      tau.hat[d] <- theta_hat
      beta.hat <- NULL
    }
  }
  ret_list <- list(
    tau.hat = tau.hat,
    beta.hat = beta.hat
  )

  return(ret_list)
}