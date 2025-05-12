# %#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
# %#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
######################################### THE CORE  ############################################
# %#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
# %#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
dgp_po <- function(n, theta.vec, gamma.vec, n.treat, is.cov = TRUE)
#-----------------------------------------------------------------------
{
  if (n.treat != length(theta.vec)) {
    stop("The number of treatments doesn't match the length of vector theta.vec.")
  }
  for (a in seq_along(theta.vec))
  {
    assign(paste("mu.", a, sep = ""), theta.vec[a])
  }
  eps.0 <- rnorm(n)
  for (a in 1:n.treat)
  {
    assign(paste("eps.", a, sep = ""), rnorm(n))
  }

  W <- sqrt(20) * (rbeta(n, 2, 2) - 1 / 2)
  x_1 <- rnorm(n, mean = 5, sd = 2)
  x_2 <- rnorm(n, mean = 2, sd = 1)

  if (is.cov == TRUE) {
    X <- data.frame(x_1, x_2)
    # X <- data.frame(x_1)
    m.0 <- gamma.vec[1] * W + gamma.vec[2] * x_1 + gamma.vec[3] * x_2
  } else {
    m.0 <- gamma.vec[1] * W
  }

  for (a in 1:n.treat)
  {
    assign(paste("m.", a, sep = ""), m.0)
  }

  Y.0 <- m.0 + eps.0

  for (a in 1:n.treat)
  {
    formula <- paste0("mu.", a, "+m.", a, "+eps.", a)
    result <- eval(parse(text = formula))
    assign(paste("Y.", a, sep = ""), result)
  }

  if (is.cov == TRUE) {
    ret.names <- c(
      paste("Y.", 0:n.treat, sep = ""),
      "W", "X", paste("m.", 0:n.treat, sep = ""),
      paste("mu.", 1:n.treat, sep = "")
    )
  } else {
    ret.names <- c(
      paste("Y.", 0:n.treat, sep = ""),
      "W", paste("m.", 0:n.treat, sep = ""),
      paste("mu.", 1:n.treat, sep = "")
    )
  }
  ret.list <- mget(ret.names)
  return(ret.list)
}

run_experiment_multi <- function(dgp_list,
                                 n.treat,
                                 k,
                                 treat_sizes,
                                 badmatch = FALSE,
                                 poutcome = FALSE) {
  # Basic checks
  stopifnot(length(treat_sizes) == (n.treat + 1))
  if (sum(treat_sizes) != k) {
    stop("Sum of treat_sizes must equal k.")
  }

  # Number of total arms is n.treat + 1, labeled 0,1,...,n.treat
  all_treat_levels <- 0:n.treat

  # -- Convert the list from dgp_po into a workable data.frame --
  #    We assume that your dgp_list has elements named "Y.0","Y.1",...,"Y.n.treat", plus W, plus X (if is.cov=TRUE).
  #    Make a data.frame that has columns for each potential outcome and your stratification variable(s).

  # First, figure out how many observations:
  n <- length(dgp_list[[paste0("Y.", 0)]]) # length of Y.0
  # Build a DF row-by-row:
  #  Potential Outcomes:
  Y_mat <- sapply(all_treat_levels, function(a) dgp_list[[paste0("Y.", a)]])
  colnames(Y_mat) <- paste0("Y.", all_treat_levels)

  # For convenience, also extract W and (optionally) X
  W <- dgp_list$W
  # If you have a data.frame X, you can splice it in; or if it's just one vector x_1, do the same.
  # Here, assume 'X' is either a single numeric or a data.frame.  We’ll just keep it as-is.
  X <- dgp_list$X

  # Combine into one data.frame:
  # We'll store W in a single column called "W".  If X is a data.frame, you may want to cbind.
  df <- data.frame(Y_mat, W = W)
  # If X is a data.frame with multiple columns, do: df <- cbind(df, X).
  # If it’s just a single vector, do:
  if (is.data.frame(X)) {
    df <- cbind(df, X)
  } else {
    # treat X as a single column
    df[["X"]] <- X
  }

  # For matching, we need a single covariate (or a score) to order or "bad‐match" by.
  # If you want to match by W, we’ll define `match_var = df[["W"]]`.
  # If you want to match by X, define `match_var = df[["X"]]`.
  # Adjust as you like.
  match_var <- df[["W"]] # or df[["X"]], etc.

  # Sort or rearrange for block creation
  if (!badmatch) {
    # Good-match approach: just sort ascending by 'match_var'
    sort_index <- order(match_var)
    df_sorted <- df[sort_index, ]
  } else {
    # "Bad matching": pair smallest with largest, etc. generalized to block size k
    # The simplest approach for pairs is to reorder: first is smallest, second is largest, third is 2nd smallest, fourth is 2nd largest, etc.
    # For block size > 2, you can do an analogous approach (though it gets trickier).
    # Here is a simple generalization that reverses half of the data row‐wise:
    o <- order(match_var) # ascending
    df_asc <- df[o, ]
    # Now we can build "bad" blocks by taking half from the front, half from the back, interleaving, etc.
    # For simplicity, we can do something like: 1st block = 1, N, 2, N-1, ...
    # The code below tries a standard “pair up extremes” approach, repeated for blocks of size k.
    # If k=2, it just does the classical pairs.  If k=3 or 4, you might want to adapt logic further.
    # This skeleton tries to keep it straightforward:
    N <- nrow(df_asc)
    # vector of indices from front/back
    half <- floor(N / 2)
    front_inds <- 1:half
    back_inds <- (N):(half + 1)
    # interleave front and back
    new_order <- as.vector(rbind(front_inds, back_inds))
    # if N is not an even multiple of k, watch out for leftover
    if (length(new_order) < N) {
      leftover <- setdiff(seq_len(N), new_order)
      new_order <- c(new_order, leftover)
    }
    df_sorted <- df_asc[new_order, ]
  }

  # Now we assign blocks of size k.  We'll loop block by block.
  N <- nrow(df_sorted)
  n.blocks <- N %/% k # integer number of full blocks
  if (n.blocks * k != N) {
    stop("N is not divisible by block size k in this simple example.")
  }

  # Prepare container for assigned treatments and block labels
  A <- rep(NA, N) # the assigned treatment
  block_id <- rep(NA, N)

  # For each block, randomly permute the treat_sizes among the k units
  # so that exactly treat_sizes[1] get treat=0, treat_sizes[2] get treat=1, etc.
  # all_treat_levels = c(0,1,...,n.treat)
  for (j in seq_len(n.blocks)) {
    # indices for block j
    these_inds <- ((j - 1) * k + 1):(j * k)

    # We create a vector of length k with exactly treat_sizes[a] copies of 'a', for a in all_treat_levels
    block_treat_vector <- unlist(mapply(function(tt, tsize) rep(tt, tsize),
      tt = all_treat_levels,
      tsize = treat_sizes
    ))
    # Randomly permute that vector so we don't always assign in the same pattern
    block_treat_vector <- sample(block_treat_vector)

    # Assign them
    A[these_inds] <- block_treat_vector
    # Record the block
    block_id[these_inds] <- j
  }

  # Now figure out the observed outcome, y_i, by picking from Y.0..Y.n.treat
  # We'll do it row by row.  For row i, the assigned treatment is A[i].  So the outcome is df_sorted[i, Y.(A[i])]
  # We can do that with indexing or by a small loop:
  y_obs <- numeric(N)
  # column names for the potential outcomes:
  pot_outcome_names <- paste0("Y.", all_treat_levels)
  for (i in seq_len(N)) {
    treat_i <- A[i]
    # which column?  "Y.0" => column 1, "Y.1" => column 2, etc.
    col_idx <- which(all_treat_levels == treat_i)
    y_obs[i] <- df_sorted[i, pot_outcome_names[col_idx]]
  }

  # Final data frame
  # Minimal columns: (y, A, W, X, block)
  out_df <- data.frame(
    y = y_obs,
    A = A,
    W = df_sorted[["W"]],
    block = block_id
  )
  # If your df_sorted includes multiple X columns (like x_1, x_2, etc.), you can add them:
  out_df <- cbind(out_df, df_sorted[, c("x_1", "x_2")])
  # if ("x_1" %in% names(df_sorted)) {
  #  out_df[["X"]] <- df_sorted[["x_1"]] # NB CORRECT FOR MULTIPLE COVARIATES!!! Working on it... No flexibility needed for dgp functions!
  # }
  # Optionally attach all potential outcomes for debugging/evaluation
  if (poutcome) {
    out_df <- cbind(
      out_df,
      df_sorted[, pot_outcome_names, drop = FALSE]
    )
  }

  # Return
  rownames(out_df) <- NULL
  return(out_df)
}

# Now generalized for the cov/nocov case with X = NULL
theta_hat_mult <- function(Y, D, X = NULL, S) {
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
          list(
            S = df$S[1],
            Y_treated = mean(df$Y[df$D == d], na.rm = TRUE),
            Y_control = mean(df$Y[df$D == 0], na.rm = TRUE),
            X_treated = list(colMeans(df[df$D == d, grep("^x_", names(df)), drop = FALSE])), # mean value among treated in each stratum
            X_control = list(colMeans(df[df$D == 0, grep("^x_", names(df)), drop = FALSE])), # mean value among control in each stratum
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

      X_mat <- as.matrix(data_full[, grepl("^x_", names(data_full))])
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
pi.hat.sreg <- function(S, D, inverse = FALSE, vector = FALSE)
#-------------------------------------------------------------------
{
  n <- length(D)
  data <- data.frame(S, D)
  counts <- data %>%
    group_by(S, D) %>%
    summarise(n = n(), .groups = "keep")
  scount <- data %>%
    group_by(S) %>%
    summarise(ns = n(), .groups = "keep")

  j <- left_join(counts, scount, by = join_by(S == S))
  j$pi_hat <- j$n / j$ns
  pi_hat_all <- j %>%
    select(c("S", "D", "pi_hat")) %>%
    spread(key = "D", value = "pi_hat")
  if (inverse) {
    n_repeat <- max(counts$D)
    ret_df <- matrix(replicate(n_repeat, pi_hat_all$"0"), nrow = nrow(pi_hat_all))
  } else {
    pi.hat.df <- select(data.frame(pi_hat_all), -c(1, 2))
    ret_df <- as.matrix(pi.hat.df)
  }
  if (vector) {
    return(as.double(c(ret_df[1, ])))
  }
  return(as.matrix(ret_df[S, ]))
}

var_hat_mult <- function(Y, D, X, S, fit) {
  # n = number of blocks
  n <- max(S)

  if (!is.null(X)) {
    # Center X and compute the augmented outcome Y_a
    X_bar <- colMeans(X)
    X_dem <- sweep(as.matrix(X), 2, X_bar)
    pi_hat_vec <- pi.hat.sreg(S, D, vector = TRUE)
    pi_hat_0 <- pi.hat.sreg(S, D, vector = TRUE, inverse = TRUE)[1]
    V <- numeric(max(D))
    for (d in 1:max(D))
    {
      # beta_hat <- theta_hat_mult(Y, D, X, S)$beta.hat[d, ]
      beta_hat <- fit$beta.hat[d, ]
      # print(beta_hat)
      # print(X_dem)
      # Y_a     <- Y - beta_hat * X_dem
      Y_a <- Y - X_dem %*% beta_hat # check carefully here and in the cluster function what is wrong with the transpose sign?

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
      V_d <- (1 / pi_hat) * v_hat_1_1 +
        (1 / pi_hat_0) * v_hat_1_0 +
        v_hat_2_11 + v_hat_2_00 -
        2 * v_hat_2_10
      V[d] <- V_d
    }
  } else{
    
  }

  return(V)
}

theta_hat_mult_cl_old <- function(Y, D, X, S, G.id, Ng) {
  tau.hat <- numeric(max(D))
  beta.hat <- numeric(max(D))

  working.df <- data.frame(Y, S, D, G.id, Ng, X)
  Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
  # print(Y.bar.g)
  cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", names(working.df)[6:ncol(working.df)])])
  cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
  data <- cl.lvl.data
  N.bar.G <- mean(data$Ng) # ??? Why is this so weird?
  # N.bar.G <- Ng
  # print(N.bar.g)

  # print(data)
  # print(head(working.df, 200))
  # print(data)
  for (d in 1:max(D))
  {
    # Compute averages for treated and control within each stratum
    agg_data <- data %>%
      group_by(S) %>%
      summarise(
        Y_treated = mean(Y.bar[D == d] * N.bar.G, na.rm = TRUE),
        Y_control = mean(Y.bar[D == 0] * N.bar.G, na.rm = TRUE),
        X_treated = mean(X[D == d], na.rm = TRUE),
        X_control = mean(X[D == 0], na.rm = TRUE),
        k = n(), # Total units in stratum (should be 2)
        l = sum(D == d), # Number of treated units (should be 1)
        q = sum(D == 0)
      ) %>%
      mutate(
        Y_diff = (1 / l) * Y_treated - (1 / q) * Y_control,
        X_diff = (1 / l) * X_treated - (1 / q) * X_control,
        k = k,
        l = l,
        q = q
      )


    data_decomp <- as.data.frame(agg_data)

    # run the linear model for covariate adjustments
    lm_model <- lm(Y_diff ~ X_diff, data = data_decomp)
    # print(summary(lm_model))
    beta_hat <- coef(lm_model)[2]
    X_bar <- mean(data$X)
    X_dem <- data$X - X_bar
    # adjusted estimator:
    theta_hat_adj <- sum((data$Y.bar * data$Ng * (data$D == d))) / sum((data$D == d) * data$Ng) -
      sum(data$Y.bar * data$Ng * (data$D == 0)) / sum((data$D == 0) * data$Ng) -
      (sum(X_dem * (data$D == d)) / sum((data$D == d) * data$Ng) -
        sum(X_dem * (data$D == 0)) / sum((data$D == 0) * data$Ng)) * beta_hat
    tau.hat[d] <- theta_hat_adj
    beta.hat[d] <- beta_hat
  }
  # return(c(tau.hat, beta.hat)) #NB! WORK ON THE WAY HOW TO EXTRACT BOTH SEPARATELY! -- DONE!
  # return(tau.hat)
  ret_list <- list(
    tau.hat = tau.hat,
    beta.hat = beta.hat
  )
  return(ret_list)
}
theta_hat_mult_cl <- function(Y, D, X, S, G.id, Ng) {
  tau.hat <- numeric(max(D))
  beta.hat <- matrix(ncol = ncol(X), nrow = max(D))

  working.df <- data.frame(Y, S, D, G.id, Ng, X)
  Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)

  cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", names(working.df)[6:ncol(working.df)])])
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
          X_treated = list(colMeans(df[df$D == d, grep("^x_", names(df)), drop = FALSE])), # mean value among treated in each stratum
          X_control = list(colMeans(df[df$D == 0, grep("^x_", names(df)), drop = FALSE])), # mean value among control in each stratum
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
    # print(summary(lm_model))
    # print(summary(lm_model))
    beta_hat <- unname(lm_model$coefficients[-1])

    X_mat <- as.matrix(data[, grepl("^x_", names(data))])

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

  ret_list <- list(
    tau.hat = tau.hat,
    beta.hat = beta.hat
  )
  return(ret_list)
}

# %#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
####################### Finally, construct a variance estimator with clusters ##################
# %#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
var_hat_mult_cl <- function(Y, D, X, S, G.id, Ng, fit) {
  n <- max(S)

  working.df <- data.frame(Y, S, D, G.id, Ng, X)
  Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
  # print(Y.bar.g)
  cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", names(working.df)[6:ncol(working.df)])])
  cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
  data <- cl.lvl.data
  N.bar.G <- mean(data$Ng) # ??? Why is this so weird?
  # n = number of blocks

  # print(n)

  # Center X and compute the augmented outcome Y_a
  # X_bar   <- mean(data$X)
  # X_dem   <- data$X - X_bar

  X_mat <- as.matrix(data[, grepl("^x_", names(data))])
  X_bar <- colMeans(X_mat)
  X_dem <- sweep(X_mat, 2, X_bar)

  pi_hat_vec <- pi.hat.creg(data$S, data$D, vector = TRUE)
  pi_hat_0 <- pi.hat.creg(data$S, data$D, vector = TRUE, inverse = TRUE)[1]


  V <- numeric(max(data$D))

  for (d in 1:max(data$D))
  {
    beta_hat <- fit$beta.hat[d, ]
    # print(X_dem)
    # print(beta_hat)
    Y_a <- (data$Ng / mean(data$Ng)) * data$Y.bar - X_dem %*% beta_hat * (1 / N.bar.G)
    # Y_a <- data$Y.bar - (X_dem %*% beta_hat) / data$Ng
    l <- sum(data$D == d) / n
    q <- sum(data$D == 0) / n
    # print(l)
    # print(q)
    pi_hat <- pi_hat_vec[d]
    # Compute Gamma_hat_1 and Gamma_hat_0
    Gamma_hat_1 <- sum(Y_a[data$D == d]) * (1 / sum(data$D == d))
    Gamma_hat_0 <- sum(Y_a[data$D == 0]) * (1 / sum(data$D == 0))

    # Precompute sums of Y_a for treated & untreated in each block
    sums_treated <- tapply(Y_a * (data$D == d), data$S, sum)
    sums_untreated <- tapply(Y_a * (data$D == 0), data$S, sum)


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
    V_d <- (1 / pi_hat) * v_hat_1_1 +
      (1 / pi_hat_0) * v_hat_1_0 +
      v_hat_2_11 + v_hat_2_00 -
      2 * v_hat_2_10
    V[d] <- V_d
  }

  return(V)
}

var_hat_mult_cl_old <- function(Y, D, X, S, G.id, Ng, tau.hat) {
  working.df <- data.frame(Y, S, D, G.id, Ng, X)
  Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
  # print(Y.bar.g)
  cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", names(working.df)[6:ncol(working.df)])])
  cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
  data <- cl.lvl.data
  N.bar.G <- mean(data$Ng) # ??? Why is this so weird?
  # n = number of blocks
  n <- max(data$S)
  # print(n)

  # Center X and compute the augmented outcome Y_a
  X_bar <- mean(data$X)
  X_dem <- data$X - X_bar
  pi_hat_vec <- pi.hat.creg(data$S, data$D, vector = TRUE)
  pi_hat_0 <- pi.hat.creg(data$S, data$D, vector = TRUE, inverse = TRUE)[1]

  V <- numeric(max(data$D))
  for (d in 1:max(data$D))
  {
    beta_hat <- tau.hat$beta.hat[d]
    Y_a <- (data$Ng / mean(data$Ng)) * data$Y.bar - beta_hat * X_dem * (1 / N.bar.G)
    l <- sum(data$D == d) / n
    q <- sum(data$D == 0) / n
    pi_hat <- pi_hat_vec[d]
    # Compute Gamma_hat_1 and Gamma_hat_0
    Gamma_hat_1 <- sum(Y_a[data$D == d]) * (1 / sum(data$D == d))
    Gamma_hat_0 <- sum(Y_a[data$D == 0]) * (1 / sum(data$D == 0))

    # Precompute sums of Y_a for treated & untreated in each block
    sums_treated <- tapply(Y_a * (data$D == d), data$S, sum)
    sums_untreated <- tapply(Y_a * (data$D == 0), data$S, sum)

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
    V_d <- (1 / pi_hat) * v_hat_1_1 +
      (1 / pi_hat_0) * v_hat_1_0 +
      v_hat_2_11 + v_hat_2_00 -
      2 * v_hat_2_10
    V[d] <- V_d
  }

  return(V)
}

#-------------------------------------------------------------------
pi.hat.creg <- function(S, D, inverse = FALSE, vector = FALSE)
#-------------------------------------------------------------------
{
  n <- length(D)
  data <- data.frame(S, D)
  counts <- data %>%
    group_by(.data$S, .data$D) %>%
    summarise(n = n(), .groups = "keep")
  scount <- data %>%
    group_by(.data$S) %>%
    summarise(ns = n(), .groups = "keep")

  j <- left_join(counts, scount, by = join_by(S == S))
  j$pi_hat <- j$n / j$ns
  pi_hat_all <- j %>%
    select(c("S", "D", "pi_hat")) %>%
    spread(key = "D", value = "pi_hat")
  if (inverse) {
    n_repeat <- max(counts$D)
    ret_df <- matrix(replicate(n_repeat, pi_hat_all$"0"), nrow = nrow(pi_hat_all))
  } else {
    pi.hat.df <- select(data.frame(pi_hat_all), -c(1, 2))
    ret_df <- as.matrix(pi.hat.df)
  }
  if (vector) {
    return(as.double(c(ret_df[1, ])))
  }
  return(as.matrix(ret_df[S, ]))
}

#-------------------------------------------------------------------
# %#   Cluster sizes generation
#-------------------------------------------------------------------
gen.cluster.sizes <- function(G, max.support)
#-------------------------------------------------------------------
{
  sample <- 10 * (rbbinom(G, max.support, alpha = 1, beta = 1) + 1)
  return(sample)
}
Nmax <- 50
max.support <- Nmax / 10 - 1
G <- 100
Ng <- gen.cluster.sizes(G, max.support)
#------------------------------------------------------------------
dgp.po.creg <- function(Ng, G, tau.vec, sigma1 = sqrt(2),
                        gamma.vec = c(0.4, 0.2, 1), n.treat)
#------------------------------------------------------------------
{
  for (a in seq_along(tau.vec))
  {
    assign(paste("mu.", a, sep = ""), tau.vec[a])
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

  for (a in 1:n.treat)
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
dgp.obs.creg <- function(baseline, n.treat, k, treat_sizes)
#-------------------------------------------------------------------
{
  # num.strata <- ncol(I.S)
  n <- baseline$G
  A <- cbind(rep(0, n))
  # Number of total arms is n.treat + 1, labeled 0,1,...,n.treat
  all_treat_levels <- 0:n.treat
  # l.seq <- num.strata / 2

  # pi.matr <- matrix(1, ncol = num.strata, nrow = n.treat)
  # pi.matr.w <- pi.matr * pi.vec

  # for (k in 1:num.strata)
  # {
  #  index <- which(I.S[, k] == 1)
  #  ns <- length(index)

  #   A[index] <- gen.treat.creg(pi.matr.w, ns, k)
  # }
  Y_mat <- sapply(all_treat_levels, function(a) baseline[[paste0("Yig.", a)]])
  colnames(Y_mat) <- paste0("Y.", all_treat_levels)

  cluster.indicator <- baseline$cl.id
  W <- baseline$Z.g.2
  X <- baseline$X

  G.seq <- seq(c(1:baseline$G))

  # Just on a cluster level cl.id, Ng, W, X
  data.short <- data.frame(
    "cl.id" = G.seq, Ng = baseline$Ng,
    W, X
  )

  # print(data.short)

  # Stratification. We use the cluster level data to first assign strata
  match_var <- data.short[["W"]] # or df[["X"]], etc.

  sort_index <- order(match_var)
  data.short_sorted <- data.short[sort_index, ]

  # Now we assign blocks of size k.  We'll loop block by block.
  N <- nrow(data.short_sorted)
  n.blocks <- N %/% k # integer number of full blocks
  if (n.blocks * k != N) {
    stop("N is not divisible by block size k in this simple example.")
  }

  # Prepare container for assigned treatments and block labels
  A <- rep(NA, N) # the assigned treatment
  block_id <- rep(NA, N)

  # For each block, randomly permute the treat_sizes among the k units
  # so that exactly treat_sizes[1] get treat=0, treat_sizes[2] get treat=1, etc.
  # all_treat_levels = c(0,1,...,n.treat)
  for (j in seq_len(n.blocks)) {
    # indices for block j
    these_inds <- ((j - 1) * k + 1):(j * k)

    # We create a vector of length k with exactly treat_sizes[a] copies of 'a', for a in all_treat_levels
    block_treat_vector <- unlist(mapply(function(tt, tsize) rep(tt, tsize),
      tt = all_treat_levels,
      tsize = treat_sizes
    ))
    # Randomly permute that vector so we don't always assign in the same pattern
    block_treat_vector <- sample(block_treat_vector)

    # Assign them
    A[these_inds] <- block_treat_vector
    # Record the block
    block_id[these_inds] <- j
  }

  data.short <- data.frame(
    A,
    S = block_id, data.short_sorted
  )
  # print(data.short)
  # print(block_id)
  # print(A)

  # print(data.short_sorted)
  # short data
  # data.short <- data.frame(
  #  "cl.id" = G.seq, A, S = strata.set$S, Ng = baseline$Ng,
  #  baseline$X
  # )

  data.long <- data.frame("cl.id" = cluster.indicator)
  # print(data.long)
  # individual level data frame without Y! Just cl.id, Ng, W, X
  merged.data <- merge(data.long, data.short, by = "cl.id")
  # print(head(merged.data, 200)) # Works! But we need to append potential outcomes!

  A <- merged.data$A
  S <- merged.data$S
  X <- merged.data[6:ncol(merged.data)]
  Ng <- merged.data$Ng

  for (a in 0:n.treat)
  {
    assign(paste("Y.", a, sep = ""), baseline[[paste("Yig.", a, sep = "")]])
  }
  # print(Y.1)
  formula <- gen.rubin.formula.creg(n.treat)
  Y.obs <- eval(parse(text = formula)) # observed outcomes for every individual
  # print(Y.obs)

  ret.list <- list(
    "Y"           = Y.obs,
    "D"           = A,
    "S"           = S,
    "W"           = baseline$Z.g.2,
    "X"           = X,
    "Ng"          = Ng,
    "G.id"        = cluster.indicator,
    "cl.lvl.data" = data.short
  )
  return(ret.list)
}
#-------------------------------------------------------------------
gen.rubin.formula.creg <- function(n.treat)
#-------------------------------------------------------------------
{
  A.values <- 0:n.treat

  formula <- "Y.obs = "

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



### Aggregate the functions to fit into the master functions res.sreg.ss, res.creg.ss ###
# Now generalized for the cov/nocov case with X = NULL
tau.hat.sreg.ss <- function(Y, D, X = NULL, S) {
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
          list(
            S = df$S[1],
            Y_treated = mean(df$Y[df$D == d], na.rm = TRUE),
            Y_control = mean(df$Y[df$D == 0], na.rm = TRUE),
            X_treated = list(colMeans(df[df$D == d, grep("^x_", names(df)), drop = FALSE])), # mean value among treated in each stratum
            X_control = list(colMeans(df[df$D == 0, grep("^x_", names(df)), drop = FALSE])), # mean value among control in each stratum
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

      X_mat <- as.matrix(data_full[, grepl("^x_", names(data_full))])
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
as.var.sreg.ss <- function(Y, D, X = NULL, S, fit = NULL) {
  # n = number of blocks
  n <- max(S)

  pi_hat_vec <- pi.hat.sreg(S, D, vector = TRUE)
  pi_hat_0 <- pi.hat.sreg(S, D, vector = TRUE, inverse = TRUE)[1]
  V <- numeric(max(D))

  if(!is.null(X))
  {
  # Center X and compute the augmented outcome Y_a
  X_bar <- colMeans(X)
  X_dem <- sweep(as.matrix(X), 2, X_bar)
  }else{
    X_dem <- 0
  }
  for (d in 1:max(D))
  {
    if(!is.null(X))
    {
      beta_hat <- fit$beta.hat[d, ]
      Y_a <- Y - X_dem %*% beta_hat # check carefully here and in the cluster function what is wrong with the transpose sign?
    }else{
      beta_hat <- 0
      Y_a <- Y
    }
    # print(beta_hat)
    # print(X_dem)
    # Y_a     <- Y - beta_hat * X_dem
    #Y_a <- Y - X_dem %*% beta_hat # check carefully here and in the cluster function what is wrong with the transpose sign?

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
    V_d <- (1 / pi_hat) * v_hat_1_1 +
      (1 / pi_hat_0) * v_hat_1_0 +
      v_hat_2_11 + v_hat_2_00 -
      2 * v_hat_2_10
    V[d] <- V_d
  }

  return(V)
}
tau.hat.creg.ss <- function(Y, D, X = NULL, S, G.id, Ng) {
  if(!is.null(X))
  {
  tau.hat <- numeric(max(D))
  beta.hat <- matrix(ncol = ncol(X), nrow = max(D))

  working.df <- data.frame(Y, S, D, G.id, Ng, X)
  Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)

  cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", names(working.df)[6:ncol(working.df)])])
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
          X_treated = list(colMeans(df[df$D == d, grep("^x_", names(df)), drop = FALSE])), # mean value among treated in each stratum
          X_control = list(colMeans(df[df$D == 0, grep("^x_", names(df)), drop = FALSE])), # mean value among control in each stratum
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
    # print(summary(lm_model))
    # print(summary(lm_model))
    beta_hat <- unname(lm_model$coefficients[-1])

    X_mat <- as.matrix(data[, grepl("^x_", names(data))])

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
  }else{
    tau.hat <- numeric(max(D))
    working.df <- data.frame(Y, S, D, G.id, Ng)
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
    beta.hat   <- NULL
  }
  }
    ret_list <- list(
    tau.hat = tau.hat,
    beta.hat = beta.hat
  )
  
  return(ret_list)
}
as.var.creg.ss  <- function(Y, D, X = NULL, S, G.id, Ng, fit = NULL) {
  n <- max(S)
  if(!is.null(X)){
    working.df <- data.frame(Y, S, D, G.id, Ng, X)
    Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
    cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", names(working.df)[6:ncol(working.df)])])
  }else{
    working.df <- data.frame(Y, S, D, G.id, Ng)
    Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
    cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng")])
  }
  cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
  data <- cl.lvl.data
  N.bar.G <- mean(data$Ng) # ??? Why is this so weird?
  # n = number of blocks
  if(!is.null(X)){
    X_mat <- as.matrix(data[, grepl("^x_", names(data))])
    X_bar <- colMeans(X_mat)
    X_dem <- sweep(X_mat, 2, X_bar)
  }
  pi_hat_vec <- pi.hat.creg(data$S, data$D, vector = TRUE)
  pi_hat_0 <- pi.hat.creg(data$S, data$D, vector = TRUE, inverse = TRUE)[1]
  V <- numeric(max(data$D))

  for (d in 1:max(data$D))
  {
  if(!is.null(X)){
    beta_hat <- fit$beta.hat[d, ]
    Y_a <- (data$Ng / mean(data$Ng)) * data$Y.bar - X_dem %*% beta_hat * (1 / N.bar.G)
  }else{
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
    sums_treated   <- tapply(Y_a * (data$D == d), data$S, sum)
    sums_untreated <- tapply(Y_a * (data$D == 0), data$S, sum)


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
    V_d <- (1 / pi_hat) * v_hat_1_1 +
        (1 / pi_hat_0) * v_hat_1_0 +
        v_hat_2_11 + v_hat_2_00 -
        2 * v_hat_2_10
    V[d] <- V_d
    }
    return(V)    
}

var_hat_mult_cl <- function(Y, D, X, S, G.id, Ng, fit) {
  n <- max(S)

  working.df <- data.frame(Y, S, D, G.id, Ng, X)
  Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
  # print(Y.bar.g)
  cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", "Ng", names(working.df)[6:ncol(working.df)])])
  cl.lvl.data <- data.frame("Y.bar" = Y.bar.g$Y, cl.lvl.data)
  data <- cl.lvl.data
  N.bar.G <- mean(data$Ng) # ??? Why is this so weird?
  # n = number of blocks

  # print(n)

  # Center X and compute the augmented outcome Y_a
  # X_bar   <- mean(data$X)
  # X_dem   <- data$X - X_bar

  X_mat <- as.matrix(data[, grepl("^x_", names(data))])
  X_bar <- colMeans(X_mat)
  X_dem <- sweep(X_mat, 2, X_bar)

  pi_hat_vec <- pi.hat.creg(data$S, data$D, vector = TRUE)
  pi_hat_0 <- pi.hat.creg(data$S, data$D, vector = TRUE, inverse = TRUE)[1]


  V <- numeric(max(data$D))

  for (d in 1:max(data$D))
  {
    beta_hat <- fit$beta.hat[d, ]
    # print(X_dem)
    # print(beta_hat)
    Y_a <- (data$Ng / mean(data$Ng)) * data$Y.bar - X_dem %*% beta_hat * (1 / N.bar.G)
    # Y_a <- data$Y.bar - (X_dem %*% beta_hat) / data$Ng
    l <- sum(data$D == d) / n
    q <- sum(data$D == 0) / n
    # print(l)
    # print(q)
    pi_hat <- pi_hat_vec[d]
    # Compute Gamma_hat_1 and Gamma_hat_0
    Gamma_hat_1 <- sum(Y_a[data$D == d]) * (1 / sum(data$D == d))
    Gamma_hat_0 <- sum(Y_a[data$D == 0]) * (1 / sum(data$D == 0))

    # Precompute sums of Y_a for treated & untreated in each block
    sums_treated <- tapply(Y_a * (data$D == d), data$S, sum)
    sums_untreated <- tapply(Y_a * (data$D == 0), data$S, sum)


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
    V_d <- (1 / pi_hat) * v_hat_1_1 +
      (1 / pi_hat_0) * v_hat_1_0 +
      v_hat_2_11 + v_hat_2_00 -
      2 * v_hat_2_10
    V[d] <- V_d
  }

  return(V)
}

