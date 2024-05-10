#-------------------------------------------------------------------
# %#    Summary method for sreg(). Provides the output table.
#-------------------------------------------------------------------
summary.sreg <- function(model)
#-------------------------------------------------------------------
{
  n <- length(model$data$Y)
  tau.hat <- as.vector(model$tau.hat)
  se.rob <- as.vector(model$se.rob)
  t.stat <- as.vector(model$t.stat)
  p.value <- as.vector(model$p.value)
  CI.left <- as.vector(model$CI.left)
  CI.right <- as.vector(model$CI.right)

  if (!is.null(model$data$x_1)) {
    cat("Saturated Model Estimation Results under CAR with linear adjustments\n")
  } else {
    cat("Saturated Model Estimation Results under CAR\n")
  }

  cat(paste("Observations:", n, "\n"))
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
    "Signif. codes:  0 `***` 0.001 `**`",
    "0.01 `*` 0.05 `.` 0.1 ` ` 1\n"
  ))
  if (any(sapply(model$ols.iter, function(x) any(is.na(x))))) {
    stop("There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
  }
}
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
    "Signif. codes:  0 `***` 0.001 `**`",
    "0.01 `*` 0.05 `.` 0.1 ` ` 1\n"
  ))
  ### Warnings: ###
  if (is.null(model$data$Ng)) {
    warning("Warning: cluster sizes have not been provided (Ng = NULL). Ng is assumed to be equal to the number of available observations in every cluster g.")
  }
  if (any(sapply(model$ols.iter, function(x) any(is.na(x))))) {
    stop("There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
  }
  if (!is.null(model$lin.adj)) {
    if (!check.cluster(data.frame("G.id" = model$data$G.id, model$lin.adj))) {
      warning("Warning: sreg cannot use individual-level covariates for covariate adjustment in cluster-randomized experiments. Any individual-level covariates have been aggregated to their cluster-level averages.")
    }
  }
}
