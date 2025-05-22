utils::globalVariables(c("CI.lower", "CI.upper", "SE", "label", "tau"))
#' Print \code{sreg} Objects
#'
#' Print the summary table of estimation results for \code{sreg} objects.
#' @param x An object of class \code{sreg}.
#' @param ... Additional arguments passed to other methods.
#' @examples
#' data <- sreg.rgen(n = 200, tau.vec = c(0.1), n.strata = 4, cluster = TRUE)
#' Y <- data$Y
#' S <- data$S
#' D <- data$D
#' X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
#' result <- sreg(Y, S, D, G.id = NULL, Ng = NULL, X)
#' print(result)
#' @method print sreg
#' @export
#' @return No return value, called for side effects.
print.sreg <- function(x, ...) {
  n <- length(x$data$Y)
  tau.hat <- as.vector(x$tau.hat)
  se.rob <- as.vector(x$se.rob)
  t.stat <- as.vector(x$t.stat)
  p.value <- as.vector(x$p.value)
  CI.left <- as.vector(x$CI.left)
  CI.right <- as.vector(x$CI.right)

  if (!is.null(x$lin.adj)) {
    cat(col_blue("Saturated Model Estimation Results under CAR with linear adjustments\n"))
  } else {
    cat(col_blue("Saturated Model Estimation Results under CAR\n"))
  }

  cat(paste(col_blue("Observations:"), n, "\n"))
  if (!is.null(x$data$G.id)) {
    G <- length(unique(x$data$G.id))
    cat(paste(col_blue("Clusters:"), G, "\n"))
  }
  cat(paste(col_blue("Number of treatments:"), max(x$data$D), "\n"))
  cat(paste(col_blue("Number of strata:"), max(x$data$S), "\n"))
  cat(paste(col_blue("Setup:"), ifelse(!is.null(x$small.strata) && x$small.strata, "small strata", "big strata"), "\n"))
  if (x$small.strata) {
    k <- length(x$data$Y) / max(x$data$S)
    cat(paste(col_blue("Strata size (k):"), k, "\n"))
  }
  cat(paste(col_blue("Standard errors:"), ifelse(!is.null(x$HC1) && x$HC1, "adjusted (HC1)", "unadjusted"), "\n"))
  cat(paste(col_blue("Treatment assignment:"), ifelse(is.null(x$data$G.id), "individual level", "cluster level"), "\n"))

  cat(paste0(
    col_blue("Covariates used in linear adjustments: "),
    paste(names(x$lin.adj), collapse = ", "), "\n"
  ))
  cat("---\n")

  cat(col_blue("Coefficients:\n"))
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
  print.data.frame(df)
  cat("---\n")
  cat(paste(
    col_cyan("Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1\n")
  ))

  if (!is.null(x$data$G.id) && is.null(x$data$Ng)) {
    warning("Warning: Cluster sizes have not been provided (Ng = NULL). Ng is assumed to be equal to the number of available observations in every cluster g.")
  }

  if (any(sapply(x$ols.iter, function(x) any(is.na(x))))) {
    stop("Error: There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
  }

  if (!is.null(x$lin.adj) && !is.null(x$data$G.id)) {
    if (!check.cluster(data.frame("G.id" = x$data$G.id, x$lin.adj))) {
      warning("Warning: sreg cannot use individual-level covariates for covariate adjustment in cluster-randomized experiments. Any individual-level covariates have been aggregated to their cluster-level averages.")
    }
  }
}

#' Plot Method for `sreg` Objects
#'
#' Visualize estimated ATEs and confidence intervals for objects of class \code{sreg}.
#'
#' @param x An object of class \code{sreg}.
#' @param treatment_labels Optional vector of treatment labels.
#' @param ... Additional arguments (not used).
#'
#' @method plot sreg
#' @export
#' @return A ggplot object (invisible).
plot.sreg <- function(x, treatment_labels = NULL, ...) {

  df <- data.frame(
    treatment = seq_along(x$tau.hat),
    tau = x$tau.hat,
    CI.lower = x$CI.left,
    CI.upper = x$CI.right,
    SE = x$se.rob
  )

  if (!is.null(treatment_labels)) {
    df$label <- factor(treatment_labels, levels = rev(treatment_labels))
  } else {
    df$label <- factor(paste("Treatment", df$treatment), levels = rev(paste("Treatment", df$treatment)))
  }

  p <- ggplot(df, aes(y = label)) +
    geom_rect(aes(
      xmin = CI.lower,
      xmax = CI.upper,
      ymin = as.numeric(label) - 0.01,
      ymax = as.numeric(label) + 0.01,
      fill = tau
    ), alpha = 0.8) +
    geom_point(aes(x = tau), shape = 23, size = 3, fill = "white", stroke = 1.2, color = "black") +
    geom_text(aes(x = tau, label = sprintf("%.2f (%.2f)", tau, SE)),
              vjust = -1.5, color = "black", size = 4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_viridis_c(option = "viridis", name = "Effect Size") +
    labs(
      x = "Treatment Effect",
      y = NULL,
      title = "Estimated ATEs with Confidence Intervals"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(hjust = 1),
      axis.title.x = element_text(face = "bold"),
      legend.position = "right"
    )

  print(p)
  invisible(p)
}
