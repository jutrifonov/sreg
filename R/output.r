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
  cat(paste(
    col_blue("Setup:"),
    if (!is.null(x$mixed.design) && x$mixed.design) {
      "mixed design (includes both small and large strata)"
    } else if (!is.null(x$small.strata) && x$small.strata) {
      "small strata"
    } else {
      "big strata"
    },
    "\n"
  ))
  if (!is.null(x$mixed.design) && x$mixed.design) {
    data_small <- x$res.small$data
    if (is.null(data_small$G.id)) {
      k <- nrow(data_small) / max(data_small$S)
    } else {
      cluster_per_stratum <- unique(data.frame(S = data_small$S, G.id = data_small$G.id))
      strata_sizes <- table(cluster_per_stratum$S)
      if (length(unique(strata_sizes)) == 1) {
        k <- unique(strata_sizes)
      } else {
        k <- paste0("varying (min=", min(strata_sizes), ", max=", max(strata_sizes), ")")
      }
    }
    cat(paste(col_blue("Strata size (k, small strata only):"), k, "\n"))
  } else if (!is.null(x$small.strata) && x$small.strata) {
    if (is.null(x$data$G.id)) {
      k <- length(x$data$Y) / max(x$data$S)
    } else {
      cluster_per_stratum <- unique(data.frame(S = x$data$S, G.id = x$data$G.id))
      strata_sizes <- table(cluster_per_stratum$S)
      if (length(unique(strata_sizes)) == 1) {
        k <- unique(strata_sizes)
      } else {
        k <- paste0("varying (min=", min(strata_sizes), ", max=", max(strata_sizes), ")")
      }
    }
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
    if (!is.null(x$mixed.design) && x$mixed.design) {
      gid_check <- x$res.small$data$G.id
    } else {
      gid_check <- x$data$G.id
    }

    if (!check.cluster(data.frame("G.id" = gid_check, x$lin.adj))) {
    }
  }
}

#' Plot Method for `sreg` Objects
#'
#' Visualize estimated ATEs and confidence intervals for objects of class \code{sreg}.
#'
#' @param x An object of class \code{sreg}.
#' @param treatment_labels Optional vector of treatment labels to display on the y-axis. If \code{NULL}, default labels like "Treatment 1", "Treatment 2", etc., are used.
#' @param title Optional plot title. Defaults to "Estimated ATEs with Confidence Intervals".
#' @param bar_fill Optional fill color(s) for the confidence interval bars. Can be \code{NULL} (default viridis scale), a single color, or a vector of two colors for a gradient.
#' @param point_shape Optional shape of the point used to mark the estimated ATE. Default is 23 (a diamond).
#' @param point_size Optional size of the point marking the ATE.
#' @param point_fill Optional fill color of the ATE point shape.
#' @param point_stroke Optional stroke (border) thickness of the ATE point shape.
#' @param point_color Optional outline color of the ATE point.
#' @param label_color Optional color of the text label displaying the estimate and standard error.
#' @param label_size Optional size of the text label displaying the estimate and standard error.
#' @param bg_color Optional background color of the plot panel. If \code{NULL}, the default theme background is used.
#' @param grid Optional logical flag. If \code{TRUE} (default), grid lines are shown; if \code{FALSE}, they are removed.
#' @param zero_line Optional logical flag. If \code{TRUE} (default), a vertical dashed line at 0 is added for reference.
#' @param y_axis_title Optional title of the y-axis. If \code{NULL}, no y-axis label is added.
#' @param x_axis_title Optional title of the x-axis. If \code{NULL}, no x-axis label is added.
#' @param ... Additional arguments passed to other methods.
#'
#' @method plot sreg
#' @export
#' @return Invisibly returns the ggplot object. Called for its side effects (i.e., generating a plot).
plot.sreg <- function(x,
                      treatment_labels = NULL,
                      title = "Estimated ATEs with Confidence Intervals",
                      bar_fill = NULL,
                      point_shape = 23,
                      point_size = 3,
                      point_fill = "white",
                      point_stroke = 1.2,
                      point_color = "black",
                      label_color = "black",
                      label_size = 4,
                      bg_color = NULL,
                      grid = TRUE,
                      zero_line = TRUE,
                      y_axis_title = NULL,
                      x_axis_title = NULL,
                      ...) {
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
    geom_point(aes(x = tau),
      shape = point_shape,
      size = point_size,
      fill = point_fill,
      stroke = point_stroke,
      color = point_color
    ) +
    geom_text(aes(x = tau, label = sprintf("%.2f (%.2f)", tau, SE)),
      vjust = -1.5, color = label_color, size = label_size
    )
  if (zero_line) {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "gray40")
  }

  # Conditional fill scale
  if (is.null(bar_fill)) {
    p <- p + scale_fill_viridis_c(option = "viridis", name = "Effect Size")
  } else if (length(bar_fill) == 1) {
    p <- p + scale_fill_gradient(low = bar_fill, high = bar_fill, name = "Effect Size")
  } else if (length(bar_fill) == 2) {
    p <- p + scale_fill_gradient(low = bar_fill[1], high = bar_fill[2], name = "Effect Size")
  } else {
    warning("bar_fill must be NULL, a single color, or a vector of two colors. Ignoring custom fill.")
    p <- p + scale_fill_viridis_c(option = "viridis", name = "Effect Size")
  }

  theme_custom <- theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(hjust = 1),
      axis.title.x = element_text(),
      legend.position = "none"
    )

  if (!grid) {
    theme_custom <- theme_custom +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  }

  if (!is.null(bg_color)) {
    theme_custom <- theme_custom +
      theme(
        panel.background = element_rect(fill = bg_color, color = NA)
      )
  }

  p <- p + theme_custom +
    labs(
      x = x_axis_title,
      y = y_axis_title,
      title = title
    )

  print(p)
  invisible(p)
}