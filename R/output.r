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
  if (is.null(x$data$G.id)) {
    n <- length(x$data$Y)
    tau.hat <- as.vector(x$tau.hat)
    se.rob <- as.vector(x$se.rob)
    t.stat <- as.vector(x$t.stat)
    p.value <- as.vector(x$p.value)
    CI.left <- as.vector(x$CI.left)
    CI.right <- as.vector(x$CI.right)

    if (!is.null(x$data$x_1)) {
      cat(col_blue("Saturated Model Estimation Results under CAR with linear adjustments\n"))
    } else {
      cat(col_blue("Saturated Model Estimation Results under CAR\n"))
    }

    cat(paste(col_blue("Observations:"), n, "\n"))
    cat(paste(
      col_blue("Number of treatments:"),
      max(x$data$D), "\n"
    ))
    cat(paste(
      col_blue("Number of strata:"),
      max(x$data$S), "\n"
    ))
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
      col_cyan(
        "Signif. codes:  0 `***` 0.001 `**`",
        "0.01 `*` 0.05 `.` 0.1 ` ` 1\n"
      )
    ))
    if (any(sapply(x$ols.iter, function(x) any(is.na(x))))) {
      stop("Error: There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
    }
  } else {
    n <- length(x$data$Y)
    G <- length(unique(x$data$G.id))
    tau.hat <- as.vector(x$tau.hat)
    se.rob <- as.vector(x$se.rob)
    t.stat <- as.vector(x$t.stat)
    p.value <- as.vector(x$p.value)
    CI.left <- as.vector(x$CI.left)
    CI.right <- as.vector(x$CI.right)

    if (!is.null(x$lin.adj)) {
      cat(col_blue("Saturated Model Estimation Results under CAR with clusters and linear adjustments\n"))
    } else {
      cat(col_blue("Saturated Model Estimation Results under CAR with clusters\n"))
    }
    cat(paste(col_blue("Observations:"), n, "\n"))
    cat(paste(col_blue("Clusters:"), G, "\n"))
    cat(paste(
      col_blue("Number of treatments:"),
      max(x$data$D), "\n"
    ))
    cat(paste(
      col_blue(
        "Number of strata:"
      ),
      max(x$data$S), "\n"
    ))
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
      col_cyan(
        "Signif. codes:  0 `***` 0.001 `**`",
        "0.01 `*` 0.05 `.` 0.1 ` ` 1\n"
      )
    ))

    if (is.null(x$data$Ng)) {
      warning("Warning: Cluster sizes have not been provided (Ng = NULL). Ng is assumed to be equal to the number of available observations in every cluster g.")
    }
    if (any(sapply(x$ols.iter, function(x) any(is.na(x))))) {
      stop("Error: There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
    }
    if (!is.null(x$lin.adj)) {
      if (!check.cluster(data.frame("G.id" = x$data$G.id, x$lin.adj))) {
        warning("Warning: sreg cannot use individual-level covariates for covariate adjustment in cluster-randomized experiments. Any individual-level covariates have been aggregated to their cluster-level averages.")
      }
    }
  }
}
