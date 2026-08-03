#' Estimate Average Treatment Effects (ATEs) and Corresponding Standard Errors
#'
#' Estimate the ATE(s) and the corresponding standard error(s) for a (collection of) treatment(s) relative to a control.
#' @import extraDistr
#' @import tidyr
#' @import dplyr
#' @import rlang
#' @import cli
#' @import ggplot2
#' @import viridis
#' @importFrom purrr map_dfr map2 map_dbl
#' @importFrom stats aggregate coef lm pnorm qnorm rbeta rnorm na.omit
#' @importFrom utils packageVersion
#'
#' @param Y a numeric \eqn{n \times 1} \code{vector/matrix/data.frame/tibble} of the observed outcomes
#' @param S a numeric \eqn{n \times 1} \code{vector/matrix/data.frame/tibble} of strata indicators indexed by \eqn{\{1, 2, 3, \ldots\}};  if \code{NULL} then the estimation is performed assuming no stratification
#' @param D a numeric \eqn{n \times 1} \code{vector/matrix/data.frame/tibble} of treatments indexed by \eqn{\{0, 1, 2, \ldots\}}, where \eqn{\code{D} = 0} denotes the control
#' @param G.id a numeric \eqn{n \times 1} \code{vector/matrix/data.frame/tibble} of cluster indicators; if \code{NULL} then estimation is performed assuming treatment is assigned at the individual level
#' @param Ng a numeric \eqn{n \times 1} \code{vector/matrix/data.frame/tibble} of cluster sizes; if \code{NULL} then \code{Ng} is assumed to be equal to the number of available observations in every cluster
#' @param X a \code{matrix/data.frame/tibble} with columns representing the covariate values for every observation; if \code{NULL} then the estimator without linear adjustments is applied. (Note: \code{sreg} cannot use individual-level covariates for covariate adjustment in cluster-randomized experiments. Any individual-level covariates will be aggregated to their cluster-level averages)
#' @param HC1 a \code{TRUE/FALSE} logical argument indicating whether the small sample correction should be applied to the variance estimator
#' @param small.strata a \code{TRUE/FALSE} logical argument indicating whether the estimators for small strata (i.e., strata with few units, such as matched pairs or n-tuples) should be used.
#' @param k an optional positive integer specifying the number of units per small stratum, or the number of clusters per small stratum in cluster-randomized designs. When \code{NULL}, mixed designs retain the automatic detection rule for matched pairs and triplets. Supply \code{k} for general k-tuple mixed designs.
#' @return An object of class \code{sreg} that is a list containing the following elements:
#' \itemize{
#' \item \code{tau.hat}: a \eqn{1 \times |\mathcal A|} \code{vector} of ATE estimates, where \eqn{|\mathcal A|} represents the number of treatments
#' \item \code{se.rob}: a \eqn{1 \times |\mathcal A|} \code{vector} of standard errors estimates, where \eqn{|\mathcal A|} represents the number of treatments
#' \item \code{t.stat}: a \eqn{1 \times |\mathcal A|} \code{vector} of \eqn{t}-statistics, where \eqn{|\mathcal A|} represents the number of treatments
#' \item \code{p.value}: a \eqn{1 \times |\mathcal A|} \code{vector} of corresponding \eqn{p}-values, where \eqn{|\mathcal A|} represents the number of treatments
#' \item \code{CI.left}: a \eqn{1 \times |\mathcal A|} \code{vector} of the left bounds of the 95\% as. confidence interval
#' \item \code{CI.right}: a \eqn{1 \times |\mathcal A|} \code{vector} of the right bounds of the 95\% as. confidence interval
#' \item \code{data}: an original data of the form \code{data.frame(Y, S, D, G.id, Ng, X)}
#' \item \code{lin.adj}: a \code{data.frame} representing the covariates that were used in implementing linear adjustments
#' \item \code{small.strata}: a \code{TRUE/FALSE} logical argument indicating whether the estimators for small strata (e.g., matched pairs or n-tuples) were used
#' \item \code{HC1}: a \code{TRUE/FALSE} logical argument indicating whether the small sample correction (HC1) was applied to the variance estimator
#' }
#' @references
#' Bugni, F. A., Canay, I. A., and Shaikh, A. M. (2018). Inference Under Covariate-Adaptive Randomization. \emph{Journal of the American Statistical Association}, 113(524), 1784–1796, \doi{10.1080/01621459.2017.1375934}.
#'
#' Bugni, F., Canay, I., Shaikh, A., and Tabord-Meehan, M. (2024+). Inference for Cluster Randomized Experiments with Non-ignorable Cluster Sizes. \emph{Forthcoming in the Journal of Political Economy: Microeconomics}, \doi{10.48550/arXiv.2204.08356}.
#'
#' Jiang, L., Linton, O. B., Tang, H., and Zhang, Y. (2023+). Improving Estimation Efficiency via Regression-Adjustment in Covariate-Adaptive Randomizations with Imperfect Compliance. \emph{Forthcoming in Review of Economics and Statistics}, \doi{10.48550/arXiv.2204.08356}.
#'
#' Bai, Y., Jiang, L., Romano, J. P., Shaikh, A. M., and Zhang, Y. (2024). Covariate adjustment in experiments with matched pairs. \emph{Journal of Econometrics}, 241(1), \doi{10.1016/j.jeconom.2024.105740}.
#'
#' Bai, Y. (2022). Optimality of Matched-Pair Designs in Randomized Controlled Trials. \emph{American Economic Review}, 112(12), \doi{10.1257/aer.20201856}.
#'
#' Bai, Y., Romano, J. P., and Shaikh, A. M. (2022). Inference in Experiments With Matched Pairs. \emph{Journal of the American Statistical Association}, 117(540), \doi{10.1080/01621459.2021.1883437}.
#'
#' Liu, J. (2024). Inference for Two-stage Experiments under Covariate-Adaptive Randomization. \doi{10.48550/arXiv.2301.09016}.
#'
#' Cytrynbaum, M. (2024). Covariate Adjustment in Stratified Experiments. \emph{Quantitative Economics}, 15(4), 971–998,  \doi{10.3982/QE2475}.
#' @author
#' \strong{Authors}:
#'
#' Juri Trifonov \email{jutrifonov@u.northwestern.edu}
#'
#' Yuehao Bai \email{yuehao.bai@usc.edu}
#'
#' Azeem Shaikh \email{amshaikh@uchicago.edu}
#'
#' Max Tabord-Meehan \email{m.tabordmeehan@utoronto.ca}
#'
#'
#' \strong{Maintainer}:
#'
#' Juri Trifonov \email{jutrifonov@u.northwestern.edu}
#'
#'
#'
#'
#' @export
#'
#' @examples
#' ### Large strata with covariate adjustment
#' set.seed(1)
#' large_data <- sreg.rgen(
#'   n = 120, tau.vec = c(0.2, 0.5), n.strata = 4,
#'   cluster = FALSE
#' )
#' fit_large <- sreg(
#'   Y = large_data$Y, S = large_data$S, D = large_data$D,
#'   X = large_data[c("x_1", "x_2")]
#' )
#' fit_large$tau.hat
#' fit_large$se.rob
#'
#' ### Uniform triplets: k is observed and may be omitted from sreg()
#' set.seed(2)
#' small_data <- sreg.rgen(
#'   n = 60, tau.vec = c(0.2, 0.5), cluster = FALSE,
#'   small.strata = TRUE, k = 3, treat.sizes = c(1, 1, 1)
#' )
#' fit_small <- sreg(
#'   Y = small_data$Y, S = small_data$S, D = small_data$D,
#'   small.strata = TRUE
#' )
#' fit_small$tau.hat
#'
#' ### Mixed design with general k-tuples: specify k in sreg()
#' set.seed(3)
#' mixed_data <- sreg.rgen(
#'   n = 60, tau.vec = 0.5, cluster = FALSE,
#'   mixed.strata = TRUE, n.small = 40, k = 4,
#'   treat.sizes = c(2, 2), n.strata = 2
#' )
#' fit_mixed <- suppressWarnings(sreg(
#'   Y = mixed_data$Y, S = mixed_data$S, D = mixed_data$D,
#'   small.strata = TRUE, k = 4
#' ))
#' fit_mixed$mixed.design
#'
#' ### Cluster-randomized design
#' set.seed(4)
#' cluster_data <- sreg.rgen(
#'   n = 24, tau.vec = 0.5, n.strata = 3, cluster = TRUE
#' )
#' fit_cluster <- sreg(
#'   Y = cluster_data$Y, S = cluster_data$S, D = cluster_data$D,
#'   G.id = cluster_data$G.id, Ng = cluster_data$Ng
#' )
#' fit_cluster$tau.hat
sreg <- function(Y, S = NULL, D, G.id = NULL, Ng = NULL, X = NULL, HC1 = TRUE,
                 small.strata = FALSE, k = NULL) {
  check.data.types(Y, S, D, G.id, Ng, X)
  check.integers(S, D, G.id, Ng)
  boolean.check(HC1)
  boolean.check.ss(small.strata)
  if (!is.null(k) &&
      (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0 ||
       k != as.integer(k))) {
    stop("k must be NULL or a positive integer.")
  }

  if (is.null(Y)) {
    stop("Error: Observed outcomes have not been provided (Y = NULL). Please provide the vector of observed outcomes.")
  }
  if (is.null(D)) {
    stop("Error: Treatments have not been provided (D = NULL). Please provide the vector of treatments.")
  }
  if (!is.null(D)) {
    check.range(D)
  }
  if (!is.null(S)) {
    check.range(S)
  }

  mixed_design <- FALSE
  if (small.strata == TRUE && !is.null(S)) {
    if (!is.null(G.id)) {
      cluster_strata <- unique(data.frame(S = S, G.id = G.id))
      strata_sizes <- table(cluster_strata$S)
    } else {
      strata_sizes <- table(S)
    }
    unique_sizes <- unique(strata_sizes)
    mixed_design <- length(unique_sizes) > 1
    if (!mixed_design && !is.null(k) && as.numeric(unique_sizes[1]) != k) {
      stop(
        paste0(
          "The supplied small-stratum size k = ", k,
          " does not match the observed stratum size of ",
          as.numeric(unique_sizes[1]), "."
        )
      )
    }
  }
  if (small.strata == FALSE && !is.null(S)) {
    if (!is.null(G.id)) {
      cluster_strata <- unique(data.frame(S = S, G.id = G.id))
      strata_sizes <- table(cluster_strata$S)
    } else {
      strata_sizes <- table(S)
    }

    unique_sizes <- unique(strata_sizes)
    mixed_design <- length(unique_sizes) > 1
    if (length(unique_sizes) == 1 && unique_sizes[1] <= 5) {
      if (!is.null(G.id)) {
        warning("Warning: All strata have the same small number of clusters (e.g., matched pairs at the cluster level), but small.strata = FALSE. Consider setting small.strata = TRUE to apply estimators designed for such designs.")
      } else {
        warning("Warning: All strata have the same small number of observations (e.g., matched pairs), but small.strata = FALSE. Consider setting small.strata = TRUE to apply estimators designed for such designs.")
      }
    }
  }

  if ("Ng" %in% names(X)) {
    names(X)[names(X) == "Ng"] <- "Ng_1"
  }
  check.df <- tibble(Y, S, D, G.id, Ng, X)

  if (any(is.na(check.df))) {
    warning("Warning: The data contains one or more NA (or NaN) values. Proceeding while ignoring these values.")
  }
  clean.df <- na.omit(check.df)

  x.ind <- max(which(colnames(clean.df) %in% c("D", "G.id", "Ng")))
  suppressWarnings({
    Y <- clean.df$Y
    S <- clean.df$S
    D <- clean.df$D
    G.id <- clean.df$G.id
    Ng <- clean.df$Ng
  })
  if ((x.ind + 1) > ncol(clean.df)) {
    X <- NULL
  } else {
    X <- clean.df[, (x.ind + 1):ncol(clean.df)]
  }

  if ("Ng_1" %in% names(X)) {
    names(X)[names(X) == "Ng_1"] <- "Ng"
  }

  if (!is.null(S)) {
    if (min(S) != 1) {
      stop(paste0("Error: The strata should be indexed by {1, 2, 3, ...}. The minimum value in the provided data is ", min(S), "."))
    }
  }
  if (!is.null(D)) {
    if (min(D) != 0) {
      stop(paste0("Error: The treatments should be indexed by {0, 1, 2, ...}, where D = 0 denotes the control. The minimum value in the provided data is ", min(D), "."))
    }
  }

  if (is.null(G.id)) {
    if (small.strata == FALSE) {
      if (!is.null(X)) {
        if (!is.null(S)) {
          dta.temp <- data.frame(S, D, X)
          if (!check.within.strata.variation(dta.temp)) {
            warning("Warning: One or more covariates do not vary within one or more strata. Proceeding with unadjusted estimator.")
            X <- NULL
          }
          if (!check.within.stratatreatment.variation(dta.temp)) {
            warning("One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE. Proceeding with the unadjusted estimator. Note: please double-check whether the design is actually a small-strata design. If so, consider setting small.strata = TRUE to use the appropriate estimator.")
            X <- NULL
          }
        }
        if (is.null(S)) {
          S.temp <- rep(1, length(D))
          dta.temp <- data.frame("S" = S.temp, D, X)
          if (!check.within.strata.variation(dta.temp)) {
            warning("Warning: One or more covariates do not vary within one or more strata. Proceeding with unadjusted estimator.")
            X <- NULL
          }
          if (!check.within.stratatreatment.variation(dta.temp)) {
            warning("Warning: One or more covariates do not vary within one or more strata-treatment combinations. Proceeding with unadjusted estimator.")
            X <- NULL
          }
        }
      }
      if (!is.null(S)) {
        data_tst <- data.frame(Y = Y, D = D, S = S)
        if (!is.null(X)) {
          data_tst <- cbind(data_tst, X)
        }
        data_tst <- design.classifier(data_tst, S = S, small.strata = small.strata, k = k)
      }
      result <- res.sreg(Y, S, D, X, HC1)
      if (!is.null(result$lin.adj)) {
        if (any(sapply(result$ols.iter, function(x) any(is.na(x))))) {
          stop("Error: There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
        }
        if (any(sapply(result$se.rob, function(x) any(is.infinite(x))))) {
          stop("Error: Variance estimate is not finite. This may be due to too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)), remove covariate adjustments, or try setting HC1 = FALSE.")
        }
      }
    } else {
      if (is.null(S)) {
        stop("Error: Strata indicator variable has not been provided (S = NULL), but small.strata = TRUE. This estimator requires stratification. Either supply a valid strata indicator S, or set small.strata = FALSE to proceed without stratification.")
      }
      if (mixed_design) {
        result <- res.sreg.mixed(Y, S, D, X, HC1, small.strata, k)
      } else {
        result <- res.sreg.ss(Y, S, D, X, HC1)
      }
      if (!is.null(result$lin.adj)) {
        if (any(sapply(result$beta.hat, function(x) any(is.na(x))))) {
          stop("Error: There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
        }
        if (any(sapply(result$se.rob, function(x) any(is.infinite(x))))) {
          stop("Error: Variance estimate is not finite. This may be due to too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)), remove covariate adjustments, or try setting HC1 = FALSE.")
        }
      }
    }
  } else {
    check.cluster.lvl(G.id, S, D, Ng)
    if (small.strata == FALSE) {
      if (!is.null(X)) {
        if (!is.null(S)) {
          dta.temp <- data.frame(S, D, X)
          cluster_df <- data.frame(G.id, S, D, X) %>%
            group_by(G.id) %>%
            summarise(across(c(S, D), ~ first(.x)),
              across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop"
            )

          if (!check.within.strata.variation(cluster_df)) {
            warning("Warning: One or more covariates do not vary within one or more strata. Proceeding with unadjusted estimator.")
            X <- NULL
          }
          if (!check.within.stratatreatment.variation(cluster_df)) {
            warning("One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE. Proceeding with the unadjusted estimator. Note: please double-check whether the design is actually a small-strata design. If so, consider setting small.strata = TRUE to use the appropriate estimator.")
            X <- NULL
          }
        }
        if (is.null(S)) {
          S.temp <- rep(1, length(D))
          dta.temp <- data.frame("S" = S.temp, D, X)
          if (!check.within.strata.variation(dta.temp)) {
            warning("Warning: One or more covariates do not vary within one or more strata. Proceeding with unadjusted estimator.")
            X <- NULL
          }
          if (!check.within.stratatreatment.variation(dta.temp)) {
            warning("Warning: One or more covariates do not vary within one or more strata-treatment combinations. Proceeding with unadjusted estimator.")
            X <- NULL
          }
        }
      }
      result <- res.creg(Y, S, D, G.id, Ng, X, HC1)

      if (is.null(Ng)) {
        warning("Warning: Cluster sizes have not been provided (Ng = NULL). Ng is assumed to be equal to the number of available observations in every cluster g.")
      }

      if (!is.null(result$lin.adj)) {
        if (any(sapply(result$ols.iter, function(x) any(is.na(x))))) {
          stop("Error: There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
        }
        if (any(sapply(result$se.rob, function(x) any(is.infinite(x))))) {
          stop("Error: Variance estimate is not finite. This may be due to too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)), remove covariate adjustments, or try setting HC1 = FALSE.")
        }
        if (!check.cluster(data.frame("G.id" = result$data$G.id, result$lin.adj))) {
          warning("Warning: sreg cannot use individual-level covariates for covariate adjustment in cluster-randomized experiments. Any individual-level covariates have been aggregated to their cluster-level averages.")
        }
      }
      if (!is.null(S)) {
        if (is.null(Ng)) {
          Ng <- rep(NA_real_, length(Y)) # placeholder to allow data.frame construction
        }

        data_tst <- data.frame(Y = Y, D = D, S = S, G.id = G.id, Ng = Ng)
        if (!is.null(X)) {
          data_tst <- cbind(data_tst, X)
        }

        data_tst <- design.classifier(data_tst, S = S, G.id = G.id, small.strata = small.strata, k = k)
      }
    } else {
      if (is.null(S)) {
        stop("Error: Strata indicator variable has not been provided (S = NULL), but small.strata = TRUE. This estimator requires stratification. Either supply a valid strata indicator S, or set small.strata = FALSE to proceed without stratification.")
      }
      if (mixed_design) {
        result <- res.creg.mixed(Y, S, D, G.id, Ng, X, HC1, small.strata, k)
      } else {
        result <- res.creg.ss(Y, S, D, G.id, Ng, X, HC1)
      }
      if (is.null(Ng)) {
        warning("Warning: Cluster sizes have not been provided (Ng = NULL). Ng is assumed to be equal to the number of available observations in every cluster g.")
      }
      if (!is.null(result$lin.adj)) {
        if (!is.null(result$lin.adj)) {
          if (isTRUE(result$mixed.design)) {
            # Use only small strata G.id for the check
            gid_check <- result$res.small$data$G.id
          } else {
            # Use full data G.id
            gid_check <- result$data$G.id
          }
          if (!check.cluster(data.frame("G.id" = gid_check, result$lin.adj))) {
            warning("Warning: sreg cannot use individual-level covariates for covariate adjustment in cluster-randomized experiments. Any individual-level covariates have been aggregated to their cluster-level averages.")
          }
        }
        if (any(sapply(result$beta.hat, function(x) any(is.na(x))))) {
          stop("Error: There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
        }
        if (any(sapply(result$se.rob, function(x) any(is.infinite(x))))) {
          stop("Error: Variance estimate is not finite. This may be due to too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)), remove covariate adjustments, or try setting HC1 = FALSE.")
        }
      }
    }
  }
  result$small.strata <- small.strata
  result$HC1 <- HC1
  return(result)
}

#' Generate a Pseudo-Random Sample under the Stratified Block Randomization Design
#'
#' The function generates the observed outcomes, treatment assignments, strata indicators, cluster indicators, cluster sizes, and covariates for estimating the treatment effect within the context of a stratified block randomization design under the covariate-adaptive randomization (CAR).
#' @param n the total number of units when \code{cluster = FALSE}, or the total number of clusters when \code{cluster = TRUE}
#' @param Nmax a maximum size of generated clusters (maximum number of observations in a cluster)
#' @param n.strata an integer specifying the number of strata
#' @param tau.vec a numeric \eqn{1 \times |\mathcal A|} \code{vector} of treatment effects, where \eqn{|\mathcal A|} represents the number of treatments
#' @param gamma.vec a numeric \eqn{1 \times 3} \code{vector} of parameters corresponding to covariates
#' @param cluster a \code{TRUE/FALSE} argument indicating whether the dgp should use a cluster-level treatment assignment or individual-level
#' @param is.cov a \code{TRUE/FALSE} argument indicating whether the dgp should include covariates or not
#' @param small.strata a \code{TRUE/FALSE} argument indicating whether the data-generating process should use a small-strata design (e.g., matched pairs, n-tuples)
#' @param k an integer specifying the number of units per stratum when \code{small.strata = TRUE}
#' @param treat.sizes a numeric \eqn{1 \times (|\mathcal A| + 1)} \code{vector} specifying the number of units assigned to each treatment within a stratum; the first element corresponds to control units (\eqn{D = 0}), the second to the first treatment (\eqn{D = 1}), and so on. When omitted for a mixed design, the \code{k} positions are allocated across arms as evenly as possible
#' @param mixed.strata a \code{TRUE/FALSE} argument indicating whether to generate both small and large strata. When \code{TRUE}, \code{small.strata} is ignored, \code{n.small} units (or clusters when \code{cluster = TRUE}) are generated in strata of size \code{k}, and the remaining units or clusters are generated in \code{n.strata} large strata
#' @param n.small the number of units (or clusters when \code{cluster = TRUE}) assigned to the small-strata component when \code{mixed.strata = TRUE}. It must be divisible by \code{k}; if \code{NULL}, the largest multiple of \code{k} not exceeding half of \code{n} is used
#' @return A \code{data.frame} containing the generated values of the following variables (with \code{n} rows for individual-level designs and one row per observation within the \code{n} generated clusters for cluster-level designs):
#' \itemize{
#' \item \code{Y}: a numeric \eqn{n \times 1} \code{vector} of observed outcomes
#' \item \code{S}: a numeric \eqn{n \times 1} \code{vector} of strata indicators
#' \item \code{D}: a numeric \eqn{n \times 1} \code{vector} of treatments indexed by \eqn{\{0, 1, 2, \ldots\}}, where \eqn{\code{D} = 0} denotes the control
#' \item \code{G.id}: a numeric \eqn{n \times 1} \code{vector} of cluster indicators
#' \item \code{X}: a \code{data.frame} with columns representing the covariate values for every observation
#' }
#' @export
#'
#' @examples
#' ### Large-strata individual-level design
#' set.seed(11)
#' large_data <- sreg.rgen(
#'   n = 80, tau.vec = 0.5, n.strata = 4, cluster = FALSE
#' )
#'
#' ### Uniform matched pairs
#' set.seed(12)
#' pair_data <- sreg.rgen(
#'   n = 40, tau.vec = 0.5, cluster = FALSE,
#'   small.strata = TRUE, k = 2, treat.sizes = c(1, 1)
#' )
#'
#' ### Mixed design with 4-tuples and large strata
#' set.seed(13)
#' mixed_data <- sreg.rgen(
#'   n = 60, tau.vec = 0.5, cluster = FALSE,
#'   mixed.strata = TRUE, n.small = 40, k = 4,
#'   treat.sizes = c(2, 2), n.strata = 2
#' )
#' sort(table(mixed_data$S))
#'
#' ### For cluster assignment, n counts clusters rather than observations
#' set.seed(14)
#' cluster_data <- sreg.rgen(
#'   n = 20, tau.vec = 0.5, n.strata = 2, cluster = TRUE
#' )
#' length(unique(cluster_data$G.id))
sreg.rgen <- function(n, Nmax = 50, n.strata = 10,
                      tau.vec = c(0), gamma.vec = c(0.4, 0.2, 1),
                      cluster = TRUE, is.cov = TRUE, small.strata = FALSE, k = 3,
                      treat.sizes = c(1, 1, 1), mixed.strata = FALSE,
                      n.small = NULL) {
  if (!is.logical(mixed.strata) || length(mixed.strata) != 1L || is.na(mixed.strata)) {
    stop("mixed.strata must be either TRUE or FALSE.")
  }

  if (mixed.strata) {
    if (missing(treat.sizes)) {
      n.arms <- length(tau.vec) + 1L
      treat.sizes <- rep(k %/% n.arms, n.arms)
      remainder <- k %% n.arms
      if (remainder > 0L) {
        treat.sizes[seq_len(remainder)] <- treat.sizes[seq_len(remainder)] + 1L
      }
    }
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0 || n != as.integer(n)) {
      stop("n must be a positive integer.")
    }
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0 || k != as.integer(k)) {
      stop("k must be a positive integer.")
    }
    if (is.null(n.small)) {
      n.small <- floor((n / 2) / k) * k
    }
    if (!is.numeric(n.small) || length(n.small) != 1L || is.na(n.small) ||
        n.small <= 0 || n.small >= n || n.small != as.integer(n.small)) {
      stop("n.small must be a positive integer smaller than n.")
    }
    if (n.small %% k != 0) {
      stop("n.small must be divisible by k.")
    }
    valid.treat.sizes <- is.numeric(treat.sizes) &&
      length(treat.sizes) == length(tau.vec) + 1L &&
      !anyNA(treat.sizes) && all(treat.sizes >= 0) &&
      all(treat.sizes == as.integer(treat.sizes)) && sum(treat.sizes) == k
    if (!valid.treat.sizes) {
      stop("treat.sizes must be a nonnegative integer vector of length length(tau.vec) + 1 that sums to k.")
    }
    n.large <- n - n.small
    if (!is.numeric(n.strata) || length(n.strata) != 1L || is.na(n.strata) ||
        n.strata <= 0 || n.strata != as.integer(n.strata)) {
      stop("n.strata must be a positive integer.")
    }
    if (n.large <= n.strata * k) {
      stop("The large-strata component must contain more than k units (or clusters) per stratum on average; decrease n.small or n.strata.")
    }
    if ((n.small / k) <= n.strata) {
      stop("The mixed design must contain more small strata than large strata; increase n.small or decrease n.strata.")
    }

    small.data <- sreg.rgen(
      n = n.small, Nmax = Nmax, n.strata = n.strata,
      tau.vec = tau.vec, gamma.vec = gamma.vec, cluster = cluster,
      is.cov = is.cov, small.strata = TRUE, k = k,
      treat.sizes = treat.sizes
    )
    large.data <- sreg.rgen(
      n = n.large, Nmax = Nmax, n.strata = n.strata,
      tau.vec = tau.vec, gamma.vec = gamma.vec, cluster = cluster,
      is.cov = is.cov, small.strata = FALSE, k = k,
      treat.sizes = treat.sizes
    )

    large.data$S <- large.data$S + max(small.data$S)
    if (cluster) {
      large.data$G.id <- large.data$G.id + max(small.data$G.id)
    }

    common.names <- intersect(names(small.data), names(large.data))
    return(rbind(small.data[, common.names, drop = FALSE],
                 large.data[, common.names, drop = FALSE]))
  }

  n.treat <- length(tau.vec)
  if (cluster == TRUE) {
    if (small.strata == TRUE) {
      G <- n
      Nmax <- Nmax
      n.treat <- length(tau.vec)
      max.support <- Nmax / 10 - 1
      Ng <- gen.cluster.sizes(G, max.support)
      data_pot <- dgp.po.creg(Ng = Ng, tau.vec = tau.vec, G = G, gamma.vec = gamma.vec, n.treat = n.treat)
      data_obs <- dgp.obs.creg.ss(baseline = data_pot, n.treat = n.treat, k = k, treat_sizes = treat.sizes)
      data_cl <- data.frame("Y" = data_obs$Y, "D" = data_obs$D, data_obs$X, "S" = data_obs$S, "Ng" = data_obs$Ng, "G.id" = data_obs$G.id)
      Y <- data_cl$Y
      D <- data_cl$D
      X <- data.frame("x_1" = data_cl$x_1, "x_2" = data_cl$x_2, "Ng" = data_cl$Ng)
      S <- data_cl$S
      G.id <- data_cl$G.id
      Ng <- data_cl$Ng
      if (is.cov) {
        data.sim <- data.frame(Y, S, D, G.id, Ng, X)
      } else {
        data.sim <- data.frame(Y, S, D, G.id, Ng)
      }
    } else {
      G <- n
      Nmax <- Nmax
      n.treat <- length(tau.vec)
      max.support <- Nmax / 10 - 1
      Ng <- gen.cluster.sizes(G, max.support)
      data.pot <- dgp.po.creg(
        Ng = Ng, tau.vec = tau.vec, G = G,
        gamma.vec = gamma.vec, n.treat = n.treat
      )
      strata <- form.strata.creg(data.pot, n.strata)
      strata.set <- data.frame(strata)
      strata.set$S <- max.col(strata.set)
      pi.vec <- rep(c(1 / (n.treat + 1)), n.treat)
      data.sim <- dgp.obs.creg(data.pot, I.S = strata, pi.vec, n.treat)
      Y <- data.sim$Y
      D <- data.sim$D
      S <- data.sim$S
      X <- data.sim$X
      Ng <- data.sim$Ng
      G.id <- data.sim$G.id
      data.sim <- data.frame(Y, S, D, G.id, Ng, X)
    }
  } else {
    if (small.strata) {
      data_raw <- dgp.po.sreg(
        n = n, theta.vec = tau.vec, gamma.vec = gamma.vec,
        n.treat = n.treat, is.cov = is.cov
      )
      data_df <- dgp.obs.sreg.ss(data_raw, n.treat = n.treat, k = k, treat_sizes = treat.sizes)
      Y <- data_df$y
      D <- data_df$A
      X <- data.frame("x_1" = data_df$x_1, "x_2" = data_df$x_2)
      S <- data_df$block
      if (is.cov) {
        data.sim <- data.frame(Y, S, D, X)
      } else {
        data.sim <- data.frame(Y, S, D)
      }
    } else {
      n.treat <- length(tau.vec)
      pot.outcomes <- dgp.po.sreg(
        n = n, tau.vec, gamma.vec = gamma.vec,
        n.treat = n.treat, is.cov = is.cov
      )
      strata <- form.strata.sreg(pot.outcomes, num.strata = n.strata)
      strata_set <- data.frame(strata)
      strata_set$S <- max.col(strata_set)
      pi.vec <- rep(c(1 / (n.treat + 1)), n.treat) # vector of target proportions (equal allocation)
      data.test <- dgp.obs.sreg(pot.outcomes,
        I.S = strata,
        pi.vec = pi.vec, n.treat = n.treat, is.cov = is.cov
      )
      Y <- as.numeric(data.test$Y)
      D <- as.numeric(data.test$D)
      S <- strata_set$S
      if (is.cov == TRUE) {
        X <- data.test$X
        data.sim <- data.frame(Y, S, D, X)
      } else {
        data.sim <- data.frame(Y, S, D)
      }
    }
  }
  return(data.sim)
}
