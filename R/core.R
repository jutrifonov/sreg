#' Estimate Average Treatment Effects (ATEs) and Corresponding Standard Errors
#'
#' Estimate the ATE(s) and the corresponding standard error(s) for a (collection of) treatment(s) relative to a control.
#' @import extraDistr
#' @import tidyr
#' @import dplyr
#' @import rlang
#' @import cli
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
#'
#'
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
#' }
#'
#' @references
#' Bugni, F. A., Canay, I. A., and Shaikh, A. M. (2018). Inference Under Covariate-Adaptive Randomization. \emph{Journal of the American Statistical Association}, 113(524), 1784–1796, \doi{10.1080/01621459.2017.1375934}.
#'
#' Bugni, F., Canay, I., Shaikh, A., and Tabord-Meehan, M. (2024+). Inference for Cluster Randomized Experiments with Non-ignorable Cluster Sizes. \emph{Forthcoming in the Journal of Political Economy: Microeconomics}, \doi{10.48550/arXiv.2204.08356}.
#'
#' Jiang, L., Linton, O. B., Tang, H., and Zhang, Y. (2023+). Improving Estimation Efficiency via Regression-Adjustment in Covariate-Adaptive Randomizations with Imperfect Compliance. \emph{Forthcoming in Review of Economics and Statistics}, \doi{10.48550/arXiv.2204.08356}.
#'
#' @author
#' \strong{Authors}:
#'
#' Juri Trifonov \email{jutrifonov@uchicago.edu}
#'
#' Yuehao Bai \email{yuehao.bai@usc.edu}
#'
#' Azeem Shaikh \email{amshaikh@uchicago.edu}
#'
#' Max Tabord-Meehan \email{maxtm@uchicago.edu}
#'
#'
#' \strong{Maintainer}:
#'
#' Juri Trifonov \email{jutrifonov@uchicago.edu}
#'
#'
#'
#'
#' @export
#'
#' @examples
#' library("sreg")
#' library("dplyr")
#' library("haven")
#' ### Example 1. Simulated Data.
#' data <- sreg.rgen(n = 1000, tau.vec = c(0), n.strata = 4, cluster = FALSE)
#' Y <- data$Y
#' S <- data$S
#' D <- data$D
#' X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
#' result <- sreg(Y, S, D, G.id = NULL, Ng = NULL, X)
#' print(result)
#' ### Example 2. Empirical Data.
#' ?AEJapp
#' data("AEJapp")
#' data <- AEJapp
#' head(data)
#' Y <- data$gradesq34
#' D <- data$treatment
#' S <- data$class_level
#' data.clean <- data.frame(Y, D, S)
#' data.clean <- data.clean %>%
#'   mutate(D = ifelse(D == 3, 0, D))
#' Y <- data.clean$Y
#' D <- data.clean$D
#' S <- data.clean$S
#' table(D = data.clean$D, S = data.clean$S)
#' result <- sreg(Y, S, D)
#' print(result)
#' pills <- data$pills_taken
#' age <- data$age_months
#' data.clean <- data.frame(Y, D, S, pills, age)
#' data.clean <- data.clean %>%
#'   mutate(D = ifelse(D == 3, 0, D))
#' Y <- data.clean$Y
#' D <- data.clean$D
#' S <- data.clean$S
#' X <- data.frame("pills" = data.clean$pills, "age" = data.clean$age)
#' result <- sreg(Y, S, D, G.id = NULL, X = X)
#' print(result)
sreg <- function(Y, S = NULL, D, G.id = NULL, Ng = NULL, X = NULL, HC1 = TRUE) {
  check.data.types(Y, S, D, G.id, Ng, X)
  check.integers(S, D, G.id, Ng)
  boolean.check(HC1)

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
    if ((x.ind + 1) >= ncol(clean.df)) {
      X <- NULL
    } else {
      X <- clean.df[, (x.ind + 1):ncol(clean.df)]
    }
  })

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
    result <- res.sreg(Y, S, D, X, HC1)
    if (any(sapply(result$ols.iter, function(x) any(is.na(x))))) {
      stop("Error: There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
    }
  } else {
    check.cluster.lvl(G.id, S, D, Ng)
    result <- res.creg(Y, S, D, G.id, Ng, X, HC1)
    if (is.null(result$data$Ng)) {
      warning("Warning: Cluster sizes have not been provided (Ng = NULL). Ng is assumed to be equal to the number of available observations in every cluster g.")
    }
    if (any(sapply(result$ols.iter, function(x) any(is.na(x))))) {
      stop("Error: There are too many covariates relative to the number of observations. Please reduce the number of covariates (k = ncol(X)) or consider estimating the model without covariate adjustments.")
    }
    if (!is.null(result$lin.adj)) {
      if (!check.cluster(data.frame("G.id" = result$data$G.id, result$lin.adj))) {
        warning("Warning: sreg cannot use individual-level covariates for covariate adjustment in cluster-randomized experiments. Any individual-level covariates have been aggregated to their cluster-level averages.")
      }
    }
  }
  return(result)
}

#' Generate a Pseudo-Random Sample under the Stratified Block Randomization Design
#'
#' The function generates the observed outcomes, treatment assignments, strata indicators, cluster indicators, cluster sizes, and covariates for estimating the treatment effect within the context of a stratified block randomization design under the covariate-adaptive randomization (CAR).
#' @param n a total number of observations in a sample
#' @param Nmax a maximum size of generated clusters (maximum number of observations in a cluster)
#' @param n.strata an integer specifying the number of strata
#' @param tau.vec a numeric \eqn{1 \times |\mathcal A|} \code{vector} of treatment effects, where \eqn{|\mathcal A|} represents the number of treatments
#' @param gamma.vec a numeric \eqn{1 \times 3} \code{vector} of parameters corresponding to covariates
#' @param cluster a \code{TRUE/FALSE} argument indicating whether the dgp should use a cluster-level treatment assignment or individual-level
#' @param is.cov a \code{TRUE/FALSE} argument indicating whether the dgp should include covariates or not
#'
#' @return An object that is a `data.frame` with \eqn{n} observations containing the generated values of the following variables:
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
#' data <- sreg.rgen(n = 1000, tau.vec = c(0), n.strata = 4, cluster = TRUE)
sreg.rgen <- function(n, Nmax = 50, n.strata,
                      tau.vec = c(0), gamma.vec = c(0.4, 0.2, 1),
                      cluster = TRUE, is.cov = TRUE) {
  if (cluster == T) {
    G <- n
    Nmax <- Nmax
    n.treat <- length(tau.vec)
    max.support <- Nmax / 10 - 1
    Ng <- gen.cluster.sizes(G, max.support)
    # Ng <- rep(Nmax, G)                                                            # uncomment and comment the previous line for a equal-size design
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
  return(data.sim)
}
