#' Estimate Average Treatment Effects (ATEs) and Corresponding Standard Errors
#'
#' This function estimates the ATE(s) and the corresponding standard error(s) for a (collection of) treatment(s) relative to a control.
#' @import extraDistr
#' @import tidyr
#' @import dplyr
#' @import rlang
#' @importFrom stats aggregate coef lm pnorm qnorm rbeta rnorm na.omit
#' @importFrom utils packageVersion
#'
#' @param Y a numeric \eqn{n \times 1} vector of the observed outcomes
#' @param S a numeric \eqn{n \times 1} vector of strata indicators;  if \code{NULL} then the estimation is performed assuming no stratification
#' @param D a numeric \eqn{n \times 1} vector of treatments indexed by \eqn{\{0, 1, 2, \ldots\}}, where \eqn{\code{D} = 0} denotes the control
#' @param G.id a numeric \eqn{n \times 1} vector of cluster indicators; if \code{NULL} then estimation is performed assuming treatment is assigned at the individual level
#' @param Ng a numeric \eqn{n \times 1} vector of cluster sizes; if \code{NULL} then \code{Ng} is assumed to be equal to the number of available observations in every cluster
#' @param X a data frame with columns representing the covariate values for every observation; if \code{NULL} then the estimator without linear adjustments is applied. (Note: sreg cannot use individual-level covariates for covariate adjustment in cluster-randomized experiments. Any individual-level covariates will be aggregated to their cluster-level averages)
#' @param HC1 a \code{TRUE/FALSE} logical argument indicating whether the small sample correction should be applied to the variance estimator
#'
#'
#' @return An object of class \code{sreg} that is a list containing the following elements:
#' \itemize{
#' \item \code{tau.hat}: a \eqn{1 \times |\mathcal A|} vector of ATE estimates, where \eqn{|\mathcal A|} represents the number of treatments
#' \item \code{se.rob}: a \eqn{1 \times |\mathcal A|} vector of standard errors estimates, where \eqn{|\mathcal A|} represents the number of treatments
#' \item \code{t.stat}: a \eqn{1 \times |\mathcal A|} vector of \eqn{t}-statistics, where \eqn{|\mathcal A|} represents the number of treatments
#' \item \code{p.value}: a \eqn{1 \times |\mathcal A|} vector of corresponding \eqn{p}-values, where \eqn{|\mathcal A|} represents the number of treatments
#' \item \code{CI.left}: a \eqn{1 \times |\mathcal A|} vector of the left bounds of the 95\% as. confidence interval
#' \item \code{CI.right}: a \eqn{1 \times |\mathcal A|} vector of the right bounds of the 95\% as. confidence interval
#' \item \code{data}: an original data of the form \code{data.frame(Y, S, D, G.id, Ng, X)}
#' \item \code{lin.adj}: a data frame representing the covariates that were used in implementing linear adjustments
#' }
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
#' ### Example 2. Data taken from Chong et al. (2016).
#' ## Data description
#' ?AEJapp
#' ## Upload the data from the package
#' data("AEJapp")
#' data <- AEJapp
#' head(data)
#' ## Replicate the empirical illustration from (Bugni et al, 2019)
#' # Prepare the data
#' Y <- data$gradesq34
#' D <- data$treatment
#' S <- data$class_level
#' data.clean <- data.frame(Y, D, S)
#' data.clean <- data.clean %>%
#'   mutate(D = ifelse(D == 3, 0, D))
#' Y <- data.clean$Y
#' D <- data.clean$D
#' S <- data.clean$S
#' # Look at the frequency table
#' table(D = data.clean$D, S = data.clean$S)
#' # Replicate the results from (Bugni et al, 2019)
#' result <- sreg::sreg(Y, S, D)
#'
#' ## Besides that, it is possible to add linear adjustments (covariates)
#' pills <- data$pills_taken
#' age <- data$age_months
#' data.clean <- data.frame(Y, D, S, pills, age)
#' data.clean <- data.clean %>%
#'   mutate(D = ifelse(D == 3, 0, D))
#' Y <- data.clean$Y
#' D <- data.clean$D
#' S <- data.clean$S
#' X <- data.frame("pills" = data.clean$pills, "age" = data.clean$age)
#' result <- sreg::sreg(Y, S, D, G.id = NULL, X = X)
sreg <- function(Y, S = NULL, D, G.id = NULL, Ng = NULL, X = NULL, HC1 = TRUE) {

  if (is.null(Y)) {
    stop("Observed outcomes have not been provided (Y = NULL). Please provide the vector of observed outcomes.")
  }
  if (is.null(D)) {
    stop("Treatments have not been provided (D = NULL). Please provide the vector of treatments.")
  }
  check.df <- tibble(Y, S, D, G.id, Ng, X) 
  if (any(is.na(check.df))) 
  {
    warning("The data contains one or more NA (or NaN) values. Proceeding while ignoring these values.")
  }
  clean.df <- na.omit(check.df)
  x.ind <- max(which(colnames(clean.df) %in% c("D", "G.id", "Ng")))
  
  suppressWarnings({
  Y <- clean.df$Y
  S <- clean.df$S
  D <- clean.df$D 
  G.id <- clean.df$G.id
  Ng <- clean.df$Ng
  if ((x.ind + 1) >= ncol(clean.df))
  {
  X <- NULL
  }else{
  X <- clean.df[, (x.ind + 1):ncol(clean.df)]
  }
  })

  if (is.null(G.id)) {
    result <- res.sreg(Y, S, D, X, HC1)
    summary.sreg(result)
  } else {

    result <- res.creg(Y, S, D, G.id, Ng, X, HC1)
    summary.creg(result)
  }
  return(result)
}

#' Generates a Pseudo-Random Sample under the Stratified Block Randomization
#'
#' The function generates the observed outcomes, treatment assignments, strata indicators, cluster indicators, cluster sizes, and covariates for estimating the treatment effect within the context of a stratified block randomization design under the covariate-adaptive randomization (CAR).
#' @param n number of observations
#' @param Nmax maximum size of clusters
#' @param n.strata number of strata
#' @param tau.vec a numeric vector of treatment effects
#' @param gamma.vec a numeric vector of parameters
#' @param cluster a \code{TRUE/FALSE} argument indicating whether the dgp should include clusters or not
#' @param is.cov a \code{TRUE/FALSE} argument indicating whether the dgp should include covariates or not
#'
#' @return a data frame containing the generated values of \code{Y}, \code{S}, \code{D}, \code{G.id}, \code{Ng}, and \code{X}
#'
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
    Ng <- gen.cluster.sizes(G, max.support) # [, 1]
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
    n.treat <- length(tau.vec) # number of treatments
    pot.outcomes <- dgp.po.sreg(
      n = n, tau.vec, gamma.vec = gamma.vec,
      n.treat = n.treat, is.cov = is.cov
    ) # generate pot. outcomes and W
    strata <- form.strata.sreg(pot.outcomes, num.strata = n.strata) # generate strata
    strata_set <- data.frame(strata) # generate strata
    strata_set$S <- max.col(strata_set) # generate strata
    pi.vec <- rep(c(1 / (n.treat + 1)), n.treat) # vector of target proportions (equal allocation)
    data.test <- dgp.obs.sreg(pot.outcomes,
      I.S = strata, # simulate observed outcomes
      pi.vec = pi.vec, n.treat = n.treat, is.cov = is.cov
    ) # simulate observed outcomes
    Y <- as.numeric(data.test$Y) # Y
    D <- as.numeric(data.test$D) # D
    S <- strata_set$S # S
    if (is.cov == TRUE) {
      X <- data.test$X
      data.sim <- data.frame(Y, S, D, X)
    } else {
      data.sim <- data.frame(Y, S, D)
    }
  }
  return(data.sim)
}
