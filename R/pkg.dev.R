#' Estimates the ATE
#'
#' @import extraDistr
#' @import tidyr
#' @import dplyr
#'
#' @param Y a numeric vector of the observed outcomes
#' @param S a numeric vector of strata indicators
#' @param D a numeric vector of treatments
#' @param G.id a numeric vector of cluster indicators
#' @param Ng a numeric vector of cluster sizes
#' @param X a data frame of covariates
#' @param Ng.cov a TRUE/FALSE argument indicating whether the Ng should be included as the only covariate in linear adjustments
#'  when X = NULL
#' @param HC1 a TRUE/FALSE argument indicating whether the small sample correction should be applied to the variance estimator.
#'
#'
#' @return a list containing the results
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
#' result <- sreg::sreg(Y, S, D, HC1 = TRUE)
#'
#' ## Besides that, it is possible to add linear adjustments (covariates)
#' x_1 <- data$pills_taken
#' x_2 <- data$age_months
#' data.clean <- data.frame(Y, D, S, x_1, x_2)
#' data.clean <- data.clean %>%
#'   mutate(D = ifelse(D == 3, 0, D))
#' Y <- data.clean$Y
#' D <- data.clean$D
#' S <- data.clean$S
#' X <- data.frame(data.clean$x_1, data.clean$x_2)
#' result <- sreg::sreg(Y, S, D, X = X, HC1 = TRUE)
sreg <- function(Y, S = NULL, D, G.id = NULL, Ng = NULL, X = NULL, Ng.cov = FALSE, HC1 = FALSE) {
  if (is.null(G.id) | is.null(Ng)) {
    result <- res.sreg(Y, S, D, X, HC1)
    summary.sreg(result)
  } else {
    result <- res.creg(Y, S, D, G.id, Ng, X, Ng.cov, HC1)
    summary.creg(result)
  }
  return(result)
}

#' Generates a pseudo-random sample for estimating ATE
#'
#' @param n number of observations
#' @param Nmax maximum size of clusters
#' @param n.strata number of strata
#' @param tau.vec a numeric vector of treatment effects
#' @param gamma.vec a numeric vector of parameters
#' @param cluster a TRUE/FALSE argument indicating whether the dgp should include clusters or not
#' @param is.cov a TRUE/FALSE argument indicating whether the dgp should include covariates or not
#'
#' @return a data frame containing the results
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
    Ng <- gen.cluster.sizes(G, max.support)#[, 1]
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
      n = n, tau.vec, gamma = gamma.vec,
      n.treat = n.treat, is.cov = is.cov
    ) # generate pot. outcomes and W
    strata <- form.strata.sreg(pot.outcomes, num.strata = n.strata) # generate strata
    strata_set <- data.frame(strata) # generate strata
    strata_set$S <- max.col(strata_set) # generate strata
    pi.vec <- rep(c(1 / (n.treat + 1)), n.treat) # vector of target proportions (equal allocation)
    # pi.vec <- c(0.1, 0.35, 0.2)
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
