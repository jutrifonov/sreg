source("~/Desktop/sreg/R/creg.R")
source("~/Desktop/sreg/R/sreg.R")



#' Estimates the ATE
#' @import extraDistr
#'
#' @param Y a numeric vector of the observed outcomes
#' @param S a numeric vector of strata indicators
#' @param D a numeric vector of treatments
#' @param G.id a numeric vector of cluster indicators
#' @param Ng a numeric vector of cluster sizes
#' @param X a data frame of covariates
#'
#' @return a list containing the results
#' @export
#'
#' @examples
#' test <- sreg(Y,S,D,X)
sreg <- function(Y,S,D,G.id = NULL, Ng = NULL, X=NULL)
{
  if (is.null(G.id) | is.null(Ng))
  {
    result <- res.sreg(Y,S,D,X)
    summary.sreg(result)
  }else{
    result <- res.creg(Y,S,D,G.id,Ng,X)
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
#'
#' @return a data frame containing the results
#' @export
#'
#' @examples
#' data <- sreg.rgen(n=1000, tau.vec = c(0), n.strata=4, cluster = T)
sreg.rgen <- function(n, Nmax=50, n.strata,
                      tau.vec = c(0), gamma.vec = c(0.4, 0.2, 1), cluster = TRUE)
{
  if (cluster == T)
  {
    G <- n
    Nmax <- 50
    n.treat <- length(tau.vec)
    max.support = Nmax/10-1
    Ng <- gen.cluster.sizes(G, max.support)[,1]
    #Ng <- rep(Nmax, G)                                                            # uncomment and comment the previous line for a equal-size design
    data.pot <- dgp.po.creg(Ng=Ng, tau.vec = (tau.vec / 0.5), G = G,
                             gamma.vec = gamma.vec, n.treat=n.treat)
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
  }else{
    n.treat <- length(tau.vec)      # number of treatments
    pot.outcomes <- dgp.po.sreg(n = n, tau.vec, gamma = gamma.vec, n.treat = n.treat)   # generate pot. outcomes and W
    strata <- form.strata.sreg(pot.outcomes, num.strata = n.strata)                     # generate strata
    strata_set <- data.frame(strata)                                                    # generate strata
    strata_set$S <- max.col(strata_set)                                                 # generate strata
    pi.vec <- rep(c(1 / (n.treat + 1)), n.treat)                                        # vector of target proportions (equal allocation)
    #pi.vec <- c(0.1, 0.35, 0.2)
    data.test <- dgp.obs.sreg(pot.outcomes,I.S = strata,                                # simulate observed outcomes
                         pi.vec = pi.vec, n.treat = n.treat)                            # simulate observed outcomes
    Y <- as.numeric(data.test$Y)                                                        # Y
    D <- as.numeric(data.test$D)                                                        # D
    S <- strata_set$S                                                                   # S
    X <- data.test$X
    data.sim <- data.frame(Y, S, D, X)
  }
  return(data.sim)
}



