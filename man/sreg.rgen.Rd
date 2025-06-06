% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{sreg.rgen}
\alias{sreg.rgen}
\title{Generate a Pseudo-Random Sample under the Stratified Block Randomization Design}
\usage{
sreg.rgen(
  n,
  Nmax = 50,
  n.strata = 10,
  tau.vec = c(0),
  gamma.vec = c(0.4, 0.2, 1),
  cluster = TRUE,
  is.cov = TRUE,
  small.strata = FALSE,
  k = 3,
  treat.sizes = c(1, 1, 1)
)
}
\arguments{
\item{n}{a total number of observations in a sample}

\item{Nmax}{a maximum size of generated clusters (maximum number of observations in a cluster)}

\item{n.strata}{an integer specifying the number of strata}

\item{tau.vec}{a numeric \eqn{1 \times |\mathcal A|} \code{vector} of treatment effects, where \eqn{|\mathcal A|} represents the number of treatments}

\item{gamma.vec}{a numeric \eqn{1 \times 3} \code{vector} of parameters corresponding to covariates}

\item{cluster}{a \code{TRUE/FALSE} argument indicating whether the dgp should use a cluster-level treatment assignment or individual-level}

\item{is.cov}{a \code{TRUE/FALSE} argument indicating whether the dgp should include covariates or not}

\item{small.strata}{a \code{TRUE/FALSE} argument indicating whether the data-generating process should use a small-strata design (e.g., matched pairs, n-tuples)}

\item{k}{an integer specifying the number of units per stratum when \code{small.strata = TRUE}}

\item{treat.sizes}{a numeric \eqn{1 \times (|\mathcal A| + 1)} \code{vector} specifying the number of units assigned to each treatment within a stratum; the first element corresponds to control units (\eqn{D = 0}), the second to the first treatment (\eqn{D = 1}), and so on}
}
\value{
An object that is a `data.frame` with \eqn{n} observations containing the generated values of the following variables:
\itemize{
\item \code{Y}: a numeric \eqn{n \times 1} \code{vector} of observed outcomes
\item \code{S}: a numeric \eqn{n \times 1} \code{vector} of strata indicators
\item \code{D}: a numeric \eqn{n \times 1} \code{vector} of treatments indexed by \eqn{\{0, 1, 2, \ldots\}}, where \eqn{\code{D} = 0} denotes the control
\item \code{G.id}: a numeric \eqn{n \times 1} \code{vector} of cluster indicators
\item \code{X}: a \code{data.frame} with columns representing the covariate values for every observation
}
}
\description{
The function generates the observed outcomes, treatment assignments, strata indicators, cluster indicators, cluster sizes, and covariates for estimating the treatment effect within the context of a stratified block randomization design under the covariate-adaptive randomization (CAR).
}
\examples{
data <- sreg.rgen(n = 1000, tau.vec = c(0), n.strata = 4, cluster = TRUE)
}
