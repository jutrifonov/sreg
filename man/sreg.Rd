% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{sreg}
\alias{sreg}
\title{Estimate Average Treatment Effects (ATEs) and Corresponding Standard Errors}
\usage{
sreg(
  Y,
  S = NULL,
  D,
  G.id = NULL,
  Ng = NULL,
  X = NULL,
  HC1 = TRUE,
  small.strata = FALSE
)
}
\arguments{
\item{Y}{a numeric \eqn{n \times 1} \code{vector/matrix/data.frame/tibble} of the observed outcomes}

\item{S}{a numeric \eqn{n \times 1} \code{vector/matrix/data.frame/tibble} of strata indicators indexed by \eqn{\{1, 2, 3, \ldots\}};  if \code{NULL} then the estimation is performed assuming no stratification}

\item{D}{a numeric \eqn{n \times 1} \code{vector/matrix/data.frame/tibble} of treatments indexed by \eqn{\{0, 1, 2, \ldots\}}, where \eqn{\code{D} = 0} denotes the control}

\item{G.id}{a numeric \eqn{n \times 1} \code{vector/matrix/data.frame/tibble} of cluster indicators; if \code{NULL} then estimation is performed assuming treatment is assigned at the individual level}

\item{Ng}{a numeric \eqn{n \times 1} \code{vector/matrix/data.frame/tibble} of cluster sizes; if \code{NULL} then \code{Ng} is assumed to be equal to the number of available observations in every cluster}

\item{X}{a \code{matrix/data.frame/tibble} with columns representing the covariate values for every observation; if \code{NULL} then the estimator without linear adjustments is applied. (Note: \code{sreg} cannot use individual-level covariates for covariate adjustment in cluster-randomized experiments. Any individual-level covariates will be aggregated to their cluster-level averages)}

\item{HC1}{a \code{TRUE/FALSE} logical argument indicating whether the small sample correction should be applied to the variance estimator}

\item{small.strata}{a \code{TRUE/FALSE} logical argument indicating whether the estimators for small strata (i.e., strata with few units, such as matched pairs or n-tuples) should be used.}
}
\value{
An object of class \code{sreg} that is a list containing the following elements:
\itemize{
\item \code{tau.hat}: a \eqn{1 \times |\mathcal A|} \code{vector} of ATE estimates, where \eqn{|\mathcal A|} represents the number of treatments
\item \code{se.rob}: a \eqn{1 \times |\mathcal A|} \code{vector} of standard errors estimates, where \eqn{|\mathcal A|} represents the number of treatments
\item \code{t.stat}: a \eqn{1 \times |\mathcal A|} \code{vector} of \eqn{t}-statistics, where \eqn{|\mathcal A|} represents the number of treatments
\item \code{p.value}: a \eqn{1 \times |\mathcal A|} \code{vector} of corresponding \eqn{p}-values, where \eqn{|\mathcal A|} represents the number of treatments
\item \code{CI.left}: a \eqn{1 \times |\mathcal A|} \code{vector} of the left bounds of the 95\% as. confidence interval
\item \code{CI.right}: a \eqn{1 \times |\mathcal A|} \code{vector} of the right bounds of the 95\% as. confidence interval
\item \code{data}: an original data of the form \code{data.frame(Y, S, D, G.id, Ng, X)}
\item \code{lin.adj}: a \code{data.frame} representing the covariates that were used in implementing linear adjustments
\item \code{small.strata}: a \code{TRUE/FALSE} logical argument indicating whether the estimators for small strata (e.g., matched pairs or n-tuples) were used
\item \code{HC1}: a \code{TRUE/FALSE} logical argument indicating whether the small sample correction (HC1) was applied to the variance estimator
}
}
\description{
Estimate the ATE(s) and the corresponding standard error(s) for a (collection of) treatment(s) relative to a control.
}
\examples{
library("sreg")
library("dplyr")
library("haven")
### Example 1. Simulated Data.
data <- sreg.rgen(n = 1000, tau.vec = c(0), n.strata = 4, cluster = FALSE)
Y <- data$Y
S <- data$S
D <- data$D
X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
result <- sreg(Y, S, D, G.id = NULL, Ng = NULL, X)
print(result)
### Example 2. Empirical Data.
?AEJapp
data("AEJapp")
data <- AEJapp
head(data)
Y <- data$gradesq34
D <- data$treatment
S <- data$class_level
data.clean <- data.frame(Y, D, S)
data.clean <- data.clean \%>\%
  mutate(D = ifelse(D == 3, 0, D))
Y <- data.clean$Y
D <- data.clean$D
S <- data.clean$S
table(D = data.clean$D, S = data.clean$S)
result <- sreg(Y, S, D)
print(result)
pills <- data$pills_taken
age <- data$age_months
data.clean <- data.frame(Y, D, S, pills, age)
data.clean <- data.clean \%>\%
  mutate(D = ifelse(D == 3, 0, D))
Y <- data.clean$Y
D <- data.clean$D
S <- data.clean$S
X <- data.frame("pills" = data.clean$pills, "age" = data.clean$age)
result <- sreg(Y, S, D, G.id = NULL, X = X)
print(result)
### Example 3. Matched Pairs (small strata).
data <- sreg.rgen(
  n = 1000, tau.vec = c(1.2), cluster = FALSE,
  small.strata = TRUE, k = 2, treat.sizes = c(1, 1)
)
Y <- data$Y
S <- data$S
D <- data$D
X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
result <- sreg(Y = Y, S = S, D = D, X = X, small.strata = TRUE)
print(result)
}
\references{
Bugni, F. A., Canay, I. A., and Shaikh, A. M. (2018). Inference Under Covariate-Adaptive Randomization. \emph{Journal of the American Statistical Association}, 113(524), 1784–1796, \doi{10.1080/01621459.2017.1375934}.

Bugni, F., Canay, I., Shaikh, A., and Tabord-Meehan, M. (2024+). Inference for Cluster Randomized Experiments with Non-ignorable Cluster Sizes. \emph{Forthcoming in the Journal of Political Economy: Microeconomics}, \doi{10.48550/arXiv.2204.08356}.

Jiang, L., Linton, O. B., Tang, H., and Zhang, Y. (2023+). Improving Estimation Efficiency via Regression-Adjustment in Covariate-Adaptive Randomizations with Imperfect Compliance. \emph{Forthcoming in Review of Economics and Statistics}, \doi{10.48550/arXiv.2204.08356}.

Bai, Y., Jiang, L., Romano, J. P., Shaikh, A. M., and Zhang, Y. (2024). Covariate adjustment in experiments with matched pairs. \emph{Journal of Econometrics}, 241(1), \doi{10.1016/j.jeconom.2024.105740}.

Liu, J. (2024). Inference for Two-stage Experiments under Covariate-Adaptive Randomization. \doi{10.48550/arXiv.2301.09016}.

Cytrynbaum, M. (2024). Covariate Adjustment in Stratified Experiments. \emph{Quantitative Economics}, 15(4), 971–998,  \doi{10.3982/QE2475}.
}
\author{
\strong{Authors}:

Juri Trifonov \email{jutrifonov@u.northwestern.edu}

Yuehao Bai \email{yuehao.bai@usc.edu}

Azeem Shaikh \email{amshaikh@uchicago.edu}

Max Tabord-Meehan \email{maxtm@uchicago.edu}


\strong{Maintainer}:

Juri Trifonov \email{jutrifonov@u.northwestern.edu}
}
