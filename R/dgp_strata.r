#-------------------------------------------------------------------
# %# Function that generates the strata from the observed covariates
#-------------------------------------------------------------------
form.strata.sreg <- function(baseline, num.strata)
#-------------------------------------------------------------------
{
  n <- length(baseline$Y.0)
  W <- baseline$W
  bounds <- seq(-2.25, 2.25, length.out = num.strata + 1)
  I.S <- matrix(0, n, num.strata)
  for (s in 1:num.strata)
  {
    I.S[, s] <- (W > bounds[s]) * (W <= bounds[s + 1])
  }

  return(I.S)
}
#---------------------------------------------------------------------
form.strata.creg <- function(baseline, num.strata)
#---------------------------------------------------------------------
{
  n <- baseline$G
  W <- baseline$Z.g.2
  bounds <- seq(min(W), max(W), length.out = num.strata + 1)
  I.S <- matrix(0, n, num.strata)
  for (s in 1:num.strata)
  {
    I.S[, s] <- (W > bounds[s]) * (W <= bounds[s + 1])
  }

  return(I.S)
}
