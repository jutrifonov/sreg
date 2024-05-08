#-------------------------------------------------------------------
# %#   Generate the formula for Y.obs (Rubin model)
#-------------------------------------------------------------------
gen.rubin.formula.sreg <- function(n.treat)
#-------------------------------------------------------------------
{
  A.values <- 0:n.treat

  formula <- "Y.obs = "

  for (a in A.values)
  {
    if (a == 0) {
      formula <- paste(formula, paste0("Y.", a, " * (A == 0)"))
    } else {
      formula <- paste(formula, paste0("Y.", a, " * (A == ", a, ")"))
    }

    if (a < n.treat) {
      formula <- paste(formula, " + ")
    }
  }
  return(formula)
}
#-------------------------------------------------------------------
gen.rubin.formula.creg <- function(n.treat)
#-------------------------------------------------------------------
{
  A.values <- 0:n.treat

  formula <- "Y.obs = "

  for (a in A.values)
  {
    if (a == 0) {
      formula <- paste(formula, paste0("Y.", a, " * (A == 0)"))
    } else {
      formula <- paste(formula, paste0("Y.", a, " * (A == ", a, ")"))
    }

    if (a < n.treat) {
      formula <- paste(formula, " + ")
    }
  }
  return(formula)
}
