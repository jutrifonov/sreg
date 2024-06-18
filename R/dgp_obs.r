#-------------------------------------------------------------------
# %#     Generates observed outcomes,
# %#     by taking as input the potential outcomes,
# %#     matrix of strata assignments, pi.vec, and
# %#     number of treatments
#----------------------------------------------------------------------
dgp.obs.sreg <- function(baseline, I.S, pi.vec, n.treat, is.cov = TRUE)
#----------------------------------------------------------------------
{
  if (n.treat != length(pi.vec)) {
    stop("The number of treatments doesn't match the length of vector pi.vec.")
  }
  num.strata <- ncol(I.S)
  n <- length(baseline$Y.0)
  A <- cbind(rep(0, n))
  l.seq <- num.strata / 2

  pi.matr <- matrix(1, ncol = num.strata, nrow = n.treat)
  pi.matr.w <- pi.matr * pi.vec

  for (k in 1:num.strata)
  {
    index <- which(I.S[, k] == 1)
    ns <- length(index)

    A[index] <- gen.treat.sreg(pi.matr.w, ns, k)
  }

  for (a in 0:n.treat)
  {
    assign(paste("Y.", a, sep = ""), baseline[[paste("Y.", a, sep = "")]])
  }
  formula <- gen.rubin.formula.sreg(n.treat)
  Y.obs <- eval(parse(text = formula))

  if (is.cov == TRUE) {
    ret.list <- list(
      "Y" = Y.obs,
      "D" = A,
      "X" = baseline$X
    )
  } else {
    ret.list <- list(
      "Y" = Y.obs,
      "D" = A
    )
  }
  return(ret.list)
}
#-------------------------------------------------------------------
dgp.obs.creg <- function(baseline, I.S, pi.vec, n.treat)
#-------------------------------------------------------------------
{
  if (n.treat != length(pi.vec)) {
    stop("The number of treatments doesn't match the length of vector pi.vec.")
  }
  num.strata <- ncol(I.S)
  n <- baseline$G
  A <- cbind(rep(0, n))
  l.seq <- num.strata / 2

  pi.matr <- matrix(1, ncol = num.strata, nrow = n.treat)
  pi.matr.w <- pi.matr * pi.vec

  for (k in 1:num.strata)
  {
    index <- which(I.S[, k] == 1)
    ns <- length(index)

    A[index] <- gen.treat.creg(pi.matr.w, ns, k)
  }
  strata.set <- data.frame(I.S)
  strata.set$S <- max.col(strata.set)
  cluster.indicator <- baseline$cl.id
  G.seq <- seq(c(1:baseline$G))
  data.short <- data.frame(
    "cl.id" = G.seq, A, S = strata.set$S, Ng = baseline$Ng,
    baseline$X
  )
  data.long <- data.frame("cl.id" = cluster.indicator)
  merged.data <- merge(data.long, data.short, by = "cl.id")
  length(merged.data$A)
  A <- merged.data$A
  S <- merged.data$S
  X <- merged.data[5:ncol(merged.data)]
  Ng <- merged.data$Ng

  for (a in 0:n.treat)
  {
    assign(paste("Y.", a, sep = ""), baseline[[paste("Yig.", a, sep = "")]])
  }
  formula <- gen.rubin.formula.creg(n.treat)
  Y.obs <- eval(parse(text = formula))

  ret.list <- list(
    "Y"           = Y.obs,
    "D"           = A,
    "S"           = S,
    "Z.2"         = baseline$Z.g.2,
    "X"           = X,
    "Ng"          = Ng,
    "G.id"        = cluster.indicator,
    "cl.lvl.data" = data.short
  )
  return(ret.list)
}
