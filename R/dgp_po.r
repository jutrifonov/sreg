#-----------------------------------------------------------------------
# %#   Potential outcomes generation
#-----------------------------------------------------------------------
dgp.po.sreg <- function(n, theta.vec, gamma.vec, n.treat, is.cov = TRUE)
#-----------------------------------------------------------------------
{
  if (n.treat != length(theta.vec)) {
    stop("The number of treatments doesn't match the length of vector theta.vec.")
  }
  for (a in seq_along(theta.vec))
  {
    assign(paste("mu.", a, sep = ""), theta.vec[a])
  }
  eps.0 <- rnorm(n)
  for (a in 1:n.treat)
  {
    assign(paste("eps.", a, sep = ""), rnorm(n))
  }

  W <- sqrt(20) * (rbeta(n, 2, 2) - 1 / 2)
  x_1 <- rnorm(n, mean = 5, sd = 2)
  x_2 <- rnorm(n, mean = 2, sd = 1)

  if (is.cov == TRUE) {
    X <- data.frame(x_1, x_2)
    m.0 <- gamma.vec[1] * W + gamma.vec[2] * x_1 + gamma.vec[3] * x_2
  } else {
    m.0 <- gamma.vec[1] * W
  }

  for (a in 1:n.treat)
  {
    assign(paste("m.", a, sep = ""), m.0)
  }

  Y.0 <- m.0 + eps.0

  for (a in 1:n.treat)
  {
    formula <- paste0("mu.", a, "+m.", a, "+eps.", a)
    result <- eval(parse(text = formula))
    assign(paste("Y.", a, sep = ""), result)
  }

  if (is.cov == TRUE) {
    ret.names <- c(
      paste("Y.", 0:n.treat, sep = ""),
      "W", "X", paste("m.", 0:n.treat, sep = ""),
      paste("mu.", 1:n.treat, sep = "")
    )
  } else {
    ret.names <- c(
      paste("Y.", 0:n.treat, sep = ""),
      "W", paste("m.", 0:n.treat, sep = ""),
      paste("mu.", 1:n.treat, sep = "")
    )
  }
  ret.list <- mget(ret.names)
  return(ret.list)
}
#------------------------------------------------------------------
dgp.po.creg <- function(Ng, G, tau.vec, sigma1 = sqrt(2),
                        gamma.vec = c(0.4, 0.2, 1), n.treat)
#------------------------------------------------------------------
{
  for (a in seq_along(tau.vec))
  {
    assign(paste("mu.", a, sep = ""), tau.vec[a])
  }

  beta.rv <- rbeta(G, 2, 2)
  Z.g.2 <- (beta.rv - 0.5) * sqrt(20)
  x_1 <- (rnorm(G, mean = 5, sd = 2) - 5) / 2
  x_2 <- (rnorm(G, mean = 2, sd = 1) - 2) / 1
  X <- data.frame(x_1, x_2)

  cluster.indicator <- rep(c(1:G), Ng)
  cl.id <- cluster.indicator
  total.sample <- length(cluster.indicator)

  for (a in 1:n.treat)
  {
    assign(paste("epsilon.ig.", a, sep = ""), rnorm(total.sample, 0, sigma1))
  }

  epsilon.ig.0 <- rnorm(total.sample, 0, 1)

  m.0 <- gamma.vec[1] * Z.g.2 + gamma.vec[2] * x_1 + gamma.vec[3] * x_2

  for (a in 1:n.treat)
  {
    assign(paste("m.", a, sep = ""), m.0)
  }

  Yig.0 <- rep(m.0, Ng) + epsilon.ig.0

  for (a in 1:n.treat)
  {
    formula <- paste0("mu.", a, "+rep(m.", a, ",Ng)", "+epsilon.ig.", a)
    result <- eval(parse(text = formula))
    assign(paste("Yig.", a, sep = ""), result)
  }

  ret.names <- c(
    paste("Yig.", 0:n.treat, sep = ""),
    "Z.g.2", "X", "G", "Ng", "cl.id", paste("m.", 0:n.treat, sep = ""),
    paste("mu.", 1:n.treat, sep = "")
  )

  ret.list <- mget(ret.names)
  return(ret.list)
}
