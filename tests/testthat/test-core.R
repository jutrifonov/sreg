test_that("simulations without clusters work", {
  set.seed(123)
  sim_data <- sreg.rgen(n = 1000, tau.vec = c(0.2, 0.5), n.strata = 10, cluster = F, is.cov = TRUE)
  Y <- sim_data$Y
  S <- sim_data$S
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)

  invisible(capture.output({
    result <- sreg::sreg(Y, S, D, G.id = NULL, Ng = NULL, X = X)
  }))

  expect_equal(round(result$tau.hat, 7), c(0.1580814, 0.4846882))
  expect_equal(round(result$se.rob, 8), c(0.07524021, 0.07616346))

  invisible(capture.output({
    result <- sreg::sreg(Y, S, D, G.id = NULL, Ng = NULL, X = NULL)
  }))

  expect_equal(round(result$tau.hat, 7), c(0.1627114, 0.4948722))
  expect_equal(round(result$se.rob, 7), c(0.1105611, 0.1124072))

  invisible(capture.output({
    result <- sreg::sreg(Y, S = NULL, D, G.id = NULL, Ng = NULL, X = X)
  }))

  expect_equal(round(result$tau.hat, 7), c(0.1578917, 0.4963735))
  expect_equal(round(result$se.rob, 8), c(0.08255663, 0.08320655))

  invisible(capture.output({
    result <- sreg::sreg(Y, S = NULL, D = D, G.id = NULL, Ng = NULL, X = NULL)
  }))

  expect_equal(round(result$tau.hat, 7), c(0.1685108, 0.5022035))
  expect_equal(round(result$se.rob, 7), c(0.1145915, 0.1161482))

  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X, HC1 = 5)
    })),
    "Error: The value of HC must be either TRUE or FALSE."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X, HC1 = "TRUE")
    })),
    "Error: The value of HC must be either TRUE or FALSE."
  )

  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(as.list(Y), S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  S[2] <- 2.5
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = NULL)
    })),
    "must contain only integer values."
  )
  S <- sim_data$S
  D[5] <- 0.5
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "must contain only integer values."
  )
  S[3] <- 5.5
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "must contain only integer values."
  )
  S <- sim_data$S
  D <- sim_data$D

  X[10:12, ] <- NaN
  Y[1:10] <- NaN
  msg <- "ignoring these values"
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = NULL, Ng = NULL, X = X, HC1 = TRUE)
    })),
    msg
  )
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  Y <- sim_data$Y
  S[12] <- NA
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = NULL, Ng = NULL, X = X, HC1 = TRUE)
    })),
    msg
  )
  S <- sim_data$S
  S[19] <- NaN
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = NULL, Ng = NULL, X = X, HC1 = TRUE)
    })),
    msg
  )
  S <- sim_data$S
  S[1] <- 0
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "Error: The strata should be indexed by"
  )
  S[10] <- -1
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "Error: The strata should be indexed by"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = NULL)
    })),
    "Error: The strata should be indexed by"
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = NULL, Ng = NULL, X = X)
    }))
  )
  S <- sim_data$S
  D[1:3] <- -1
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "Error: The treatments should be indexed by"
  )
  D[4] <- -2
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "Error: The treatments should be indexed by"
  )
  D <- sim_data$D
  D[1] <- -1
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "Error: The treatments should be indexed by"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = NULL)
    })),
    "Error: The treatments should be indexed by"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = NULL, Ng = NULL, X = NULL)
    })),
    "Error: The treatments should be indexed by"
  )
  D <- sim_data$D
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    }))
  )
})

test_that("simulations with clusters work", {
  set.seed(123) # fix the random seed
  sim_data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50, is.cov = TRUE
  )
  Y <- sim_data$Y
  S <- sim_data$S
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)

  G.id <- sim_data$G.id
  Ng <- sim_data$Ng

  invisible(capture.output({
    result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
  }))

  expect_equal(round(result$tau.hat, 7), c(0.3350243, 0.8488310))
  expect_equal(round(result$se.rob, 7), c(0.1529870, 0.2183116))

  invisible(capture.output({
    result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = NULL, HC1 = TRUE)
  }))

  expect_equal(round(result$tau.hat, 7), c(0.3003475, 1.0875317))
  expect_equal(round(result$se.rob, 7), c(0.3126555, 0.3266127))

  invisible(capture.output({
    result <- sreg::sreg(Y, S = NULL, D, G.id, Ng, X = X, HC1 = TRUE)
  }))

  expect_equal(round(result$tau.hat, 7), c(0.1687191, 0.8792778))
  expect_equal(round(result$se.rob, 7), c(0.1641397, 0.1636365))

  invisible(capture.output({
    result <- sreg::sreg(Y, S = NULL, D, G.id, Ng, X = NULL, HC1 = TRUE)
  }))

  expect_equal(round(result$tau.hat, 7), c(0.2613412, 1.0528130))
  expect_equal(round(result$se.rob, 7), c(0.3047530, 0.3084314))

  invisible(capture.output({
    result <- sreg::sreg(Y, S = S, D, G.id = NULL, Ng = NULL, X = X)
  }))

  expect_equal(round(result$tau.hat, 7), c(0.1548014, 0.6955543))
  expect_equal(round(result$se.rob, 8), c(0.06201006, 0.05871032))

  X <- data.frame("Ng" = sim_data$Ng, "x_1" = sim_data$x_1, "x_2" = sim_data$x_2)

  invisible(capture.output({
    result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.3041215, 0.7184981))
  expect_equal(round(result$se.rob, 7), c(0.1381628, 0.1343521))


  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X, HC1 = 5)
    })),
    "Error: The value of HC must be either TRUE or FALSE."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X, HC1 = "TRUE")
    })),
    "Error: The value of HC must be either TRUE or FALSE."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = NULL, X = X, HC1 = "TRUE")
    })),
    "Error: The value of HC must be either TRUE or FALSE."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X, HC1 = "TRUE")
    })),
    "Error: The value of HC must be either TRUE or FALSE."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y = NULL, S = S, D, G.id = G.id, Ng = NULL, X = X)
    })),
    "Observed outcomes have not been provided"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = NULL, G.id = NULL, Ng = Ng, X = X)
    })),
    "Treatments have not been provided"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = as.list(X))
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = as.character(S), D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(as.character(Y), S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = as.character(D), G.id = G.id, Ng = Ng, X = X)
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = as.character(X))
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = as.list(S), D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(as.list(Y), S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = as.list(D), G.id = G.id, Ng = Ng, X = X)
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = as.list(G.id), Ng = Ng, X = X)
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = as.list(Ng), X = X)
    })),
    "variable has a different type than matrix, numeric vector, or data frame."
  )
  S <- sim_data$S
  S[1] <- 0
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "Error: The strata should be indexed by"
  )
  S[10] <- -1
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "Error: The strata should be indexed by"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = NULL, X = X)
    })),
    "Error: The strata should be indexed by"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = NULL, X = NULL)
    })),
    "Error: The strata should be indexed by"
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = G.id, Ng = Ng, X = X)
    }))
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = G.id, Ng = Ng, X = NULL)
    }))
  )
  S <- sim_data$S
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(as.matrix(Y), S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    }))
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(as.numeric(Y), S = as.numeric(S), D = as.integer(D), G.id = as.integer(G.id), Ng = as.integer(Ng), X = X)
    }))
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(tibble(Y), S = as.numeric(S), D = as.vector(D), G.id = as.integer(G.id), Ng = as.matrix(Ng), X = tibble(X))
    }))
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(tibble(Y), S = as.data.frame(S), D = as.data.frame(D), G.id = as.data.frame(G.id), Ng = as.data.frame(Ng), X = as.data.frame(X))
    }))
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(as.data.frame(Y), S = as.data.frame(S), D = as.data.frame(D), G.id = as.data.frame(G.id), Ng = as.data.frame(Ng), X = as.data.frame(X))
    }))
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(as.matrix(Y), S = as.matrix(S), D = as.vector(D), G.id = as.matrix(G.id), Ng = as.matrix(Ng), X = X)
    }))
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(tibble(Y), S = tibble(S), D = tibble(D), G.id = tibble(G.id), Ng = tibble(Ng), X = tibble(X))
    }))
  )
  S[2] <- 2.5
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "must contain only integer values."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = NULL, X = X)
    })),
    "must contain only integer values."
  )
  S <- sim_data$S
  D[5] <- 1.5
  D[7] <- 1.5
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "must contain only integer values."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = NULL, X = X)
    })),
    "must contain only integer values."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = G.id, Ng = NULL, X = NULL)
    })),
    "must contain only integer values."
  )
  S[2] <- 2.5
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = NULL)
    })),
    "must contain only integer values."
  )
  S <- sim_data$S
  D <- sim_data$D
  Ng[2] <- 37.64
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "must contain only integer values."
  )
  Ng <- sim_data$Ng
  G.id[29] <- 29.23480
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "must contain only integer values."
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = G.id, Ng = Ng, X = NULL)
    })),
    "must contain only integer values."
  )
  set.seed(123) # fix the random seed
  sim_data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50
  )
  Y <- sim_data$Y
  S <- sim_data$S
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  G.id <- sim_data$G.id
  Ng <- sim_data$Ng
  D[1:40] <- -1
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "Error: The treatments should be indexed by"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = NULL)
    })),
    "Error: The treatments should be indexed by"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "Error: The treatments should be indexed by"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = G.id, Ng = Ng, X = NULL)
    })),
    "Error: The treatments should be indexed by"
  )
  D <- sim_data$D
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    }))
  )
})

test_that("One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE", {
  set.seed(123)
  sim_data <- sreg.rgen(n = 25, tau.vec = c(0.2, 0.5), n.strata = 3, cluster = F, is.cov = TRUE)
  Y <- sim_data$Y
  S <- sim_data$S
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)

  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S, D, G.id = NULL, Ng = NULL, X = X)
    })),
    "One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE"
  )

  set.seed(123)
  sim_data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 10, cluster = T, Nmax = 50
  )
  Y <- sim_data$Y
  S <- sim_data$S
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  G.id <- sim_data$G.id
  Ng <- sim_data$Ng

  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    "One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE."
  )
  set.seed(123)
  sim_data <- sreg.rgen(
    n = 50, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50
  )
  Y <- sim_data$Y
  S <- sim_data$S
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  G.id <- sim_data$G.id
  Ng <- sim_data$Ng

  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    "There are too many covariates relative to the number of observations."
  )
})

test_that("individual level X warning works", {
  set.seed(123)
  sim_data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50
  )
  Y <- sim_data$Y
  S <- sim_data$S
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  G.id <- sim_data$G.id
  Ng <- sim_data$Ng
  X[1, 1] <- 2.3894

  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    "sreg cannot use individual-level covariates"
  )
})

test_that("no cluster sizes warning works", {
  set.seed(123)
  sim_data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50
  )
  Y <- sim_data$Y
  S <- sim_data$S
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  G.id <- sim_data$G.id
  Ng <- sim_data$Ng

  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = NULL, X = X, HC1 = TRUE)
    })),
    "Cluster sizes have not been provided"
  )
})

test_that("data contains one or more NA (or NaN) values warning works", {
  set.seed(123) # fix the random seed
  # Generate a pseudo-random sample with clusters and two treatments = c(0.2, 0.8)

  sim_data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50
  )
  Y <- sim_data$Y
  S <- sim_data$S
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  Y[1:10] <- NA
  G.id <- sim_data$G.id
  Ng <- sim_data$Ng

  msg <- "ignoring these values"
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  Y <- sim_data$Y
  X[1, 1] <- NA
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  X[1:5, ] <- NA
  X[10:12, ] <- NA
  Y[1:10] <- NaN
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  X[10:12, ] <- NaN
  Y <- sim_data$Y
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  G.id[100:105] <- NA
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  G.id <- sim_data$G.id
  Ng[4:5] <- NaN
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  Ng <- sim_data$Ng
  D[24:25] <- NA
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  D <- sim_data$D
  S[23] <- NaN
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
})

test_that("skipped values in range of S/D works", {
  set.seed(123)
  sim_data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.5),
    n.strata = 5, cluster = F, is.cov = TRUE
  )
  Y <- sim_data$Y
  S <- sim_data$S
  S[S == 4] <- 1
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "There are skipped values in the range"
  )
  set.seed(123)
  sim_data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.5),
    n.strata = 5, cluster = T, is.cov = TRUE
  )
  Y <- sim_data$Y
  S <- sim_data$S
  S[S == 3] <- 4
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  G.id <- sim_data$G.id
  Ng <- sim_data$Ng
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "There are skipped values in the range"
  )
  S <- sim_data$S
  G.id[41:50] <- 3
  S[41:50] <- 3
  Ng[41:50] <- 30
  D[41:50] <- 2
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "sreg cannot use individual-level covariates for covariate adjustment"
  )
})

test_that("non cluster-level error for S, D, Ng works", {
  set.seed(123)
  sim_data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.5),
    n.strata = 5, cluster = T, is.cov = TRUE
  )
  Y <- sim_data$Y
  S <- sim_data$S
  D <- sim_data$D
  X <- data.frame("x_1" = sim_data$x_1, "x_2" = sim_data$x_2)
  G.id <- sim_data$G.id
  Ng <- sim_data$Ng
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    }))
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = G.id, Ng = Ng, X = X)
    }))
  )

  S[41] <- 3
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "The values for S, D, and Ng must be consistent within each cluster"
  )
  S <- sim_data$S
  D[13:20] <- 1
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "The values for S, D, and Ng must be consistent within each cluster"
  )
  D <- sim_data$D
  Ng[300] <- 30
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "The values for S, D, and Ng must be consistent within each cluster"
  )
  Ng <- sim_data$Ng
  S[41] <- 3
  D[13:20] <- 1
  Ng[300] <- 30
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "The values for S, D, and Ng must be consistent within each cluster"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "The values for S, D, and Ng must be consistent within each cluster"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S, D = D, G.id = G.id, Ng = NULL, X = X)
    })),
    "The values for S, D, and Ng must be consistent within each cluster"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = NULL, D = D, G.id = G.id, Ng = NULL, X = X)
    })),
    "The values for S, D, and Ng must be consistent within each cluster"
  )
})


test_that("empirical example works", {
  library(haven)
  data("AEJapp")
  sim_data <- AEJapp
  Y <- sim_data$gradesq34
  D <- sim_data$treatment
  S <- sim_data$class_level
  pills <- sim_data$pills_taken
  age <- sim_data$age_months
  data.clean <- data.frame(Y, D, S, pills, age)
  data.clean <- data.clean %>%
    mutate(D = ifelse(D == 3, 0, D))
  Y <- data.clean$Y
  D <- data.clean$D
  S <- data.clean$S
  X <- data.frame("pills" = data.clean$pills, "age" = data.clean$age)

  invisible(capture.output({
    result <- sreg::sreg(Y, S, D)
  }))

  expect_equal(round(result$tau.hat, 8), c(-0.05112971, 0.40903373))
  expect_equal(round(result$se.rob, 7), c(0.2064541, 0.2065146))

  invisible(capture.output({
    result <- sreg::sreg(Y, S, D, X = X)
  }))

  expect_equal(round(result$tau.hat, 8), c(-0.02861589, 0.34608688))
  expect_equal(round(result$se.rob, 7), c(0.1796427, 0.1836229))

  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y = NULL, S, D)
    })),
    "Observed outcomes have not been provided"
  )
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S, D = NULL)
    })),
    "Treatments have not been provided"
  )
})
test_that("dgp.po warning work", {
  set.seed(123)
  expect_error(
    invisible(capture.output({
      result <- dgp.po.sreg(n = 100, theta.vec = c(0, 0.5), n.treat = 3, gamma.vec = c(0.4, 0.2, 1))
    })),
    "The number of treatments doesn't match the length of vector theta.vec."
  )
})

# SREG 2.0 testing enviroment #
### no clusters
test_that("data: small strata, option: small strata", {
  # sim_data: small strata, option: small strata
  set.seed(123)
  tau.vec <- c(0.2, 0.8)
  n.treat <- length(tau.vec)
  n_1 <- 300
  sim_data <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 4, cluster = FALSE, small.strata = TRUE, treat.sizes = c(1, 1, 1), k = 3)
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, HC1 = TRUE, small.strata = TRUE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.4045371, 0.7833705))
  expect_equal(round(result$se.rob, 7), c(0.2215056, 0.1960352))
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, HC1 = FALSE, small.strata = TRUE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.4045371, 0.7833705))
  expect_equal(round(result$se.rob, 7), c(0.2215056, 0.1960352))
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = TRUE, small.strata = TRUE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.3073525, 0.8224690))
  expect_equal(round(result$se.rob, 7), c(0.1496169, 0.1450199))
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = TRUE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.3073525, 0.8224690))
  expect_equal(round(result$se.rob, 7), c(0.1472941, 0.1427597))
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1), HC1 = TRUE, small.strata = TRUE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.4156365, 0.8045963))
  expect_equal(round(result$se.rob, 7), c(0.2160681, 0.1947304))
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1), HC1 = FALSE, small.strata = TRUE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.4156365, 0.8045963))
  expect_equal(round(result$se.rob, 7), c(0.2137443, 0.1926630))

  expect_error(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, D = sim_data$D, S = NULL, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = TRUE)
    })),
    "Strata indicator variable has not been provided (S = NULL), but small.strata = TRUE. This estimator requires stratification. Either supply a valid strata indicator S, or set small.strata = FALSE to proceed without stratification",
    fixed = TRUE
  )
  expect_error(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, D = NULL, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = TRUE)
    })),
    "Treatments have not been provided (D = NULL). Please provide the vector of treatments.",
    fixed = TRUE
  )
  expect_error(
    invisible(capture.output({
      result <- sreg(Y = NULL, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = TRUE)
    })),
    "Observed outcomes have not been provided (Y = NULL). Please provide the vector of observed outcomes.",
    fixed = TRUE
  )
})

test_that("data: big strata, option: big strata", {
  # sim_data: small strata, option: small strata
  set.seed(123)
  tau.vec <- c(0.2, 0.9, 1.5)
  n.treat <- length(tau.vec)
  n_1 <- 1000
  sim_data <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 4, cluster = FALSE, small.strata = FALSE, treat.sizes = c(1, 1, 1), k = 3)
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, HC1 = TRUE, small.strata = FALSE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.4237099, 1.0070806, 1.4392785))
  expect_equal(round(result$se.rob, 7), c(0.1319973, 0.1300569, 0.1325378))
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, HC1 = FALSE, small.strata = FALSE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.4237099, 1.0070806, 1.4392785))
  expect_equal(round(result$se.rob, 7), c(0.1309383, 0.1290127, 0.1314774))
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = TRUE, small.strata = FALSE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.2114763, 0.8411237, 1.4103687))
  expect_equal(round(result$se.rob, 7), c(0.0887169, 0.0886237, 0.0888706))
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = FALSE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.2114763, 0.8411237, 1.4103687))
  expect_equal(round(result$se.rob, 7), c(0.0880103, 0.0879149, 0.0881632))
  invisible(capture.output({
    result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1), HC1 = TRUE, small.strata = FALSE)
  }))
  expect_equal(round(result$tau.hat, 7), c(0.4104883, 0.9986809, 1.4669068))
  expect_equal(round(result$se.rob, 7), c(0.1278097, 0.1237450, 0.1269005))
})

test_that("data: small strata, option: big strata", {
  # sim_data: small strata, option: big strata
  set.seed(123)
  tau.vec <- c(0.2, 0.8)
  n.treat <- length(tau.vec)
  n_1 <- 900
  sim_data <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 4, cluster = FALSE, small.strata = TRUE, treat.sizes = c(1, 1, 1), k = 3)

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, HC1 = TRUE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(-0.0866568, 0.6284132))
  expect_equal(round(result$se.rob, 7), c(0.0676520, 0.0689729))
  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, HC1 = TRUE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("All strata have the same small number of observations", warnings, fixed = TRUE)))
      expect_true(any(grepl("At least 25% of strata are small", warnings, fixed = TRUE)))
      expect_true(any(grepl("HC1 adjustment unstable or undefined", warnings, fixed = TRUE)))
    },
    regexp = NA
  )

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = TRUE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(-0.0866568, 0.6284132))
  expect_equal(round(result$se.rob, 7), c(0.0676520, 0.0689729))

  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = TRUE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("All strata have the same small number of observations", warnings, fixed = TRUE)))
      expect_true(any(grepl("One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("At least 25% of strata are small", warnings, fixed = TRUE)))
      expect_true(any(grepl("HC1 adjustment unstable or undefined", warnings, fixed = TRUE)))
    },
    regexp = NA
  )

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(-0.0866568, 0.6284132))
  expect_equal(round(result$se.rob, 7), c(0.0676520, 0.0689729))

  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("All strata have the same small number of observations", warnings, fixed = TRUE)))
      expect_true(any(grepl("One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("At least 25% of strata are small", warnings, fixed = TRUE)))
    },
    regexp = NA
  )

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1), HC1 = TRUE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(-0.0866568, 0.6284132))
  expect_equal(round(result$se.rob, 7), c(0.0676520, 0.0689729))

  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1), HC1 = TRUE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("All strata have the same small number of observations", warnings, fixed = TRUE)))
      expect_true(any(grepl("One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("At least 25% of strata are small", warnings, fixed = TRUE)))
      expect_true(any(grepl("HC1 adjustment unstable or undefined", warnings, fixed = TRUE)))
    },
    regexp = NA
  )
})

test_that("data: big strata, option: small strata", {
  set.seed(123)
  tau.vec <- c(0.2, 0.9)
  n.treat <- length(tau.vec)
  n_1 <- 1000
  sim_data <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 20, cluster = FALSE, small.strata = FALSE, treat.sizes = c(1, 1), k = 2)
  expect_error(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = TRUE)
    })),
    "Invalid input: Either all strata are large or too few strata qualify as 'small' to proceed with small.strata = TRUE. Please set small.strata = FALSE.",
    fixed = TRUE
  )
  expect_error(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = TRUE, small.strata = TRUE)
    })),
    "Invalid input: Either all strata are large or too few strata qualify as 'small' to proceed with small.strata = TRUE. Please set small.strata = FALSE.",
    fixed = TRUE
  )
  expect_error(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, D = sim_data$D, S = sim_data$S, X = NULL, HC1 = TRUE, small.strata = TRUE)
    })),
    "Invalid input: Either all strata are large or too few strata qualify as 'small' to proceed with small.strata = TRUE. Please set small.strata = FALSE.",
    fixed = TRUE
  )
  expect_error(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, D = sim_data$D, S = NULL, X = NULL, HC1 = TRUE, small.strata = TRUE)
    })),
    "Strata indicator variable has not been provided (S = NULL), but small.strata = TRUE.",
    fixed = TRUE
  )
  expect_error(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, D = sim_data$D, S = NULL, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = TRUE, small.strata = TRUE)
    })),
    "Strata indicator variable has not been provided (S = NULL), but small.strata = TRUE.",
    fixed = TRUE
  )
})

test_that("data: mixed design, option: small strata", {
  set.seed(123)
  tau.vec <- c(0.2, 0.9)
  n.treat <- length(tau.vec)
  n_1 <- 3000
  n_2 <- 50
  data_s <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 4, cluster = FALSE, small.strata = TRUE, treat.sizes = c(1, 1, 1), k = 3)


  data_b <- sreg.rgen(n = n_2, tau.vec = tau.vec, n.strata = 2, cluster = FALSE, small.strata = FALSE, treat.sizes = c(1, 1), k = 2)
  # Step 1: Get the max stratum ID in data_s
  max_id <- max(data_s$S)

  # Step 2: Get unique strata in data_b and assign new IDs
  unique_b_strata <- sort(unique(data_b$S))
  num_b_strata <- length(unique_b_strata)

  # Create a named mapping from old to new stratum IDs
  new_ids <- seq(max_id + 1, max_id + num_b_strata)
  stratum_map <- setNames(new_ids, unique_b_strata)

  # Step 3: Relabel data_b$S
  data_b$S <- stratum_map[as.character(data_b$S)]

  sim_data <- rbind(data_s, data_b)

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = TRUE, small.strata = TRUE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1679006, 0.8306431))
  expect_equal(round(result$se.rob, 7), c(0.0449655, 0.0451653))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = NULL, HC1 = TRUE, small.strata = TRUE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1475503, 0.8076297))
  expect_equal(round(result$se.rob, 7), c(0.0649087, 0.0660474))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = TRUE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1679006, 0.8306431))
  expect_equal(round(result$se.rob, 7), c(0.0448600, 0.0450745))

  expect_warning(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = TRUE)
    })),
    "Mixed design detected: at least 25% of strata are small. Weighted estimators will be used.",
    fixed = TRUE
  )
  expect_warning(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = NULL, HC1 = FALSE, small.strata = TRUE)
    })),
    "Mixed design detected: at least 25% of strata are small. Weighted estimators will be used.",
    fixed = TRUE
  )
  expect_warning(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = NULL, HC1 = TRUE, small.strata = TRUE)
    })),
    "Mixed design detected: at least 25% of strata are small. Weighted estimators will be used.",
    fixed = TRUE
  )
  expect_error(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, S = NULL, D = sim_data$D, X = data.frame(sim_data$x_1), HC1 = TRUE, small.strata = TRUE)
    })),
    "Strata indicator variable has not been provided (S = NULL), but small.strata = TRUE.",
    fixed = TRUE
  )
  expect_error(
    invisible(capture.output({
      result <- sreg(Y = sim_data$Y, S = sim_data$S, D = NULL, X = data.frame(sim_data$x_1), HC1 = TRUE, small.strata = TRUE)
    })),
    "Treatments have not been provided (D = NULL). Please provide the vector of treatments.",
    fixed = TRUE
  )
})

test_that("data: mixed design, option: big strata", {
  set.seed(123)
  tau.vec <- c(0.2, 0.9)
  n.treat <- length(tau.vec)
  n_1 <- 3000
  n_2 <- 50
  data_s <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 4, cluster = FALSE, small.strata = TRUE, treat.sizes = c(1, 1, 1), k = 3)

  data_b <- sreg.rgen(n = n_2, tau.vec = tau.vec, n.strata = 2, cluster = FALSE, small.strata = FALSE, treat.sizes = c(1, 1), k = 2)
  # Step 1: Get the max stratum ID in data_s
  max_id <- max(data_s$S)

  # Step 2: Get unique strata in data_b and assign new IDs
  unique_b_strata <- sort(unique(data_b$S))
  num_b_strata <- length(unique_b_strata)

  # Create a named mapping from old to new stratum IDs
  new_ids <- seq(max_id + 1, max_id + num_b_strata)
  stratum_map <- setNames(new_ids, unique_b_strata)

  # Step 3: Relabel data_b$S
  data_b$S <- stratum_map[as.character(data_b$S)]

  sim_data <- rbind(data_s, data_b)

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1475503, 0.8076297))
  expect_equal(round(result$se.rob, 7), c(0.0376357, 0.0382750))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = TRUE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1475503, 0.8076297))
  expect_equal(round(result$se.rob, 7), c(0.0792277, 0.0798430))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = NULL, HC1 = FALSE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1475503, 0.8076297))
  expect_equal(round(result$se.rob, 7), c(0.0376357, 0.0382750))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = NULL, HC1 = TRUE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1475503, 0.8076297))
  expect_equal(round(result$se.rob, 7), c(0.0792277, 0.0798430))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, S = NULL, D = sim_data$D, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = TRUE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1644669, 0.8312292))
  expect_equal(round(result$se.rob, 7), c(0.0483432, 0.0477374))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = sim_data$Y, S = NULL, D = sim_data$D, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1644669, 0.8312292))
  expect_equal(round(result$se.rob, 7), c(0.0483194, 0.0477140))


  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = FALSE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("At least 25% of strata are small", warnings, fixed = TRUE)))
    },
    regexp = NA
  )

  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = data.frame(sim_data$x_1, sim_data$x_2), HC1 = TRUE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("At least 25% of strata are small", warnings, fixed = TRUE)))
    },
    regexp = NA
  )

  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = NULL, HC1 = FALSE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("At least 25% of strata are small", warnings, fixed = TRUE)))
    },
    regexp = NA
  )

  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = NULL, HC1 = TRUE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("At least 25% of strata are small", warnings, fixed = TRUE)))
    },
    regexp = NA
  )

  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = sim_data$Y, S = sim_data$S, D = sim_data$D, X = NULL, HC1 = FALSE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("At least 25% of strata are small", warnings, fixed = TRUE)))
    },
    regexp = NA
  )
})

### WITH clusters
test_that("data: small strata, option: small strata", {
  set.seed(123)
  tau.vec <- c(0.2, 0.8)
  n.treat <- length(tau.vec)
  n_1 <- 900

  data_sim <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 4, cluster = TRUE, small.strata = TRUE, treat.sizes = c(1, 1, 1), k = 3)

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, G.id = data_sim$G.id, Ng = data_sim$Ng, HC1 = FALSE, small.strata = TRUE)
  expect_equal(round(result$tau.hat, 7), c(0.1321043, 0.8705829))
  expect_equal(round(result$se.rob, 7), c(0.1015669, 0.1052238))

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, G.id = data_sim$G.id, Ng = data_sim$Ng, HC1 = TRUE, small.strata = TRUE)
  expect_equal(round(result$tau.hat, 7), c(0.1321043, 0.8705829))
  expect_equal(round(result$se.rob, 7), c(0.1015669, 0.1052238))

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), G.id = data_sim$G.id, Ng = data_sim$Ng, HC1 = TRUE, small.strata = TRUE)
  expect_equal(round(result$tau.hat, 7), c(0.1399561, 0.8402642))
  expect_equal(round(result$se.rob, 7), c(0.0455854, 0.0522454))

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), G.id = data_sim$G.id, Ng = data_sim$Ng, HC1 = FALSE, small.strata = TRUE)
  expect_equal(round(result$tau.hat, 7), c(0.1399561, 0.8402642))
  expect_equal(round(result$se.rob, 7), c(0.0453538, 0.0519852))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), G.id = data_sim$G.id, Ng = NULL, HC1 = FALSE, small.strata = TRUE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1399561, 0.8402642))
  expect_equal(round(result$se.rob, 7), c(0.0453538, 0.0519852))

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2, data_sim$Ng), G.id = data_sim$G.id, Ng = data_sim$Ng, HC1 = TRUE, small.strata = TRUE)
  expect_equal(round(result$tau.hat, 7), c(0.1400713, 0.8409700))
  expect_equal(round(result$se.rob, 7), c(0.0456655, 0.0525556))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$Ng), G.id = data_sim$G.id, Ng = NULL, HC1 = TRUE, small.strata = TRUE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1254917, 0.8627286))
  expect_equal(round(result$se.rob, 7), c(0.1016849, 0.1038816))

  expect_error(
    invisible(capture.output({
      result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), G.id = NULL, Ng = data_sim$Ng, HC1 = TRUE, small.strata = TRUE)
    })),
    "Either all strata are large or too few strata qualify as 'small' to proceed with small.strata = TRUE. Please set small.strata = FALSE.",
    fixed = TRUE
  )

  expect_error(
    invisible(capture.output({
      result <- sreg(Y = data_sim$Y, D = data_sim$D, S = NULL, X = data.frame(data_sim$x_1, data_sim$x_2), G.id = data_sim$G.id, Ng = data_sim$Ng, HC1 = TRUE, small.strata = TRUE)
    })),
    "Strata indicator variable has not been provided (S = NULL), but small.strata = TRUE. This estimator requires stratification.",
    fixed = TRUE
  )
  expect_error(
    invisible(capture.output({
      result <- sreg(Y = data_sim$Y, D = NULL, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), G.id = data_sim$G.id, Ng = data_sim$Ng, HC1 = TRUE, small.strata = TRUE)
    })),
    "Treatments have not been provided (D = NULL). Please provide the vector of treatments.",
    fixed = TRUE
  )
})

test_that("data: big strata, option: big strata", {
  set.seed(123)
  tau.vec <- c(0.2, 0.9, 1.5)
  n.treat <- length(tau.vec)
  n_1 <- 1000

  data_sim <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 4, cluster = TRUE, small.strata = FALSE, treat.sizes = c(1, 1, 1), k = 3)

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, HC1 = TRUE, G.id = data_sim$G.id, Ng = data_sim$Ng, small.strata = FALSE)
  expect_equal(round(result$tau.hat, 7), c(0.1366736, 1.0663200, 1.5402160))
  expect_equal(round(result$se.rob, 7), c(0.1076841, 0.1073968, 0.1061023))

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, HC1 = FALSE, G.id = data_sim$G.id, Ng = data_sim$Ng, small.strata = FALSE)
  expect_equal(round(result$tau.hat, 7), c(0.1366736, 1.0663200, 1.5402160))
  expect_equal(round(result$se.rob, 7), c(0.1068247, 0.1065375, 0.1052502))

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), HC1 = TRUE, G.id = data_sim$G.id, Ng = data_sim$Ng, small.strata = FALSE)
  expect_equal(round(result$tau.hat, 7), c(0.0695722, 0.8947400, 1.3867602))
  expect_equal(round(result$se.rob, 7), c(0.0525225, 0.0545726, 0.0594523))

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), HC1 = FALSE, G.id = data_sim$G.id, Ng = data_sim$Ng, small.strata = FALSE)
  expect_equal(round(result$tau.hat, 7), c(0.0695722, 0.8947400, 1.3867602))
  expect_equal(round(result$se.rob, 7), c(0.0521126, 0.0541451, 0.0589783))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), HC1 = FALSE, G.id = data_sim$G.id, Ng = NULL, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.0695722, 0.8947400, 1.3867602))
  expect_equal(round(result$se.rob, 7), c(0.0521126, 0.0541451, 0.0589783))

  set.seed(123)
  tau.vec <- c(0.2, 0.9, 1.5)
  n.treat <- length(tau.vec)
  n_1 <- 100

  data_sim <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 4, cluster = TRUE, small.strata = FALSE, treat.sizes = c(1, 1, 1), k = 3)

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), HC1 = FALSE, G.id = NULL, Ng = data_sim$Ng, small.strata = FALSE)

  expect_equal(round(result$tau.hat, 7), c(0.2199062, 0.6091449, 1.3856801))
  expect_equal(round(result$se.rob, 7), c(0.0682736, 0.0689731, 0.0652965))

  set.seed(123)
  tau.vec <- c(0.2, 0.9, 1.5)
  n.treat <- length(tau.vec)
  n_1 <- 1000
  data_sim <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 4, cluster = TRUE, small.strata = FALSE, treat.sizes = c(1, 1, 1), k = 3)

  expect_warning(
    invisible(capture.output({
      result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), HC1 = FALSE, G.id = data_sim$G.id, Ng = NULL, small.strata = FALSE)
    })),
    "Cluster sizes have not been provided (Ng = NULL).",
    fixed = TRUE
  )

  result <- sreg(Y = data_sim$Y, D = data_sim$D, S = NULL, X = data.frame(data_sim$x_1, data_sim$x_2), HC1 = FALSE, G.id = data_sim$G.id, Ng = data_sim$Ng, small.strata = FALSE)
  expect_equal(round(result$tau.hat, 7), c(0.0751089, 0.9047198, 1.4042165))
  expect_equal(round(result$se.rob, 7), c(0.0624560, 0.0635305, 0.0675028))

  expect_error(
    invisible(capture.output({
      result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, X = data.frame(data_sim$x_1, data_sim$x_2), HC1 = FALSE, G.id = data_sim$G.id, Ng = data_sim$Ng, small.strata = "abcd")
    })),
    "The value of small.strata must be either TRUE or FALSE. A non-boolean value was provided.",
    fixed = TRUE
  )
})

test_that("data: small strata, option: big strata", {
  set.seed(123)
  tau.vec <- c(0.2, 0.8)
  n.treat <- length(tau.vec)
  n_1 <- 900

  data_sim <- sreg.rgen(n = n_1, tau.vec = tau.vec, n.strata = 4, cluster = TRUE, small.strata = TRUE, treat.sizes = c(1, 1, 1), k = 3)

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, G.id = data_sim$G.id, Ng = data_sim$Ng, HC1 = FALSE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1399429, 0.8776288))
  expect_equal(round(result$se.rob, 7), c(0.0571977, 0.0592139))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, G.id = data_sim$G.id, Ng = data_sim$Ng, HC1 = TRUE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1399429, 0.8776288))
  expect_equal(round(result$se.rob, 7), c(0.0571977, 0.0592139))

  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, G.id = data_sim$G.id, Ng = data_sim$Ng, HC1 = TRUE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("All strata have the same small number of clusters (e.g., matched pairs at the cluster level), but small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("HC1 adjustment unstable or undefined due to degenerate strata-treatment structure; reverting to unadjusted estimator.", warnings, fixed = TRUE)))
      expect_true(any(grepl("At least 25% of strata are small, but small.strata = FALSE.", warnings, fixed = TRUE)))
    },
    regexp = NA
  )

  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, G.id = data_sim$G.id, X = data.frame(data_sim$x_1, data_sim$x_2), Ng = data_sim$Ng, HC1 = TRUE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("All strata have the same small number of clusters (e.g., matched pairs at the cluster level), but small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("HC1 adjustment unstable or undefined due to degenerate strata-treatment structure; reverting to unadjusted estimator.", warnings, fixed = TRUE)))
      expect_true(any(grepl("At least 25% of strata are small, but small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE.", warnings, fixed = TRUE)))
    },
    regexp = NA
  )

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, G.id = data_sim$G.id, Ng = data_sim$Ng, X = data.frame(data_sim$x_1, data_sim$x_2), HC1 = TRUE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1399429, 0.8776288))
  expect_equal(round(result$se.rob, 7), c(0.0571977, 0.0592139))

  invisible(
    suppressWarnings(
      capture.output({
        result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, G.id = data_sim$G.id, Ng = NULL, X = data.frame(data_sim$x_1, data_sim$x_2, data_sim$Ng), HC1 = TRUE, small.strata = FALSE)
      })
    )
  )
  expect_equal(round(result$tau.hat, 7), c(0.1399429, 0.8776288))
  expect_equal(round(result$se.rob, 7), c(0.0571977, 0.0592139))

  expect_warning(
    {
      warnings <- character()

      withCallingHandlers(
        {
          invisible(capture.output({
            result <- sreg(Y = data_sim$Y, D = data_sim$D, S = data_sim$S, G.id = data_sim$G.id, Ng = NULL, X = data.frame(data_sim$x_1, data_sim$x_2, data_sim$Ng), HC1 = TRUE, small.strata = FALSE)
          }))
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      expect_true(any(grepl("All strata have the same small number of clusters (e.g., matched pairs at the cluster level), but small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("HC1 adjustment unstable or undefined due to degenerate strata-treatment structure; reverting to unadjusted estimator.", warnings, fixed = TRUE)))
      expect_true(any(grepl("At least 25% of strata are small, but small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("One or more covariates do not vary within one or more stratum-treatment combinations while small.strata = FALSE.", warnings, fixed = TRUE)))
      expect_true(any(grepl("Cluster sizes have not been provided (Ng = NULL).", warnings, fixed = TRUE)))

    },
    regexp = NA
  )
})
