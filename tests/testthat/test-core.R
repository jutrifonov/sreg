test_that("simulations without clusters work", {
  set.seed(123)
  data <- sreg.rgen(n = 1000, tau.vec = c(0.2, 0.5), n.strata = 10, cluster = F, is.cov = TRUE)
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)

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
  S <- data$S
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
  S <- data$S
  D <- data$D

  X[10:12, ] <- NaN
  Y[1:10] <- NaN
  msg <- "ignoring these values"
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = NULL, Ng = NULL, X = X, HC1 = TRUE)
    })),
    msg
  )
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  Y <- data$Y
  S[12] <- NA
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = NULL, Ng = NULL, X = X, HC1 = TRUE)
    })),
    msg
  )
  S <- data$S
  S[19] <- NaN
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = NULL, Ng = NULL, X = X, HC1 = TRUE)
    })),
    msg
  )
  S <- data$S
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
  S <- data$S
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
  D <- data$D
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
  D <- data$D
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    }))
  )
})

test_that("simulations with clusters work", {
  set.seed(123) # fix the random seed
  data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50, is.cov = TRUE
  )
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng

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

  X <- data.frame("Ng" = data$Ng, "x_1" = data$x_1, "x_2" = data$x_2)

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
  S <- data$S
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
  S <- data$S
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(as.matrix(Y), S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    }))
  )
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(as.numeric(Y), S = as.numeric(S), D = as.integer(D), G.id = as.integer(G.id), Ng = as.integer(Ng), X = as.matrix(X))
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
      result <- sreg::sreg(as.matrix(Y), S = as.matrix(S), D = as.vector(D), G.id = as.matrix(G.id), Ng = as.matrix(Ng), X = as.matrix(X))
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
  S <- data$S
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
  S <- data$S
  D <- data$D
  Ng[2] <- 37.64
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "must contain only integer values."
  )
  Ng <- data$Ng
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
  data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50
  )
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng
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
  D <- data$D
  expect_silent(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    }))
  )
})

test_that("degrees of freedom error works", {
  set.seed(123)
  data <- sreg.rgen(n = 25, tau.vec = c(0.2, 0.5), n.strata = 3, cluster = F, is.cov = TRUE)
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)

  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S, D, G.id = NULL, Ng = NULL, X = X)
    })),
    "too many covariates relative to the number of observations"
  )

  set.seed(123)
  data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 10, cluster = T, Nmax = 50
  )
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng

  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    "too many covariates relative to the number of observations"
  )
})

test_that("individual level X warning works", {
  set.seed(123)
  data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50
  )
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng
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
  data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50
  )
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng

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

  data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.8),
    n.strata = 4, cluster = T, Nmax = 50
  )
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  Y[1:10] <- NA
  G.id <- data$G.id
  Ng <- data$Ng

  msg <- "ignoring these values"
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  Y <- data$Y
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
  Y <- data$Y
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id[100:105] <- NA
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  G.id <- data$G.id
  Ng[4:5] <- NaN
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  Ng <- data$Ng
  D[24:25] <- NA
  expect_warning(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    })),
    msg
  )
  D <- data$D
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
  data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.5),
    n.strata = 5, cluster = F, is.cov = TRUE
  )
  Y <- data$Y
  S <- data$S
  S[S == 4] <- 1
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = NULL, Ng = NULL, X = X)
    })),
    "There are skipped values in the range"
  )
  set.seed(123)
  data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.5),
    n.strata = 5, cluster = T, is.cov = TRUE
  )
  Y <- data$Y
  S <- data$S
  S[S == 3] <- 4
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "There are skipped values in the range"
  )
  S <- data$S
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
  data <- sreg.rgen(
    n = 100, tau.vec = c(0.2, 0.5),
    n.strata = 5, cluster = T, is.cov = TRUE
  )
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng
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
  S <- data$S
  D[13:20] <- 1
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "The values for S, D, and Ng must be consistent within each cluster"
  )
  D <- data$D
  Ng[300] <- 30
  expect_error(
    invisible(capture.output({
      result <- sreg::sreg(Y, S = S, D = D, G.id = G.id, Ng = Ng, X = X)
    })),
    "The values for S, D, and Ng must be consistent within each cluster"
  )
  Ng <- data$Ng
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
  data <- AEJapp
  Y <- data$gradesq34
  D <- data$treatment
  S <- data$class_level
  pills <- data$pills_taken
  age <- data$age_months
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
