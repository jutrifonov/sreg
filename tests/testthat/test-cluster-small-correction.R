test_that("small-strata cluster point estimator uses expanded outcomes and a common denominator", {
  set.seed(20260812)
  dat <- sreg.rgen(
    n = 400, tau.vec = c(0.5, 0.8), cluster = TRUE, is.cov = FALSE,
    small.strata = TRUE, k = 4, treat.sizes = c(2, 1, 1)
  )

  fit <- tau.hat.creg.ss(
    dat$Y, dat$D, X = NULL, dat$S, dat$G.id, dat$Ng
  )
  cluster_data <- .creg_ss_cluster_data(
    dat$Y, dat$S, dat$D, dat$G.id, dat$Ng
  )
  expected <- vapply(1:2, function(d) {
    (mean(cluster_data$T[cluster_data$D == d]) -
      mean(cluster_data$T[cluster_data$D == 0])) / mean(cluster_data$Ng)
  }, numeric(1))

  expect_equal(fit$tau.hat, expected)

  order <- sample.int(length(dat$Y))
  shuffled <- tau.hat.creg.ss(
    dat$Y[order], dat$D[order], X = NULL, dat$S[order],
    dat$G.id[order], dat$Ng[order]
  )
  expect_equal(shuffled$tau.hat, fit$tau.hat)
})

test_that("small-strata cluster inference supports multiple arms", {
  set.seed(20260813)
  dat <- sreg.rgen(
    n = 400, tau.vec = c(0.5, 0.8), cluster = TRUE, is.cov = TRUE,
    small.strata = TRUE, k = 4, treat.sizes = c(2, 1, 1)
  )
  X <- data.frame(x_1 = dat$x_1, x_2 = dat$x_2)

  adjusted <- sreg(
    Y = dat$Y, S = dat$S, D = dat$D, G.id = dat$G.id, Ng = dat$Ng,
    X = X, HC1 = TRUE, small.strata = TRUE
  )
  unadjusted <- sreg(
    Y = dat$Y, S = dat$S, D = dat$D, G.id = dat$G.id, Ng = dat$Ng,
    X = NULL, HC1 = TRUE, small.strata = TRUE
  )

  expect_length(adjusted$tau.hat, 2)
  expect_length(adjusted$se.rob, 2)
  expect_true(all(is.finite(adjusted$tau.hat)))
  expect_true(all(is.finite(adjusted$se.rob)))
  expect_true(all(adjusted$se.rob > 0))
  expect_true(all(is.finite(unadjusted$tau.hat)))
  expect_true(all(is.finite(unadjusted$se.rob)))
  expect_true(all(unadjusted$se.rob > 0))
})

test_that("binary small-strata variance is the multi-arm formula with two arms", {
  set.seed(20260814)
  dat <- sreg.rgen(
    n = 400, tau.vec = 0.8, cluster = TRUE, is.cov = TRUE,
    small.strata = TRUE, k = 2, treat.sizes = c(1, 1)
  )
  X <- data.frame(x_1 = dat$x_1, x_2 = dat$x_2)
  fit <- tau.hat.creg.ss(dat$Y, dat$D, X, dat$S, dat$G.id, dat$Ng)
  variance <- as.var.creg.ss(
    dat$Y, dat$D, X, dat$S, dat$G.id, dat$Ng, fit, HC1 = TRUE
  )

  expect_length(variance, 1)
  expect_true(is.finite(variance))
  expect_gt(variance, 0)
})

test_that("small-strata cluster adjustment works when cluster sizes are inferred", {
  set.seed(20260815)
  dat <- sreg.rgen(
    n = 400, tau.vec = c(0.5, 0.8), cluster = TRUE, is.cov = TRUE,
    small.strata = TRUE, k = 4, treat.sizes = c(2, 1, 1)
  )
  X <- data.frame(x_1 = dat$x_1, x_2 = dat$x_2)

  supplied <- sreg(
    dat$Y, dat$S, dat$D, dat$G.id, dat$Ng, X,
    small.strata = TRUE
  )
  inferred <- suppressWarnings(sreg(
    dat$Y, dat$S, dat$D, dat$G.id, Ng = NULL, X = X,
    small.strata = TRUE
  ))

  expect_equal(inferred$tau.hat, supplied$tau.hat)
  expect_equal(inferred$se.rob, supplied$se.rob)

  one_covariate <- suppressWarnings(sreg(
    dat$Y, dat$S, dat$D, dat$G.id, Ng = NULL,
    X = data.frame(cluster_size = dat$Ng), small.strata = TRUE
  ))
  expect_true(all(is.finite(one_covariate$tau.hat)))
  expect_true(all(is.finite(one_covariate$se.rob)))

  named_ng <- sreg(
    dat$Y, dat$S, dat$D, dat$G.id, dat$Ng,
    X = data.frame(x_1 = dat$x_1, Ng = dat$Ng), small.strata = TRUE
  )
  expect_true(all(is.finite(named_ng$tau.hat)))
  expect_true(all(is.finite(named_ng$se.rob)))
})
