test_that("sreg estimates individual-level mixed designs with 4-tuples", {
  set.seed(4101)
  sim <- sreg.rgen(
    n = 120, tau.vec = 0.5, cluster = FALSE,
    mixed.strata = TRUE, n.small = 80, k = 4,
    treat.sizes = c(2, 2), n.strata = 4
  )

  expect_warning(
    fit <- sreg(
      Y = sim$Y, S = sim$S, D = sim$D,
      small.strata = TRUE, k = 4
    ),
    "same size (k = 4), which is used as the small-stratum size",
    fixed = TRUE
  )
  expect_true(fit$mixed.design)
  expect_true(all(table(fit$res.small$data$S) == 4))
  expect_equal(nrow(fit$res.small$data), 80)
})

test_that("sreg estimates cluster-level mixed designs with 4-tuples", {
  set.seed(4102)
  sim <- sreg.rgen(
    n = 80, tau.vec = 0.5, cluster = TRUE,
    mixed.strata = TRUE, n.small = 48, k = 4,
    treat.sizes = c(2, 2), n.strata = 3
  )

  expect_warning(
    fit <- sreg(
      Y = sim$Y, S = sim$S, D = sim$D,
      G.id = sim$G.id, Ng = sim$Ng,
      small.strata = TRUE, k = 4
    ),
    "same size (k = 4), which is used as the small-stratum size",
    fixed = TRUE
  )
  small_cluster_strata <- unique(fit$res.small$data[c("S", "G.id")])
  expect_true(fit$mixed.design)
  expect_true(all(table(small_cluster_strata$S) == 4))
  expect_equal(length(unique(small_cluster_strata$G.id)), 48)
})

test_that("general k must be supplied when automatic detection cannot identify it", {
  set.seed(4103)
  sim <- sreg.rgen(
    n = 120, tau.vec = 0.5, cluster = FALSE,
    mixed.strata = TRUE, n.small = 80, k = 4,
    treat.sizes = c(2, 2), n.strata = 4
  )

  expect_error(
    sreg(Y = sim$Y, S = sim$S, D = sim$D, small.strata = TRUE),
    "too few strata qualify as 'small'",
    fixed = TRUE
  )
  expect_error(
    sreg(Y = sim$Y, S = sim$S, D = sim$D, small.strata = TRUE, k = 5),
    "too few strata qualify as 'small'",
    fixed = TRUE
  )
})

test_that("explicit k works for uniform k-tuple designs and is validated", {
  set.seed(4104)
  sim <- sreg.rgen(
    n = 80, tau.vec = 0.5, cluster = FALSE,
    small.strata = TRUE, k = 4, treat.sizes = c(2, 2)
  )

  expect_silent(
    fit <- sreg(
      Y = sim$Y, S = sim$S, D = sim$D,
      small.strata = TRUE, k = 4
    )
  )
  expect_false(isTRUE(fit$mixed.design))
  expect_error(
    sreg(Y = sim$Y, S = sim$S, D = sim$D, small.strata = TRUE, k = 3),
    "does not match the observed stratum size",
    fixed = TRUE
  )
  expect_error(
    sreg(Y = sim$Y, S = sim$S, D = sim$D, small.strata = TRUE, k = 0),
    "k must be NULL or a positive integer",
    fixed = TRUE
  )
})
