test_that("large-strata custom DGP defaults preserve generated data", {
  set.seed(2026)
  baseline <- sreg.rgen(
    n = 300, n.strata = 5, tau.vec = c(0.5, 0.8),
    cluster = FALSE, small.strata = FALSE
  )
  set.seed(2026)
  explicit_null <- sreg.rgen(
    n = 300, n.strata = 5, tau.vec = c(0.5, 0.8),
    cluster = FALSE, small.strata = FALSE,
    allocation.probs = NULL, stratum.effects = NULL,
    treatment.effects.by.stratum = NULL
  )
  expect_identical(baseline, explicit_null)
})

test_that("large-strata generator accepts stratum-specific allocations", {
  allocation <- rbind(
    c(0.60, 0.20, 0.20),
    c(0.50, 0.30, 0.20),
    c(0.34, 0.33, 0.33),
    c(0.20, 0.30, 0.50),
    c(0.20, 0.20, 0.60)
  )
  set.seed(30)
  generated <- sreg.rgen(
    n = 3000, n.strata = 5, tau.vec = c(0.5, 0.8),
    cluster = FALSE, small.strata = FALSE,
    allocation.probs = allocation
  )

  counts <- table(generated$S, generated$D)
  for (stratum in seq_len(5)) {
    stratum_n <- sum(counts[stratum, ])
    expect_equal(
      as.integer(counts[stratum, c("1", "2")]),
      floor(allocation[stratum, c(2, 3)] * stratum_n)
    )
  }
})

test_that("large-strata generator applies custom outcome effects", {
  stratum_effects <- c(-2, -1, 0, 1, 2)
  treatment_effects <- cbind(
    c(-1.5, -0.5, 0.5, 1.5, 2.5),
    c(2.4, 1.6, 0.8, 0, -0.8)
  )

  set.seed(91)
  baseline <- sreg.rgen(
    n = 600, n.strata = 5, tau.vec = c(0.5, 0.8),
    cluster = FALSE, small.strata = FALSE
  )
  set.seed(91)
  customized <- sreg.rgen(
    n = 600, n.strata = 5, tau.vec = c(0.5, 0.8),
    cluster = FALSE, small.strata = FALSE,
    stratum.effects = stratum_effects,
    treatment.effects.by.stratum = treatment_effects
  )

  expect_identical(customized$S, baseline$S)
  expect_identical(customized$D, baseline$D)
  expected_shift <- stratum_effects[baseline$S]
  treated <- baseline$D > 0
  expected_shift[treated] <- expected_shift[treated] +
    treatment_effects[cbind(baseline$S[treated], baseline$D[treated])] -
    c(0.5, 0.8)[baseline$D[treated]]
  expect_equal(customized$Y - baseline$Y, expected_shift, tolerance = 1e-12)
})

test_that("custom large-strata arguments are validated", {
  expect_error(
    sreg.rgen(
      n = 100, n.strata = 5, tau.vec = c(0.5, 0.8), cluster = FALSE,
      allocation.probs = matrix(1 / 3, nrow = 4, ncol = 3)
    ),
    "allocation.probs"
  )
  expect_error(
    sreg.rgen(
      n = 100, n.strata = 5, tau.vec = c(0.5, 0.8), cluster = FALSE,
      stratum.effects = 1:4
    ),
    "stratum.effects"
  )
  expect_error(
    sreg.rgen(
      n = 100, n.strata = 5, tau.vec = c(0.5, 0.8), cluster = TRUE,
      stratum.effects = rep(0, 5)
    ),
    "currently supported only"
  )
})
