test_that("is.cov controls covariate columns in large-strata cluster designs", {
  set.seed(20260822)
  without_covariates <- sreg.rgen(
    n = 20, n.strata = 2, tau.vec = 0.5,
    cluster = TRUE, is.cov = FALSE
  )

  set.seed(20260822)
  with_covariates <- sreg.rgen(
    n = 20, n.strata = 2, tau.vec = 0.5,
    cluster = TRUE, is.cov = TRUE
  )

  expect_false(any(c("x_1", "x_2") %in% names(without_covariates)))
  expect_true(all(c("x_1", "x_2") %in% names(with_covariates)))
  expect_true(all(c("Y", "S", "D", "G.id", "Ng") %in%
    names(without_covariates)))
})
