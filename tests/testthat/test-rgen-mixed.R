test_that("sreg.rgen generates mixed individual-level designs", {
  set.seed(2026)
  sim <- sreg.rgen(
    n = 120, tau.vec = c(0.2, 0.8), n.strata = 4,
    cluster = FALSE, mixed.strata = TRUE, n.small = 90,
    k = 3, treat.sizes = c(1, 1, 1)
  )

  sizes <- as.integer(table(sim$S))
  expect_equal(nrow(sim), 120)
  expect_true(any(sizes == 3))
  expect_true(any(sizes > 3))
  expect_equal(sum(sizes == 3), 30)
  expect_equal(sort(unique(sim$D)), 0:2)

  fit <- suppressWarnings(sreg(
    Y = sim$Y, S = sim$S, D = sim$D,
    X = NULL, small.strata = TRUE
  ))
  expect_true(fit$mixed.design)
})

test_that("sreg.rgen generates mixed cluster-level designs", {
  set.seed(2026)
  sim <- sreg.rgen(
    n = 120, tau.vec = c(0.2, 0.8), n.strata = 4,
    cluster = TRUE, mixed.strata = TRUE, n.small = 90,
    k = 3, treat.sizes = c(1, 1, 1)
  )

  cluster.data <- unique(sim[c("S", "D", "G.id")])
  sizes <- as.integer(table(cluster.data$S))
  expect_equal(length(unique(sim$G.id)), 120)
  expect_equal(nrow(cluster.data), 120)
  expect_false(anyDuplicated(sim$G.id[match(unique(sim$G.id), sim$G.id)]) > 0)
  expect_true(any(sizes == 3))
  expect_true(any(sizes > 3))
  expect_equal(sum(sizes == 3), 30)

  fit <- suppressWarnings(sreg(
    Y = sim$Y, S = sim$S, D = sim$D, G.id = sim$G.id, Ng = sim$Ng,
    X = NULL, small.strata = TRUE
  ))
  expect_true(fit$mixed.design)
})

test_that("mixed sreg.rgen validates its component sizes", {
  expect_error(
    sreg.rgen(n = 120, mixed.strata = TRUE, n.small = 91),
    "n.small must be divisible by k"
  )
  expect_error(
    sreg.rgen(n = 30, n.strata = 10, tau.vec = c(0.2, 0.8),
              mixed.strata = TRUE, n.small = 15),
    "large-strata component"
  )
  expect_error(
    sreg.rgen(n = 120, tau.vec = c(0.2, 0.8), mixed.strata = TRUE,
              n.small = 90, treat.sizes = c(1, 1)),
    "treat.sizes must be a nonnegative integer vector"
  )
})

test_that("existing sreg.rgen calls retain their behavior", {
  set.seed(11)
  old <- sreg.rgen(n = 30, cluster = FALSE, small.strata = FALSE)
  set.seed(11)
  explicit <- sreg.rgen(
    n = 30, cluster = FALSE, small.strata = FALSE,
    mixed.strata = FALSE
  )
  expect_identical(old, explicit)
})

test_that("mixed sreg.rgen derives an allocation when treat.sizes is omitted", {
  set.seed(7)
  sim <- sreg.rgen(
    n = 60, n.strata = 3, cluster = FALSE,
    mixed.strata = TRUE, n.small = 42, k = 3
  )
  small <- sim[as.integer(table(sim$S))[match(sim$S, as.integer(names(table(sim$S))))] == 3, ]
  allocations <- table(small$S, small$D)
  expect_true(all(rowSums(allocations) == 3))
  expect_true(all(allocations[, "0"] == 2))
  expect_true(all(allocations[, "1"] == 1))
})
