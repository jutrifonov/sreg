mixed_cluster_fixture <- function() {
  set.seed(8126)
  sreg.rgen(
    n = 600, tau.vec = c(0.4, 0.9), n.strata = 4,
    cluster = TRUE, mixed.strata = TRUE, n.small = 360,
    k = 3, treat.sizes = c(1, 1, 1)
  )
}

test_that("mixed cluster weights use supplied cluster population sizes", {
  sim <- mixed_cluster_fixture()

  cluster_strata <- unique(sim[c("G.id", "S")])
  stratum_sizes <- table(cluster_strata$S)
  small_strata <- names(stratum_sizes)[stratum_sizes == 3]
  is_small <- sim$S %in% as.numeric(small_strata)

  # Make represented population sizes differ from the numbers sampled.
  Ng_population <- sim$Ng
  Ng_population[is_small] <- 2L * Ng_population[is_small]

  fit <- suppressWarnings(sreg(
    Y = sim$Y, S = sim$S, D = sim$D,
    G.id = sim$G.id, Ng = Ng_population,
    small.strata = TRUE
  ))

  cluster_sizes <- unique(fit$data[c("G.id", "stratum_type", "Ng")])
  N_small <- sum(cluster_sizes$Ng[cluster_sizes$stratum_type == "small"])
  N_big <- sum(cluster_sizes$Ng[cluster_sizes$stratum_type == "big"])
  expected <- (N_small * fit$res.small$tau.hat +
    N_big * fit$res.big$tau.hat) / (N_small + N_big)

  expect_equal(fit$tau.hat, expected)
  expect_false(isTRUE(all.equal(N_small, nrow(fit$res.small$data))))
})

test_that("mixed cluster weights infer cluster sizes from available observations", {
  sim <- mixed_cluster_fixture()

  fit <- suppressWarnings(sreg(
    Y = sim$Y, S = sim$S, D = sim$D,
    G.id = sim$G.id, Ng = NULL,
    small.strata = TRUE
  ))

  cluster_sizes <- unique(fit$data[c("G.id", "stratum_type", "Ng")])
  N_small <- sum(cluster_sizes$Ng[cluster_sizes$stratum_type == "small"])
  N_big <- sum(cluster_sizes$Ng[cluster_sizes$stratum_type == "big"])
  expected <- (N_small * fit$res.small$tau.hat +
    N_big * fit$res.big$tau.hat) / (N_small + N_big)

  expect_equal(fit$tau.hat, expected)
  expect_equal(N_small, nrow(fit$res.small$data))
  expect_equal(N_big, nrow(fit$res.big$data))
})

test_that("mixed cluster adjustment is applied to both components", {
  sim <- mixed_cluster_fixture()

  fit <- suppressWarnings(sreg(
    Y = sim$Y, S = sim$S, D = sim$D,
    G.id = sim$G.id, Ng = sim$Ng,
    X = sim[c("x_1", "x_2")], small.strata = TRUE
  ))

  expect_false(is.null(fit$res.small$lin.adj))
  expect_false(is.null(fit$res.big$lin.adj))
  expect_equal(ncol(fit$res.small$lin.adj), 2L)
  expect_equal(ncol(fit$res.big$lin.adj), 2L)
})

test_that("mixed cluster variance uses cluster-level component-share variation", {
  sim <- mixed_cluster_fixture()
  Ng_population <- sim$Ng + (sim$G.id %% 5L)

  fit <- suppressWarnings(sreg(
    Y = sim$Y, S = sim$S, D = sim$D,
    G.id = sim$G.id, Ng = Ng_population,
    small.strata = TRUE, HC1 = FALSE
  ))

  clusters <- unique(fit$data[c("G.id", "stratum_type", "Ng")])
  G_total <- nrow(clusters)
  N_small <- sum(clusters$Ng[clusters$stratum_type == "small"])
  N_big <- sum(clusters$Ng[clusters$stratum_type == "big"])
  p_small <- N_small / (N_small + N_big)
  p_big <- N_big / (N_small + N_big)
  N_bar <- (N_small + N_big) / G_total
  V_p <- mean(
    clusters$Ng^2 *
      (as.numeric(clusters$stratum_type == "big") - p_big)^2
  ) / N_bar^2

  expected_variance <- p_small^2 * fit$res.small$se.rob^2 +
    p_big^2 * fit$res.big$se.rob^2 +
    (V_p / G_total) *
      (fit$res.small$tau.hat - fit$res.big$tau.hat)^2

  expect_equal(fit$se.rob^2, expected_variance)
})
