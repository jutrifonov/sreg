mixed_adjustment_data <- function(cluster = FALSE) {
  set.seed(if (cluster) 9202 else 9201)
  k <- 3L
  n_small_strata <- 30L
  G_small <- n_small_strata * k
  G_big <- 6L
  G <- G_small + G_big

  S_cluster <- c(
    rep(seq_len(n_small_strata), each = k),
    rep(n_small_strata + 1L, G_big)
  )
  D_cluster <- c(
    rep(0:2, n_small_strata),
    rep(0:2, each = 2L)
  )
  x1 <- rnorm(G)
  x2 <- rnorm(G)
  Y_cluster <- 1 + 0.5 * x1 - 0.3 * x2 + 0.4 * D_cluster + rnorm(G)

  if (!cluster) {
    return(data.frame(
      Y = Y_cluster, S = S_cluster, D = D_cluster,
      x1 = x1, x2 = x2
    ))
  }

  members <- 3L
  rows <- rep(seq_len(G), each = members)
  data.frame(
    Y = Y_cluster[rows] + rnorm(length(rows), sd = 0.2),
    S = S_cluster[rows], D = D_cluster[rows], G.id = rows,
    Ng = rep(members, length(rows)), x1 = x1[rows], x2 = x2[rows]
  )
}

test_that("individual mixed adjustment is applied to both components", {
  set.seed(9203)
  sim <- sreg.rgen(
    n = 600, tau.vec = c(0.4, 0.9), n.strata = 4,
    cluster = FALSE, mixed.strata = TRUE, n.small = 360,
    k = 3, treat.sizes = c(1, 1, 1)
  )

  fit <- suppressWarnings(sreg(
    Y = sim$Y, S = sim$S, D = sim$D,
    X = sim[c("x_1", "x_2")], small.strata = TRUE
  ))

  expect_false(is.null(fit$res.small$lin.adj))
  expect_false(is.null(fit$res.big$lin.adj))
  expect_equal(ncol(fit$res.big$lin.adj), 2L)
})

test_that("individual mixed adjustment reports unidentified large regressions", {
  dat <- mixed_adjustment_data(cluster = FALSE)

  expect_error(
    suppressWarnings(sreg(
      Y = dat$Y, S = dat$S, D = dat$D, X = dat[c("x1", "x2")],
      small.strata = TRUE, k = 3
    )),
    "The large-strata component of the mixed design cannot support the requested covariate adjustment",
    fixed = TRUE
  )
})

test_that("cluster mixed adjustment reports unidentified large regressions", {
  dat <- mixed_adjustment_data(cluster = TRUE)

  expect_error(
    suppressWarnings(sreg(
      Y = dat$Y, S = dat$S, D = dat$D, G.id = dat$G.id, Ng = dat$Ng,
      X = dat[c("x1", "x2")], small.strata = TRUE, k = 3
    )),
    "The large-strata component of the mixed design cannot support the requested covariate adjustment",
    fixed = TRUE
  )
})

test_that("unadjusted mixed estimation remains available after adjustment failure", {
  individual <- mixed_adjustment_data(cluster = FALSE)
  cluster <- mixed_adjustment_data(cluster = TRUE)

  expect_silent(suppressWarnings(sreg(
    Y = individual$Y, S = individual$S, D = individual$D,
    X = NULL, small.strata = TRUE, k = 3
  )))
  expect_silent(suppressWarnings(sreg(
    Y = cluster$Y, S = cluster$S, D = cluster$D,
    G.id = cluster$G.id, Ng = cluster$Ng,
    X = NULL, small.strata = TRUE, k = 3
  )))
})
