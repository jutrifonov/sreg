test_that("cluster large-strata variance includes clusters in other treatment arms", {
  S <- rep(1:2, each = 9)
  D <- rep(rep(0:2, each = 3), 2)
  A <- ifelse(D == 1, 1, ifelse(D == 0, 0, -999999))
  Ng <- c(2, 4, 7, 3, 6, 8, 5, 9, 11,
    3, 5, 8, 4, 7, 10, 6, 9, 13)
  T <- 1.5 + 0.4 * Ng + 0.8 * (D == 1) + 1.3 * (D == 2) +
    0.5 * (S == 2)
  Y.bar.g <- T / Ng
  tau <- 0.35
  pi <- rep(1 / 3, length(D))
  mu.0 <- 0.2 + 0.15 * Ng + 0.1 * (S == 2)
  mu.d <- 0.6 + 0.25 * Ng + 0.2 * (S == 2)

  data.list <- list(data.frame(
    S = S, D = D, A = A, I = as.numeric(D %in% c(0, 1)), Ng = Ng
  ))
  fit <- list(
    tau.hat = tau,
    Y.bar.g = Y.bar.g,
    Ng = Ng,
    mu.hat = list(cbind(mu.0, mu.d)),
    pi.hat = list(pi),
    pi.hat.0 = pi,
    data.list = data.list
  )

  manual.se <- function(mu0, mud) {
    M <- mud - mu0
    Xi.tilde.d <- M + (T - mud) / pi
    Xi.tilde.0 <- M - (T - mu0) / pi
    N.bar <- ave(Ng, S, FUN = mean)
    correction <- tau * (Ng - N.bar)

    d.mean <- ave(Xi.tilde.d[D == 1], S[D == 1], FUN = mean)
    d.mean <- setNames(d.mean[!duplicated(S[D == 1])], unique(S[D == 1]))
    zero.mean <- ave(Xi.tilde.0[D == 0], S[D == 0], FUN = mean)
    zero.mean <- setNames(zero.mean[!duplicated(S[D == 0])], unique(S[D == 0]))

    Xi.d <- Xi.tilde.d - d.mean[as.character(S)] - correction
    Xi.0 <- Xi.tilde.0 - zero.mean[as.character(S)] - correction
    Xi.other <- M - ave(M, S, FUN = mean) - correction

    T.d <- ave(T[D == 1], S[D == 1], FUN = mean)
    T.d <- setNames(T.d[!duplicated(S[D == 1])], unique(S[D == 1]))
    T.0 <- ave(T[D == 0], S[D == 0], FUN = mean)
    T.0 <- setNames(T.0[!duplicated(S[D == 0])], unique(S[D == 0]))
    Xi.2 <- T.d[as.character(S)] - T.0[as.character(S)] - tau * N.bar

    within <- mean(
      as.numeric(D == 1) * Xi.d^2 +
        as.numeric(D == 0) * Xi.0^2 +
        as.numeric(!D %in% c(0, 1)) * Xi.other^2
    )
    variance <- (within + mean(Xi.2^2)) / mean(Ng)^2

    legacy.within <- mean(
      as.numeric(D == 1) * Xi.d^2 + as.numeric(D == 0) * Xi.0^2
    )
    legacy.variance <- (legacy.within + mean(Xi.2^2)) / mean(Ng)^2

    c(corrected = sqrt(variance / length(D)),
      legacy = sqrt(legacy.variance / length(D)))
  }

  adjusted <- manual.se(mu.0, mu.d)
  expect_equal(
    sreg:::as.var.creg(model = TRUE, fit = fit, HC1 = FALSE),
    unname(adjusted["corrected"])
  )
  expect_gt(adjusted["corrected"], adjusted["legacy"])

  unadjusted <- manual.se(rep(0, length(D)), rep(0, length(D)))
  expect_equal(
    sreg:::as.var.creg(model = NULL, fit = fit, HC1 = FALSE),
    unname(unadjusted["corrected"])
  )
  expect_gt(unadjusted["corrected"], unadjusted["legacy"])
})
