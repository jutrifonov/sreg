cluster_order_fixture <- function(mixed) {
  strata <- c(rep(1:2, each = 18), if (mixed) rep(3:14, each = 3))
  g <- seq_along(strata)
  sizes <- 2 + g %% 4
  clusters <- data.frame(
    G.id = g, S = strata, D = rep(0:2, length.out = length(g)),
    Ng = 2 * sizes + g %% 3, x1 = sin(g * 0.7), x2 = cos(g * 0.3)
  )
  d <- clusters[rep(g, sizes), ]
  d$Y <- 0.3 * d$D + 0.2 * d$S + 0.4 * d$x1 +
    sin(d$G.id * 1.1) + 0.05 * cos(seq_len(nrow(d)))
  rownames(d) <- NULL
  d
}

for (mixed in c(FALSE, TRUE)) {
  local({
    use_mixed <- mixed
    test_that(paste(if (use_mixed) "mixed" else "large", "cluster inference is invariant to observation order"), {
      d <- cluster_order_fixture(use_mixed)
      set.seed(90217)
      orders <- list(rev(seq_len(nrow(d))), sample.int(nrow(d)))
      fields <- c("tau.hat", "se.rob", "t.stat", "p.value", "CI.left", "CI.right")
      for (adjust in c(FALSE, TRUE)) for (hc1 in c(FALSE, TRUE)) {
        for (supplied_sizes in c(FALSE, TRUE)) {
          fit <- function(z) suppressWarnings(sreg(
            Y = z$Y, S = z$S, D = z$D, G.id = z$G.id,
            Ng = if (supplied_sizes) z$Ng else NULL,
            X = if (adjust) z[c("x1", "x2")] else NULL,
            HC1 = hc1, small.strata = use_mixed,
            k = if (use_mixed) 3 else NULL
          ))
          reference <- fit(d)
          for (idx in orders) {
            reordered <- fit(d[idx, ])
            for (field in fields) expect_equal(
              reordered[[field]], reference[[field]], tolerance = 1e-10,
              info = paste(field, "adjust", adjust, "HC1", hc1, "Ng", supplied_sizes)
            )
          }
        }
      }
    })
  })
}
