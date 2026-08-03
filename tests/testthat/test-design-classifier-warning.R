test_that("mixed-design warning reports the detected individual-level k", {
  make_data <- function(sizes) {
    data.frame(
      Y = seq_len(sum(sizes)),
      D = rep(0:1, length.out = sum(sizes)),
      S = rep(seq_along(sizes), sizes)
    )
  }

  expect_warning(
    classified <- design.classifier(
      make_data(c(2, 4, 5, 6)), S = S, small.strata = TRUE
    ),
    paste0(
      "Mixed design detected: at least 25% of all strata have the same size ",
      "(k = 2), which is used as the small-stratum size. ",
      "Weighted estimators will be used."
    ),
    fixed = TRUE
  )
  expect_equal(unique(classified$stratum_type[classified$S == 1]), "small")
  expect_true(all(classified$stratum_type[classified$S != 1] == "big"))
})

test_that("mixed-design warning detects k from cluster counts", {
  cluster_sizes <- c(3, 3, 3, 3, 8)
  cluster_data <- data.frame(
    Y = seq_len(sum(cluster_sizes)),
    D = rep(0:1, length.out = sum(cluster_sizes)),
    S = rep(seq_along(cluster_sizes), cluster_sizes),
    G.id = seq_len(sum(cluster_sizes))
  )

  expect_warning(
    classified <- design.classifier(
      cluster_data, S = S, G.id = G.id, small.strata = TRUE
    ),
    "same size (k = 3), which is used as the small-stratum size",
    fixed = TRUE
  )
  expect_equal(sort(unique(classified$S[classified$stratum_type == "small"])), 1:4)
  expect_equal(unique(classified$S[classified$stratum_type == "big"]), 5)
})

test_that("fewer than 25 percent at one size does not produce a mixed warning", {
  sizes <- c(2, 4, 5, 6, 7)
  data <- data.frame(
    Y = seq_len(sum(sizes)),
    D = rep(0:1, length.out = sum(sizes)),
    S = rep(seq_along(sizes), sizes)
  )

  expect_error(
    design.classifier(data, S = S, small.strata = TRUE),
    "too few strata qualify as 'small'",
    fixed = TRUE
  )
})
