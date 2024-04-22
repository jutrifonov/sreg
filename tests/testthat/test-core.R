test_that("simulations without clusters work", {
  set.seed(123) # fix the random seed
  data <- sreg.rgen(n = 1000, tau.vec = c(0.2, 0.5), n.strata = 10, cluster = F, is.cov = TRUE)
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)

  invisible(capture.output({result <- sreg::sreg(Y, S, D, G.id = NULL, Ng = NULL, X = X)}))

  expect_equal(round(result$tau, 7), c(0.1580814, 0.4846882))
  expect_equal(round(result$se, 8), c(0.07524021, 0.07616346))

  invisible(capture.output({result <- sreg::sreg(Y, S, D, G.id = NULL, Ng = NULL, X = NULL)}))

  expect_equal(round(result$tau, 7), c(0.1627114, 0.4948722))
  expect_equal(round(result$se, 7), c(0.1105611, 0.1124072))

  invisible(capture.output({result <- sreg::sreg(Y, S = NULL, D, G.id = NULL, Ng = NULL, X = X)}))

  expect_equal(round(result$tau, 7), c(0.1578917, 0.4963735))
  expect_equal(round(result$se, 8), c(0.08255663, 0.08320655))

  invisible(capture.output({result <- sreg::sreg(Y, S = NULL, D = D, G.id = NULL, Ng = NULL, X = NULL)}))

  expect_equal(round(result$tau, 7), c(0.1685108, 0.5022035))
  expect_equal(round(result$se, 7), c(0.1145915, 0.1161482))
})

test_that("simulations with clusters work", {
  set.seed(123) # fix the random seed
  data <- sreg.rgen(
  n = 100, tau.vec = c(0.2, 0.8),
  n.strata = 4, cluster = T, Nmax = 50
  )
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng

  invisible(capture.output({result <- sreg:: sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)}))

  expect_equal(round(result$tau.hat, 7), c(0.3350243, 0.8488310))
  expect_equal(round(result$se.rob, 7), c(0.1529870, 0.2183116))

  invisible(capture.output({result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = NULL, HC1 = TRUE)}))

  expect_equal(round(result$tau.hat, 7), c(0.3003475, 1.0875317))
  expect_equal(round(result$se.rob, 7), c(0.3126555, 0.3266127))

  invisible(capture.output({result <- sreg::sreg(Y, S = NULL, D, G.id, Ng, X = X, HC1 = TRUE)}))

  expect_equal(round(result$tau.hat, 7), c(0.1687191, 0.8792778))
  expect_equal(round(result$se.rob, 7), c(0.1641397, 0.1636365))

  invisible(capture.output({result <- sreg::sreg(Y, S = NULL, D, G.id, Ng, X = NULL, HC1 = TRUE)}))

  expect_equal(round(result$tau.hat, 7), c(0.2613412, 1.0528130))
  expect_equal(round(result$se.rob, 7), c(0.3047530, 0.3084314))

  invisible(capture.output({result <- sreg::sreg(Y, S = S, D, G.id = NULL, Ng = NULL, X = X)}))

  expect_equal(round(result$tau.hat, 7), c(0.1548014, 0.6955543))
  expect_equal(round(result$se.rob, 8), c(0.06201006, 0.05871032))

  invisible(capture.output({result <- sreg::sreg(Y, S = S, D, G.id = NULL, Ng = NULL, X = X)}))

  expect_error(invisible(capture.output({result <- sreg::sreg(Y = NULL, S = S, D, G.id = G.id, Ng = NULL, X = X)})), 
                  "Observed outcomes have not been provided")
  expect_error(invisible(capture.output({result <- sreg::sreg(Y, S = S, D = NULL, G.id = NULL, Ng = Ng, X = X)})), 
                  "Treatments have not been provided")

})

test_that("degrees of freedom warning works", {
  set.seed(123)
  data <- sreg.rgen(n = 25, tau.vec = c(0.2, 0.5), n.strata = 3, cluster = F, is.cov = TRUE)
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)

  expect_warning(invisible(capture.output({result <- sreg::sreg(Y, S, D, G.id = NULL, Ng = NULL, X = X)})), 
                            "too many covariates relative to the number of observations")

  set.seed(123)
  data <- sreg.rgen(
  n = 100, tau.vec = c(0.2, 0.8),
  n.strata = 10, cluster = T, Nmax = 50)
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng   

  expect_warning(invisible(capture.output({result <- sreg:: sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)})), 
                            "too many covariates relative to the number of observations")
                  
})

test_that("individual level X warning works", {
  set.seed(123)
  data <- sreg.rgen(
  n = 100, tau.vec = c(0.2, 0.8),
  n.strata = 4, cluster = T, Nmax = 50)
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng
  X[1, 1] <- 2.3894

  expect_warning(invisible(capture.output({result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)})), 
                            "sreg cannot use individual-level covariates")                  
})

test_that("no cluster sizes warning works", {
  set.seed(123)
  data <- sreg.rgen(
  n = 100, tau.vec = c(0.2, 0.8),
  n.strata = 4, cluster = T, Nmax = 50)
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  G.id <- data$G.id
  Ng <- data$Ng

  expect_warning(invisible(capture.output({result <- sreg::sreg(Y, S = S, D, G.id = G.id, Ng = NULL, X = X, HC1 = TRUE)})), 
                            "cluster sizes have not been provided")                  
})

test_that("empirical example works", {
  library(haven)
  data("AEJapp")
  data <- AEJapp
  Y <- data$gradesq34
  D <- data$treatment
  S <- data$class_level
  pills <- data$pills_taken
  age <- data$age_months
  data.clean <- data.frame(Y, D, S, pills, age)
  data.clean <- data.clean %>%
  mutate(D = ifelse(D == 3, 0, D))
  Y <- data.clean$Y
  D <- data.clean$D
  S <- data.clean$S
  X <- data.frame("pills" = data.clean$pills, "age" = data.clean$age)

  invisible(capture.output({result <- sreg::sreg(Y, S, D)}))

  expect_equal(round(result$tau, 8), c(-0.05112971, 0.40903373))
  expect_equal(round(result$se, 7), c(0.2064541, 0.2065146))

  invisible(capture.output({result <- sreg::sreg(Y, S, D, X = X)}))

  expect_equal(round(result$tau, 8), c(-0.02861589, 0.34608688))
  expect_equal(round(result$se, 7), c(0.1796427, 0.1836229))

  expect_error(invisible(capture.output({result <- sreg::sreg(Y = NULL, S, D)})), 
                  "Observed outcomes have not been provided")
  expect_error(invisible(capture.output({result <- sreg::sreg(Y, S, D = NULL)})), 
                  "Treatments have not been provided")
})
