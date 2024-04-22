test_that("empirical example works", {

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
})
