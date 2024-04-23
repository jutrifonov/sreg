#-------------------------------------------------------------------
# %#     Function that implements \hat{\sigma^2} --
# %#     i.e. the variance estimator
#-------------------------------------------------------------------
as.var.sreg <- function(Y, S, D, X = NULL, model = NULL, tau, HC1)
#-------------------------------------------------------------------
{
  var.vec <- numeric(max(D))
  n.vec <- numeric(max(D))

  if (!is.null(X)) {
    for (d in 1:max(D))
    {
      data <- data.frame(Y, S, D, X)
      data$pi <- pi.hat.sreg(S, D)[, d]
      data$pi.0 <- pi.hat.sreg(S, D, inverse = T)[, 1]
      n <- length(Y)
      data$A <- ifelse(D == d, 1, ifelse(D == 0, 0, -999999))
      data$I <- as.numeric(data$A != -999999)

      mu.hat.d <- lin.adj.sreg(d, data$S, data[4:(4 + ncol(X) - 1)], model)
      mu.hat.0 <- lin.adj.sreg(0, data$S, data[4:(4 + ncol(X) - 1)], model)

      Xi.tilde.1 <- (mu.hat.d - mu.hat.0) +
        (data$Y - mu.hat.d) / data$pi

      Xi.tilde.0 <- (mu.hat.d - mu.hat.0) -
        (data$Y - mu.hat.0) / data$pi.0

      data <- data.frame(data, Xi.tilde.1, Xi.tilde.0, Y.tau.D = data$Y - tau[d] * data$A * data$I)

      count.Xi.1 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.1 = mean(.data$Xi.tilde.1)) %>%
        filter(.data$A != -999999)
      count.Xi.0 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.0 = mean(.data$Xi.tilde.0)) %>%
        filter(.data$A != -999999)
      count.Y <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Y.tau = mean(.data$Y.tau.D)) %>%
        filter(.data$A != -999999)

      j <- left_join(count.Xi.1, count.Xi.0, by = join_by("S" == "S", "A" == "A")) %>% left_join(count.Y, by = join_by("S" == "S", "A" == "A"))

      Xi.tilde.1.all <- j %>%
        select(c("S", "A", "Xi.mean.1")) %>%
        spread(key = "A", value = "Xi.mean.1")
      Xi.tilde.0.all <- j %>%
        select(c("S", "A", "Xi.mean.0")) %>%
        spread(key = "A", value = "Xi.mean.0")
      Y.tau.D.all <- j %>%
        select(c("S", "A", "Y.tau")) %>%
        spread(key = "A", value = "Y.tau")

      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.tau.D.mean <- as.matrix(select(data.frame(Y.tau.D.all), -1))

      Xi.1.mean <- Xi.tilde.1.mean[S, 2]
      Xi.0.mean <- Xi.tilde.0.mean[S, 1]
      Y.tau.D.1.mean <- Y.tau.D.mean[S, 2]
      Y.tau.D.0.mean <- Y.tau.D.mean[S, 1]

      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean
      Xi.hat.2 <- Y.tau.D.1.mean - Y.tau.D.0.mean

      sigma.hat.sq <- mean(data$I * (data$A * (Xi.hat.1)^2 + (1 - data$A) * (Xi.hat.0)^2) + Xi.hat.2^2)
      if (HC1 == TRUE) {
        var.vec[d] <- (mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2))) * (n / (n - (max(S) + max(D) * max(S)))) +
          mean(Xi.hat.2^2)
      } else {
        var.vec[d] <- sigma.hat.sq
      }
      n.vec[d] <- n
    }
  } else {
    for (d in 1:max(D))
    {
      data <- data.frame(Y, S, D)
      data$pi <- pi.hat.sreg(S, D)[, d]
      data$pi.0 <- pi.hat.sreg(S, D, inverse = T)[, 1]
      n <- length(Y)
      data$A <- ifelse(D == d, 1, ifelse(D == 0, 0, -999999))
      data$I <- as.numeric(data$A != -999999)

      mu.hat.d <- 0
      mu.hat.0 <- 0

      Xi.tilde.1 <- (mu.hat.d - mu.hat.0) +
        (data$Y - mu.hat.d) / data$pi

      Xi.tilde.0 <- (mu.hat.d - mu.hat.0) -
        (data$Y - mu.hat.0) / data$pi.0

      data <- data.frame(data, Xi.tilde.1, Xi.tilde.0, Y.tau.D = data$Y - tau[d] * data$A * data$I)

      count.Xi.1 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.1 = mean(.data$Xi.tilde.1)) %>%
        filter(.data$A != -999999)
      count.Xi.0 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.0 = mean(.data$Xi.tilde.0)) %>%
        filter(.data$A != -999999)
      count.Y <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Y.tau = mean(.data$Y.tau.D)) %>%
        filter(.data$A != -999999)

      j <- left_join(count.Xi.1, count.Xi.0, by = join_by("S" == "S", "A" == "A")) %>% left_join(count.Y, by = join_by("S" == "S", "A" == "A"))

      Xi.tilde.1.all <- j %>%
        select(c("S", "A", "Xi.mean.1")) %>%
        spread(key = "A", value = "Xi.mean.1")
      Xi.tilde.0.all <- j %>%
        select(c("S", "A", "Xi.mean.0")) %>%
        spread(key = "A", value = "Xi.mean.0")
      Y.tau.D.all <- j %>%
        select(c("S", "A", "Y.tau")) %>%
        spread(key = "A", value = "Y.tau")

      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.tau.D.mean <- as.matrix(select(data.frame(Y.tau.D.all), -1))

      Xi.1.mean <- Xi.tilde.1.mean[S, 2]
      Xi.0.mean <- Xi.tilde.0.mean[S, 1]
      Y.tau.D.1.mean <- Y.tau.D.mean[S, 2]
      Y.tau.D.0.mean <- Y.tau.D.mean[S, 1]

      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean
      Xi.hat.2 <- Y.tau.D.1.mean - Y.tau.D.0.mean

      sigma.hat.sq <- mean(data$I * (data$A * (Xi.hat.1)^2 + (1 - data$A) * (Xi.hat.0)^2) + Xi.hat.2^2)
      if (HC1 == TRUE) {
        var.vec[d] <- (mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2))) * (n / (n - (max(S) + max(D) * max(S)))) +
          mean(Xi.hat.2^2)
      } else {
        var.vec[d] <- sigma.hat.sq
      }
      n.vec[d] <- n
    }
  }

  se.vec <- sqrt(var.vec / n.vec)
  return(se.vec)
}
#-------------------------------------------------------------------
as.var.creg <- function(model = NULL, fit, HC1)
#-------------------------------------------------------------------
{
  var.vec <- numeric(length(fit$tau.hat))
  n.vec <- numeric(length(fit$tau.hat))

  if (!is.null(model)) {
    for (d in 1:length(fit$tau.hat))
    {
      Y.bar.g <- fit$Y.bar.g
      Ng <- fit$Ng
      mu.hat.0 <- fit$mu.hat[[d]][, 1]
      mu.hat.d <- fit$mu.hat[[d]][, 2]
      tau.est <- fit$tau.hat
      pi.hat <- fit$pi.hat[[d]]
      pi.hat.0 <- fit$pi.hat.0
      data <- fit$data.list[[d]]
      n <- length(Y.bar.g)

      Xi.tilde.1 <- (mu.hat.d - mu.hat.0) +
        (Ng * Y.bar.g - mu.hat.d) / pi.hat

      Xi.tilde.0 <- (mu.hat.d - mu.hat.0) -
        (Ng * Y.bar.g - mu.hat.0) / pi.hat.0

      data <- data.frame(data, Xi.tilde.1, Xi.tilde.0, Y.Ng = Y.bar.g * Ng)

      count.Xi.1 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.1 = mean(.data$Xi.tilde.1)) %>%
        filter(.data$A != -999999)
      count.Xi.0 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.0 = mean(.data$Xi.tilde.0)) %>%
        filter(.data$A != -999999)
      count.Y <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Y.bar = mean(.data$Y.Ng)) %>%
        filter(.data$A != -999999)
      count.Ng <- data %>%
        group_by(.data$S) %>%
        summarise(Ng.bar = mean(.data$Ng))

      j <- left_join(count.Xi.1, count.Xi.0, by = join_by("S" == "S", "A" == "A")) %>%
        left_join(count.Y, by = join_by("S" == "S", "A" == "A")) %>%
        left_join(count.Ng, by = join_by("S" == "S"))

      Xi.tilde.1.all <- j %>%
        select(c("S", "A", "Xi.mean.1")) %>%
        spread(key = "A", value = "Xi.mean.1")
      Xi.tilde.0.all <- j %>%
        select(c("S", "A", "Xi.mean.0")) %>%
        spread(key = "A", value = "Xi.mean.0")
      Y.Ng.all <- j %>%
        select(c("S", "A", "Y.bar")) %>%
        spread(key = "A", value = "Y.bar")
      Ng.bar.all <- j %>%
        select(c("S", "A", "Ng.bar")) %>%
        spread(key = "A", value = "Ng.bar")

      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.Ng.mean <- as.matrix(select(data.frame(Y.Ng.all), -1))
      Ng.bar.mean <- as.matrix(select(data.frame(Ng.bar.all), -1))

      Xi.1.mean <- Xi.tilde.1.mean[data$S, 2]
      Xi.0.mean <- Xi.tilde.0.mean[data$S, 1]
      Y.g.bar.cl.1 <- Y.Ng.mean[data$S, 2]
      Y.g.bar.cl.0 <- Y.Ng.mean[data$S, 1]
      N.g.bar.cl <- Ng.bar.mean[data$S, 1]

      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.2 <- Y.g.bar.cl.1 - Y.g.bar.cl.0 - tau.est[d] * N.g.bar.cl

      sigma.hat.sq <- mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2) + Xi.hat.2^2) / (mean(Ng))^2

      if (HC1 == TRUE) {
        var.vec[d] <- ((mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2))) * (n / (n - (max(data$S) + max(data$D) * max(data$S)))) +
          mean(Xi.hat.2^2)) / (mean(Ng))^2
      } else {
        var.vec[d] <- sigma.hat.sq
      }
      n.vec[d] <- n
    }
  } else {
    for (d in 1:length(fit$tau.hat))
    {
      Y.bar.g <- fit$Y.bar.g
      Ng <- fit$Ng
      tau.est <- fit$tau.hat
      pi.hat <- fit$pi.hat[[d]]
      pi.hat.0 <- fit$pi.hat.0
      data <- fit$data.list[[d]]
      n <- length(Y.bar.g)

      mu.hat.0 <- 0
      mu.hat.d <- 0

      Xi.tilde.1 <- (mu.hat.d - mu.hat.0) +
        (Ng * Y.bar.g - mu.hat.d) / pi.hat

      Xi.tilde.0 <- (mu.hat.d - mu.hat.0) -
        (Ng * Y.bar.g - mu.hat.0) / pi.hat.0

      data <- data.frame(data, Xi.tilde.1, Xi.tilde.0, Y.Ng = Y.bar.g * Ng)

      count.Xi.1 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.1 = mean(.data$Xi.tilde.1)) %>%
        filter(.data$A != -999999)
      count.Xi.0 <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Xi.mean.0 = mean(.data$Xi.tilde.0)) %>%
        filter(.data$A != -999999)
      count.Y <- data %>%
        group_by(.data$S, .data$A) %>%
        summarise(Y.bar = mean(.data$Y.Ng)) %>%
        filter(.data$A != -999999)
      count.Ng <- data %>%
        group_by(.data$S) %>%
        summarise(Ng.bar = mean(.data$Ng))

      j <- left_join(count.Xi.1, count.Xi.0, by = join_by("S" == "S", "A" == "A")) %>%
        left_join(count.Y, by = join_by("S" == "S", "A" == "A")) %>%
        left_join(count.Ng, by = join_by("S" == "S"))

      Xi.tilde.1.all <- j %>%
        select(c("S", "A", "Xi.mean.1")) %>%
        spread(key = "A", value = "Xi.mean.1")
      Xi.tilde.0.all <- j %>%
        select(c("S", "A", "Xi.mean.0")) %>%
        spread(key = "A", value = "Xi.mean.0")
      Y.Ng.all <- j %>%
        select(c("S", "A", "Y.bar")) %>%
        spread(key = "A", value = "Y.bar")
      Ng.bar.all <- j %>%
        select(c("S", "A", "Ng.bar")) %>%
        spread(key = "A", value = "Ng.bar")

      Xi.tilde.1.mean <- as.matrix(select(data.frame(Xi.tilde.1.all), -1))
      Xi.tilde.0.mean <- as.matrix(select(data.frame(Xi.tilde.0.all), -1))
      Y.Ng.mean <- as.matrix(select(data.frame(Y.Ng.all), -1))
      Ng.bar.mean <- as.matrix(select(data.frame(Ng.bar.all), -1))

      Xi.1.mean <- Xi.tilde.1.mean[data$S, 2]
      Xi.0.mean <- Xi.tilde.0.mean[data$S, 1]
      Y.g.bar.cl.1 <- Y.Ng.mean[data$S, 2]
      Y.g.bar.cl.0 <- Y.Ng.mean[data$S, 1]
      N.g.bar.cl <- Ng.bar.mean[data$S, 1]

      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean - tau.est[d] * (Ng - N.g.bar.cl)
      Xi.hat.2 <- Y.g.bar.cl.1 - Y.g.bar.cl.0 - tau.est[d] * N.g.bar.cl

      if (HC1 == TRUE) {
        var.vec[d] <- ((mean(data$I * (data$A * Xi.hat.1^2 + (1 - data$A) * Xi.hat.0^2))) * (n / (n - (max(data$S) + max(data$D) * max(data$S)))) +
          mean(Xi.hat.2^2)) / mean(Ng)^2
      } else {
        sigma.hat.sq <- mean(data$I * (data$A * (Xi.hat.1)^2 + (1 - data$A) * (Xi.hat.0)^2) + Xi.hat.2^2) / (mean(Ng))^2
        var.vec[d] <- sigma.hat.sq
      }
      n.vec[d] <- n
    }
  }
  se.vec <- sqrt(var.vec / n.vec)
  return(se.vec)
}
