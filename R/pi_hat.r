#-------------------------------------------------------------------
# %#     Function that implements the calculation of \hat{\pi} --
# %#     i.e., calculates the proportions assigned to treatments
#-------------------------------------------------------------------
pi.hat.sreg <- function(S, D, inverse = FALSE)
#-------------------------------------------------------------------
{
  n <- length(D)
  data <- data.frame(S, D)
  counts <- data %>%
    group_by(S, D) %>%
    summarise(n = n())
  scount <- data %>%
    group_by(S) %>%
    summarise(ns = n())

  j <- left_join(counts, scount, by = join_by(S == S))
  j$pi_hat <- j$n / j$ns
  pi_hat_all <- j %>%
    select(c("S", "D", "pi_hat")) %>%
    spread(key = "D", value = "pi_hat")
  if (inverse) {
    n_repeat <- max(counts$D)
    ret_df <- matrix(replicate(n_repeat, pi_hat_all$"0"), nrow = nrow(pi_hat_all))
  } else {
    pi.hat.df <- select(data.frame(pi_hat_all), -c(1, 2))
    ret_df <- as.matrix(pi.hat.df)
  }
  return(as.matrix(ret_df[S, ]))
}
#-------------------------------------------------------------------
pi.hat.creg <- function(S, D, inverse = FALSE)
#-------------------------------------------------------------------
{
  n <- length(D)
  data <- data.frame(S, D)
  counts <- data %>%
    group_by(.data$S, .data$D) %>%
    summarise(n = n())
  scount <- data %>%
    group_by(.data$S) %>%
    summarise(ns = n())

  j <- left_join(counts, scount, by = join_by(S == S))
  j$pi_hat <- j$n / j$ns
  pi_hat_all <- j %>%
    select(c("S", "D", "pi_hat")) %>%
    spread(key = "D", value = "pi_hat")
  if (inverse) {
    n_repeat <- max(counts$D)
    ret_df <- matrix(replicate(n_repeat, pi_hat_all$"0"), nrow = nrow(pi_hat_all))
  } else {
    pi.hat.df <- select(data.frame(pi_hat_all), -c(1, 2))
    ret_df <- as.matrix(pi.hat.df)
  }
  return(as.matrix(ret_df[S, ]))
}
