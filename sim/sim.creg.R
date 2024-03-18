# %##%##%##%###%##%##%##%###%##%##%##%###%##%
### This R file provides the template   ####
###     for Monte-Carlo simulations     ####
### under CAR with multiple treatments  ####
# %##%##%##%###%##%##%##%###%##%##%##%###%##%
####      The code is developed by      ####
####      @Juri Trifonov, UChicago      ####
####            Supervisors:            ####
####      @Azeem Shaikh, UChicago       ####
####    @Max Tabord-Meehan, UChicago    ####
# %##%##%##%###%##%##%##%###%##%##%##%###%##%
# %##%##%##%##
# %# v.1.2.0 #%#
# %##%##%##%##
#-------------------------------------------------------------------
# install.packages(c(
#  "compiler",
#  "extraDistr",
#  "VGAM",
#  "Matrix",
#  "parallel",
#  "progress"), dependencies = TRUE)
install.packages(devtools) # install devtools
install_github("yurytrifonov/sreg") # install sreg
library(devtools) # install devtools
library(sreg)
library(sandwich)
library(lubridate)
library(compiler)
library(extraDistr)
library(VGAM)
library(Matrix)
library(parallel)
library(progress)
library(purrr)
library(SimDesign)
library(pbapply)

###############
rm(list = ls())
###############
options(scipen = 999) # disable scientific notation
###############

# Set up a cluster with available CPU cores
num_cores <- detectCores()
cl <- makeCluster(num_cores, outfile = "")

# Upload the libraries to the cluster
clusterEvalQ(cl, {
  library(sandwich)
  library(lubridate)
  library(compiler)
  library(extraDistr)
  library(Matrix)
  library(progress)
  library(parallel)
  library(sreg)
  library(pbapply)
})

# The main function for the Lapply loop
# Function that performs simulations and takes as input
# Only the number of simulation, sim.id
sim.func <- function(sim.id) {
  #output <- capture.output({
  seed <- 1000 + sim.id
  set.seed(seed)

  n <- 200
  tau.vec <- c(0.8, 0.4)
  n.treat <- length(tau.vec)
  n.strata <- 2
  data <- sreg.rgen(n = n, Nmax = 50, n.strata = n.strata, tau.vec = tau.vec, cluster = T, is.cov = TRUE)
  Y <- data$Y
  S <- data$S
  D <- data$D
  #X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
  X <- data.frame("Ng" = data$Ng, "x_1" = data$x_1, "x_2" = data$x_2)
  Ng <- data$Ng
  G.id <- data$G.id
  result <- tryCatch(
    {
      #sreg(Y, S = S, D, G.id, Ng, X = X)
      sreg(Y, S = NULL, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
    },
    error = function(e) { # tryCatch to avoid errors that stop the execution
      # Print the error message if an error occurs
      cat("Simulation", sim.id, "encountered an error:", conditionMessage(e), "\n")
      # Return a default value or NULL when an error occurs
      NA
    }
  )

  # if condition for NA cases
  if (anyNA(result) == TRUE) {
    tau <- NA
    se <- NA
    tstat <- NA
    CI.left <- NA
    CI.right <- NA
    ci.hit <- NA
    results <- list(
      tau = tau,
      se = se,
      tstat = tstat,
      CI.left = CI.left,
      CI.right = CI.right,
      ci.hit = ci.hit
    )
  } else {
    # else condition for Non-NA cases
    tau <- result$tau.hat
    se <- result$se.rob
    tstat <- result$t.stat
    CI.left <- result$CI.left
    CI.right <- result$CI.right
    ci.hit <- as.numeric(tau.vec >= result$CI.left &
      tau.vec <= result$CI.right)
    results <- list(
      tau = tau,
      se = se,
      tstat = tstat,
      CI.left = CI.left,
      CI.right = CI.right,
      ci.hit = ci.hit
    )
  }
  message(paste("Simulation", sim.id, "is completed succesfully"))
  
  return(results)
  #})
  #return(output)
}

# Parallelize the simulations and store the results
simres <- parLapply(cl, 1:1000, sim.func)
simres <- pblapply(1:10000, sim.func, cl=cl)
# mb <- microbenchmark(parLapply(cl, 1:100, sim.func), times = 1)
save(simres, file = "/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/v.1.2.5/creg.cov (all 100k iter)/1000.RData")
###################
# Close the cluster
stopCluster(cl)
# Close the cluster
###################

# Extract parameters of interest from the results
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
rowMeans(tau)
apply(tau, 1, sd)
rowMeans(se)
rowMeans(ci.hit[, 1:5000])

mean(tau)
sd(tau)
mean(se)
mean(ci.hit)