#%##%##%##%###%##%##%##%###%##%##%##%###%##%
### This R file provides the template   ####
###     for Monte-Carlo simulations     ####
### under CAR with multiple treatments  ####
#%##%##%##%###%##%##%##%###%##%##%##%###%##%
####      The code is developed by      ####
####      @Juri Trifonov, UChicago      ####
####            Supervisors:            ####
####      @Azeem Shaikh, UChicago       ####
####    @Max Tabord-Meehan, UChicago    ####
#%##%##%##%###%##%##%##%###%##%##%##%###%##%
#%##%##%##%##
#%# v.4.1 #%#
#%##%##%##%##
#-------------------------------------------------------------------
#install.packages(c(
#  "compiler",
#  "extraDistr",
#  "VGAM",
#  "Matrix",
#  "parallel",
#  "progress"), dependencies = TRUE)
install.packages(devtools)                                                      # install devtools
library(devtools)                                                               # install devtools
install_github("yurytrifonov/sreg")                                             # install sreg
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
#%##%##%##%###%##%##%##%###%##%##%##%###%##%#%##%##%##%###%##%##%##%##
#%##%##%##%###%##%##%##%###%##%##%##%###%##%#%##%##%##%###%##%##%##%##
#         Please, provide the path to the corresponding source
#                  file with functions on your PC
#                    ↓↓↓↓↓↓↓↓↓↓↓HERE↓↓↓↓↓↓↓↓↓↓↓
rm(list = ls())
#%##%##%##%###%##%##%##%###%##%##%##%###%##%#%##%##%##%###%##%##%##%##
#%##%##%##%###%##%##%##%###%##%##%##%###%##%#%##%##%##%###%##%##%##%##

options(scipen = 999)        # disable scientific notation

# Set up a cluster with available CPU cores
num_cores <- detectCores()
cl <- makeCluster(num_cores,outfile = "")

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
})

# The main function for the Lapply loop
# Function that performs simulations and takes as input
# Only the number of simulation, sim.id
sim.func <- function(sim.id)
{
  seed <- 1000 + sim.id
  set.seed(seed)

  n = 1000
  tau.vec <- c(0.4)
  n.treat <- length(tau.vec)
  n.strata <- 4
  data <- sreg.rgen(n=n, Nmax = 50, n.strata = n.strata, tau.vec = tau.vec, cluster = T, is.cov = TRUE)
  Y <- data$Y
  S <- data$S
  D <- data$D
  X <- data.frame("x_1"= data$x_1, "x_2" = data$x_2)
  Ng <- data$Ng
  G.id <- data$G.id


  # Estimate the ATE, s.e., etc.
  #test <- sreg(Y,S,D,G.id = NULL, Ng = NULL, X = NULL)
  #test$tau.hat

  #model <- lm.iter(Y,D,S,G.id,Ng,X, exp.option =T) # change for exp.option = T if the equal-size design
  #fit <- tau.hat(Y,D,S,G.id,Ng,X,model, exp.option = T)
  result <- tryCatch({sreg(Y, S, D, G.id, Ng, X=NULL, Ng.cov = T, HC1 = FALSE)}, error = function(e) { # tryCatch to avoid errors that stop the execution
    # Print the error message if an error occurs
    cat("Simulation", sim.id, "encountered an error:", conditionMessage(e), "\n")
    # Return a default value or NULL when an error occurs
    NA
  })

  # if condition for NA cases
  if (anyNA(result) == TRUE)
  {
    tau      <- NA
    se       <- NA
    tstat    <- NA
    CI.left  <- NA
    CI.right <- NA
    ci.hit   <- NA
    results <- list(tau = tau,
                    se = se,
                    tstat = tstat,
                    CI.left = CI.left,
                    CI.right = CI.right,
                    ci.hit = ci.hit)
  } else{
    # else condition for Non-NA cases
    tau      <- result$tau.hat
    se       <- result$se.rob
    tstat    <- result$t.stat
    CI.left  <- result$CI.left
    CI.right <- result$CI.right
    ci.hit   <- as.numeric(tau.vec >= result$CI.left &
                             tau.vec <= result$CI.right)
    results <- list(tau = tau,
                    se = se,
                    tstat = tstat,
                    CI.left = CI.left,
                    CI.right = CI.right,
                    ci.hit = ci.hit)}
  message(paste("Simulation", sim.id, "completed succesfully"))

  return(results)

}

# Parallelize the simulations and store the results
simres <- parLapply(cl, 1:5000, sim.func)
#mb <- microbenchmark(parLapply(cl, 1:100, sim.func), times = 1)
save(simres, file = "/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/creg/2000.RData")
###################
# Close the cluster
stopCluster(cl)
# Close the cluster
###################

# Extract parameters of interest from the results
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
mean(tau)
sd(tau)
mean(se)
mean(ci.hit)

rowMeans(tau)
apply(tau, 1, sd)
rowMeans(se)
rowMeans(ci.hit)



G = 1000
Nmax=500;
tau.vec <- c(0)
n.treat <- length(tau.vec)
max.support = Nmax/10-1;
gamma.vec <-c(0.4, 0.2, 1)
n.strata <- 4

Ng <- gen.cluster.sizes(G, max.support)[,1]
#Ng <- rep(Nmax, G)                                                            # uncomment and comment the previous line for a equal-size design
data.pot <- gen.data.pot(Ng=Ng, tau.vec = (tau.vec / 0.5), G = G,
                         gamma.vec = gamma.vec, n.treat=n.treat)
strata <- form.strata(data.pot, n.strata)
strata.set <- data.frame(strata)
strata.set$S <- max.col(strata.set)
pi.vec <- rep(c(1 / (n.treat + 1)), n.treat)
data.sim <- dgp.obs(data.pot, I.S = strata, pi.vec, n.treat)

finale <- data.frame('Y'= data.sim$Y, 'D' = data.sim$D)
Y <- data.sim$Y
D <- data.sim$D
S <- data.sim$S
X <- data.sim$X
Ng <- data.sim$Ng
G.id <- data.sim$G.id

model <- lm.iter(Y,D,S,G.id,Ng,X)
df <- data.frame(Y,D,S,G.id,Ng,X)
lin.adj(1, model = model, data = df)

est <- tau.hat(Y,D,S,G.id,Ng,X,model)
est$tau.hat
as.var(model, est)
#model <- lm.iter(Y,D,S,G.id,Ng,X, exp.option =T) # change for exp.option = T if the equal-size design
#fit <- tau.hat(Y,D,S,G.id,Ng,X,model, exp.option = T)
result <- tryCatch({sreg(Y,D,S,G.id,Ng,X, exp.option = F)}, error = function(e) { # tryCatch to avoid errors that stop the execution
  # Print the error message if an error occurs
  cat("Simulation", sim.id, "encountered an error:", conditionMessage(e), "\n")
  # Return a default value or NULL when an error occurs
  NA
})


set.seed(123)
n = 1000
tau.vec <- c(0.8)
n.treat <- length(tau.vec)
n.strata <- 4
data <- sreg.rgen(n=n, tau.vec=tau.vec, n.strata=n.strata, cluster = F, is.cov = TRUE)
data <- sreg.rgen(n=n, Nmax = 50, n.strata = 4, tau.vec = tau.vec, cluster = T, is.cov = TRUE)
Y <- data$Y
S <- data$S
D <- data$D
X <- data.frame("x_1"= data$x_1, "x_2" = data$x_2)
Ng <- data$Ng
G.id <- data$G.id

test <- sreg(Y, S, D, G.id, Ng, X, Ng.cov = T, HC1 = FALSE)
