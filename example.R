install.packages(devtools) # install devtools
install_github("yurytrifonov/sreg") # install sreg
library(devtools) # install devtools
library(sreg) # install sreg
packageDescription("sreg") # package description
# ?sreg # R documentation for sreg()
# ?sreg.rgen # R documentation for sreg.rgen()
# %#%#%#%#%#%#%#%#%#%
# %#%#%#%#%#%#%#%#%#%

set.seed(123) # fix the random seed

# Generate a pseudo-random sample without clusters and with only one treatment = 0
data <- sreg.rgen(n = 20, tau.vec = c(0.8), n.strata = 5, cluster = F, is.cov = TRUE)
Y <- data$Y
S <- data$S
D <- data$D
X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
# Estimate the ATE, s.e., etc.
test <- sreg(Y, S, D, G.id = NULL, Ng = NULL, X = X)

set.seed(123) # fix the random seed
# Generate a pseudo-random sample with clusters and two treatments = c(0.2, 0.8)
data <- sreg.rgen(
  n = 1000, tau.vec = c(0.2, 0.8),
  n.strata = 4, cluster = T, Nmax = 50
)
Y <- data$Y
S <- data$S
D <- data$D
X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
head(X)
Ng <- data$Ng
#Ng[1:30] <- 40
X[1, 1] <- 2.3894
X[1, 2] <- 5.2090398
X[2, 1] <- 4.2049587
X[2, 2] <- 9.01023948
G.id <- data$G.id
save <- res.creg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X, HC1 = TRUE)
save$tau.hat
save$se.rob
summary.creg(save)

#X <- data.frame("Ng" = data$Ng, "x_1" = data$x_1, "x_2" = data$x_2)
# X <- data.frame("Ng" = data$Ng)

data <- data.frame(G.id, X)
head(data)
# Estimate the ATE, s.e., etc.
test <- sreg(Y, S = S, D, G.id = G.id, Ng = Ng, X = X)
test <- sreg(Y, S = NULL, D, G.id, Ng, X = NULL, HC1 = FALSE)

# Problem with not enough df
set.seed(123) # fix the random seed
# Generate a pseudo-random sample without clusters and with only one treatment = 0
data <- sreg.rgen(n = 100, tau.vec = c(0.8, 0), n.strata = 6, cluster = F, is.cov = TRUE)
Y <- data$Y
S <- data$S
D <- data$D
X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
# Estimate the ATE, s.e., etc.
test <- sreg(Y, S, D, G.id = NULL, Ng = NULL, X = X)
