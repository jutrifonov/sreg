install.packages(devtools) # install devtools
library(devtools) # install devtools
install_github("yurytrifonov/sreg") # install sreg
library(sreg) # install sreg
packageDescription("sreg") # package description
?sreg # R documentation for sreg()
?sreg.rgen # R documentation for sreg.rgen()
# %#%#%#%#%#%#%#%#%#%
# %#%#%#%#%#%#%#%#%#%

set.seed(123) # fix the random seed

# Generate a pseudo-random sample without clusters and with only one treatment = 0
data <- sreg.rgen(n = 250, tau.vec = c(0.2, 0), n.strata = 4, cluster = F)
Y <- data$Y
S <- data$S
D <- data$D
X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)

# Estimate the ATE, s.e., etc.
test <- sreg(Y, S, D, G.id = NULL, Ng = NULL, X)

# Generate a pseudo-random sample with clusters and two treatments = c(0.2, 0.8)
data <- sreg.rgen(
  n = 1000, tau.vec = c(0.2, 0.8),
  n.strata = 4, cluster = T, Nmax = 50
)
Y <- data$Y
S <- data$S
D <- data$D
X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
Ng <- data$Ng
G.id <- data$G.id

# Estimate the ATE, s.e., etc.
test <- sreg(Y, S, D, G.id, Ng, X = NULL, Ng.cov = T)
test <- sreg(Y, S, D, G.id, Ng, X = X, Ng.cov = T, HC1 = F)
