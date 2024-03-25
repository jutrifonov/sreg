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
data <- sreg.rgen(n = 10000, tau.vec = c(0, 0.2), n.strata = 4, cluster = F, is.cov = TRUE)
Y <- data$Y
S <- data$S
D <- data$D
X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
# Estimate the ATE, s.e., etc.
test <- sreg(Y, S, D, G.id = NULL, Ng = NULL, X = X)

set.seed(123) # fix the random seed

# Generate a pseudo-random sample with clusters and two treatments = c(0.2, 0.8)
data <- sreg.rgen(
  n = 100, tau.vec = c(0.2, 0.8),
  n.strata = 4, cluster = T, Nmax = 50
)
Y <- data$Y
S <- data$S
D <- data$D
X <- data.frame("x_1" = data$x_1, "x_2" = data$x_2)
#X <- data.frame("Ng" = data$Ng, "x_1" = data$x_1, "x_2" = data$x_2)
# X <- data.frame("Ng" = data$Ng)
Ng <- data$Ng
Ng[1:30] <- 40
G.id <- data$G.id

# Estimate the ATE, s.e., etc.
test <- sreg(Y, S = S, D, G.id = G.id, Ng = NULL, X = X)
test <- sreg(Y, S = NULL, D, G.id, Ng, X = NULL, HC1 = FALSE)
