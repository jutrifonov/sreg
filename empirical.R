install.packages("devtools") # install sreg
install.packages("haven")
install.packages("dplyr")
library(devtools)
install_github("yurytrifonov/sreg") # install sreg

library(devtools) # load devtools
library(sreg) # load sreg
library(dplyr) # load dplyr (to transform data)
library(haven)
packageDescription("sreg") # package description
?sreg # R documentation for sreg()
?sreg.rgen # R documentation for sreg.rgen()
?AEJapp # R documentation for AEJapp dataset
# %#%#%#%#%#%#%#%#%#%
# %#%#%#%#%#%#%#%#%#%
data("AEJapp") # upload the Chong et al. (2016) data from the package
data <- AEJapp # rename for convenience
head(data)
# %#%#%#%#%#%#%#%#%#%
# Replicate the empirical illustration
# from (Bugni et al, 2019)
# %#%#%#%#%#%#%#%#%#%
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
table(D = data.clean$D, S = data.clean$S)
# %#%#%#%#%#%#%#%#%#%
result <- sreg::sreg(Y, S, D, HC1 = T)
# %#%#%#%#%#%#%#%#%#%
result <- sreg::sreg(Y = Y, S = S, D = D, X = X)
