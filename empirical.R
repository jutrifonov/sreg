install.packages("devtools")                                                    # install sreg
install.packages("haven")
install.packages("dplyr")
install_github("yurytrifonov/sreg")                                             # install sreg

library(devtools)                                                               # load devtools
library(sreg)                                                                   # load sreg
library(haven)                                                                  # load haven (to read dta files)
library(dplyr)                                                                  # load dplyr (to transform data)
packageDescription("sreg")                                                      # package description
?sreg                                                                           # R documentation for sreg()
?sreg.rgen                                                                      # R documentation for sreg.rgen()
?AEJapp                                                                         # R documentation for AEJapp dataset
#%#%#%#%#%#%#%#%#%#%
#%#%#%#%#%#%#%#%#%#%
data("AEJapp")                                                                  # upload the Chong et al. (2016) data from the package
data <- AEJapp                                                                  # rename for convenience
head(data)
View(data)
#%#%#%#%#%#%#%#%#%#%
# Replicate the empirical illustration
# from (Bugni et al, 2019)
#%#%#%#%#%#%#%#%#%#%
Y <- data$gradesq34
D <- data$treatment
S <- data$class_level
data.clean <- data.frame(Y,D,S)
data.clean <- data.clean %>%
  mutate(D = ifelse(D == 3, 0, D))
Y = data.clean$Y
D = data.clean$D
S = data.clean$S
table(D = data.clean$D, S = data.clean$S)
#%#%#%#%#%#%#%#%#%#%
result <- sreg::sreg(Y,S,D, HC1 = T)
#%#%#%#%#%#%#%#%#%#%









################################################################################
##v.1.4
#result <- sreg(Y,D,S)                                                            # estimate the ATEs
#summary.sreg(result)
#
#table(D = data.clean$D, S = data.clean$S)
#source('/Users/trifonovjuri/Desktop/sreg.source/sreg.func_v.2.0.R')
#result <- sreg(Y,S,D)
#result$vcov.rob
#result <- sreg::sreg(Y,S,D)
#summary.sreg(result)
#summary.sreg(result)
#model <- lm.iter.sreg(Y = data.clean$Y, S = data.clean$S, D = data.clean$D)
#tau.hat.sreg(Y = data.clean$Y, S = data.clean$S, D = data.clean$D)
#
#tau.hat(Y = data.clean$Y, S = data.clean$S, D = data.clean$D)
#tau.hat.sreg.alt(Y = data.clean$Y, S = data.clean$S, D = data.clean$D)
#tau <- tau.hat.sreg(Y = data.clean$Y, S = data.clean$S, D = data.clean$D)
#length(Y)
#tau <- c(-0.051, 0.409 )
#as.var.sreg(Y = data.clean$Y, S = data.clean$S, D = data.clean$D, tau = tau)
#
#
#### Checking the variance estimator
#source("~/Desktop/sreg.source/sreg.func_v.2.0.R")
#result <- sreg(Y = data.clean$Y, S = data.clean$S, D = data.clean$D)
#se.corr <- result$se.rob
#se.raw <- result$se.rob
#
#result$se.rob
#se.corr / se.raw
#summary.sreg(result)
#
#
#
#
#
#data.bcs <- read.csv("/Users/trifonovjuri/Downloads/Application/data_case1.csv")
#Y<-data.bcs[,1];
#D<-data.bcs[,2];
#n<-length(Y);
#I.S<-matrix(unlist(data.bcs[,3:7]),n,5)
#strata.set <-as.data.frame(I.S)
#strata.set$S <- max.col(strata.set)
#S <- strata.set$S
#res <- sreg::sreg(Y,S,D)










