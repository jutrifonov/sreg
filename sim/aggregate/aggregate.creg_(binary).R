library(ggplot2)

load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/creg/100.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.100 <- mean(tau)
se.num.100 <- sd(tau)
se.analyt.100 <- mean(se)
ci.hit.100 <- mean(ci.hit)


load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/creg/250.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.250 <- mean(tau)
se.num.250 <- sd(tau)
se.analyt.250 <- mean(se)
ci.hit.250 <- mean(ci.hit)

load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/creg/500.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.500 <- mean(tau)
se.num.500 <- sd(tau)
se.analyt.500 <- mean(se)
ci.hit.500 <- mean(ci.hit)

load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/creg/750.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.750 <- mean(tau)
se.num.750 <- sd(tau)
se.analyt.750 <- mean(se)
ci.hit.750 <- mean(ci.hit)


load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/creg/1000.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.1000 <- mean(tau)
se.num.1000 <- sd(tau)
se.analyt.1000 <- mean(se)
ci.hit.1000 <- mean(ci.hit)

load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/creg/1500.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.1500 <- mean(tau)
se.num.1500 <- sd(tau)
se.analyt.1500 <- mean(se)
ci.hit.1500 <- mean(ci.hit)

load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/creg/2000.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.2000 <- mean(tau)
se.num.2000 <- sd(tau)
se.analyt.2000 <- mean(se)
ci.hit.2000 <- mean(ci.hit)

tau.vec <- c(tau.100, tau.250, tau.500, tau.750,
            tau.1000, tau.1500, tau.2000)
se.num.vec <- c(se.num.100, se.num.250, se.num.500, se.num.750,
                se.num.1000, se.num.1500, se.num.2000)
se.analyt.vec <- c(se.analyt.100, se.analyt.250, se.analyt.500, se.analyt.750, se.analyt.1000, se.analyt.1500, se.analyt.2000)
ci.hit.vec <- c(ci.hit.100, ci.hit.250, ci.hit.500, ci.hit.750, ci.hit.1000, ci.hit.1500, ci.hit.2000)
sample.set <- c(100, 250, 500, 750, 1000, 1500, 2000)

graph.data <- data.frame('s.size' = sample.set,
                         'tau' = tau.vec,
                         'se.num' = se.num.vec,
                         'se.analyt' = se.analyt.vec,
                         'se.diff' = abs(se.num.vec - se.analyt.vec),
                         'ci.hit' = ci.hit.vec)


ggplot(graph.data) +
  geom_line(aes(x = s.size, y = se.diff, color = "SE"), size = 1) +
  geom_point(aes(x = s.size, y = se.diff, color = "SE"), size = 3) +
  labs(x = "Sample Size", y = "SE Difference") +
  ggtitle("Graph of Sample Size vs. SE Difference") +
  theme_minimal()

ggplot(graph.data) +
  geom_line(aes(x = s.size, y = se.num, color = "SE Num"), size = 1) +
  geom_point(aes(x = s.size, y = se.num, color = "SE Num"), size = 3) +
  geom_line(aes(x = s.size, y = se.analyt, color = "SE Analyt"), size = 1) +
  geom_point(aes(x = s.size, y = se.analyt, color = "SE Analyt"), size = 3) +
  labs(x = "Sample Size", y = "SE") +
  theme_minimal()

ggplot(graph.data) +
  geom_line(aes(x = s.size, y = ci.hit, color = "hit"), size = 1) +
  geom_point(aes(x = s.size, y = ci.hit, color = "hit"), size = 3) +
  labs(x = "Sample Size", y = "Coverage %") + scale_y_continuous(limits = c(0.6, 1))



