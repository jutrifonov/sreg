library(ggplot2)

load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/sreg/adj/New version/100.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.100 <- rowMeans(tau)
se.num.100 <- apply(tau, 1, sd)
se.analyt.100 <- rowMeans(se)
ci.hit.100 <- rowMeans(ci.hit)


load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/sreg/adj/New version/250.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.250 <- rowMeans(tau)
se.num.250 <- apply(tau, 1, sd)
se.analyt.250 <- rowMeans(se)
ci.hit.250 <- rowMeans(ci.hit)

load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/sreg/adj/New version/500.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.500 <- rowMeans(tau)
se.num.500 <- apply(tau, 1, sd)
se.analyt.500 <- rowMeans(se)
ci.hit.500 <- rowMeans(ci.hit)

load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/sreg/adj/New version/750.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.750 <- rowMeans(tau)
se.num.750 <- apply(tau, 1, sd)
se.analyt.750 <- rowMeans(se)
ci.hit.750 <- rowMeans(ci.hit)


load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/sreg/adj/New version/1000.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.1000 <- rowMeans(tau)
se.num.1000 <- apply(tau, 1, sd)
se.analyt.1000 <- rowMeans(se)
ci.hit.1000 <- rowMeans(ci.hit)

#load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/sreg/adj/New version/2000.RData')
#tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
#se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
#ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
#tau.1500 <- rowMeans(tau)
#se.num.1500 <- apply(tau, 1, sd)
#se.analyt.1500 <- rowMeans(se)
#ci.hit.1500 <- rowMeans(ci.hit)

load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/sreg/adj/New version/2000.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.2000 <- rowMeans(tau)
se.num.2000 <- apply(tau, 1, sd)
se.analyt.2000 <- rowMeans(se)
ci.hit.2000 <- rowMeans(ci.hit)

load('/Users/trifonovjuri/Desktop/sreg.source/mc.files/res/sreg/adj/New version/3000.RData')
tau <- na.omit(as.matrix(sapply(simres, function(simres) simres$tau)))
se <- na.omit(as.matrix(sapply(simres, function(simres) simres$se)))
ci.hit <- na.omit(as.matrix(sapply(simres, function(simres) simres$ci.hit)))
tau.3000 <- rowMeans(tau)
se.num.3000 <- apply(tau, 1, sd)
se.analyt.3000 <- rowMeans(se)
ci.hit.3000 <- rowMeans(ci.hit)

tau_1.vec <- c(tau.100[1], tau.250[1], tau.500[1], tau.750[1],
            tau.1000[1], tau.2000[1], tau.3000[1])
tau_2.vec <- c(tau.100[2], tau.250[2], tau.500[2], tau.750[2],
               tau.1000[2], tau.2000[2], tau.3000[2])
se.num_1.vec <- c(se.num.100[1], se.num.250[1], se.num.500[1], se.num.750[1],
                se.num.1000[1], se.num.2000[1], se.num.3000[1])
se.num_2.vec <- c(se.num.100[2], se.num.250[2], se.num.500[2], se.num.750[2],
                  se.num.1000[2], se.num.2000[2], se.num.3000[2])
se.analyt_1.vec <- c(se.analyt.100[1], se.analyt.250[1], se.analyt.500[1], se.analyt.750[1], se.analyt.1000[1], se.analyt.2000[1], se.analyt.3000[1])
se.analyt_2.vec <- c(se.analyt.100[2], se.analyt.250[2], se.analyt.500[2], se.analyt.750[2], se.analyt.1000[2], se.analyt.2000[2], se.analyt.3000[2])
ci.hit_1.vec <- c(ci.hit.100[1], ci.hit.250[1], ci.hit.500[1], ci.hit.750[1], ci.hit.1000[1], ci.hit.2000[1], ci.hit.3000[1])
ci.hit_2.vec <- c(ci.hit.100[2], ci.hit.250[2], ci.hit.500[2], ci.hit.750[2], ci.hit.1000[2],  ci.hit.2000[2], ci.hit.3000[2])
sample.set <- c(100, 250, 500, 750, 1000, 2000, 3000)

graph.data <- data.frame('s.size' = sample.set,
                         'tau_1' = tau_1.vec, 'tau_2' = tau_2.vec,
                         'se.num_1' = se.num_1.vec, 'se.num_2' = se.num_2.vec,
                         'se.analyt_1' = se.analyt_1.vec, 'se.analyt_2' = se.analyt_2.vec,
                         'se.diff_1' = abs(se.num_1.vec - se.analyt_1.vec), 'se.diff_2' = abs(se.num_2.vec - se.analyt_2.vec),
                         'ci.hit_1' = ci.hit_1.vec,'ci.hit_2' = ci.hit_2.vec)


ggplot(graph.data) +
  geom_line(aes(x = s.size, y = se.diff_1, color = "SE1"), size = 1) +
  geom_point(aes(x = s.size, y = se.diff_1, color = "SE1"), size = 3) +
  geom_line(aes(x = s.size, y = se.diff_2, color = "SE2"), size = 1) +
  geom_point(aes(x = s.size, y = se.diff_2, color = "SE2"), size = 3) +
  labs(x = "Sample Size", y = "SE Difference") +
  ggtitle("Graph of Sample Size vs. SE Difference") +
  scale_color_manual(values = c("SE1" = "blue", "SE2" = "red")) +
  theme_minimal()

ggplot(graph.data) +
  geom_line(aes(x = s.size, y = se.num_1, color = "SE Num_1"), size = 1, linetype = "dashed") +
  geom_point(aes(x = s.size, y = se.num_1, color = "SE Num_1"), size = 3) +
  geom_line(aes(x = s.size, y = se.analyt_1, color = "SE Analyt_1"), size = 1) +
  geom_point(aes(x = s.size, y = se.analyt_1, color = "SE Analyt_1"), size = 3) +
  geom_line(aes(x = s.size, y = se.num_2, color = "SE Num_2"), size = 1, linetype = "dashed") +
  geom_point(aes(x = s.size, y = se.num_2, color = "SE Num_2"), size = 3) +
  geom_line(aes(x = s.size, y = se.analyt_2, color = "SE Analyt_2"), size = 1) +
  geom_point(aes(x = s.size, y = se.analyt_2, color = "SE Analyt_2"), size = 3) +
  labs(x = "Sample Size", y = "SE") +
  scale_color_manual(values = c("SE Num_1" = "blue", "SE Analyt_1" = "red", "SE Num_2" = "grey", "SE Analyt_2" = "orange")) +
  theme_minimal()

ggplot(graph.data) +
  geom_line(aes(x = s.size, y = se.num_2, color = "SE Num_2"), size = 1, linetype = "dashed") +
  geom_point(aes(x = s.size, y = se.num_2, color = "SE Num_2"), size = 3) +
  geom_line(aes(x = s.size, y = se.analyt_2, color = "SE Analyt_2"), size = 1) +
  geom_point(aes(x = s.size, y = se.analyt_2, color = "SE Analyt_2"), size = 3) +
  labs(x = "Sample Size", y = "SE") +
  scale_color_manual(values = c("SE Num_2" = "grey", "SE Analyt_2" = "orange")) +
  theme_minimal()

ggplot(graph.data) +
  geom_line(aes(x = s.size, y = ci.hit_1, color = "hit_1"), size = 1) +
  geom_point(aes(x = s.size, y = ci.hit_1, color = "hit_1"), size = 3) +
  geom_line(aes(x = s.size, y = ci.hit_2, color = "hit_2"), size = 1) +
  geom_point(aes(x = s.size, y = ci.hit_2, color = "hit_2"), size = 3) +
  labs(x = "Sample Size", y = "Coverage %") + scale_y_continuous(limits = c(0.6, 1))



