#%##%##%##%###%##%##%##%###%##%##%##%###%##%
### This R file provides the collection ####
### of functions to estimate the ATE    ####
### under CAR with multiple treatments  ####
###        & cluster-level treatment    ####
#%##%##%##%###%##%##%##%###%##%##%##%###%##%
####      The code is developed by      ####
####      @Juri Trifonov, UChicago      ####
####            Supervisors:            ####
####      @Azeem Shaikh, UChicago       ####
####    @Max Tabord-Meehan, UChicago    ####
#%##%##%##%###%##%##%##%###%##%##%##%###%##%
#%##%##%##%##
#%# v.4.6 #%#
#%##%##%##%##
#-------------------------------------------------------------------
#%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%
#%##%##%##%##%##%#      I. ATE  estimator     #%##%##%##%##%##%##%##
#%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%
#-------------------------------------------------------------------
#------------------------------------------------------------------
#%# (1) Auxiliary function providing the appropriate data.frame
#%#     for the subsequent iterative OLS estimation. Takes into account
#%#     the number of observations and creates indicators.
#%source function for theta.est.str()
# NB CHANGE FOR MULTIVARIATE!
#-------------------------------------------------------------------
filter.ols.creg <- function(Y,S,D,G.id,Ng,X,s,d)
  #-------------------------------------------------------------------
{
  working.df <- data.frame(Y,S,D,G.id,Ng,X)
  Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
  cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", 'Ng')]) # created data on a cluster level for estimating pi.hat(s)
  X.unique <- unique(working.df[, 6:ncol(working.df)])
  cl.lvl.data <- data.frame(cl.lvl.data, 'Y.bar' = Y.bar.g$Y, X.unique)
  data <- cl.lvl.data
  
  keep.s <- s
  keep.d <- d
  filtered.data <- data[data$D %in% keep.d & data$S %in% keep.s, ]
  
  data.ols <- filtered.data
  return(data.ols)
}
#-------------------------------------------------------------------
lm.iter.creg <- function(Y,S,D,G.id,Ng,X,exp.option =FALSE)
  #-------------------------------------------------------------------
{ 
 
  theta.list <- rep(list(matrix(NA, ncol = ncol(X) + 1, nrow = max(S))), (max(D) + 1))
  if (exp.option == TRUE)
  {
    theta.list <- rep(list(matrix(NA, ncol = ncol(X), nrow = max(S))), (max(D) + 1))
  }
  
  for (d in 0:max(D))
  {
    for (s in 1:max(S))
    {
      data.filtered <- filter.ols.creg(Y,S,D,G.id,Ng,X,s,d)
      data.X <- data.filtered[, 6:(6 + ncol(X)-1)]
      data.filtered.adj <- data.frame(Y.bar.Ng = data.filtered$Y.bar * data.filtered$Ng, Ng = data.filtered$Ng, data.X)
      if (exp.option == TRUE)
      {
        data.filtered.adj <- data.frame(Y.bar.Ng = data.filtered$Y.bar * data.filtered$Ng, data.X)
      }
      result <- lm(Y.bar.Ng ~ ., data = data.filtered.adj)
      
      if (exp.option == TRUE)
      {
        theta.list[[d+1]][s, ] <- coef(result)[2:(2 + ncol(X)-1)]
      }
      if (exp.option == FALSE)
      {
        theta.list[[d+1]][s, ] <- coef(result)[2:(2 + ncol(X))]
      }
    }
  }
  list.rtrn <- theta.list
  return(list.rtrn)
}
#-------------------------------------------------------------------
lin.adj.creg <- function(a,data,model, exp.option = FALSE)
  #-------------------------------------------------------------------
{
  working.df <- data
  Y.bar.g <- aggregate(Y ~ G.id, working.df, mean)
  cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", 'Ng')]) # created data on a cluster level for estimating pi.hat(s)
  X.unique <- unique(working.df[, 6:ncol(working.df)])
  cl.lvl.data <- data.frame(cl.lvl.data, 'Y.bar' = Y.bar.g$Y, X.unique)
  data <- cl.lvl.data
  X.data <- cl.lvl.data[, 6:ncol(cl.lvl.data)]
  
  theta.mtrx <- model[[a+1]]
  theta.vec.matched <- theta.mtrx[data$S, ]
  if (exp.option == FALSE)
  {
    Ng.hat <- theta.vec.matched[, 1] * data$Ng
    X.hat <- diag(as.matrix(X.data) %*% t(theta.vec.matched[, -1]))
    mu.hat <- Ng.hat + X.hat 
  }else{
    X.hat <- diag(as.matrix(X.data) %*% t(theta.vec.matched))
    mu.hat <-X.hat 
  }
  
  return(mu.hat)
}
#-------------------------------------------------------------------
pi.hat.creg <- function(cl.lvl.data)
  #-------------------------------------------------------------------
{
  S <- cl.lvl.data$S
  D <- cl.lvl.data$D
  n <- length(S)
  data <- data.frame(S,D)
  
  pi.hat.mtrx <- matrix(NA, nrow = n, ncol = max(D))
  for (d in 1:max(D))
  {
    for (i in 1:n)
    {
      n.1.s <-  length(data[data$D %in% d & data$S %in% data$S[i], 2])
      n.0.s <-  length(data[data$D %in% 0 & data$S %in% data$S[i], 2])
      n.s <- n.1.s + n.0.s
      pi.hat.mtrx[i,d] <- n.1.s / n.s
    }
  }
  return(pi.hat.mtrx)
}
#-------------------------------------------------------------------
tau.hat.creg <- function(Y,S,D,G.id,Ng,X=NULL,model=NULL, exp.option = FALSE)
  #-------------------------------------------------------------------
{
  tau.hat.vec <- rep(NA, max(D))
  Y.bar.g.list <- rep(list(NA),max(D))
  mu.hat.list <- rep(list(NA),max(D))
  pi.hat.list <- rep(list(NA),max(D))
  data.bin.list <- rep(list(NA),max(D))
  
  if(!is.null(X))
  {
    working.df <- data.frame(Y,S,D,G.id,Ng,X)
   
    cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", 'Ng')]) # created data on a cluster level for estimating pi.hat(s)
    
    for (d in 1:max(D))
    {
      data.Y.bar <- working.df[working.df$D %in% c(d,0), ]
      Y.bar.g <- aggregate(Y ~ G.id, data.Y.bar, mean)
      Y.bar.g.list[[d]] <- Y.bar.g$Y
      
      data <- cl.lvl.data
      data$pi <- pi.hat.creg(cl.lvl.data)[, d]
      data.bin <- data[data$D %in% c(d,0), ]
      data.bin$A <- as.numeric(data.bin$D != 0)
      data.bin.mu <- working.df[working.df$D %in% c(d,0), ]
      pi.hat.list[[d]] <- data.bin$pi
      data.bin.list[[d]] <- data.bin
      
      
      
      mu.hat.d <- lin.adj.creg(d, data = data.bin.mu, model, exp.option = exp.option)
      
      mu.hat.0 <- lin.adj.creg(0, data = data.bin.mu, model, exp.option = exp.option)
      
      Xi.g <- ((data.bin$A * (Y.bar.g$Y * data.bin$Ng - mu.hat.d)) / data.bin$pi) - (((1 - data.bin$A) * (Y.bar.g$Y * data.bin$Ng - mu.hat.0)) / (1 - data.bin$pi)) + mu.hat.d - mu.hat.0
      
      mu.hat.list[[d]] <- as.matrix(cbind(mu.hat.0,mu.hat.d), ncol = 2)
      
      Ng.ind <- data.bin$Ng
      tau.hat <- sum(Xi.g) / sum(Ng.ind)
      
      tau.hat.vec[d] <- tau.hat
    }
    rtrn.list <- list('tau.hat' = tau.hat.vec,
                      'mu.hat' = mu.hat.list,
                      'pi.hat' = pi.hat.list,
                      'Y.bar.g' = Y.bar.g.list,
                      'data.bin' = data.bin.list)
  }else{
    working.df <- data.frame(Y,S,D,G.id,Ng)
    cl.lvl.data <- unique(working.df[, c("G.id", "D", "S", 'Ng')]) # created data on a cluster level for estimating pi.hat(s)
    
    for (d in 1:max(D))
    {
      data.Y.bar <- working.df[working.df$D %in% c(d,0), ]
      Y.bar.g <- aggregate(Y ~ G.id, data.Y.bar, mean)
      Y.bar.g.list[[d]] <- Y.bar.g$Y
      
      data <- cl.lvl.data
      data$pi <- pi.hat.creg(cl.lvl.data)[, d]
      data.bin <- data[data$D %in% c(d,0), ]
      data.bin$A <- as.numeric(data.bin$D != 0)
      data.bin.mu <- working.df[working.df$D %in% c(d,0), ]
      pi.hat.list[[d]] <- data.bin$pi
      data.bin.list[[d]] <- data.bin
      
      mu.hat.d <- 0
      mu.hat.0 <- 0
      
      Xi.g <- ((data.bin$A * (Y.bar.g$Y * data.bin$Ng - mu.hat.d)) / data.bin$pi) - 
              (((1 - data.bin$A) * (Y.bar.g$Y * data.bin$Ng - mu.hat.0)) / (1 - data.bin$pi)) + 
              mu.hat.d - mu.hat.0
      
      mu.hat.list[[d]] <- as.matrix(cbind(mu.hat.0,mu.hat.d), ncol = 2)
      
      Ng.ind <- data.bin$Ng
      tau.hat <- sum(Xi.g) / sum(Ng.ind)  
      
      tau.hat.vec[d] <- tau.hat
    }
    rtrn.list <- list('tau.hat' = tau.hat.vec,
                      'mu.hat' = mu.hat.list,
                      'pi.hat' = pi.hat.list,
                      'Y.bar.g' = Y.bar.g.list,
                      'data.bin' = data.bin.list)
  }
  return(rtrn.list)
}

#-------------------------------------------------------------------
#%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%
#%##%##%##%##%##%#    II. Variance  estimator     #%##%##%##%##%##%#
#%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# Variance Estimator
#-------------------------------------------------------------------
as.var.creg <- function(model=NULL, fit)
{
  var.vec <- rep(NA, length(fit$tau.hat))
  n.vec <- rep(NA, length(fit$tau.hat))
  
  if(!is.null(X))
  {
    for (d in 1:length(fit$tau.hat))
    {
      Y.bar.g <- fit$Y.bar.g[[d]]
      mu.hat.0 <- fit$mu.hat[[d]][,1]
      mu.hat.d <- fit$mu.hat[[d]][,2]
      pi.hat <- fit$pi.hat[[d]]
      tau.est <- fit$tau.hat
      
      data.filter <- fit$data.bin[[d]]
      
      n.d <- length(data.filter$G.id)

      Xi.tilde.1 <- (1 - (1/pi.hat)) * mu.hat.d - mu.hat.0 + 
        (data.filter$Ng * Y.bar.g / pi.hat) - tau.est[d] * data.filter$Ng
      
      Xi.tilde.0 <- ((1 / (1 - pi.hat)) - 1) * mu.hat.0 + mu.hat.d -
        (data.filter$Ng * Y.bar.g / (1 - pi.hat)) - tau.est[d] * data.filter$Ng
      
      data.bin <- data.frame(data.filter, Xi.tilde.1, Xi.tilde.0, Y.tau.D = Y.bar.g * data.filter$Ng - tau.est[d] * data.filter$Ng)
  
      n.d <- length(data.bin$G.id)
      Ng.d <-data.bin$Ng
      
      Xi.1.mean <- rep(NA,n.d)
      Xi.0.mean <- rep(NA,n.d)
      Y.g.bar.cl.1 <- rep(NA,n.d)
      Y.g.bar.cl.0 <- rep(NA,n.d)
      Y.g.mean.1 <- rep(NA,n.d)
      Y.g.mean.0 <- rep(NA,n.d)
      
      for (i in 1:n.d)
      {
        Xi.1.mean[i] <- mean(data.bin[data.bin$A %in% 1 & data.bin$S %in% data.bin$S[i], ]$Xi.tilde.1)
        Xi.0.mean[i] <- mean(data.bin[data.bin$A %in% 0 & data.bin$S %in% data.bin$S[i], ]$Xi.tilde.0)
        Y.g.bar.cl.1[i] <- mean(data.bin[data.bin$A %in% 1 & data.bin$S %in% data.bin$S[i], ]$Y.tau.D)
        Y.g.bar.cl.0[i] <- mean(data.bin[data.bin$A %in% 0 & data.bin$S %in% data.bin$S[i], ]$Y.tau.D)
      }
      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean
      Xi.hat.2 <- Y.g.bar.cl.1 - Y.g.bar.cl.0
      
      sigma.hat.sq <-  mean((data.bin$A * (Xi.hat.1)^2 + (1 - data.bin$A) * (Xi.hat.0)^2 + (Xi.hat.2)^2)) /  (mean(Ng.d))^2
      
      var.vec[d] <- sigma.hat.sq
      n.vec[d]   <- n.d
    }
  }else{
    for (d in 1:length(fit$tau.hat))
    {
      Y.bar.g <- fit$Y.bar.g[[d]]
      pi.hat <- fit$pi.hat[[d]]
      tau.est <- fit$tau.hat
      
      data.filter <- fit$data.bin[[d]]
      
      n.d <- length(data.filter$G.id)
      
      mu.hat.0 <- 0
      mu.hat.d <- 0

      Xi.tilde.1 <- (1 - (1/pi.hat)) * mu.hat.d - mu.hat.0 + 
        (data.filter$Ng * Y.bar.g / pi.hat) - tau.est[d] * data.filter$Ng
      
      Xi.tilde.0 <- ((1 / (1 - pi.hat)) - 1) * mu.hat.0 + mu.hat.d -
        (data.filter$Ng * Y.bar.g / (1 - pi.hat)) - tau.est[d] * data.filter$Ng
      
      data.bin <- data.frame(data.filter, Xi.tilde.1, Xi.tilde.0, Y.tau.D = Y.bar.g * data.filter$Ng - tau.est[d] * data.filter$Ng)
      
      n.d <- length(data.bin$G.id)
      Ng.d <-data.bin$Ng
      
      Xi.1.mean <- rep(NA,n.d)
      Xi.0.mean <- rep(NA,n.d)
      Y.g.bar.cl.1 <- rep(NA,n.d)
      Y.g.bar.cl.0 <- rep(NA,n.d)
      Y.g.mean.1 <- rep(NA,n.d)
      Y.g.mean.0 <- rep(NA,n.d)
      
      for (i in 1:n.d)
      {
        Xi.1.mean[i] <- mean(data.bin[data.bin$A %in% 1 & data.bin$S %in% data.bin$S[i], ]$Xi.tilde.1)
        Xi.0.mean[i] <- mean(data.bin[data.bin$A %in% 0 & data.bin$S %in% data.bin$S[i], ]$Xi.tilde.0)
        Y.g.bar.cl.1[i] <- mean(data.bin[data.bin$A %in% 1 & data.bin$S %in% data.bin$S[i], ]$Y.tau.D)
        Y.g.bar.cl.0[i] <- mean(data.bin[data.bin$A %in% 0 & data.bin$S %in% data.bin$S[i], ]$Y.tau.D)
      }
      
      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean
      Xi.hat.2 <- Y.g.bar.cl.1 - Y.g.bar.cl.0
      
      sigma.hat.sq <-  mean((data.bin$A * (Xi.hat.1)^2 + (1 - data.bin$A) * (Xi.hat.0)^2 + (Xi.hat.2)^2)) /  (mean(Ng.d))^2
      
      var.vec[d] <- sigma.hat.sq
      n.vec[d]   <- n.d
    }
    
  }
  se.vec <- sqrt(var.vec/n.vec)
  return(se.vec)
}

#-------------------------------------------------------------------
#%# (10) The core function. It provides estimates of ATE, their s.e., 
#%#     calculates t-stats and corresponding p-values
#-------------------------------------------------------------------
res.creg <- function(Y,S,D,G.id,Ng,X,model, exp.option = FALSE)
  #-------------------------------------------------------------------
{
  n <- length(Y)
  if(!is.null(X))
  {
    if (exp.option == FALSE)
    {
      model <- lm.iter.creg(Y,S,D,G.id,Ng,X, exp.option = FALSE)
      fit <- tau.hat.creg(Y,S,D,G.id,Ng,X,model, exp.option = FALSE)
    }else{
      model <- lm.iter.creg(Y,S,D,G.id,Ng,X, exp.option = TRUE)
      fit <- tau.hat.creg(Y,S,D,G.id,Ng,X,model, exp.option = TRUE)
    }
    tau.est <- fit$tau.hat
    se.rob <- as.var.creg(model,fit)
    t.stat <- tau.est / se.rob
    p.value <- 2 * pmin(pnorm(t.stat), 1 - pnorm(t.stat))
    CI.left <- tau.est - qnorm(0.975) * se.rob
    CI.right <- tau.est + qnorm(0.975) * se.rob
    res.list <- list('tau.hat' = tau.est,
                     'se.rob'  = se.rob,
                     't.stat'  = t.stat,
                     'p.value' = p.value,
                     'as.CI'   = c(CI.left, CI.right),
                     'CI.left' = CI.left,
                     'CI.right'= CI.right,
                     'data' = data.frame(Y,S,D,G.id,Ng,X))
  }else{
    if (exp.option == FALSE)
    {
      fit <- tau.hat.creg(Y,S,D,G.id,Ng,X=NULL,model=NULL, exp.option = FALSE)
    }else{
      fit <- tau.hat.creg(Y,S,D,G.id,Ng,X=NULL,model=NULL, exp.option = TRUE)
    }
    tau.est <- fit$tau.hat
    se.rob <- as.var.creg(model=NULL,fit)
    t.stat <- tau.est / se.rob
    p.value <- 2 * pmin(pnorm(t.stat), 1 - pnorm(t.stat))
    CI.left <- tau.est - qnorm(0.975) * se.rob
    CI.right <- tau.est + qnorm(0.975) * se.rob
    res.list <- list('tau.hat' = tau.est,
                     'se.rob'  = se.rob,
                     't.stat'  = t.stat,
                     'p.value' = p.value,
                     'as.CI'   = c(CI.left, CI.right),
                     'CI.left' = CI.left,
                     'CI.right'= CI.right,
                     'data' = data.frame(Y,S,D,G.id,Ng))
  }
  return(res.list)
}

#-------------------------------------------------------------------
#%# (11) Summary method for sreg(). Provide the output table.
#-------------------------------------------------------------------
summary.creg <- function(model)
  #-------------------------------------------------------------------
{
  n         <- length(model$data$Y)
  G         <- length(unique(model$data$G.id))
  tau.hat   <- as.vector(model$tau.hat)
  se.rob    <- as.vector(model$se.rob)
  t.stat    <- as.vector(model$t.stat)
  p.value   <- as.vector(model$p.value)
  CI.left   <- as.vector(model$CI.left)
  CI.right  <- as.vector(model$CI.right)
  
  if(!is.null(model$data$x_1))
  {
    cat("Saturated Model Estimation Results under CAR with clusters and linear adjustments\n")                  
  }else{
    cat("Saturated Model Estimation Results under CAR with clusters\n")    
  }
  cat(paste("Observations:", n, "\n"))
  cat(paste("Clusters:", G, "\n"))  
  cat(paste("Number of treatments:",                          
            max(model$data$D), "\n"))  
  cat(paste("Number of strata:",                           
            max(model$data$S), "\n"))  
  
  cat("---\n")                                           
  
  cat("Coefficients:\n")                               
  
  m <- length(tau.hat)  
  
  stars <- rep("", m)                                    
  stars[p.value <= 0.001] <- "***"
  stars[(p.value > 0.001) & (p.value < 0.01)] <- "**"
  stars[(p.value > 0.01) & (p.value <= 0.05)] <- "*"
  stars[(p.value > 0.05) & (p.value <= 0.1)] <- "."
  
  df <- data.frame("Tau" = tau.hat,                  
                   "As.se" = se.rob,
                   "T-stat" = t.stat,
                   "P-value" = p.value,
                   "CI.left" = CI.left,
                   "CI.right" = CI.right,
                   "Significance" = stars)
  is.df.num.col <- sapply(df, is.numeric)                
  df[, is.df.num.col] <- round(df[, is.df.num.col], 5)   
  print(df)                                              
  cat("---\n")                                           
  cat(paste("Signif. codes:  0 ‘***’ 0.001 ‘**’",       
            "0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n"))                      
}

#-------------------------------------------------------------------------------
#%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%#%##%##%##%#
#%##%##%##%##%##%## III. DGP functions for simulations #%##%##%##%##%##%##%##%#%
#%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%#%##%##%##%#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------
#%# (1) Cluster sizes generation
#-------------------------------------------------------------------
gen.cluster.sizes <-function(G,max.support,s=1.5)
#-------------------------------------------------------------------
{
  #-------------------------------------------------------------------
  # X-plain: Draws G clusters with a maximum sample size Nmax given Z
  #-------------------------------------------------------------------
  # INPUTS: - G: the number of clusters
  #         - max.support: n-1 in the beta binomial | Nmax = 10*(max.support +1)  
  #         - s is the parameter in the zeta distribution (s-1 in wiki notation)
  #------------------------------------------------------------------
  # RETURNS: 4 designs for sample sizes for G clusters. 
  #          The output is a (G x 4) matrix 
  #------------------------------------------------------------------
  # Setup sample sizes per cluster
  
  # Design 1: Beta binomial uniform
  sample.1 = 10*(rbbinom(G,max.support, alpha = 1, beta = 1)+1);
  sample.2 = 10*(rbbinom(G,max.support, alpha = 0.4, beta = 0.4)+1);
  sample.3 = 10*(rbbinom(G,max.support, alpha = 10, beta = 100)+1);
  #sample.4 = rzeta(G,s)*10;
  
  samples = matrix(c(sample.1,sample.2,sample.3),G,3)
  
  # Return the data frame
  return(samples);
}
#-------------------------------------------------------------------
#%# (2) Potential outcomes generation
#-------------------------------------------------------------------
dgp.po.creg <- function(Ng, G, tau.vec, sigma1=sqrt(2), 
                         gamma.vec = c(0.4, 0.2, 1, 0.1, 0.8), n.treat)
  #------------------------------------------------------------------
{
  x_1 <- rbeta(G, 2, 4)
  x_2 <- rbeta(G, 2, 4)
  X <- data.frame(x_1, x_2)
  beta.rv <- rbeta(G,2,2)
  Z.g.2 <- x_1
  mu.0 <- 10 * (x_1 - (1/3)) + 6 * (x_2 - (1/3)) + 2
  mu.1 <- mu.0
  
  cluster.indicator = rep(c(1:G),Ng);
  cl.id <-  cluster.indicator
  total.sample = length(cluster.indicator)
  
  epsilon.ig.0 = rnorm(total.sample, 0, sigma1);
  epsilon.ig.1 = rnorm(total.sample, 0, sigma1);
  
  Yig.0 <- rep(mu.0, Ng) + 2 * epsilon.ig.0
  Yig.1 <- rep(mu.1, Ng) + 2 * epsilon.ig.1
  
  ret.names <- c(paste("Yig.", 0:n.treat, sep = ""),
                 "X", "G", "Ng", "cl.id", "Z.g.2", paste("mu.", 0:n.treat, sep = ""))
  
  ret.list <- mget(ret.names)
  return (ret.list)
}

#-------------------------------------------------------------------
#%# (3) Random Treatment Assignment
#%source function for dgp.obs()
#-------------------------------------------------------------------
gen.treat.creg <- function(pi.matr.w, ns, k) 
#-------------------------------------------------------------------
{
  rows <- nrow(pi.matr.w)
  code.elements <- character(rows + 1)
  
  for (i in 1:rows) 
  {
    code.elements[i] <- paste0("rep(", i,
                               ", floor(pi.matr.w[", i, ",", k, "]*ns))")
  }
  
  code.elements[rows + 1] <- paste0("rep(0, (ns - ",
                                    paste0("floor(pi.matr.w[", 1:rows,
                                           ",",  k, "]*ns)",
                                           collapse = " - "), "))")
  
  code <- paste(code.elements, collapse = ", ")
  
  result <- eval(parse(text = paste("sample(c(", code, "))")))
  
  return(result)
}

#-------------------------------------------------------------------
#%# (4) Generate the formula for Y.obs (Rubin model)
#%source function for dgp.obs()
#-------------------------------------------------------------------
gen.rubin.formula.creg <- function(n.treat) 
#-------------------------------------------------------------------
{
  # Create a sequence of A values from 0 to max.A
  A.values <- 0:n.treat
  
  # Initialize an empty formula string
  formula <- "Y.obs = "
  
  # Generate the formula dynamically with indicators
  for (a in A.values) 
  {
    if (a == 0) 
    {
      formula <- paste(formula, paste0("Y.", a, " * (A == 0)"))
    } else 
    {
      formula <- paste(formula, paste0("Y.", a, " * (A == ", a, ")"))
    }
    
    if (a < n.treat) 
    {
      formula <- paste(formula, " + ")
    }
  }
  return(formula)
}

#-------------------------------------------------------------------
#%# (5) Generate observed outcomes,
#%#     by taking as input the potential outcomes,
#%#     matrix of strata assignments, pi.vec, and 
#%#     number of treatments
#-------------------------------------------------------------------
dgp.obs.creg <- function(baseline, I.S, pi.vec, n.treat)
  #-------------------------------------------------------------------
{
  if (n.treat != length(pi.vec))
  {
    return(print("The number of treatments doesn't
                 match the length of vector pi.vec"))
  }
  num.strata <- ncol(I.S)
  n  <- baseline$G
  A  <- cbind(rep(0,n))  # treatment Assignment
  l.seq<-num.strata/2
  
  pi.matr <- matrix(1, ncol = num.strata, nrow = n.treat)
  #pi.vec <- rep(c(1 / (n.treat + 1)), n.treat)
  pi.matr.w <- pi.matr * pi.vec
  
  for (k in 1:num.strata)
  {
    index   <- which(I.S[,k]==1)
    ns      <- length(index)
    
    # pick a random permutation of elements in \mathbb{A} and 0 
    A[index]<- gen.treat.creg(pi.matr.w, ns, k)
  }
  strata.set <- data.frame(I.S)
  strata.set$S <- max.col(strata.set)
  cluster.indicator <- baseline$cl.id 
  G.seq <- seq(c(1:baseline$G))
  data.short <- data.frame('cl.id'=G.seq, A, S = strata.set$S, Ng = baseline$Ng, 
                           baseline$X)
  data.long <- data.frame('cl.id'= cluster.indicator)
  merged.data <- merge(data.long, data.short, by = "cl.id")
  length(merged.data$A)
  A <- merged.data$A
  S <- merged.data$S
  X <- merged.data[5:ncol(merged.data)]
  Ng <- merged.data$Ng
  #now we need to generate observed outcomes via Rubin model
  for (a in 0:n.treat) 
  {
    assign(paste("Y.", a, sep = ""), baseline[[paste("Yig.", a, sep = "")]])
  }
  formula <- gen.rubin.formula.creg(n.treat)
  Y.obs <- eval(parse(text = formula))
  
  ret.list <- list('Y' = Y.obs,
                   'D' = A,
                   'S' = S,
                   'Z.2' = baseline$Z.g.2,
                   'X' = X,
                   'Ng' = Ng,
                   'G.id' = cluster.indicator, 
                   'cl.lvl.data' = data.short)
  return(ret.list)
}

#---------------------------------------------------------------------
#%# (6) Function taken from Ivan's website
#%# to generate the strata from the observed covariates
#%# NB! Works only if we form strata from one W.
#---------------------------------------------------------------------
form.strata.creg <- function(baseline,num.strata)
  #---------------------------------------------------------------------
{
  #-------------------------------------------------------------------
  # X-plain: Generates strata indicators from covariates (W).
  #         - W: covariates (must be scalar)
  #         - num.strata: the number of strata to be formed
  #         - model: model in baseline.
  #-------------------------------------------------------------------
  n <- baseline$G
  W<-baseline$Z.g.2;
  bounds<-seq(min(W),max(W),length.out = num.strata+1); # careful with bounds
  I.S   <- matrix(0,n,num.strata);
  for (s in 1:num.strata)
  {
    I.S[,s]<- (W>bounds[s])*(W<=bounds[s+1]);
  }
  
  return(I.S)
}




