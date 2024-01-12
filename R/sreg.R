#%##%##%##%###%##%##%##%###%##%##%##%###%##%
### This R file provides the collection ####
### of functions to estimate the ATE    ####
### under CAR with multiple treatments  ####
#%##%##%##%###%##%##%##%###%##%##%##%###%##%
####      The code is developed by      ####
####      @Juri Trifonov, UChicago      ####
####            Supervisors:            ####
####      @Azeem Shaikh, UChicago       ####
####    @Max Tabord-Meehan, UChicago    ####
#%##%##%##%###%##%##%##%###%##%##%##%###%##%
#%##%##%##%##
#%# v.3.2 #%#
#%##%##%##%##
#-------------------------------------------------------------------
#%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%
#%##%##%##%##%##%#      I. ATEs estimation    #%##%##%##%##%##%##%##
#%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%##%
#-------------------------------------------------------------------
#-------------------------------------------------------------------
### (1) Function for automatically creating indicators for S or D
#%source function for full.satur()
#-------------------------------------------------------------------
indicator.sreg <- function(variable, type)
  #-------------------------------------------------------------------
{
  if (min(variable) == 0)
  {
    strg <- matrix(NA, nrow = length(variable), ncol = max(variable))            # matrix for storing I{D=d} or I{S=s}
    
    for (j in (min(variable) + 1) : max(variable)) 
    {
      strg[,j] <- as.numeric(variable == j)
    }
    ind.data <- as.data.frame(strg)
    if (type == 's')
    {
      col.names <- paste0("s_", min(variable) : max(variable))
    }
    if (type == 'd')
    {
      col.names <- paste0("d_", (min(variable) + 1) : max(variable))
    }
  } else{
    strg <- matrix(NA, nrow = length(variable), ncol = max(variable))
    for (j in min(variable):max(variable)) 
    {
      strg[,j] <- as.numeric(variable == j)
    }
    ind.data <- as.data.frame(strg)
    if (type == 's')
    {
      col.names <- paste0("s_", min(variable) : max(variable))
    }
    if (type == 'd')
    {
      col.names <- paste0("d_", min(variable) : max(variable))
    }
  }
  colnames(ind.data) <- col.names
  return(ind.data)
}


#-------------------------------------------------------------------
#%# (4) Auxiliary function providing the appropriate data.frame
#%#     for the subsequent iterative OLS estimation. Takes into account
#%#     the number of observations and creates indicators.
#%source function for theta.est.str()
# NB CHANGE FOR MULTIVARIATE!
#-------------------------------------------------------------------
filter.ols.sreg <- function(Y,S,D,X,s,d)
#-------------------------------------------------------------------
{
  X <- as.matrix(X)
  data <- data.frame(Y,S,D,X)
  keep.s <- s
  keep.d <- d
  filtered.data <- data[data$D %in% keep.d & data$S %in% keep.s, ]
  
  data.ols <- filtered.data
  return(data.ols)
}


#-------------------------------------------------------------------
#%# (4) Function that implements the OLS estimation of the 
#%#     fully-saturated regression via lm() with the generated
#%#     appropriate formula via full.satur()
#%source function for theta.est.str()
#-------------------------------------------------------------------
lm.iter.sreg <- function(Y,S,D,X)
#-------------------------------------------------------------------
{ 
  theta.list <- rep(list(matrix(NA, ncol = ncol(X), nrow = max(S))), (max(D) + 1))
 
  for (d in 0:max(D))
  {
    for (s in 1:max(S))
    {
      data.filtered <- filter.ols.sreg(Y,S,D,X,s,d)
      data.X <- data.filtered[, 4:(4 + ncol(X)-1)]
      data.filtered.adj <- data.frame(Y = data.filtered$Y, data.X)
      result <- lm(Y ~ ., data = data.filtered.adj)
      theta.list[[d+1]][s, ] <- coef(result)[2:(2 + ncol(X)-1)]
    }
  }
  list.rtrn <- theta.list
  return(list.rtrn)
}


#-------------------------------------------------------------------
lin.adj.sreg <- function(a, S, X, model)
#-------------------------------------------------------------------
{
  data <- data.frame(S,X)
  theta.mtrx <- model[[a+1]]
  theta.vec.matched <- theta.mtrx[data$S, ]
  mu.hat <- diag(as.matrix(X) %*% t(theta.vec.matched)) 
  return(mu.hat)
}


#-------------------------------------------------------------------
pi.hat.sreg <- function(S,D)
#-------------------------------------------------------------------
{
  n <- length(D)
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
tau.hat.sreg <- function(Y,S,D,X=NULL, model=NULL)
#-------------------------------------------------------------------
{
  tau.hat <- rep(NA, max(D))
  for (d in 1:max(D))
  {
    if(!is.null(X))
    {
    data <- data.frame(Y,S,D,X)
    data$pi <- pi.hat.sreg(S,D)[, d]
    data.bin <- data[data$D %in% c(d,0), ]
    data.bin$A <- as.numeric(data.bin$D != 0)
    
    mu.hat.d <- lin.adj.sreg(d, data.bin$S, data.bin[4:(4+ncol(X)-1)], model)
    
    mu.hat.0 <- lin.adj.sreg(0, data.bin$S, data.bin[4:(4+ncol(X)-1)], model)
    
    Ksi.vec <- ((data.bin$A * (data.bin$Y - mu.hat.d)) / data.bin$pi) - 
      (((1 - data.bin$A) * (data.bin$Y - mu.hat.0)) / (1 - data.bin$pi)) + 
      mu.hat.d - mu.hat.0
    
    tau.hat[d] <- mean(Ksi.vec)
    }else{
      data <- data.frame(Y,S,D)
      data$pi <- pi.hat.sreg(S,D)[, d]
      data.bin <- data[data$D %in% c(d,0), ]
      data.bin$A <- as.numeric(data.bin$D != 0)
      
      mu.hat.d <- 0
      mu.hat.0 <- 0
      
      Ksi.vec <- ((data.bin$A * (data.bin$Y - mu.hat.d)) / data.bin$pi) - 
        (((1 - data.bin$A) * (data.bin$Y - mu.hat.0)) / (1 - data.bin$pi)) + 
        mu.hat.d - mu.hat.0
      
      tau.hat[d] <- mean(Ksi.vec)
    }
    
  }
  return(tau.hat)
}


#-------------------------------------------------------------------
# Variance Estimator
#-------------------------------------------------------------------
 as.var.sreg <- function(Y,S,D,X=NULL, model=NULL, tau)
{
  var.vec <- rep(NA, max(D))
  n.vec <- rep(NA, max(D))
  
  if(!is.null(X))
  {
    for (d in 1:max(D))
    {
      data <- data.frame(Y,S,D,X)
      data$pi <- pi.hat.sreg(S,D)[, d]  
      data.bin <- data[data$D %in% c(d,0), ]
      data.bin$A <- as.numeric(data.bin$D != 0)
      n.d <- length(data.bin$Y)
      
      mu.hat.d <- lin.adj.sreg(d, data.bin$S, data.bin[4:(4+ncol(X)-1)], model)
      mu.hat.0 <- lin.adj.sreg(0, data.bin$S, data.bin[4:(4+ncol(X)-1)], model)
      data.bin <- data.frame(data.bin, mu.hat.d, mu.hat.0)
      
      Xi.tilde.1 <- (1 - (1/data.bin$pi)) * mu.hat.d - mu.hat.0 + 
        (data.bin$Y / data.bin$pi) 
      
      Xi.tilde.0 <- ((1 / (1 - data.bin$pi)) - 1) * mu.hat.0 + mu.hat.d -
        (data.bin$Y / (1 - data.bin$pi)) 
      
      data.bin <- data.frame(data.bin, Xi.tilde.1, Xi.tilde.0, Y.tau.D = data.bin$Y)
      
      Xi.1.mean <- rep(NA,n.d)
      Xi.0.mean <- rep(NA,n.d)
      Y.tau.D.1.mean <- rep(NA,n.d)
      Y.tau.D.0.mean <- rep(NA,n.d)
      
      for (i in 1:n.d)
      {
        Xi.1.mean[i] <- mean(data.bin[data.bin$A %in% 1 & data.bin$S %in% data.bin$S[i], ]$Xi.tilde.1)
        Xi.0.mean[i] <- mean(data.bin[data.bin$A %in% 0 & data.bin$S %in% data.bin$S[i], ]$Xi.tilde.0)
        Y.tau.D.1.mean[i] <- mean(data.bin[data.bin$A %in% 1 & data.bin$S %in% data.bin$S[i], ]$Y.tau.D)
        Y.tau.D.0.mean[i] <- mean(data.bin[data.bin$A %in% 0 & data.bin$S %in% data.bin$S[i], ]$Y.tau.D)
      }
      
      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean
      Xi.hat.2 <- Y.tau.D.1.mean - Y.tau.D.0.mean
      
      
      sigma.hat.sq <-  mean(data.bin$A * (Xi.hat.1)^2 + (1 - data.bin$A) * (Xi.hat.0)^2 + (Xi.hat.2)^2)
      
      var.vec[d] <- sigma.hat.sq
      n.vec[d]   <- n.d
      
    }
  }else{
    for (d in 1:max(D))
    {
      data <- data.frame(Y,S,D)
      data$pi <- pi.hat.sreg(S,D)[, d]  
      data.bin <- data[data$D %in% c(d,0), ]
      data.bin$A <- as.numeric(data.bin$D != 0)
      n.d <- length(data.bin$Y)
      
      mu.hat.d <- 0
      mu.hat.0 <- 0
      data.bin <- data.frame(data.bin, mu.hat.d, mu.hat.0)
      
      Xi.tilde.1 <- (1 - (1/data.bin$pi)) * mu.hat.d - mu.hat.0 + 
        (data.bin$Y / data.bin$pi) 
      
      Xi.tilde.0 <- ((1 / (1 - data.bin$pi)) - 1) * mu.hat.0 + mu.hat.d -
        (data.bin$Y / (1 - data.bin$pi)) 
      
      data.bin <- data.frame(data.bin, Xi.tilde.1, Xi.tilde.0, Y.tau.D = data.bin$Y) 
      
      Xi.1.mean <- rep(NA,n.d)
      Xi.0.mean <- rep(NA,n.d)
      Y.tau.D.1.mean <- rep(NA,n.d)
      Y.tau.D.0.mean <- rep(NA,n.d)
      
      for (i in 1:n.d)
      {
        Xi.1.mean[i] <- mean(data.bin[data.bin$A %in% 1 & data.bin$S %in% data.bin$S[i], ]$Xi.tilde.1)
        Xi.0.mean[i] <- mean(data.bin[data.bin$A %in% 0 & data.bin$S %in% data.bin$S[i], ]$Xi.tilde.0)
        Y.tau.D.1.mean[i] <- mean(data.bin[data.bin$A %in% 1 & data.bin$S %in% data.bin$S[i], ]$Y.tau.D)
        Y.tau.D.0.mean[i] <- mean(data.bin[data.bin$A %in% 0 & data.bin$S %in% data.bin$S[i], ]$Y.tau.D)
      }
      
      Xi.hat.1 <- Xi.tilde.1 - Xi.1.mean
      Xi.hat.0 <- Xi.tilde.0 - Xi.0.mean
      Xi.hat.2 <- Y.tau.D.1.mean - Y.tau.D.0.mean
      
      
      sigma.hat.sq <-  mean(data.bin$A * (Xi.hat.1)^2 + (1 - data.bin$A) * (Xi.hat.0)^2 + (Xi.hat.2)^2)
      
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
res.sreg <- function(Y,S,D,X)
  #-------------------------------------------------------------------
{
  n <- length(Y)
  if(!is.null(X))
  {
    model <- lm.iter.sreg(Y,S,D,X)
    tau.est <- tau.hat.sreg(Y,S,D,X,model)
    se.rob <- as.var.sreg(Y,S,D,X,model,tau.est)
    
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
                     'data' = data.frame(Y,S,D,X))
  }else{
    tau.est <- tau.hat.sreg(Y,S,D,X=NULL,model=NULL)
    se.rob <- as.var.sreg(Y,S,D,X=NULL,model=NULL,tau.est)
    
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
                     'data' = data.frame(Y,S,D))
    
  }
  
  return(res.list)
}

#-------------------------------------------------------------------
#%# (11) Summary method for sreg(). Provide the output table.
#-------------------------------------------------------------------
summary.sreg <- function(model)
  #-------------------------------------------------------------------
{
  n         <- length(model$data$Y)                     
  tau.hat   <- as.vector(model$tau.hat)
  se.rob    <- as.vector(model$se.rob)
  t.stat    <- as.vector(model$t.stat)
  p.value   <- as.vector(model$p.value)
  CI.left   <- as.vector(model$CI.left)
  CI.right  <- as.vector(model$CI.right)
  
  if(!is.null(model$data$x_1))
  {
  cat("Saturated Model Estimation Results under CAR with linear adjustments\n")                  
  }else{
    cat("Saturated Model Estimation Results under CAR\n")                  
  }
  
  cat(paste("Observations:", n, "\n"))                  
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
#%# (1) Potential outcomes generation
#-------------------------------------------------------------------
dgp.po.sreg <- function(n, theta.vec, gamma.vec, n.treat)
  #-------------------------------------------------------------------
{
  if (n.treat !=length(theta.vec))
  {
    return(print("The number of treatments doesn't
                 match the length of vector theta.vec"))
  }
  for (a in seq_along(theta.vec)) 
  {
    assign(paste("mu.", a, sep = ""), theta.vec[a])                              #create mu.a, where a \in \mathbb{A}
  }                             
  eps.0 <- rnorm(n)                                                              # create epsilons for every treatment a \in \mathbb{A}
  for (a in 1:n.treat)                              
  {                             
    assign(paste("eps.", a, sep = ""), rnorm(n))                             
  }                             
  
  W <- sqrt(20)*(rbeta(n,2,2)-1/2)       
  x_1 <- rnorm(n, mean = 5, sd = 2)
  x_2 <- rnorm(n, mean = 2, sd = 1)
  #x_3 <- rnorm(n, mean = 1, sd = 3)
  #x_4 <- rnorm(n, mean = 10, sd = 5)
  X <- data.frame(x_1, x_2)#, x_3, x_4)
  m.0 <- gamma.vec[1] * W + gamma.vec[2] * x_1 + gamma.vec[3] * x_2 #+ gamma.vec[4] * x_3 + gamma.vec[5] * x_4 
  #m.0 <- gamma.vec[1] * W
  for (a in 1:n.treat)                                                           # create m() functions m.a for every treatment a \in \mathbb{A}
  {
    assign(paste("m.", a, sep = ""), m.0)
  }
  
  Y.0 <- m.0 + eps.0
  
  for (a in 1:n.treat)
  {
    formula <- paste0("mu.", a, "+m.", a, "+eps.", a)
    result <- eval(parse(text = formula))
    assign(paste("Y.", a, sep = ""), result)
  }
  
  ret.names <- c(paste("Y.", 0:n.treat, sep = ""),
                 "W", "X", paste("m.", 0:n.treat, sep = ""),
                 paste("mu.", 1:n.treat, sep = ""))
  ret.list <- mget(ret.names)
  return (ret.list)
}
#-------------------------------------------------------------------
#%# (2) Random Treatment Assignment
#%source function for dgp.obs()
#-------------------------------------------------------------------
gen.treat.sreg <- function(pi.matr.w, ns, k) 
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
#%# (3) Generate the formula for Y.obs (Rubin model)
#%source function for dgp.obs()
#-------------------------------------------------------------------
gen.rubin.formula.sreg <- function(n.treat) {
  #-------------------------------------------------------------------
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
#%# (4) Generate observed outcomes,
#%#     by taking as input the potential outcomes,
#%#     matrix of strata assignments, pi.vec, and 
#%#     number of treatments
#-------------------------------------------------------------------
dgp.obs.sreg <- function(baseline, I.S, pi.vec, n.treat)
  #-------------------------------------------------------------------
{
  if (n.treat != length(pi.vec))
  {
    return(print("The number of treatments doesn't
                 match the length of vector pi.vec"))
  }
  num.strata <- ncol(I.S)
  n  <- length(baseline$Y.0)
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
    A[index]<- gen.treat.sreg(pi.matr.w, ns, k)
  }
  
  #now we need to generate observed outcomes via Rubin model
  for (a in 0:n.treat) 
  {
    assign(paste("Y.", a, sep = ""), baseline[[paste("Y.", a, sep = "")]])
  }
  formula <- gen.rubin.formula.sreg(n.treat)
  Y.obs <- eval(parse(text = formula))
  
  ret.list <- list('Y' = Y.obs,
                   'D' = A,
                   'X' = baseline$X)
  return(ret.list)
}
#-------------------------------------------------------------------
#%# Function taken from Ivan's website
#%# to generate the strata from the observed covariates
#%# NB! Works only if we form strata from one W.
#-------------------------------------------------------------------
form.strata.sreg <- function(baseline,num.strata)
  #-------------------------------------------------------------------
{
  #-------------------------------------------------------------------
  # X-plain: Generates strata indicators from covariates (W).
  #         - W: covariates (must be scalar)
  #         - num.strata: the number of strata to be formed
  #         - model: model in baseline.
  #-------------------------------------------------------------------
  n <- length(baseline$Y.0)
  W<-baseline$W;
  bounds<-seq(-2.25,2.25,length.out = num.strata+1); # careful with bounds
  I.S   <- matrix(0,n,num.strata);
  for (s in 1:num.strata)
  {
    I.S[,s]<- (W>bounds[s])*(W<=bounds[s+1]);
  }
  
  return(I.S)
}











