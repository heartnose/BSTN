rm(list = ls())
source("./functions_tensor.R") # load required functions

# BSTN (Bayesian Skewed Tensor Normal)
#-------------------------------------------------#
#   Apply BSTN to the simulated data set          #
#-------------------------------------------------#

load('simulated data.Rda') # load simulated data 
vecy <- data[[1]]; Y <- data[[2]]; X <- data[[3]];

BSTN_matrix <- function(Y,X,vecy, n.burn = 10, n.save = 100, thin = 1){
  
  #Inputs:
  #y = teeth x site x nsubs responses
  #X = p x nsubs matrix of subject-level covariates 
  #vecy = (teeth x site) x nsubs matrix
  #n.burn = number to discard as burn-in
  #n.save = number of MCMC samples
  #thin = number of iterations to skip before collecting the next posterior sample
  
  t = dim(Y)[1]; s = dim(Y)[2]; p = dim(X)[1]; n = dim(X)[2];
  
  # Store MCMC out
  
  B.est.save <- array(NA, c(p, t*s, n.save))
  W.save <-array(NA, c(t*s, n, n.save))
  rho.save <- matrix(NA, 2, n.save)
  sigma.sq.save <- matrix(NA, 1, n.save)
  lam.est.save <- matrix(NA, 1, n.save)
  
  #initial values
  
  BOLS<-solve(X%*%t(X))%*%X%*%t(vecy) # vectorized OLS estimator
  
  B.est <- BOLS  # B_{(3)}
  
  W <- matrix(1, t*s,n) # vecotrize W = |Z_{2i}|
  
  rho <- matrix(c(0.5,0.5), 2, 1)
  
  sigma.sq <- 1.5
  
  lam.est <- 1.5
  
  source("./functions_tensor.R")
  
  #load required packages
  library(expm)   # calculate square root of a matrix 
  library(Matrix)
  library(matrixcalc) 
  library(LaplacesDemon) 
  library(MASS) # sample from multivariate normal distribution
  library(msm)  # sample from truncated normal distribution
  library(truncnorm) # samaple from truncated univariate multivariate normal distribution
  library(abind)
  library(magrittr)
  library(doParallel)
  library(TruncatedNormal)
  registerDoParallel(cores=2)  
  
  #---------------------------------------------------------------------------------
  
  # Update rho1 & rho2 (M-H algorithm) using closed form of equicorrelation assumption for R1 & R2
  
  rho.update <- function(t,s,rho,vecy,X,B.est,lam.est,W,sigma.sq){
    
    rho.prop = c(rbeta(1,2,2), rbeta(1,2,2)) #(rho1.prop, rho2.prop)
    
    det.curr = (sigma.sq^(t*s))*((((1-rho[2])^(s-1))*(1 + (s-1)*rho[2]))^t)*((((1-rho[1])^(t-1))*(1 + (t-1)*rho[1]))^s)
    det.prop = (sigma.sq^(t*s))*((((1-rho.prop[2])^(s-1))*(1 + (s-1)*rho.prop[2]))^t)*((((1-rho.prop[1])^(t-1))*(1 + (t-1)*rho.prop[1]))^s)
    
    inv.curr = (1/(sigma.sq))*kron(((1/(1-rho[2]))*(diag(s) - (rho[2]/ (1 + (s-1)*rho[2]))*matrix(1,s,s))), ((1/(1-rho[1]))*(diag(t) - (rho[1]/ (1 + (t-1)*rho[1]))*matrix(1,t,t)))) # (1/sigma^2)*{R_{\rho_2}^{-1} \otimes R_{\rho_1}^{-1}}
    inv.prop = (1/(sigma.sq))*kron(((1/(1-rho.prop[2]))*(diag(s) - (rho.prop[2]/ (1 + (s-1)*rho.prop[2]))*matrix(1,s,s))), ((1/(1-rho.prop[1]))*(diag(t) - (rho.prop[1]/ (1 + (t-1)*rho.prop[1]))*matrix(1,t,t))))
    
    inv.R.curr = kron(((1/(1-rho[2]))*(diag(s) - (rho[2]/ (1 + (s-1)*rho[2]))*matrix(1,s,s))), ((1/(1-rho[1]))*(diag(t) - (rho[1]/ (1 + (t-1)*rho[1]))*matrix(1,t,t))))
    inv.R.prop = kron(((1/(1-rho.prop[2]))*(diag(s) - (rho.prop[2]/ (1 + (s-1)*rho.prop[2]))*matrix(1,s,s))), ((1/(1-rho.prop[1]))*(diag(t) - (rho.prop[1]/ (1 + (t-1)*rho.prop[1]))*matrix(1,t,t))))
    
    logdens.curr = -0.5*n*log(det.curr) - foreach(l = 1:n,.combine='+') %dopar% {0.5*((t(vecy)[l,] - t(X)[l,]%*%B.est - lam.est*t(W[,l]))%*%inv.curr%*%(vecy[,l] - t(B.est)%*%X[,l] - lam.est*W[,l]))}
    logdens.prop = -0.5*n*log(det.prop) - foreach(l = 1:n,.combine='+') %dopar% {0.5*((t(vecy)[l,] - t(X)[l,]%*%B.est - lam.est*t(W[,l]))%*%inv.prop%*%(vecy[,l] - t(B.est)%*%X[,l] - lam.est*W[,l]))}
    
    logratio = logdens.prop - logdens.curr
    if(log(runif(1)) > logratio) {rho = rho} else {rho = rho.prop}
    if(log(runif(1)) > logratio) {inv.Sigma = inv.curr} else {inv.Sigma = inv.prop}
    if(log(runif(1)) > logratio) {inv.R = inv.R.curr} else {inv.R = inv.R.prop}
    
    return(list(rho, inv.Sigma, inv.R))
    
  }
  
  #---------------------------------------------------------------------------------
  # Update sigma.sq
  
  sigma.sq.update <- function (t,s,n,vecy,X,B.est,lam.est,W,inv.R){
    
    g1 = 2; g2 = 2; #prior distribution for inv.sigma.sq ~ Ga(2,2)
    Sww <- foreach(l = 1:n,.combine='+') %dopar% {t(W[,l])%*%W[,l]}
    S <- foreach(l = 1:n,.combine='+') %dopar% {((t(vecy)[l,] - t(X)[l,]%*%B.est - lam.est*t(W[,l]))%*%inv.R%*%(vecy[,l] - t(B.est)%*%X[,l] - lam.est*W[,l]))}
    inv.sigma.sq <- rgamma(1, g1 + n*t*s, 0.5*S + 0.5*Sww + g2)
    sigma.sq <- 1/inv.sigma.sq
    
    return(sigma.sq)
  }
  
  #---------------------------------------------------------------------------------
  # Update skewness parameter lam.est using full-conditional distribution (set b = 10)
  
  lam.est.update <- function (W,inv.Sigma,B.est,X,n,vecy){
    
    b = 1; # prior distribution for lambda ~ N(0,b)
    Swinvw <- foreach(l = 1:n,.combine='+') %dopar% {t(W[,l])%*%inv.Sigma%*%W[,l]}
    A.lam <- (b^2)*Swinvw + 1
    B.lam <- foreach(l = 1:n,.combine='+') %dopar% {(t(vecy)[l,] - t(X)[l,]%*%B.est)%*%inv.Sigma%*%W[,l]*(b^2) + 1} # skewness \pi(\lambda) ~ N(1,1)
    lam.est <- rnorm(1, mean = B.lam/A.lam, sd = b/(sqrt(A.lam)))
    
    return(lam.est)
    
  }
  
  #---------------------------------------------------------------------------------
  # Update W = abs(Z_2) : Note that the each component is sampled from univariate truncated normal distribution
  
  W.update <- function(t,s,n, lam.est, inv.Sigma, B.est, vecy, X, W){
    
    for (jk in 1:(t*s)){
      for (N in 1:n){
        D <- ((lam.est^2)*inv.Sigma + diag(t*s))[jk,jk] 
        E <- ((lam.est*inv.Sigma)%*%(vecy - t(B.est)%*%X))[jk,N] - sum(((lam.est^2)*inv.Sigma)[jk,-jk]*W[-jk,N])
        W[jk,N] = rtruncnorm(1, a=0, b=Inf, mean = E/D, sd = sqrt(1/D)) 
      }
    }
    
    return(W)  
    
  }
  
  
  #---------------------------------------------------------------------------------
  # Update B.est
  B.est.update <- function(X,vecy,inv.Sigma,W,lam.est, t,s,p){
    
    Sxy <- X%*%t(vecy)%*%inv.Sigma
    Sxw <- X%*%t(W)%*%inv.Sigma*lam.est
    Sxx <- X%*%t(X)
    A.inv <- solve(kron(Sxx,inv.Sigma) + 0.001*diag(t*s*p)) 
    vecB <- mvrnorm(1, mu = A.inv%*%(as.vector(t(Sxy)) - as.vector(t(Sxw))), Sigma = A.inv, tol = 1e-3) # \pi(\beta) ~ N_{TSp}(vec(1), \textbf{I}) 
    B.est <- mat(array(vecB, dim = c(t,s,p)),3) # reshape to B_{(3)} : p by ts matrix
    
    return(B.est)
    
  }
  
  begin_sampler <- proc.time()[3]
  
  # MCMC iterations 
  for (i in 1:(n.burn + n.save*thin)) { #}
    
    rho.result = rho.update(t,s,rho,vecy,X,B.est,lam.est,W,sigma.sq)
    rho = rho.result[[1]]; inv.Sigma = rho.result[[2]]; inv.R = rho.result[[3]];  
    
    sigma.sq = sigma.sq.update(t,s,n,vecy,X,B.est,lam.est,W,inv.R)   
    lam.est = lam.est.update(W,inv.Sigma,B.est,X,n,vecy)  
    W = W.update(t,s,n, lam.est, inv.Sigma, B.est, vecy, X, W) 
    B.est = B.est.update(X,vecy,inv.Sigma,W,lam.est, t,s,p)
    
    
    ##Keep track of MCMC output:
    if(i > n.burn & (i - n.burn)%%thin==0){
      ii = (i - n.burn)/thin
      
      rho.save[,ii] <- rho
      sigma.sq.save[,ii] <- sigma.sq
      lam.est.save[,ii] <- lam.est
      W.save[,,ii] <- W
      B.est.save[,,ii] <- B.est
      
    }
    
    if(i %% ceiling((n.save)/10) == 0){
      cat(
        paste0(
          "##### ",
          Sys.time(),
          " Iteration # ",i," of ", (n.save*thin + n.burn),
          " #####\n",
          "##### Time elapsed: ",proc.time()[3] - begin_sampler, " seconds #####\n",
          "##### Time per each iteration: ",(proc.time()[3] - begin_sampler)/i, " seconds #####\n"
        )
      )
    }
    #print(c("Iteration",i, round(rho,4), round(sigma.sq,4), round(lam.est,4)))     
    
  } # end sampler
  
  result = list(rho.save, sigma.sq.save, lam.est.save, W.save, B.est.save)
  names(result) = c('rho','sigma.sq','lam.est','W','B.est')
  
  return(result)
  
} # end BSTN function

result <- BSTN_matrix(Y,X,vecy, n.burn = 100, n.save = 1000, thin = 5)