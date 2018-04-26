rm(list = ls())

library(readr)
library(MASS)
library(ggplot2)
library(grid)
library(reshape2)
library(scales)
library(QZ)
library(gmm)
library(mgcv)
library(boot)

source("basis function.R")
source("param estimation.R") 
source("vfn estimation(coef).R")
source("vfn estimation(basis).R")

set.seed(666)
param_hat <- value.Fn <- c()
CFL_data <- as.matrix(read_csv("CFL_QE_2013.csv"))

# Configuration
Experiment <- matrix(c(50,20,0.8,5.2,0.7,0.9,0.8,1.2),4,2)

#for(K in 4:nrow(Experiment)){
K=4
  # initial values
  beta <- 0.98
  theta <- as.numeric(Experiment[K,1])
  rho <- as.numeric(Experiment[K,2])
  
  t0 <- c(beta,theta,rho)
  #
  #initial {a} of value function  
  m <- 3
  kk <- m*m
  itmax <- 100
  a <- runif(kk+2,0.6,0.9)
  tet0 <- c(vc_0=a[1],alpha0=a[2],alpha=a[3:(kk+2)])
  
  #Loop
  # Estimation
  for(i in 1:itmax){
    
    # Estimation of value function
    sdf.est <- optim(tet0, vfn.est, method = "Nelder-Mead")
    vfn_CMT_temp1 <- vfn.est(tet0)
    #
    #Extract coefficients
    a.hat <- as.vector(sdf.est$par)
    vfn_CMT_temp2 <- vfn.est(a.hat)
    #
    #stop the loop if the GMM objective is small enough
    if(abs(vfn_CMT_temp1 - vfn_CMT_temp2) > 0.0001){
      tet0 <- a.hat
    } 
    #form estimated value
    vfn0 <- vfn.hat(tet0)
    #
    #Estimate parameters
    #
    param.temp1 <- optim(t0, param.est, method = "Nelder-Mead")
    param_CMT_temp1 <- param.est(t0)
    #
    param.temp2 <- as.vector(param.temp1$par)
    param_CMT_temp2 <- param.est(param.temp2)
    #
    #stop the loop if the GMM objective is small enough
    if(abs(vfn_CMT_temp1 - param_CMT_temp2) < 0.0001){
      param.result <- param.temp2
      cat("Iteration stoped because the GMM objective is converged. \n")
      remove(t0,tet0)
      break
    }
    #
    est.diff <- abs(param.temp2 - t0)
    if(all(est.diff < 0.01)==TRUE){
      param.result <- param.temp2
      cat("Convergence succeed!")
      break
    }  else {
      t0 <- param.temp2
    }
    #
    if (i%%10 == 0){
      cat("Current Iteration : ", i,"/Itmax \n")
    }
    if(i >= itmax){
      cat("Reached the itmax, convergence failed! \n")
    }
  }
  
  value.Fn <- cbind(value.Fn, vfn0)
  
  param_hat <- rbind(param_hat, param.result)
  
  #remove(vfn_CMT_temp1,vfn_CMT_temp2,param_CMT_temp1,param_CMT_temp2,param.result,vfn0,tet0,t0)
  #remove(param.temp1,param.temp2,a, a.hat,est.diff,sdf.est)
  
  cat("Completed combination: ", Experiment[K,],"\n")

#}

mean(vfn0)
sd(vfn0)
arima(vfn0,order = c(1,0,0))
