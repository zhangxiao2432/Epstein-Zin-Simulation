# Run this file to perform simulations 
#
rm(list=ls())
#
# load packages and codes
library("MASS")
library("mgcv")
library("tsDyn")
library("QZ")
library("ggplot2")
library("gmm")
library("vars")
library("tables")
library("xtable")
library("gridExtra")
library("grid")
source("Tauchen.R")
source("data_sim.R")
source("GMM_CMT.R")

# Configuration
Experiment <- matrix(c(0.80,0.80,1.35,1.35,0.8,5.20,1.35,5.20),4,2)

N <- 500  # number of simulations 
TIME <- 90 # number of observations
S1 <- 4   # number of states of each variable
S2 <- 4   

NVAR <- 2 # number of state variable
NLAG <- 1 # number of lag
NVAL <- matrix(c(S1,S2), NVAR, 1)
B <- matrix(c(0.004,0.021), NVAR, 1) #mean
AUTO <- matrix(c(0.117,0.017,0.414,-0.161), NVAR, NVAR*NLAG)     # autoregression
SIGMA <- matrix(c(0.01400,0.00177,0.00177,0.00120), NVAR, NVAR)  # covariance of epsilon

# Initialize
TS_cons <- TS_divids <- Rsp_fit_coef <- Rw_fit_coef <- lamda_fit_coef <- matrix(NA,N,5) 
param_hat_series1 <- param_hat_series2 <- param_hat_series3 <- RRT_mean <- matrix(NA,N,4)
param_hat_series1_A <- param_hat_series2_A <- param_hat_series3_A <- matrix(NA,N,4)
Rsp_fit_sdr <- Rsp_fit_tvalue <- Rw_fit_sdr <- Rw_fit_tvalue <- lamda_fit_sdr <- lamda_fit_tvalue <- matrix(NA,N,4)
RRT_sd <- TS_sder_cons <- TS_sder_divids <- matrix(NA,N,3)
AR_cons_coef <- AR_cons_sd <- AR_cons_tratio <- matrix(NA,N,6)
AR_Rsp_coef <- AR_Rsp_sd <- AR_Rsp_tratio <- AR_Rf_coef <- AR_Rf_sd <- AR_Rf_tratio <- matrix(NA,N,6)
GT_RRT <- GGAR_Rsp <- GGAR_Rf <- GGfit_Rsp <- GGfit_Rw <- GGfit_lamda <- GGAR <- GFIT <- c()
param_hat_series <- param_hat_series_A <- c()

# Use Tauchen's method to approximate a markov chain
Output <-Tauchen(NVAR,NLAG,NVAL,B,AUTO,SIGMA)
# 
set.seed(666)

for(K in 1:nrow(Experiment)){
  
  # Parameters
  beta  <- 0.98 #CDF
  alpha <- Experiment[K,1] #RRA
  rho <- Experiment[K,2] #EIS
  gamma <- round((1-alpha)/(1-rho),digits = 2)
  
  for(i in 1:N){
    
    # Simulate data set
    t0 <- Parameters <- c(beta,alpha,rho)
    #t0 <- c(beta,alpha,rho)
    #Parameters <- c(0.98,0.8,0.8)
   
    
    data1 <- data.sim(Output,Parameters,TIME)
    
    # Creat instruments
    cons_lead <- data1[-1,1,]
    cons_lag <- data1[-TIME,1]
    Rsp_lead <- data1[-1,5]
    Rsp_lag <- data1[-TIME,5]
    Rw_lead <- data1[-1,3]
    Rw_lag <- data1[-TIME,3]
    Rf_lead <- data1[-1,4]
    Rf_lag <- data1[-TIME,4]
    data2 <- data.frame(cbind(cons_lead,cons_lag,Rf_lead,Rf_lag,Rw_lead,Rw_lag,Rsp_lead,Rsp_lag))
    colnames(data2) <- c("Cons_lead","Cons_lag","Rf_lead","Rf_lag","Rw_lead","Rw_lag","Rsp_lead","Rsp_lag")
    
    # series1 
    series1 <- data.frame(cbind(matrix(1,TIME-1,1),data2[,c(2)]))
    P1 <- tensor.prod.model.matrix(list(series1,series1))
    
    GMT1 <- function(t0){
      
      #parameters
      beta <- t0[1]
      alpha <- t0[2]
      rho <- t0[3]
      gamma <- (1-alpha)/(1-rho)
      
      condi_mt <- as.matrix(beta^gamma*cons_lead^(-rho*gamma)*Rw_lead^(gamma-1)*Rsp_lead-1)
      gmat.temp <- tensor.prod.model.matrix(list(P1,condi_mt))
      gmat <- gmat.temp[,-c(3)]
      gbar <- as.matrix(colMeans(gmat))
      w <- crossprod(gmat)/NROW(gmat) 
      #GMT <- crossprod(gbar,solve(w,gbar))
      GMT <- crossprod(gbar)
      
      return(GMT)
    }
    
    # series2 
    series2 <- data.frame(cbind(matrix(1,TIME-1,1),data2[,c(2,4)]))
    P2 <- tensor.prod.model.matrix(list(series2,series2))
    
    GMT2 <- function(t0){
      
      #parameters
      beta <- t0[1]
      alpha <- t0[2]
      rho <- t0[3]
      gamma <- (1-alpha)/(1-rho)
      
      condi_mt <- as.matrix(beta^gamma*cons_lead^(-rho*gamma)*Rw_lead^(gamma-1)*Rsp_lead-1)
      gmat.temp <- tensor.prod.model.matrix(list(P2,condi_mt))
      gmat <- gmat.temp[,-c(4,7,8)]
      #GMT <- norm(t(gmat)%*%gmat,"F")
      gbar <- as.matrix(colMeans(gmat))
      w <- crossprod(gmat)/NROW(gmat) 
      #GMT <- crossprod(gbar,solve(w,gbar))
      GMT <- crossprod(gbar)
      
      return(GMT)
    }
    
    # series3 
    series3 <- data.frame(cbind(matrix(1,TIME-1,1),data2[,c(2,4,8)]))
    P3 <- tensor.prod.model.matrix(list(series3,series3))
    
    GMT3 <- function(t0){
      
      #parameters
      beta <- t0[1]
      alpha <- t0[2]
      rho <- t0[3]
      gamma <- (1-alpha)/(1-rho)
      
      condi_mt <- as.matrix(beta^gamma*cons_lead^(-rho*gamma)*Rw_lead^(gamma-1)*Rsp_lead-1)
      gmat.temp <- tensor.prod.model.matrix(list(P3,condi_mt))
      gmat <- gmat.temp[,-c(5,9,10,13,14,15)]
      #GMT <- norm(t(gmat)%*%gmat,"F")
      gbar <- as.matrix(colMeans(gmat))
      w <- crossprod(gmat)/NROW(gmat) 
      #GMT <- crossprod(gbar,solve(w,gbar))
      GMT <- crossprod(gbar)
      
      return(GMT)
    }
    
    # estimation
    param.est_series1 <- gmm(GMM_CMT1, x=as.matrix(data2),t0,type="twoStep",optfct = "nlminb",wmatrix="ident",vcov ="iid")
    param.est_series2 <- gmm(GMM_CMT2, x=as.matrix(data2),t0,type="twoStep",optfct = "nlminb",wmatrix="ident",vcov ="iid")
    param.est_series3 <- gmm(GMM_CMT3, x=as.matrix(data2),t0,type="twoStep",optfct = "nlminb",wmatrix="ident",vcov ="iid")
    
    # Alternative estimation
    param.est_series1_A <- optim(t0,GMT1, method = "Nelder-Mead")
    param.est_series2_A <- optim(t0,GMT2, method = "Nelder-Mead")
    param.est_series3_A <- optim(t0,GMT3, method = "Nelder-Mead")
    
    #Table1
    #Output_Consgroweh & Divgrowth
    TS <- VAR(data1[,c(1,2)],p=1)
    TS_coef <- coef(TS)
    TS_coef_cons <- t(TS_coef$cons_t)[1,]
    TS_coef_divids <- t(TS_coef$divid_t)[1,]
    TS_cov_cons <- summary(TS)$covres[1,]
    TS_cov_divids <- summary(TS)$covres[2,c(2,1)]
    
    # cons_lag/divid_lag/const/Var/Cov
    TS_cons[i,] <- c(TS_coef_cons,TS_cov_cons)
    TS_divids[i,] <- c(TS_coef_divids,TS_cov_divids)
    TS_sder_cons[i,] <- t(TS_coef$cons_t)[2,]
    TS_sder_divids[i,] <- t(TS_coef$divid_t)[2,]
    remove(TS,TS_coef,TS_coef_cons,TS_coef_divids,TS_cov_cons,TS_cov_divids)
    
    # Output_return
    RRT_mean[i,] <- matrix(c(mean(data1[,3])-1,mean(data1[,5])-1,mean(data1[,4])-1,mean(data1[,5])-mean(data1[,4])),1,4)
    RRT_sd[i,] <- matrix(c(sd(data1[,3]),sd(data1[,5]),sd(data1[,4])),1,3)
    
    # Table2
    # Estimation of instruments autoregression
    AR_cons <- arima(data1[,1],order = c(5,0,0))
    AR_Rsp <- arima(data1[,5],order = c(5,0,0))
    AR_Rf <- arima(data1[,4],order = c(5,0,0))
    
    # Consumption growth
    #ar1/ar2/.../ar5/constant
    AR_cons_coef[i,] <- AR_cons$coef
    AR_cons_sd[i,] <- sqrt(diag(AR_cons$var.coef))
    AR_cons_tratio[i,] <- AR_cons$coef/sqrt(diag(AR_cons$var.coef))
    
    # Rsp
    AR_Rsp_coef[i,] <- AR_Rsp$coef
    AR_Rsp_sd[i,] <- sqrt(diag(AR_Rsp$var.coef))
    AR_Rsp_tratio[i,] <- AR_Rsp$coef/sqrt(diag(AR_Rsp$var.coef))
    
    # Rf
    AR_Rf_coef[i,] <- AR_Rf$coef
    AR_Rf_sd[i,] <- sqrt(diag(AR_Rf$var.coef))
    AR_Rf_tratio[i,] <- AR_Rf$coef/sqrt(diag(AR_Rf$var.coef))
    
    # Regress returns on instruments
    Rsp_fit <- lm(Rsp_lead ~ Rsp_lag + Cons_lag + Rf_lead, data=data2)
    Rw_fit <- lm(Rw_lead ~ Rsp_lag + Cons_lag + Rf_lead, data=data2)
    lamda_fit <- lm(cons_lead ~ Rsp_lag + Cons_lag + Rf_lead, data=data2)
    
    # constant/Rsp_lag/Cons_lag/Rf_lag/adj_Rsquare
    Rsp_fit_coef[i,] <- c(t(summary(Rsp_fit)$coefficients)[1,],summary(Rsp_fit)$adj.r.squared)
    Rsp_fit_sdr[i,] <- t(summary(Rsp_fit)$coefficients)[2,]
    Rsp_fit_tvalue[i,] <- t(summary(Rsp_fit)$coefficients)[3,]
    
    Rw_fit_coef[i,] <- c(t(summary(Rw_fit)$coefficients)[1,],summary(Rw_fit)$adj.r.squared)
    Rw_fit_sdr[i,] <- t(summary(Rw_fit)$coefficients)[2,]
    Rw_fit_tvalue[i,] <- t(summary(Rw_fit)$coefficients)[3,]
    
    lamda_fit_coef[i,] <- c(t(summary(lamda_fit)$coefficients)[1,],summary(lamda_fit)$adj.r.squared)
    lamda_fit_sdr[i,] <- t(summary(lamda_fit)$coefficients)[2,]
    lamda_fit_tvalue[i,] <- t(summary(lamda_fit)$coefficients)[3,]
    
    # Table3
    # Output GMM estimation
    # Series 1
    beta_hat_series1 <- as.numeric(param.est_series1$coefficients[1])
    alpha_hat_series1 <- as.numeric(param.est_series1$coefficients[2])
    rho_hat_series1 <- as.numeric(param.est_series1$coefficients[3])
    gamma_hat_series1 <- as.numeric((1-alpha_hat_series1)/(1-rho_hat_series1))
    
    if (beta_hat_series1 > 0.9 && beta_hat_series1 <1 ){
      param_hat_series1[i,] <- c(beta_hat_series1,alpha_hat_series1,rho_hat_series1,gamma_hat_series1)
    } else {param_hat_series1[i,] <- c(NA,NA,NA,NA)
    }
    remove(beta_hat_series1,alpha_hat_series1,rho_hat_series1,gamma_hat_series1)
    
    # Series 2
    beta_hat_series2 <- as.numeric(param.est_series2$coefficients[1])
    alpha_hat_series2 <- as.numeric(param.est_series2$coefficients[2])
    rho_hat_series2 <- as.numeric(param.est_series2$coefficients[3])
    gamma_hat_series2 <- as.numeric((1-alpha_hat_series2)/(1-rho_hat_series2))
    
    if (beta_hat_series2 > 0.9 && beta_hat_series2 < 1 ){
      param_hat_series2[i,] <- c(beta_hat_series2,alpha_hat_series2,rho_hat_series2,gamma_hat_series2)
    } else {param_hat_series1[i,] <- c(NA,NA,NA,NA)
    }
    remove(beta_hat_series2,alpha_hat_series2,rho_hat_series2,gamma_hat_series2)
    
    # Series 3
    beta_hat_series3 <- as.numeric(param.est_series3$coefficients[1])
    alpha_hat_series3 <- as.numeric(param.est_series3$coefficients[2])
    rho_hat_series3 <- as.numeric(param.est_series3$coefficients[3])
    gamma_hat_series3 <- as.numeric((1-alpha_hat_series3)/(1-rho_hat_series3))
    
    if (beta_hat_series3 > 0.9 && beta_hat_series3 < 1 ){
      param_hat_series3[i,] <- c(beta_hat_series3,alpha_hat_series3,rho_hat_series3,gamma_hat_series3)
    } else {param_hat_series1[i,] <- c(NA,NA,NA,NA)
    }
    remove(beta_hat_series3,alpha_hat_series3,rho_hat_series3,gamma_hat_series3)
    
    # Output GMM estimation alternative method
    # Series 1
    beta_hat_series1_A <- as.numeric(param.est_series1_A$par[1])
    alpha_hat_series1_A <- as.numeric(param.est_series1_A$par[2])
    rho_hat_series1_A <- as.numeric(param.est_series1_A$par[3])
    gamma_hat_series1_A <- as.numeric((1-alpha_hat_series1_A)/(1-rho_hat_series1_A))
    
    if (beta_hat_series1_A > 0.9 && beta_hat_series1_A <1 ){
      param_hat_series1_A[i,] <- c(beta_hat_series1_A,alpha_hat_series1_A,rho_hat_series1_A,gamma_hat_series1_A)
    } else {param_hat_series1_A[i,] <- c(NA,NA,NA,NA)
    }
    remove(beta_hat_series1_A,alpha_hat_series1_A,rho_hat_series1_A,gamma_hat_series1_A)
    
    # Series 2
    beta_hat_series2_A <- as.numeric(param.est_series2_A$par[1])
    alpha_hat_series2_A <- as.numeric(param.est_series2_A$par[2])
    rho_hat_series2_A <- as.numeric(param.est_series2_A$par[3])
    gamma_hat_series2_A <- as.numeric((1-alpha_hat_series2_A)/(1-rho_hat_series2_A))
    
    if (beta_hat_series2_A > 0.9 && beta_hat_series2_A <1 ){
      param_hat_series2_A[i,] <- c(beta_hat_series2_A,alpha_hat_series2_A,rho_hat_series2_A,gamma_hat_series2_A)
    } else {param_hat_series2_A[i,] <- c(NA,NA,NA,NA)
    }
    remove(beta_hat_series2_A,alpha_hat_series2_A,rho_hat_series2_A,gamma_hat_series2_A)
    
    # Series 3
    beta_hat_series3_A <- as.numeric(param.est_series3_A$par[1])
    alpha_hat_series3_A <- as.numeric(param.est_series3_A$par[2])
    rho_hat_series3_A <- as.numeric(param.est_series3_A$par[3])
    gamma_hat_series3_A <- as.numeric((1-alpha_hat_series3_A)/(1-rho_hat_series3_A))
    
    if (beta_hat_series3_A > 0.9 && beta_hat_series3_A <1 ){
      param_hat_series3_A[i,] <- c(beta_hat_series3_A,alpha_hat_series3_A,rho_hat_series3_A,gamma_hat_series3_A)
    } else {param_hat_series3_A[i,] <- c(NA,NA,NA,NA)
    }
    remove(beta_hat_series3_A,alpha_hat_series3_A,rho_hat_series3_A,gamma_hat_series3_A)
    
    if (i%%50 == 0){
      cat("Currently completed : ", i,"/500 \n")
    }
  }
  
  # Generate table1 panel A
  TS_cons_mean <- colMeans(TS_cons[,c(3,2,1,4,5)])
  TS_cons_sd <- c(colMeans(TS_sder_cons),NA,NA)
  TS_divid_mean <- colMeans(TS_divids[,c(3,2,1,4,5)])
  TS_divid_sd <- c(colMeans(TS_sder_divids),NA,NA)
  
  Variables1 <- matrix(c("Dividends","","Consumption",""),4,1)
  Variables2 <- matrix(c("Estimator","Std.error","Estimator","Std.error"),4,1)
  TS_CB_temp <- rbind(TS_divid_mean,TS_divid_sd,TS_cons_mean,TS_cons_sd)
  TS_CB <- data.frame(as.matrix(cbind(Variables1,Variables2,round(TS_CB_temp,digits = 3))))
  rownames(TS_CB) <-NULL
  colnames(TS_CB) <- c("Variable","","Intercept","Divid_lag","Cons_lag","Var","Cov")
  
  # Remove used datasets
  remove(TS_cons_mean,TS_cons_sd,TS_divid_mean,TS_divid_sd,Variables1,Variables2,TS_CB_temp)
  
  # Generate table1 panel B
  Variables1 <- matrix(c(alpha,NA,rho,NA),2,2)
  Variables2 <- matrix(c("Estimator","Std.error"),2,1)
  G_RRT_mean <- colMeans(RRT_mean)
  G_RRT_sd <- c(colMeans(RRT_sd),NA)
  G_RRT_temp <- round(rbind(G_RRT_mean,G_RRT_sd),digits = 3)
  G_RRT <- as.matrix(cbind(Variables1,Variables2,G_RRT_temp))
  
  GT_RRT <- rbind(GT_RRT,G_RRT)
  colnames(GT_RRT) <- c("Alpha","Rho","","Rw","Rsp","Rf","Equity premium")
  rownames(GT_RRT) <- NULL
  remove(Variables1,Variables2,G_RRT_mean,G_RRT_sd,G_RRT_temp,G_RRT)
  
  # Generate table2 panel A
  # Consumption growth
  Variables1 <- matrix(c("Cons_growth","",""),3,1)
  Variables2 <- matrix(c("Estimator","Std.error","T-value"),3,1)
  AR_cons_coef_mean <- colMeans(AR_cons_coef)
  AR_cons_sd_mean <- colMeans(AR_cons_sd)
  AR_cons_tratio_mean <- colMeans(AR_cons_tratio)
  GAR_cons_temp <- round(rbind(AR_cons_coef_mean,AR_cons_sd_mean,AR_cons_tratio_mean),digits = 3)
  GAR_cons <- as.matrix(cbind(Variables1,Variables2,GAR_cons_temp[,c(6,1,2,3,4,5)]))
  colnames(GAR_cons) <- c("Variable","","Intercept","ar1","ar2","ar3","ar4","ar5")
  rownames(GAR_cons) <- NULL
  remove(Variables1,Variables2,AR_cons_coef_mean,AR_cons_sd_mean,AR_cons_tratio_mean,GAR_cons_temp)
  
  #Rsp
  Variables1 <- matrix(c("Rsp","",""),3,1)
  Variables2 <- matrix(c("Estimator","Std.error","T-value"),3,1)
  AR_Rsp_coef_mean <- colMeans(AR_Rsp_coef)
  AR_Rsp_sd_mean <- colMeans(AR_Rsp_sd)
  AR_Rsp_tratio_mean <- colMeans(AR_Rsp_tratio)
  GAR_Rsp_temp <- round(rbind(AR_Rsp_coef_mean,AR_Rsp_sd_mean,AR_Rsp_tratio_mean),digits = 3)
  GAR_Rsp <- as.matrix(cbind(Variables1,Variables2,GAR_Rsp_temp[,c(6,1,2,3,4,5)]))
  colnames(GAR_Rsp) <- c("Variable","","Intercept","ar1","ar2","ar3","ar4","ar5")
  rownames(GAR_Rsp) <- NULL
  remove(Variables1,Variables2,AR_Rsp_coef_mean,AR_Rsp_sd_mean,AR_Rsp_tratio_mean,GAR_Rsp_temp)
  
  GGAR_Rsp <- rbind(GGAR_Rsp,GAR_Rsp)
  colnames(GGAR_Rsp) <- c("Variable","","Intercept","ar1","ar2","ar3","ar4","ar5")
  rownames(GGAR_Rsp) <- NULL
  
  #Rf
  Variables1 <- matrix(c("Rf","",""),3,1)
  Variables2 <- matrix(c("Estimator","Std.error","T-value"),3,1)
  AR_Rf_coef_mean <- colMeans(AR_Rf_coef)
  AR_Rf_sd_mean <- colMeans(AR_Rf_sd)
  AR_Rf_tratio_mean <- colMeans(AR_Rf_tratio)
  GAR_Rf_temp <- round(rbind(AR_Rf_coef_mean,AR_Rf_sd_mean,AR_Rf_tratio_mean),digits = 3)
  GAR_Rf <- as.matrix(cbind(Variables1,Variables2,GAR_Rf_temp[,c(6,1,2,3,4,5)]))
  colnames(GAR_Rf) <- c("Variable","","Intercept","ar1","ar2","ar3","ar4","ar5")
  rownames(GAR_Rf) <- NULL
  remove(Variables1,Variables2,AR_Rf_coef_mean,AR_Rf_sd_mean,AR_Rf_tratio_mean,GAR_Rf_temp)
  
  GGAR_Rf <- rbind(GGAR_Rf,GAR_Rf)
  colnames(GGAR_Rf) <- c("Variable","","Intercept","ar1","ar2","ar3","ar4","ar5")
  rownames(GGAR_Rf) <- NULL
  
  title <- matrix(c(beta,alpha,rho,NA,NA,NA,NA,NA),1,8)
  GGAR <- rbind(GGAR,title,GAR_cons,GAR_Rsp,GAR_Rf)
  colnames(GGAR) <- c("Variable","","Intercept","ar1","ar2","ar3","ar4","ar5")
  rownames(GGAR) <- NULL
  remove(GAR_cons,GAR_Rsp,GAR_Rf,title)
  
  # Generate table2 panel B
  # Rsp
  Variables1 <- matrix(c("Rsp","",""),3,1)
  Variables2 <- matrix(c("Estimator","Std.error","T-value"),3,1)
  Rsp_fit_coef_mean <- colMeans(Rsp_fit_coef)
  Rsp_fit_sdr_mean <- c(colMeans(Rsp_fit_sdr),NA)
  Rsp_fit_tvalue_mean <- c(colMeans(Rsp_fit_tvalue),NA)
  Gfit_Rsp_temp <- round(rbind(Rsp_fit_coef_mean,Rsp_fit_sdr_mean,Rsp_fit_tvalue_mean),digits = 3)
  Gfit_Rsp <- as.matrix(cbind(Variables1,Variables2,Gfit_Rsp_temp))
  colnames(Gfit_Rsp) <- c("Variable","","Intercept","Rsp_lag","Cons_lag","Rf_lag","Adj_Rsquare")
  rownames(Gfit_Rsp) <- NULL
  remove(Variables1,Variables2,Rsp_fit_coef_mean,Rsp_fit_sdr_mean,Rsp_fit_tvalue_mean,Gfit_Rsp_temp)
  
  GGfit_Rsp <- rbind(GGfit_Rsp,Gfit_Rsp)
  colnames(Gfit_Rsp) <- c("Variable","","Intercept","Rsp_lag","Cons_lag","Rf_lag","Adj_Rsquare")
  rownames(Gfit_Rsp) <- NULL
  
  # Rw
  Variables1 <- matrix(c("Rw","",""),3,1)
  Variables2 <- matrix(c("Estimator","Std.error","T-value"),3,1)
  Rw_fit_coef_mean <- colMeans(Rw_fit_coef)
  Rw_fit_sdr_mean <- c(colMeans(Rw_fit_sdr),NA)
  Rw_fit_tvalue_mean <- c(colMeans(Rw_fit_tvalue),NA)
  Gfit_Rw_temp <- round(rbind(Rw_fit_coef_mean,Rw_fit_sdr_mean,Rw_fit_tvalue_mean),digits = 3)
  Gfit_Rw <- as.matrix(cbind(Variables1,Variables2,Gfit_Rw_temp))
  colnames(Gfit_Rw) <- c("Variable","","Intercept","Rsp_lag","Cons_lag","Rf_lag","Adj_Rsquare")
  rownames(Gfit_Rw) <- NULL
  remove(Variables1,Variables2,Rw_fit_coef_mean,Rw_fit_sdr_mean,Rw_fit_tvalue_mean,Gfit_Rw_temp)
  
  GGfit_Rw <- rbind(GGfit_Rw,Gfit_Rw)
  colnames(Gfit_Rw) <- c("Variable","","Intercept","Rsp_lag","Cons_lag","Rf_lag","Adj_Rsquare")
  rownames(Gfit_Rw) <- NULL
  
  # lamda
  Variables1 <- matrix(c("lamda","",""),3,1)
  Variables2 <- matrix(c("Estimator","Std.error","T-value"),3,1)
  lamda_fit_coef_mean <- colMeans(lamda_fit_coef)
  lamda_fit_sdr_mean <- c(colMeans(lamda_fit_sdr),NA)
  lamda_fit_tvalue_mean <- c(colMeans(lamda_fit_tvalue),NA)
  Gfit_lamda_temp <- round(rbind(lamda_fit_coef_mean,lamda_fit_sdr_mean,lamda_fit_tvalue_mean),digits = 3)
  Gfit_lamda <- as.matrix(cbind(Variables1,Variables2,Gfit_lamda_temp))
  colnames(Gfit_lamda) <- c("Variable","","Intercept","Rsp_lag","Cons_lag","Rf_lag","Adj_Rsquare")
  rownames(Gfit_lamda) <- NULL
  remove(Variables1,Variables2,lamda_fit_coef_mean,lamda_fit_sdr_mean,lamda_fit_tvalue_mean,Gfit_lamda_temp)
  
  GGfit_lamda <- rbind(GGfit_lamda,Gfit_lamda)
  colnames(Gfit_lamda) <- c("Variable","","Intercept","Rsp_lag","Cons_lag","Rf_lag","Adj_Rsquare")
  rownames(Gfit_lamda) <- NULL
  
  title <- matrix(c(beta,alpha,rho,NA,NA,NA,NA,NA),1,7)
  GFIT <- rbind(GFIT,title,Gfit_lamda,Gfit_Rsp,Gfit_Rw)
  colnames(Gfit_lamda) <- c("Variable","","Intercept","Rsp_lag","Cons_lag","Rf_lag","Adj_Rsquare")
  rownames(Gfit_lamda) <- NULL
  remove(title,Gfit_lamda,Gfit_Rsp,Gfit_Rw)
  
  
  # Generate table3
  # series1
  Variables1 <- matrix(c("Series1","",""),3,1)
  Variables2 <- matrix(c("Mean","Median","std.dev"),3,1)
  param_hat_series1_mean <- c(mean(param_hat_series1[,1],na.rm = TRUE),mean(param_hat_series1[,2],na.rm = TRUE),
                              mean(param_hat_series1[,3],na.rm = TRUE),mean(param_hat_series1[,4],na.rm = TRUE))
  param_hat_series1_median <- c(median(param_hat_series1[,1],na.rm = TRUE),median(param_hat_series1[,2],na.rm = TRUE),
                                median(param_hat_series1[,3],na.rm = TRUE),median(param_hat_series1[,4],na.rm = TRUE))
  param_hat_series1_sd <- c(sd(param_hat_series1[,1],na.rm = TRUE),sd(param_hat_series1[,2],na.rm = TRUE),
                                sd(param_hat_series1[,3],na.rm = TRUE),sd(param_hat_series1[,4],na.rm = TRUE))
  param_hat_series1_temp1 <- round(rbind(param_hat_series1_mean,param_hat_series1_median,
                                         param_hat_series1_sd),digits = 3)
  param_hat_series1_temp2 <- as.matrix(cbind(Variables1,Variables2,param_hat_series1_temp1))
  colnames(param_hat_series1_temp2) <- c("Series","","Beta","Alpha","Rho","Gamma") 
  rownames(param_hat_series1_temp2) <- NULL
  remove(Variables1,Variables2,param_hat_series1_mean,param_hat_series1_median,param_hat_series1_sd,
         param_hat_series1_temp1)
  
  # series2
  Variables1 <- matrix(c("Series2","",""),3,1)
  Variables2 <- matrix(c("Mean","Median","std.dev"),3,1)
  param_hat_series2_mean <- c(mean(param_hat_series2[,1],na.rm = TRUE),mean(param_hat_series2[,2],na.rm = TRUE),
                              mean(param_hat_series2[,3],na.rm = TRUE),mean(param_hat_series2[,4],na.rm = TRUE))
  param_hat_series2_median <- c(median(param_hat_series2[,1],na.rm = TRUE),median(param_hat_series2[,2],na.rm = TRUE),
                                median(param_hat_series2[,3],na.rm = TRUE),median(param_hat_series2[,4],na.rm = TRUE))
  param_hat_series2_sd <- c(sd(param_hat_series2[,1],na.rm = TRUE),sd(param_hat_series2[,2],na.rm = TRUE),
                            sd(param_hat_series2[,3],na.rm = TRUE),sd(param_hat_series2[,4],na.rm = TRUE))
  param_hat_series2_temp1 <- round(rbind(param_hat_series2_mean,param_hat_series2_median,
                                         param_hat_series2_sd),digits = 3)
  param_hat_series2_temp2 <- as.matrix(cbind(Variables1,Variables2,param_hat_series2_temp1))
  colnames(param_hat_series2_temp2) <- c("Series","","Beta","Alpha","Rho","Gamma") 
  rownames(param_hat_series2_temp2) <- NULL
  remove(Variables1,Variables2,param_hat_series2_mean,param_hat_series2_median,param_hat_series2_sd,
         param_hat_series2_temp1)
  
  # series3
  Variables1 <- matrix(c("Series3","",""),3,1)
  Variables2 <- matrix(c("Mean","Median","std.dev"),3,1)
  param_hat_series3_mean <- c(mean(param_hat_series3[,1],na.rm = TRUE),mean(param_hat_series3[,2],na.rm = TRUE),
                              mean(param_hat_series3[,3],na.rm = TRUE),mean(param_hat_series3[,4],na.rm = TRUE))
  param_hat_series3_median <- c(median(param_hat_series3[,1],na.rm = TRUE),median(param_hat_series3[,2],na.rm = TRUE),
                                median(param_hat_series3[,3],na.rm = TRUE),median(param_hat_series3[,4],na.rm = TRUE))
  param_hat_series3_sd <- c(sd(param_hat_series3[,1],na.rm = TRUE),sd(param_hat_series3[,2],na.rm = TRUE),
                            sd(param_hat_series3[,3],na.rm = TRUE),sd(param_hat_series3[,4],na.rm = TRUE))
  param_hat_series3_temp1 <- round(rbind(param_hat_series3_mean,param_hat_series3_median,
                                         param_hat_series3_sd),digits = 3)
  param_hat_series3_temp2 <- as.matrix(cbind(Variables1,Variables2,param_hat_series3_temp1))
  colnames(param_hat_series3_temp2) <- c("Series","","Beta","Alpha","Rho","Gamma") 
  rownames(param_hat_series3_temp2) <- NULL
  remove(Variables1,Variables2,param_hat_series3_mean,param_hat_series3_median,param_hat_series3_sd,
         param_hat_series3_temp1)
  
  title <- matrix(c(NA,NA,beta,alpha,rho,gamma),1,6)
  param_hat_series <- rbind(param_hat_series,title,param_hat_series1_temp2,param_hat_series2_temp2,
                            param_hat_series3_temp2)
  colnames(param_hat_series) <- c("Series","","Beta","Alpha","Rho","Gamma") 
  rownames(param_hat_series) <- NULL
  remove(param_hat_series1_temp2,param_hat_series2_temp2,param_hat_series3_temp2)
  
  # Alternative method
  # series1
  Variables1 <- matrix(c("Series1","",""),3,1)
  Variables2 <- matrix(c("Mean","Median","std.dev"),3,1)
  param_hat_series1_A_mean <- c(mean(param_hat_series1_A[,1],na.rm = TRUE),mean(param_hat_series1_A[,2],na.rm = TRUE),
                                mean(param_hat_series1_A[,3],na.rm = TRUE),mean(param_hat_series1_A[,4],na.rm = TRUE))
  param_hat_series1_A_median <- c(median(param_hat_series1_A[,1],na.rm = TRUE),median(param_hat_series1_A[,2],na.rm = TRUE),
                                median(param_hat_series1_A[,3],na.rm = TRUE),median(param_hat_series1_A[,4],na.rm = TRUE))
  param_hat_series1_A_sd <- c(sd(param_hat_series1_A[,1],na.rm = TRUE),sd(param_hat_series1_A[,2],na.rm = TRUE),
                            sd(param_hat_series1_A[,3],na.rm = TRUE),sd(param_hat_series1_A[,4],na.rm = TRUE))
  param_hat_series1_A_temp1 <- round(rbind(param_hat_series1_A_mean,param_hat_series1_A_median,
                                         param_hat_series1_A_sd),digits = 3)
  param_hat_series1_A_temp2 <- as.matrix(cbind(Variables1,Variables2,param_hat_series1_A_temp1))
  colnames(param_hat_series1_A_temp2) <- c("Series","","Beta","Alpha","Rho","Gamma") 
  rownames(param_hat_series1_A_temp2) <- NULL
  remove(Variables1,Variables2,param_hat_series1_A_mean,param_hat_series1_A_median,param_hat_series1_A_sd,
         param_hat_series1_A_temp1)
  
  # series2
  Variables1 <- matrix(c("Series2","",""),3,1)
  Variables2 <- matrix(c("Mean","Median","std.dev"),3,1)
  param_hat_series2_A_mean <- c(mean(param_hat_series2_A[,1],na.rm = TRUE),mean(param_hat_series2_A[,2],na.rm = TRUE),
                                mean(param_hat_series2_A[,3],na.rm = TRUE),mean(param_hat_series2_A[,4],na.rm = TRUE))
  param_hat_series2_A_median <- c(median(param_hat_series2_A[,1],na.rm = TRUE),median(param_hat_series2_A[,2],na.rm = TRUE),
                                median(param_hat_series2_A[,3],na.rm = TRUE),median(param_hat_series2_A[,4],na.rm = TRUE))
  param_hat_series2_A_sd <- c(sd(param_hat_series2_A[,1],na.rm = TRUE),sd(param_hat_series2_A[,2],na.rm = TRUE),
                            sd(param_hat_series2_A[,3],na.rm = TRUE),sd(param_hat_series2_A[,4],na.rm = TRUE))
  param_hat_series2_A_temp1 <- round(rbind(param_hat_series2_A_mean,param_hat_series2_A_median,
                                         param_hat_series2_A_sd),digits = 3)
  param_hat_series2_A_temp2 <- as.matrix(cbind(Variables1,Variables2,param_hat_series2_A_temp1))
  colnames(param_hat_series2_A_temp2) <- c("Series","","Beta","Alpha","Rho","Gamma") 
  rownames(param_hat_series2_A_temp2) <- NULL
  remove(Variables1,Variables2,param_hat_series2_A_mean,param_hat_series2_A_median,param_hat_series2_A_sd,
         param_hat_series2_A_temp1)
  
  # series3
  Variables1 <- matrix(c("Series3","",""),3,1)
  Variables2 <- matrix(c("Mean","Median","std.dev"),3,1)
  param_hat_series3_A_mean <- c(mean(param_hat_series3_A[,1],na.rm = TRUE),mean(param_hat_series3_A[,2],na.rm = TRUE),
                                mean(param_hat_series3_A[,3],na.rm = TRUE),mean(param_hat_series3_A[,4],na.rm = TRUE))
  param_hat_series3_A_median <- c(median(param_hat_series3_A[,1],na.rm = TRUE),median(param_hat_series3_A[,2],na.rm = TRUE),
                                median(param_hat_series3_A[,3],na.rm = TRUE),median(param_hat_series3_A[,4],na.rm = TRUE))
  param_hat_series3_A_sd <- c(sd(param_hat_series3_A[,1],na.rm = TRUE),sd(param_hat_series3_A[,2],na.rm = TRUE),
                            sd(param_hat_series3_A[,3],na.rm = TRUE),sd(param_hat_series3_A[,4],na.rm = TRUE))
  param_hat_series3_A_temp1 <- round(rbind(param_hat_series3_A_mean,param_hat_series3_A_median,
                                         param_hat_series3_A_sd),digits = 3)
  param_hat_series3_A_temp2 <- as.matrix(cbind(Variables1,Variables2,param_hat_series3_A_temp1))
  colnames(param_hat_series3_A_temp2) <- c("Series","","Beta","Alpha","Rho","Gamma") 
  rownames(param_hat_series3_A_temp2) <- NULL
  remove(Variables1,Variables2,param_hat_series3_A_mean,param_hat_series3_A_median,param_hat_series3_A_sd,
         param_hat_series3_A_temp1)
  
  title <- matrix(c(NA,NA,beta,alpha,rho,gamma),1,6)
  param_hat_series_A <- rbind(param_hat_series_A,title,param_hat_series1_A_temp2,param_hat_series2_A_temp2,
                            param_hat_series3_A_temp2)
  colnames(param_hat_series_A) <- c("Series","","Beta","Alpha","Rho","Gamma") 
  rownames(param_hat_series_A) <- NULL
  remove(param_hat_series1_A_temp2,param_hat_series2_A_temp2,param_hat_series3_A_temp2,title)
  
  cat("Completed combination: ", Experiment[K,],"\n")
  
}

# Remove used datasets
remove(TS_cons,TS_divids,TS_sder_cons,TS_sder_divids)
#
remove(AR_cons_coef,AR_cons_sd,AR_cons_tratio)
remove(AR_Rsp_coef,AR_Rsp_sd,AR_Rsp_tratio)
remove(AR_Rf_coef,AR_Rf_sd,AR_Rf_tratio)
remove(GGAR_Rf,GGAR_Rsp)
remove(AR_Rf,AR_Rsp,AR_cons)
#
remove(Rsp_fit_coef,Rsp_fit_sdr,Rsp_fit_tvalue)
remove(Rw_fit_coef,Rw_fit_sdr,Rw_fit_tvalue)
remove(lamda_fit_coef,lamda_fit_sdr,lamda_fit_tvalue)
remove(lamda_fit,Rsp_fit,Rw_fit)
remove(GGfit_lamda,GGfit_Rsp,GGfit_Rw)
#
remove(param.est_series1,param.est_series2,param.est_series3)
remove(param.est_series1_A,param.est_series2_A,param.est_series3_A)
remove(param_hat_series1,param_hat_series2,param_hat_series3)
remove(param_hat_series1_A,param_hat_series2_A,param_hat_series3_A)
remove(P1,P2,P3)
remove(cons_lead,cons_lag,Rf_lead,Rf_lag,Rw_lead,Rw_lag,Rsp_lead,Rsp_lag)
#
# Print results
# Table 1
# A
TS_CB <- as.matrix(TS_CB)
TS_CB[is.na(TS_CB)] <- " "
GT_RRT[is.na(GT_RRT)] <- " "
T1A_title <- textGrob("Table 1A",just="left",gp=gpar(fontsize=10))
T1A <- tableGrob(as.matrix(TS_CB), cols = colnames(TS_CB),theme = ttheme_default())
T1B <- tableGrob(GT_RRT, cols = colnames(GT_RRT),theme = ttheme_default())

h <- grobHeight(T1A)
w <- grobWidth(T1A)
T1A_title <- textGrob("Table 1A", y=unit(0.5,"npc") + h,
                  vjust=0, gp=gpar(fontsize=15,fontface="italic"))
G_T1A <- gTree(children=gList(T1A, T1A_title))
remove(h,w)

# B
h <- grobHeight(T1B)
w <- grobWidth(T1B)
T1B_title <- textGrob("Table 1B", y=unit(0.5,"npc") + h,
                      vjust=0, gp=gpar(fontsize=15,fontface="italic"))
G_T1B <- gTree(children=gList(T1B, T1B_title))
Tab1 <- grid.arrange(G_T1A,G_T1B,nrow=2)
remove(h,w,T1A,T1A_title,T1B,T1B_title,G_T1A,G_T1B)

# Table 2
# A
GGAR[is.na(GGAR)] <- " "
T2A_title <- textGrob("Table 2A",just="left",gp=gpar(fontsize=10))
T2A <- tableGrob(GGAR, cols = colnames(GGAR),theme = ttheme_default())

h <- grobHeight(T2A)
w <- grobWidth(T2A)
T2A_title <- textGrob("Table 2A", y=unit(0.5,"npc") + h,
                      vjust=0, gp=gpar(fontsize=15,fontface="italic"))
G_T2A <- gTree(children=gList(T2A, T2A_title))
remove(h,w,T2A,T2A_title)

# B
GFIT[is.na(GFIT)] <- " "
T2B_title <- textGrob("Table 2B",just="left",gp=gpar(fontsize=10))
T2B <- tableGrob(GFIT, cols = colnames(GFIT),theme = ttheme_default())

h <- grobHeight(T2B)
w <- grobWidth(T2B)
T2B_title <- textGrob("Table 2B", y=unit(0.5,"npc") + h,
                      vjust=0, gp=gpar(fontsize=15,fontface="italic"))
G_T2B <- gTree(children=gList(T2B, T2B_title))
remove(h,w,T2B,T2B_title)

# Table 3
# A
param_hat_series[is.na(param_hat_series)] <- " "
T3A_title <- textGrob("Table 3A",just="left",gp=gpar(fontsize=10))
T3A <- tableGrob(param_hat_series, cols = colnames(param_hat_series),theme = ttheme_default())

h <- grobHeight(T3A)
w <- grobWidth(T3A)
T3A_title <- textGrob("Table 3A", y=unit(0.5,"npc") + h,
                      vjust=0, gp=gpar(fontsize=15,fontface="italic"))
G_T3A <- gTree(children=gList(T3A, T3A_title))
remove(h,w,T3A,T3A_title)

# B
param_hat_series_A[is.na(param_hat_series_A)] <- " "
T3B_title <- textGrob("Table 3B",just="left",gp=gpar(fontsize=10))
T3B <- tableGrob(param_hat_series_A, cols = colnames(param_hat_series_A),theme = ttheme_default())

h <- grobHeight(T3B)
w <- grobWidth(T3B)
T3B_title <- textGrob("Table 3B", y=unit(0.5,"npc") + h,
                      vjust=0, gp=gpar(fontsize=15,fontface="italic"))
G_T3B <- gTree(children=gList(T3B, T3B_title))
remove(h,w,T3B,T3B_title)

###
pdf("data_output.pdf", height=12)
grid.draw(Tab1)
grid.newpage()
grid.draw(G_T2A)
grid.newpage()
grid.draw(G_T2B)
grid.newpage()
grid.draw(G_T3A)
grid.newpage()
grid.draw(G_T3B)
dev.off()

#remove(Tab1,G_T2A,G_T2B,G_T3A,G_T3B)

#Generate LaTex code
#Table 1A
xtable(TS_CB)
#############################################################

#Table 1B
xtable(GT_RRT)
#############################################################

#Table 2
table2 <- cbind(param_hat_series,param_hat_series_A)
xtable(table2[,-(7:8)])
#############################################################

#Table A1
xtable(GFIT)
#############################################################

#Table A2
xtable(GGAR)
#############################################################

