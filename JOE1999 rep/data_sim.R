#Function to generate dataset

library("rootSolve")

data.sim <- function(Tauchen,parameters,TT){
  
  #parameters
  beta <- parameters[1]
  alpha <- parameters[2]
  rho <- parameters[3]
  gamma <- (1-alpha)/(1-rho)
  
  #Simulate consumption growth and dividend growth time series.
  C <- cbind(B,AUTO)
  data1 <- VAR.sim(C,n=TIME,lag = 1,varcov = SIGMA)
  colnames(data1) <- c("Xi_t","Lamda_t")
  divid_t <- as.matrix(exp(data1[,1]))
  cons_t <- as.matrix(exp(data1[,2]))
  
  #Generate state value for consumption growth(g) and dividend growth(d).
  stat_value_d <- exp(Output$YMAT[,1])
  stat_value_g <- exp(Output$YMAT[,2])
  
  #Generate transition probability matrix.
  trans_prob <- Output$PIOUT
  
  #Generate stationary probability distribution.
  stat_prob <- Output$PISTA
  
  #R_w
  #Solve the price/dividends ratio v_star
  equation1 <- function(v){
    
    f1 <- f2 <- f3 <- matrix(0,S1*S2,1)
    for (i in 1:(S1*S2)) {
      for (j in 1:(S1*S2)) {
        f1[j] <- trans_prob[i,j]*(beta^gamma)*(stat_value_g[j])^(1-alpha)*(1+v[j])^gamma
      }
      f2[i] <- sum(f1)
      f3[i] <- f2[i]-v[i]^gamma
    }
    c(eq = f3)
  }
  v_star <- multiroot(f = equation1, start = c(1:(S1*S2)))$root
  
  R_w <- matrix(0,(S1*S2),(S1*S2))
  for (i in 1:(S1*S2)) {
    for(j in 1:(S1*S2)){
      R_w[i,j] <- (1+v_star[j])*stat_value_g[j]/v_star[i]
    }
  }
  
  #R_f
  R_f <- matrix(0,(S1*S2),1)
  R_f_temp <- matrix(0,(S1*S2),1)
  for (i in 1:(S1*S2)) {
    for (j in 1:(S1*S2)) {
      R_f_temp[j] <- trans_prob[i,j]*(stat_value_g[j])^(-alpha)*((1+v_star[j])/v_star[i])^((rho-alpha)/(1-rho))
    }
    R_f[i] <- 1/((beta^gamma)*sum(R_f_temp))
  }
  remove("R_f_temp")
  
  #R_sp
  #Solve the price/dividends ratio vsp_star
  equation2 <- function(v){
    
    f1 <- f2 <- f3 <- matrix(0,S1*S2,1)
    for (i in 1:(S1*S2)) {
      for (j in 1:(S1*S2)) {
        f1[j] <- trans_prob[i,j]*(beta^gamma)*(stat_value_g[j])^(1-alpha)*(1+v[j])^gamma
      }
      f2[i] <- sum(f1)
      f3[i] <- f2[i]-v[i]^gamma
    }
    c(eq = f3)
  }
  vsp_star <- multiroot(f = equation2, start = c(1:(S1*S2)))$root
  
  R_sp <- matrix(0,(S1*S2),(S1*S2))
  for (i in 1:(S1*S2)) {
    for(j in 1:(S1*S2)){
      R_sp[i,j] <- (1+vsp_star[j])/vsp_star[i]*stat_value_d[j]
    }
  }
  
  # Generate time series from Markov chain
  # Return on aggregate wealth
  Rw_t <- function(TT){
    RR <- matrix(0, TT, 1)
    for (i in 1:TT){
      u <- runif(1,min = 0, max = 1)
      P <- PP <- 0
      j <- jj <- 1
      while (P < u) {
        P <- P + stat_prob[j]
        j <- j + 1
      }
      s <- j-1
      uu <- runif(1,min = 0, max = 1)
      while (PP < uu) {
        PP <- PP + trans_prob[s,jj]
        jj <- jj + 1
      }
      ss <- jj-1
      RR[i] <- R_w[s,ss] 
    }
    return(RR)
  }
  
  # Return on risk-free assets
  Rf_t <- function(TT){
    RR <- matrix(0, TT, 1)
    for (i in 1:TT){
      u <- runif(1,min = 0, max = 1)
      P <- PP <- 0
      j <- jj <- 1
      while (P < u) {
        P <- P + stat_prob[j]
        j <- j + 1
      }
      s <- j-1
      RR[i] <- R_f[s] 
    }
    return(RR)
  }
  
  # Return on S&P500 assets
  Rsp_t <- function(TT){
    RR <- matrix(0, TT, 1)
    for (i in 1:TT){
      u <- runif(1,min = 0, max = 1)
      P <- PP <- 0
      j <- jj <- 1
      while (P < u) {
        P <- P + stat_prob[j]
        j <- j + 1
      }
      s <- j-1
      uu <- runif(1,min = 0, max = 1)
      while (PP < uu) {
        PP <- PP + trans_prob[s,jj]
        jj <- jj + 1
      }
      ss <- jj-1
      RR[i] <- R_sp[s,ss] 
    }
    return(RR)
  }
  
  Rw_t <- Rw_t(TIME)
  Rf_t <- Rf_t(TIME)
  Rsp_t <- Rsp_t(TIME)
  
  data <- data.frame(cbind(cons_t, divid_t, Rw_t, Rf_t, Rsp_t))
  colnames(data) <- c("cons_t", "divid_t", "Rw_t", "Rf_t", "Rsp_t")
  
  return(data)
}