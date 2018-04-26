# Function to generate GMM moment condition

# series1 
GMM_CMT1 <- function(t0,data){
  
  series1 <- data.frame(cbind(matrix(1,TIME-1,1),data[,c(2)]))
  P1 <- tensor.prod.model.matrix(list(series1,series1))
  
  #parameters
  beta <- t0[1]
  alpha <- t0[2]
  rho <- t0[3]
  gamma <- (1-alpha)/(1-rho)
  
  condi_mt <- as.matrix(beta^gamma*cons_lead^(-rho*gamma)*Rw_lead^(gamma-1)*Rsp_lead-1)
  
  gmat.temp <- tensor.prod.model.matrix(list(P1,condi_mt))
  gmat <- gmat.temp[,-c(3)]
  
  return(gmat)
}

# series2 
GMM_CMT2 <- function(t0,data){
  
  series2 <- data.frame(cbind(matrix(1,TIME-1,1),data[,c(2,4)]))
  P2 <- tensor.prod.model.matrix(list(series2,series2))
  
  #parameters
  beta <- t0[1]
  alpha <- t0[2]
  rho <- t0[3]
  gamma <- (1-alpha)/(1-rho)
  
  condi_mt <- as.matrix(beta^gamma*cons_lead^(-rho*gamma)*Rw_lead^(gamma-1)*Rsp_lead-1)
  
  gmat.temp <- tensor.prod.model.matrix(list(P2,condi_mt))
  gmat <- gmat.temp[,-c(4,7,8)]
  
  return(gmat)
}

# series3 
GMM_CMT3 <- function(t0,data){
  
  series3 <- data.frame(cbind(matrix(1,TIME-1,1),data[,c(2,4,8)]))
  P3 <- tensor.prod.model.matrix(list(series3,series3))
  
  #parameters
  beta <- t0[1]
  alpha <- t0[2]
  rho <- t0[3]
  gamma <- (1-alpha)/(1-rho)
  
  condi_mt <- as.matrix(beta^gamma*cons_lead^(-rho*gamma)*Rw_lead^(gamma-1)*Rsp_lead-1)
  
  gmat.temp <- tensor.prod.model.matrix(list(P3,condi_mt))
  gmat <- gmat.temp[,-c(5,9,10,13,14,15)]
  
  return(gmat)
}
