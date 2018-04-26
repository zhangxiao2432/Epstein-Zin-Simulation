# value function estimation
vfn.est <- function(tet){

  x <- CFL_data
  
  # Preference Parameter
  beta <- t0[1]
  theta <- t0[2]
  rho <- t0[3]
  
  cons <- x[, 2]  # consumption growth
  asset <- x[, 7:12] # asset returns
  
  lead_cons <-  cons[-1] # remove the first obsevation (1952Q1)
  lead_asset <- asset[-1,]  # remove the first obsevation (1952Q1)
  
  # Cardinal Spline for consumption
  source("basis function.R")
  cons_1 <- Cardinal_B((cons-min(cons))/(max(cons)-min(cons))*5+1-1+0.5,3)
  cons_2 <- Cardinal_B((cons-min(cons))/(max(cons)-min(cons))*5+1-2,3)
  cons_3 <- Cardinal_B((cons-min(cons))/(max(cons)-min(cons))*5+1-3-0.5,3)
  bs_cons <- as.matrix(cbind(cons_1, cons_2, cons_3))
  
  kk=m*m # number of additive term
  vc_0 <- tet[1]
  
  bs_vc0 <- matrix(c(Cardinal_B(vc_0+0.5, 3),Cardinal_B(vc_0-0.5, 3), Cardinal_B(vc_0-1.5, 3)), 1,3)  
  vc <- matrix(0, nrow(x),1)
  bs_vc <- bs_vc0
  
  tempt1 <-tensor.prod.model.matrix(list(t(as.matrix(bs_cons[1,])), as.matrix(bs_vc0)))
  tempt2 <-tet[2]+as.vector(tet[3:(kk+2)])%*%t(tempt1)

  vc[1] <- tempt2
  remove("tempt1", "tempt2")
  
  for (i in 2:nrow(x)){
    bs_vci <- matrix(c(Cardinal_B(vc[i-1]+0.5, 3),Cardinal_B(vc[i-1]-0.5, 3), Cardinal_B(vc[i-1]-1.5, 3)), 1,3)  
    tempt3 <- tensor.prod.model.matrix(list(t(as.matrix(bs_cons[i,])), as.matrix(bs_vci)))
    tempt4 <- tet[2]+as.vector(tet[3: (kk+2)])%*%t(tempt3)
    
    vc[i] <- tempt4
    bs_vc <- rbind(bs_vc, bs_vci)
    remove("tempt3","tempt4")
  }
  
  bs_vcn <- matrix(c(Cardinal_B(vc[nrow(x)]-1+1.5, 3),Cardinal_B(vc[nrow(x)]-2+1.5, 3), Cardinal_B(vc[nrow(x)]-3+1.5, 3)), 1,3) 
  bs_vc <- rbind(bs_vc, bs_vcn)
  bs_vc <- bs_vc[-1,]
  
  lead_vc <- vc[-1]
  lag_vc <- vc[-length(vc)]
 
  # for rho=1, the risk factor takes a limit form
  if (rho >1-1e-4 & rho < 1+1e-4){
    risk_fc <- lead_vc*lead_cons/(lag_vc)^(1/beta)
  } else {risk_fc <- lead_vc*lead_cons/(1/beta*((lag_vc)^(1-rho)-1+beta)^(1/(1-rho)))
  }
  
  condi_mt <- as.matrix(beta*(lead_cons)^(-rho)*risk_fc^(rho-theta)*lead_asset-1) # conditional moment
  
  # instrument for the conditional moments
  iv <- x[-213,2:5]  # instruments variables for step 1 -- [consgrowth cay rrel spex] starting from 1952Q1 --w_t
  tempt5<- tensor.prod.model.matrix(list(iv,iv))
  tempt6 <- cbind(1,2,3, 4, 7, 11, 12, 5, 9, 10, 13, 14, 15) # remove consgrowth-cross terms, RREL-cross terms,  and repeated cross products
  tempt7 <- tempt5[,-tempt6]
  P <- cbind(1, iv, tempt7)    # P is an instrument matrix of 15 variables.
  remove("tempt5","tempt6","tempt7")
  
  # P <- cbind(1, iv)
  gmat <- tensor.prod.model.matrix(list(P,condi_mt))
  
  gbar <- as.matrix(colMeans(gmat))
  w <- crossprod(gmat)/NROW(gmat) 
  #GMT <- crossprod(gbar,solve(w,gbar))
  GMT <- crossprod(gbar)
  
  return(GMT) 
}
