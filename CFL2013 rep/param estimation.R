# estimate parameters
param.est <- function(delta){

  x <- CFL_data
  #
  beta <- delta[1]
  theta <- delta[2]
  rho <- delta[3]
  #
  vc <- vfn0
  #
  cons <- x[, 2]  # consumption growth
  asset <- x[, 7:12] # asset returns
  
  lead_cons <-  cons[-1] # remove the first obsevation (1952Q1)
  lead_asset <- asset[-1,]  # remove the first obsevation (1952Q1)
  #
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