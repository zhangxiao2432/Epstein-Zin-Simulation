vfn.hat<- function(tet){
  
  x <- CFL_data
  
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
    tempt1 <- tensor.prod.model.matrix(list(t(as.matrix(bs_cons[i,])), as.matrix(bs_vci)))
    tempt2 <- tet[2]+as.vector(tet[3: (2+kk)])%*%t(tempt1)
    
    vc[i] <- tempt2
    bs_vc <- rbind(bs_vc, bs_vci)
    remove("tempt1", "tempt2")
  }
  
  return(vc) 
}