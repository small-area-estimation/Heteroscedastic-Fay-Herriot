###############################################################################
###############################################################################
###
###                      Analytical MSE HFH              
###                                                                                                        
###


### Authors: Lola Esteban-Lefler - Agustín Pérez Martín - Esteban Cabello

library(pracma)

mse0.HFH <- function(X,W,eta,H.inv,V.inv,sigmaud2.hat, sigmaed2){
  
  D <- nrow(X)
  tX <- t(X)
  tW <- t(W)
  Q <- solve(tX%*%V.inv%*%X)
  
  
  g1 <- sigmaud2.hat*sigmaed2/(sigmaud2.hat+sigmaed2)
  g2 <- avar.sigmaud2 <- c()
  for(d in 1:D){
    g2[d] <- sigmaed2[d]^2/(sigmaud2.hat[d]+sigmaed2[d])^2*(X[d,]%*%Q%*%matrix(tX[,d]))
    avar.sigmaud2[d] <- exp(2*W[d,]%*%eta)*(exp(2*W[d,]%*%H.inv%*%matrix(tW[,d]))-exp(W[d,]%*%H.inv%*%matrix(tW[,d])))
  }

  g3 <- sigmaed2^2/(sigmaud2.hat+sigmaed2)^3*avar.sigmaud2
  
  mse0 <- g1+g2+2*g3
  
  return(list(mse = mse0, g1g2 = g1+g2, g3 = 2*g3))
  
}
