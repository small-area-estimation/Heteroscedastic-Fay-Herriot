###############################################################################
###############################################################################
###
###                          REML HFH algorithm                        
###                                                                           
###                              May 2025                                
###


### Authors: Lola Esteban - Agustin Perez - Esteban Cabello




REML.HFH.BoxCox <- function(X, W, yd, D, sigmaed2, eta, lambda, MAXITER = 100, precision = 10^-4) {
  
  
  X <- cbind(1, X)
  tX <- t(X)
  W <- cbind(1, W)
  p <- ncol(X)
  q <- ncol(W)
  
  eta.ini <- rep(0,q)
  eta.ini[1:length(eta)] <- eta
  
  eta.hat.prev <- matrix(eta.ini, nrow=q)
  #eta.hat.prev <- matrix(c(0, 0), nrow=q)
  if(lambda == 0){
    for(ITER in 1:MAXITER){
      
      sigmaud2 <- as.vector(exp(W%*%eta.hat.prev))
      
      # Creation of variance matrix of the target variable (diagonal)
      V <- diag(sigmaud2 + sigmaed2)
      # cat(det(V),"\n")
      V.inv <- solve(V)
      
      Q <- solve(tX%*%V.inv%*%X)
      K <- V.inv %*% X %*% Q %*% tX %*% V.inv
      P <- V.inv - K
      
      Vi <- Pi <- Si <- Ki <- list()
      for(i.q in 1:q){
        Vi[[i.q]] <- diag(sigmaud2*W[,i.q])
        Pi[[i.q]] <- -V.inv%*%Vi[[i.q]]%*%V.inv + V.inv%*%Vi[[i.q]]%*%K - K%*%Vi[[i.q]]%*%K + K%*%Vi[[i.q]]%*%V.inv
        Si[[i.q]] <- 0.5*( -sum(diag(V.inv%*%Vi[[i.q]])) + sum(diag(K%*%Vi[[i.q]])) - t(yd)%*%Pi[[i.q]]%*%yd )
        Ki[[i.q]] <- -V.inv%*%Vi[[i.q]]%*%K + K%*%Vi[[i.q]]%*%K - K%*%Vi[[i.q]]%*%V.inv
      }
      
      # The i.q loop is repeated because all the calculated Vi values are needed for use in the Sij values.
      Vij <- S1ij <- S2ij <- S3ij <- list()
      P1ij <- P2ij <- P3ij <- list()
      Hij <- list()
      for(i.q in 1:q){
        Vij[[i.q]] <- S1ij[[i.q]] <- S2ij[[i.q]] <- S3ij[[i.q]] <- list()
        P1ij[[i.q]] <- P2ij[[i.q]] <- P3ij[[i.q]] <- list()
        Hij[[i.q]] <- list()
        for(j.q in 1:q){
          Vij[[i.q]][[j.q]] <- diag(sigmaud2*W[,i.q]*W[,j.q])
          S1ij[[i.q]][[j.q]] <- -sum(diag(V.inv%*%Vi[[j.q]]%*%V.inv%*%Vi[[i.q]])) + sum(diag(V.inv%*%Vij[[i.q]][[j.q]]))
          S2ij[[i.q]][[j.q]] <- sum(diag(Ki[[j.q]]%*%Vi[[i.q]])) + sum(diag(K%*%Vij[[i.q]][[j.q]]))
          
          P1ij[[i.q]][[j.q]] <- V.inv%*%Vij[[i.q]][[j.q]]%*%V.inv - 2*V.inv%*%Vi[[j.q]]%*%V.inv%*%Vi[[i.q]]%*%V.inv
          P2ij[[i.q]][[j.q]] <- -V.inv%*%Vi[[j.q]]%*%V.inv%*%Vi[[i.q]]%*%K + V.inv%*%Vij[[i.q]][[j.q]]%*%K + V.inv%*%Vi[[i.q]]%*%Ki[[j.q]]
          P3ij[[i.q]][[j.q]] <- Ki[[j.q]]%*%Vi[[i.q]]%*%K + K%*%Vij[[i.q]][[j.q]]%*%K + K%*%Vi[[i.q]]%*%Ki[[j.q]]
          S3ij[[i.q]][[j.q]] <- t(yd)%*%(-P1ij[[i.q]][[j.q]] + 2*P2ij[[i.q]][[j.q]] - P3ij[[i.q]][[j.q]])%*%yd
          
          Hij[[i.q]][[j.q]] <- 0.5*(-S1ij[[i.q]][[j.q]] + S2ij[[i.q]][[j.q]] -S3ij[[i.q]][[j.q]])
        }
      }
      
      Hij <- matrix(unlist(Hij), ncol=q)
      Hij.inv <- solve(Hij)
      
      dif <- -Hij.inv%*%matrix(unlist(Si), ncol=1)
      
      eta.hat <- eta.hat.prev + dif
      eta.hat.prev <- eta.hat
      
      # Stop criterion
      if(identical(as.numeric(abs(dif)<precision),rep(1,q))){
        sigmaud2.hat <- as.vector(exp(W%*%eta.hat))
        
        V.hat <- diag(as.vector(sigmaud2.hat + sigmaed2))
        V.hat.inv <- solve(V.hat)
        Q.hat <- solve(tX%*%V.hat.inv%*%X)
        
        
        beta.hat <- Q.hat%*%tX%*%V.hat.inv%*%yd
        
        tau.hat <- matrix(c(as.vector(beta.hat), as.vector(eta.hat)), ncol=1)
        
        
        mud.hat <- (as.vector(sigmaud2.hat)/(as.vector(sigmaud2.hat) + sigmaed2))*yd + (sigmaed2/(as.vector(sigmaud2.hat) + sigmaed2))*X%*%beta.hat
        
        det.tXX <- log(det(tX%*%X))
        det.V <- log(det(V.hat))
        det.Q <- log(det(Q.hat))
        det.tyPy <- log(det(t(yd)%*%P%*%yd))
        
        log.like <- -0.5*(D-p)*log(2*pi) + 0.5*det.tXX - 0.5*det.V - 0.5*det.Q - 0.5*det.tyPy
        
        break
      }
      
      
    } 
  }
  
  
  else{
    for(ITER in 1:MAXITER){
      
      sigmaud2 <- as.vector((lambda*W%*%eta.hat.prev+1)^(1/lambda))
      
      # Creation of variance matrix of the target variable (diagonal)
      V <- diag(sigmaud2 + sigmaed2)
      # cat(det(V),"\n")
      V.inv <- solve(V)
      
      Q <- solve(tX%*%V.inv%*%X)
      K <- V.inv %*% X %*% Q %*% tX %*% V.inv
      P <- V.inv - K
      
      Vi <- Pi <- Si <- Ki <- list()
      for(i.q in 1:q){
        Vi[[i.q]] <- diag(sigmaud2^(1-lambda)*W[,i.q])
        Pi[[i.q]] <- -V.inv%*%Vi[[i.q]]%*%V.inv + V.inv%*%Vi[[i.q]]%*%K - K%*%Vi[[i.q]]%*%K + K%*%Vi[[i.q]]%*%V.inv
        Si[[i.q]] <- 0.5*( -sum(diag(V.inv%*%Vi[[i.q]])) + sum(diag(K%*%Vi[[i.q]])) - t(yd)%*%Pi[[i.q]]%*%yd )
        Ki[[i.q]] <- -V.inv%*%Vi[[i.q]]%*%K + K%*%Vi[[i.q]]%*%K - K%*%Vi[[i.q]]%*%V.inv
      }
      
      # The i.q loop is repeated because all the calculated Vi values are needed for use in the Sij values.
      Vij <- S1ij <- S2ij <- S3ij <- list()
      P1ij <- P2ij <- P3ij <- list()
      Hij <- list()
      for(i.q in 1:q){
        Vij[[i.q]] <- S1ij[[i.q]] <- S2ij[[i.q]] <- S3ij[[i.q]] <- list()
        P1ij[[i.q]] <- P2ij[[i.q]] <- P3ij[[i.q]] <- list()
        Hij[[i.q]] <- list()
        for(j.q in 1:q){
          Vij[[i.q]][[j.q]] <- diag((1-lambda)*(sigmaud2^(1-2*lambda)*W[,i.q]*W[,j.q]))
          S1ij[[i.q]][[j.q]] <- -sum(diag(V.inv%*%Vi[[j.q]]%*%V.inv%*%Vi[[i.q]])) + sum(diag(V.inv%*%Vij[[i.q]][[j.q]]))
          S2ij[[i.q]][[j.q]] <- sum(diag(Ki[[j.q]]%*%Vi[[i.q]])) + sum(diag(K%*%Vij[[i.q]][[j.q]]))
          
          P1ij[[i.q]][[j.q]] <- V.inv%*%Vij[[i.q]][[j.q]]%*%V.inv - 2*V.inv%*%Vi[[j.q]]%*%V.inv%*%Vi[[i.q]]%*%V.inv
          P2ij[[i.q]][[j.q]] <- -V.inv%*%Vi[[j.q]]%*%V.inv%*%Vi[[i.q]]%*%K + V.inv%*%Vij[[i.q]][[j.q]]%*%K + V.inv%*%Vi[[i.q]]%*%Ki[[j.q]]
          P3ij[[i.q]][[j.q]] <- Ki[[j.q]]%*%Vi[[i.q]]%*%K + K%*%Vij[[i.q]][[j.q]]%*%K + K%*%Vi[[i.q]]%*%Ki[[j.q]]
          S3ij[[i.q]][[j.q]] <- t(yd)%*%(-P1ij[[i.q]][[j.q]] + 2*P2ij[[i.q]][[j.q]] - P3ij[[i.q]][[j.q]])%*%yd
          
          Hij[[i.q]][[j.q]] <- 0.5*(-S1ij[[i.q]][[j.q]] + S2ij[[i.q]][[j.q]] -S3ij[[i.q]][[j.q]])
        }
      }
      
      Hij <- matrix(unlist(Hij), ncol=q)
      Hij.inv <- solve(Hij)
      
      dif <- -Hij.inv%*%matrix(unlist(Si), ncol=1)
      
      eta.hat <- eta.hat.prev + dif
      eta.hat.prev <- eta.hat
      
      # Stop criterion
      if(identical(as.numeric(abs(dif)<precision),rep(1,q))){
        sigmaud2.hat <- as.vector((lambda*W%*%eta.hat+1)^(1/lambda))
        
        V.hat <- diag(as.vector(sigmaud2.hat + sigmaed2))
        V.hat.inv <- solve(V.hat)
        Q.hat <- solve(tX%*%V.hat.inv%*%X)
        
        
        beta.hat <- Q.hat%*%tX%*%V.hat.inv%*%yd
        
        tau.hat <- matrix(c(as.vector(beta.hat), as.vector(eta.hat)), ncol=1)
        
        
        mud.hat <- (as.vector(sigmaud2.hat)/(as.vector(sigmaud2.hat) + sigmaed2))*yd + (sigmaed2/(as.vector(sigmaud2.hat) + sigmaed2))*X%*%beta.hat
        
        det.tXX <- log(det(tX%*%X))
        det.V <- log(det(V.hat))
        det.Q <- log(det(Q.hat))
        det.tyPy <- log(det(t(yd)%*%P%*%yd))
        
        log.like <- -0.5*(D-p)*log(2*pi) + 0.5*det.tXX - 0.5*det.V - 0.5*det.Q - 0.5*det.tyPy
        
        break
      }
      
      
      
    } 
    
  } 
  
  
  return(list(param = tau.hat, eblups = mud.hat, ITER = ITER, 
              Hij.inv = Hij.inv, lambda = lambda, loglike = log.like))
  
  
}
