###############################################################################
###############################################################################
###
###                          Algoritmo REML HFHmodel                          
###                                                                           
###                              Mayo 2025                                
###


### Autor: Lola Esteban - Agustin Perez - Esteban Cabello


### Editado por:
### con fecha: 




REML.HFH.BoxCox <- function(X, W, yd, D, sigmaed2, eta, lambda, MAXITER = 100, precision = 10^-4) {
  
  
  X <- cbind(1, X)
  tX <- t(X)
  W <- cbind(1, W)
  p <- ncol(X)
  q <- ncol(W)
  
  eta.ini <- rep(0,q)
  eta.ini[1:length(eta)] <- eta
  
  eta.gorro.prev <- matrix(eta.ini, nrow=q)
  #eta.gorro.prev <- matrix(c(0, 0), nrow=q)
  if(lambda == 0){
    for(ITER in 1:MAXITER){
      
      sigmaud2 <- as.vector(exp(W%*%eta.gorro.prev))
      
      # Creación matriz de varianzas (diagonal)
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
      
      # Se repite el bucle i.q debido a que se necesitan todos los Vi calculados para usarlos en los Sij
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
      
      eta.gorro <- eta.gorro.prev + dif
      eta.gorro.prev <- eta.gorro
      
      # Criterio de parada
      if(identical(as.numeric(abs(dif)<precision),rep(1,q))){
        sigmaud2.gorro <- as.vector(exp(W%*%eta.gorro))
        
        V.gorro <- diag(as.vector(sigmaud2.gorro + sigmaed2))
        V.gorro.inv <- solve(V.gorro)
        Q.gorro <- solve(tX%*%V.gorro.inv%*%X)
        
        
        beta.gorro <- Q.gorro%*%tX%*%V.gorro.inv%*%yd
        
        tau.gorro <- matrix(c(as.vector(beta.gorro), as.vector(eta.gorro)), ncol=1)
        
        
        mud.gorro <- (as.vector(sigmaud2.gorro)/(as.vector(sigmaud2.gorro) + sigmaed2))*yd + (sigmaed2/(as.vector(sigmaud2.gorro) + sigmaed2))*X%*%beta.gorro
        
        det.tXX <- log(det(tX%*%X))
        det.V <- log(det(V.gorro))
        det.Q <- log(det(Q.gorro))
        det.tyPy <- log(det(t(yd)%*%P%*%yd))
        
        log.like <- -0.5*(D-p)*log(2*pi) + 0.5*det.tXX - 0.5*det.V - 0.5*det.Q - 0.5*det.tyPy
        
        break
      }
      
      
    } 
  }
  
  
  else{
    for(ITER in 1:MAXITER){
      
      sigmaud2 <- as.vector((lambda*W%*%eta.gorro.prev+1)^(1/lambda))
      
      # Creación matriz de varianzas (diagonal)
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
      
      # Se repite el bucle i.q debido a que se necesitan todos los Vi calculados para usarlos en los Sij
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
      
      eta.gorro <- eta.gorro.prev + dif
      eta.gorro.prev <- eta.gorro
    
      # Criterio de parada
      if(identical(as.numeric(abs(dif)<precision),rep(1,q))){
        sigmaud2.gorro <- as.vector((lambda*W%*%eta.gorro+1)^(1/lambda))
        
        V.gorro <- diag(as.vector(sigmaud2.gorro + sigmaed2))
        V.gorro.inv <- solve(V.gorro)
        Q.gorro <- solve(tX%*%V.gorro.inv%*%X)
        
        
        beta.gorro <- Q.gorro%*%tX%*%V.gorro.inv%*%yd
        
        tau.gorro <- matrix(c(as.vector(beta.gorro), as.vector(eta.gorro)), ncol=1)
        
        
        mud.gorro <- (as.vector(sigmaud2.gorro)/(as.vector(sigmaud2.gorro) + sigmaed2))*yd + (sigmaed2/(as.vector(sigmaud2.gorro) + sigmaed2))*X%*%beta.gorro
        
        det.tXX <- log(det(tX%*%X))
        det.V <- log(det(V.gorro))
        det.Q <- log(det(Q.gorro))
        det.tyPy <- log(det(t(yd)%*%P%*%yd))
        
        log.like <- -0.5*(D-p)*log(2*pi) + 0.5*det.tXX - 0.5*det.V - 0.5*det.Q - 0.5*det.tyPy
        
        break
      }
    
      
      
    } 
    
  } 
  
  
  return(list(param = tau.gorro, eblups = mud.gorro, ITER = ITER, Hij.inv = Hij.inv, lambda = lambda, loglike = log.like))
  
  
}
