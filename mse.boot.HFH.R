###############################################################################
###############################################################################
###
###                   Bootstrap-based MSE HFH
##                                                         
###

### Authors: María Dolores Esteban Lefler, Agustín Pérez Martín, Esteban Cabello García



mse.boot.HFH <- function(X, W, yd, D, sigmaed2, tau.hat, MAXITER = 100, precision = 10^-4, BB=50){


  # Initial values parameters 

  difmud.b <- difmudFH.b <- difyd.b <- 0 
  mud.hat.b  <- matrix(0, ncol = BB, nrow = D)

  X.1 <- cbind(1, X)
  beta.hat <- tau.hat[1:4]
  eta.hat <- tau.hat[5:6]
  W.1 <- cbind(1, W)
  sigmaud2.hat <- as.vector(exp(W.1%*%eta.hat))
  
  
  b <- BadTot_b <- BadTot2_b <- 0
  
  # Start bootstrap loop (b)
  while(b<BB){
    b <- b+1
    
    cat("Bootstrap replicate", b, "\n")
    
    # Simulating e_d and u_d
    
    ed.b <- rnorm(D, 0, sqrt(sigmaed2))
    ud.b <- rnorm(D, 0, sqrt(sigmaud2.hat))
    
    
    # Simulating y_d and mu_d
    mud.b <- as.matrix(X.1%*%beta.hat + ud.b, ncol=1)
    
    yd.b <- mudb + ed.b
    

    # Bootstrap fitting
    fit.b <- try(REML.HFH.BoxCox(X, W, yd.b, D, sigmaed2, eta.hat, lambda = 0, MAXITER = 100, precision = 10^-4), TRUE)

    
    if(class(fit.b)=="try-error"){
      # write.table(data.frame(class(fit),D,b), file=paste0("WARNING_BOOTSTRAP B=", BB, ".txt"), append=TRUE, col.names=FALSE)
      cat("\t Muestra_bootstrap", b, " rechazada por try-error\n")
      b <- b-1
      BadTot_b <- BadTot_b + 1
    }
    else{
      if(fit.b[[3]] < 100) {  	# MAXITER <100 iterations (convergence)
        mud.hat.ast <- fit.b[[2]]
        eta.hat.ast <- as.vector(fit.b[[1]][3:4,1])
        sigmaud2.hat.ast <- as.vector(exp(W.1%*%eta.hat.ast))
        
        mudFH.hat.ast <- fit.FH.b[[1]]
        
        
        difmud.b <- difmud.b + (mud.hat.ast - mud.b)^2        
        difyd.b <- difyd.b + (yd.b -mud.b)^2

      }
      else {
        b <- b-1
        BadTot2_b <- BadTot2_b + 1
        # write.table(data.frame(class(fit), D, boot=b, BadTot2_b), file=paste0("WARNING_BOOTSTRAP B=", BB,".txt"), append=TRUE, col.names=FALSE)
        cat("\t Bootstrap_Sample", b, " rejected by MAXITER\n")
      }
      cat( "\n")
    }
  }# Loop BB
  
  mse_mud_ast <- difmud.b/BB    
  mse_yd_ast <- difyd.b/BB
  
  return(list(data.frame(EBLUP = mse_mud_ast, Direct = mse_yd_ast), 
              BadTot_b, BadTot2_b))
}
