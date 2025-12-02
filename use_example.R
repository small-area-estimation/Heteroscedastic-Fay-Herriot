###############################################################################
###############################################################################
###
###                   HFH-BC(0) model - Use example           
###                                                            
###                                                    
###
### Authors: Esteban Cabello, Lola Esteban and Agustin Perez
### e-mail: ecabello@umh.es


library(sae)

data <- read.csv("data.csv", header = T, sep = ",")


##### Initial values eta (standard FH model)

FH <- eblupFH(yd ~ X, sigmaed2, method = "REML", MAXITER = 100, PRECISION = 0.0001, 
              B = 0, data = data)

eta.FH <- log(FH$fit$refvar)


sigmaed2 <- data$sigmaed2
yd <- data$yd
X <- model.matrix(yd~X, data = data)

V.inv.FH <- diag(1/(FH$fit$refvar+sigmaed2)) ### V^{-1} matrix
Vu <- diag(rep(FH$fit$refvar,nrow(data))) ### V_u matrix
random.eff.FH <- Vu%*%V.inv.FH%*%(yd-X%*%FH$fit$estcoef$beta) ## Area-level random effects


## Table for Regression Parameters (Beta's)
alpha <- 0.05
se.b <-  FH$fit$estcoef$std.error
lim.inf <- FH$fit$estcoef$beta - qnorm(1-alpha/2)*se.b
lim.sup <- FH$fit$estcoef$beta + qnorm(1-alpha/2)*se.b
pv <- FH$fit$estcoef$pvalue
coefBeta <- cbind(FH$fit$estcoef$beta, se.b, lim.inf, lim.sup, pv)
colnames(coefBeta) <- c("Beta", "std.error", "lim.inf", "lim.sup",
                        "p.value")
coefBeta


resid.FH.std <- (yd - FH$eblup)/sqrt(sigmaed2) ## Standardized residuals (sresiduals)

mseFH <- mseFH(yd ~ X, sigmaed2, method = "REML", MAXITER = 100, PRECISION = 0.0001, 
               B = 0, data = data)

loglike.FH <- mseFH$est$fit$goodness[1]
AIC.FH <- -2*loglike.FH + 2*(2+1) ## 1 explanatory variable + Intercept + 1 variance component



###### HFH-BC(0) model 

source("REMLHFH.BoxCox.R")

X <- model.matrix(yd ~ X, data)
W <- model.matrix(yd ~ W, data) 
D <- nrow(data)
eta <- c(eta.FH, 0)

model.HFH <- REML.HFH.BoxCox(X[,-1], W[,-1], yd, D, sigmaed2, eta, lambda = 0, MAXITER = 100, precision = 10^-4)

n.params.HFH <- length(model.HFH[[1]])

beta.fitted <- model.HFH[[1]][1:2]
eta.fitted <- model.HFH[[1]][3:4]
eblups.HFH <- model.HFH[[2]]
solve.H.eta <- -model.HFH[[4]]
V.inv.FH <- diag(1/(FH$fit$refvar+sigmaed2)) ### V^{-1} matrix

loglike.HFH <- model.HFH[[6]]
AIC.HFH <- -2*loglike.HFH + 2*n.params.HFH

loglike.HFH > loglike.FH
AIC.HFH < AIC.FH

########## Table for eta's ########
alpha <- 0.05
se.eta <- sqrt(diag(solve.H.eta))
t.val <- eta.fitted/se.eta
pv <- 2 * pnorm(as.vector(abs(t.val)), lower.tail = FALSE)
lim.inf <- eta.fitted - qnorm(1-alpha/2)*se.eta 
lim.sup <- eta.fitted + qnorm(1-alpha/2)*se.eta 
coefEta <- cbind(eta.fitted, se.eta, t.val, lim.inf,lim.sup,pv)

colnames(coefEta) = c("Eta's", "std.error", "t.statistics", "lim.inf", "lim.sup",
                      "p.value")
rownames(coefEta) <- colnames(W)
coefEta

## Table for Beta
sigmaud2.hat <- as.vector(exp(W%*%eta.fitted))
V.hat <- diag(sigmaud2.hat + sigmaed2)
solve.V.hat <- solve(V.hat)

se.b <-  sqrt(diag(solve(t(X)%*%solve(V.hat)%*%X)))
t.val <- beta.fitted/se.b
pv <- 2 * pnorm(as.vector(abs(t.val)), lower.tail = FALSE)
lim.inf <- beta.fitted - qnorm(1-alpha/2)*se.b
lim.sup <- beta.fitted + qnorm(1-alpha/2)*se.b
coefBeta <- cbind(beta.fitted, se.b, t.val, lim.inf, lim.sup, pv)
colnames(coefBeta) <- c("Beta", "std.error", "t.statistics", "lim.inf", "lim.sup",
                        "p.value")
coefBeta


### Random effects

ud <- sigmaud2.hat/(sigmaud2.hat+sigmaed2)*(yd-X%*%beta.fitted)


resid.HFH.std <- (yd - eblups.HFH)/sqrt(sigmaed2) ## Standardized residuals (sresiduals) HFH model


###### QQ-plots

par(mfrow = c(1,2))

qqnorm(resid.HFH.std, axes = F, cex.main=2, cex.lab = 1.6, main = "")
qqline(resid.HFH.std, col = "red")
title(main=expression(paste("Normal Q-Q Plot HFH-BC(0) sresiduals ")), cex.main=2)
axis(side = 1, col = "grey", cex.axis = 1.6)
axis(side = 2, col = "grey", cex.axis = 1.6)

qqnorm(resid.FH.std, axes = F, cex.main=2, cex.lab = 1.6, main = "")
qqline(resid.FH.std, col = "red")
title(main=expression(paste("Normal Q-Q Plot FH sresiduals")), cex.main=2)
axis(side = 1, col = "grey", cex.axis = 1.6)
axis(side = 2, col = "grey", cex.axis = 1.6)


###### Analytical MSE

mse.ana <- mse0.HFH(X, W, eta.fitted, solve.H.eta, solve.V.hat, sigmaud2.hat, sigmaed2)

##### Bootstrap-based MSE
do <- F

if(do == T){
  
  source("boot.mse.HFH.R")
  tau.gorro <-  c(beta.fitted, eta.fitted)
  
  set.seed(03202)
  mse.boot.HFH <- BOOT.HFH(X[,-1], W[,-1], yd, D, sigmaed2, tau.gorro, data, MAXITER = 100, precision = 10^-4, BB=500)
  rmse.boot <- sqrt(mse.boot.HFH[[1]])

}

