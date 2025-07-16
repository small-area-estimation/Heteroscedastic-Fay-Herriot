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


##### Initial values eta (univariant FH model)

FH <- eblupFH(yd ~ X, sigmaed2, method = "REML", MAXITER = 100, PRECISION = 0.0001, 
              B = 0, data = data)

eta.FH <- log(FH$fit$refvar)


sigmaed2 <- data$sigmaed2
yd <- data$yd
X <- model.matrix(yd~X, data = data)

V.inv.FH <- diag(1/(FH$fit$refvar+sigmaed2)) ### V^{-1} matrix
Vu <- diag(rep(FH$fit$refvar,nrow(data))) ### V_u matrix
random.eff.FH <- Vu%*%V.inv.FH%*%(yd-X%*%FH$fit$estcoef$beta) ## area-level random effects


## Table for Regression Parameters (Beta's)
se.b <-  FH$fit$estcoef$std.error
lim.inf <- FH$fit$estcoef$beta - qnorm(1-0.05/2)*se.b
lim.sup <- FH$fit$estcoef$beta + qnorm(1-0.05/2)*se.b
pv <- FH$fit$estcoef$pvalue
coefBeta <- cbind(FH$fit$estcoef$beta, se.b, lim.inf, lim.sup, pv)
colnames(coefBeta) <- c("Beta", "std.error", "lim.inf", "lim.sup",
                        "p.value")
coefBeta


resid.FH.std <- (yd - FH$eblup)/sqrt(sigmaed2)
hist(resid.FH.std, breaks = 15, freq = F, main = "Histogram", xlab = "sresiduals", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
lines(density(resid.FH.std), col = "red")

boxplot(resid.FH.std, cex.lab = 1.6, cex.axis = 1.8, axes = F)
title(main=expression(paste("Boxplot of sresiduals (FH)")), cex.main=2)
axis(2, col = "gray70", cex.axis=1.5)


shapiro.test(resid.FH.std)
nortest::lillie.test(resid.FH.std)


library(ggpubr)
ggqqplot(data = data.frame(resid.FH.std), x = "resid.FH.std",
         title = "Normal Q-Q Plot FH sresiduals",
         color= c("#0073C2FF"),
         ylab = "sresiduals",
         ggtheme = theme_classic()) 

qqnorm(resid.FH.std, axes = F, cex.main=2, cex.lab = 1.6, main = "")
qqline(resid.FH.std, col = "red")
title(main=expression(paste("Normal Q-Q Plot FH sresiduals")), cex.main=2)
axis(side = 1, col = "grey", cex.axis = 1.6)
axis(side = 2, col = "grey", cex.axis = 1.6)




### Por aqui




mseFH <- mseFH(Y1 ~ EDU5 + R5000 + NAC2, varY1, method = "REML", MAXITER = 100, PRECISION = 0.0001, 
               B = 0, data = data)

loglike.FH <- mseFH$est$fit$goodness[1]
AIC.FH <- -2*loglike.FH + 2*(4+1)
rmseFH.analytic <- sqrt(mseFH$mse)


###### HFH-BC(0) model 

source("REMLHFH.BoxCox.R")

X <- model.matrix(Y1 ~ EDU5 + R5000 + NAC2, data)
W <- model.matrix(Y1 ~ EDU3, data) 
D <- nrow(data)
eta <- c(eta.FH,0)

model.0 <- REML.HFH.BoxCox(X[,-1], W[,-1], yd, D, sigmaed2, eta, lambda = 0, MAXITER = 100, precision = 10^-4)

n.params.0 <- length(model.0[[1]])

beta.fitted <- model.0[[1]][1:4]
eta.fitted <- model.0[[1]][5:6]
eblups.0 <- model.0[[2]]
solve.H.eta <- -model.0[[4]]

loglike.HFH <- model.0[[6]]
AIC.HFH <- -2*loglike.HFH + 2*n.params.0

loglike.HFH > loglike.FH
AIC.HFH < AIC.FH

########## Table for eta's ########
alpha <- 0.05
se.eta <- sqrt(diag(solve.H.eta))
t.val <- eta.fitted/se.eta
pv <- 2 * pnorm(as.vector(abs(t.val)), lower.tail = FALSE)
lim.inf <- eta.fitted - qnorm(1-0.05/2)*se.eta 
lim.sup <- eta.fitted + qnorm(1-0.05/2)*se.eta 
coefEta <- cbind(eta.fitted, se.eta, t.val, lim.inf,lim.sup,pv)

colnames(coefEta) = c("Eta's", "std.error", "t.statistics", "lim.inf","lim.sup",
                      "p.value")
rownames(coefEta) <- colnames(W)
coefEta

## Table for Beta
sigmaud2.hat <- as.vector(exp(W%*%eta.fitted))
V.hat <- diag(sigmaud2.hat + sigmaed2)

se.b <-  sqrt(diag(solve(t(X)%*%solve(V.hat)%*%X)))
t.val <- beta.fitted/se.b
pv <- 2 * pnorm(as.vector(abs(t.val)), lower.tail = FALSE)
lim.inf <- beta.fitted - qnorm(1-0.05/2)*se.b
lim.sup <- beta.fitted + qnorm(1-0.05/2)*se.b
coefBeta <- cbind(beta.fitted, se.b, t.val, lim.inf, lim.sup, pv)
colnames(coefBeta) <- c("Beta", "std.error", "t.statistics", "lim.inf", "lim.sup",
                        "p.value")
coefBeta


xtable(coefEta,digits = 4)
xtable(coefBeta,digits = 4)

### Random effects

ud <- sigmaud2.hat/(sigmaud2.hat+sigmaed2)*(yd-X%*%beta.fitted)

###### Model validation


resid.std.0 <- (yd-eblups.0)/sqrt(sigmaed2)
hist(resid.std.0, breaks = 15, freq = F, main = "Histogram", xlab = "sresiduals", cex.main = 1.8, cex.lab = 1.6, cex.axis = 1.8)
lines(density(resid.std.0), col = "red")

boxplot(resid.std.0, cex.lab = 1.6, cex.axis = 1.8, axes = F)
title(main=expression(paste("Boxplot of sresiduals (HFH)")), cex.main=2)
axis(2, at = seq(-2,2,by = 1), col = "gray70", cex.axis=1.5)

# shapiro.test(resid.st)
# nortest::lillie.test(resid.st)


library(ggpubr)
ggqqplot(data = data.frame(resid.std.0), x = "resid.std.0",
         title = "Normal Q-Q Plot",
         color= c("#0073C2FF"),
         ylab = "sresiduals",
         ggtheme = theme_classic()
) + theme(text = element_text(size=18), plot.title = element_text(hjust = 0.5)) 





qqnorm(resid.std.0, axes = F, cex.main=2, cex.lab = 1.6, main = "")
qqline(resid.std.0, col = "red")
title(main=expression(paste("Normal Q-Q Plot HFH-BC(0) sresiduals ")), cex.main=2)
axis(side = 1, col = "grey", cex.axis = 1.6)
axis(side = 2, col = "grey", cex.axis = 1.6)

shapiro.test(resid.std.0)
shapiro.test(resid.FH.std)
