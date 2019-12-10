rm(list=ls()); par(mfrow=c(1, 1), pch=20)
library(mvtnorm); library(ROCR); library(gtools); library(ROCR)
setwd('/home/robin/Bureau/RECHERCHE/EXPRESSION/IFCAM/IFCAM-genomics/Programs/Integration/Tests')
source('../Functions/F_KerFdr.R')
source('../Functions/F_JointPvalueFromObserved.R')

# Dims
p = 1e4; alpha = .05; Q = 5

# Z configurations
Hconfig = as.matrix(expand.grid(lapply(1:Q, function(q){c(0, 1)})))
Hnumber = 2^Q

# Parms: class = 00, 01, 10, 11 -> K = 4
piq = rdirichlet(Q, c(8, 2))
pi = sapply(1:Hnumber, function(h){prod((piq[,1]**(1-Hconfig[h, ]))*(piq[,2]**Hconfig[h, ]))})
pi[Hnumber] = 1e-2
pi = rdirichlet(1, p*pi)
delta = sqrt(3*(1:Q))
rho0= 0.2; rho1 = .5

# Simulation of the test statistics
K = length(pi)
mu = Hconfig * (rep(1, Hnumber)%o%delta)
Sigma0 = matrix(rho0, Q, Q); diag(Sigma0) = 1
Sigma1 = matrix(rho1, Q, Q); diag(Sigma1) = 1
Sigma = lapply(1:Hnumber, function(h){
   Res <- Sigma0
   Res[which(Hconfig[h, ]==1), which(Hconfig[h, ]==1)] <- Sigma1[which(Hconfig[h, ]==1), which(Hconfig[h, ]==1)]
   return(Res)
})
Z = as.vector(t(rmultinom(p, 1, pi))%*%(1:K))
Y = matrix(0, p, Q)
sapply(1:K, function(k){Y[which(Z==k), ] <<- rmvnorm(sum(Z==k), mean=mu[k, ], sigma=Sigma[[k]])})
P = pnorm(Y, lower.tail=FALSE)

###############################################################################
# Union-intersection test
pmax = apply(P, 1, max)
# boxplot(pmax ~ Z); hist(pmax, breaks=sqrt(p))
pmaxBH = p.adjust(pmax, method='fdr')

###############################################################################
# Marginal mixtures
Ker = lapply(1:Q, function(q){F_KerFdr(P[, q], rescaleKernWidth=TRUE)})
crossProb = sapply(1:Hnumber, function(h){
   Res = sapply(1:Q, function(q){Hconfig[h, q]*Ker[[q]]$tau + (1-Hconfig[h, q])*(1-Ker[[q]]$tau)})
   mean(apply(Res, 1, prod))
})
F_Ker2Cdf <- function(P, Tau){
   orderP = order(P); rankP = rank(P) 
   Cdf = cumsum(Tau[orderP]) / sum(Tau)
   return(Cdf[rankP])
}
Cdf = lapply(1:Q, function(q){F_Ker2Cdf(P[, q], Ker[[q]]$tau)})
# trueCdf = lapply(1:Q, function(q){1 - pnorm((qnorm(1 - P[, q])-delta[q]))}) 
# par(mfrow=c(2, 2)); sapply(1:Q, function(q){plot(P[, q], Cdf[[q]], log='y'); points(P[, q], trueCdf[[q]], col=2)})

# ###############################################################################
# # Old version : naive estimate of the joint H1 cdf
# CdfProduct = sapply(1:(Hnumber-1), function(h){
#    cdf <- rep(1, p)
#    sapply(1:Q, function(q){cdf <<- cdf * ((1-Hconfig[h, q])*pmax + Hconfig[h, q]*Cdf[[q]])})
#    return(cdf)
# })

###############################################################################
# Estimated alternative cdf for each Hconfig
CdfH1 = CdfProduct = matrix(0, p, (Hnumber-1))
sapply(1:(Hnumber-1), function(h){
   cat(h, '')
   if(sum((Hconfig[h, ]==1)) > 0){
      pmax.h = apply(P[, which(Hconfig[h, ]==1), drop=FALSE], 1, max)
      tau.h = rep(1, p)
      sapply(which(Hconfig[h, ]==1), function(q){tau.h <<- tau.h * Ker[[q]]$tau})
      CdfH1[, h] <<- F_Ker2Cdf(pmax.h, tau.h)
      CdfProduct[, h] <<- CdfH1[, h] * (pmax^sum(Hconfig[h, ]==0))
   }else{CdfProduct[, h] <<- pmax^Q}
})

###############################################################################
# Alternative union p-value
PvalUnion = CdfProduct%*%crossProb[-Hnumber] / (1 - crossProb[Hnumber])
boxplot(PvalUnion ~ Z)   

KKer = F_KerFdr(PvalUnion, p0=1-crossProb[Hconfig], rescaleKernWidth=TRUE, plotting=TRUE)
orderTau = order(1-KKer$tau); rankTau = rank(1-KKer$tau) 
fdr = cumsum(1-KKer$tau[orderTau]) / (1:p)
KerUnion =  fdr[rankTau]

par(mfrow=c(2, 2), pch=20, mex=.5)
hist(PvalUnion, breaks=sqrt(p))

plot(performance(prediction(PvalUnion, (Z==Z^Q)), "tpr", "fpr"))
plot(performance(prediction(pmax, (Z==Z^Q)), "tpr", "fpr"), add=TRUE, col=2)
abline(0, 1)

plot(performance(prediction(PvalUnion, (Z==Z^Q)), "prec", "rec"))
plot(performance(prediction(pmax, (Z==Z^Q)), "prec", "rec"), add=TRUE, col=2)

table((Z==K), (pmaxBH<alpha))
table((Z==K), (p.adjust(PvalUnion, method='fdr')<alpha))
table((Z==K), (KerUnion<alpha))

