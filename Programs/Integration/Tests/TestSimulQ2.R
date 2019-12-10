rm(list=ls()); par(mfrow=c(1, 1), pch=20)
library(mvtnorm); library(ROCR)
setwd('/home/robin/Bureau/RECHERCHE/EXPRESSION/IFCAM/IFCAM-genomics/Programs/Integration/Tests')
source('../Functions/F_KerFdr.R')
source('../Functions/F_JointPvalueFromObserved.R')


# Dims
p = 1e4; alpha = .05

# Parms: class = 00, 01, 10, 11 -> K = 4
pi = c(.8, .08, .08, .04); 
delta = 2*c(1, 2)
rho0= .3; rho1 = .5

# Simulation of the test statistics
K = length(pi)
mu = matrix(c(0, 0, delta[1], 0, 0, delta[2], delta[1], delta[2]), K, 2, byrow=TRUE)
Sigma00 = matrix(rho0, 2, 2); diag(Sigma00) = 1
Sigma11 = matrix(rho1, 2, 2); diag(Sigma11) = 1
Sigma = list(Sigma00, Sigma00, Sigma00, Sigma11)
Z = as.vector(t(rmultinom(p, 1, pi))%*%(1:K))
Y = matrix(0, p, 2)
sapply(1:K, function(k){Y[which(Z==k), ] <<- rmvnorm(sum(Z==k), mean=mu[k, ], sigma=Sigma[[k]])})
# plot(Y[, 1], Y[, 2], col=Z); abline(h=0, v=0)

# Pvalues
P = pnorm(Y, lower.tail=FALSE)
# plot(P[, 1], P[, 2], col=Z); abline(h=0, v=0)
# hist(P[, 1], breaks=sqrt(p)); hist(P[, 2], breaks=sqrt(p))

###############################################################################
# Union-intersection
pmax = apply(P, 1, max)
boxplot(pmax ~ Z)
hist(pmax[which(Z < 4)], breaks=sqrt(p))
hist(pmax[which(Z == 4)], breaks=sqrt(p))
hist(pmax, breaks=sqrt(p))

# Multiple testing on Pmax
pmaxBH = p.adjust(pmax, method='fdr')
boxplot(pmaxBH ~ Z, title='FDR corrected Pmax')
abline(h = max(pmax[which(pmaxBH<alpha)]), col=2)

###############################################################################
# Marginal mixtures
Ker1 = F_KerFdr(P[, 1], rescaleKernWidth=TRUE)
Ker2 = F_KerFdr(P[, 2], rescaleKernWidth=TRUE)
crossProb = t(cbind(1-Ker1$tau, Ker1$tau))%*%cbind(1-Ker2$tau, Ker2$tau)/p

# Union testing
testTab = data.frame(stat = qnorm(apply(P, 1, max), lower.tail=FALSE), tau = 1 - (1-Ker1$tau)*(1-Ker2$tau))
res = F_JointPvalueFromObserved(testTab)
hist(res$pval, breaks=sqrt(p))
# hist(res$pval[which(Z < 4)], breaks=sqrt(p))
# hist(res$pval[which(Z == 4)], breaks=sqrt(p))
boxplot(res$pval ~ Z)

###############################################################################
# Multiple testing for our union pvalue
# BH
res$pBH = p.adjust(res$pval, method='fdr')
boxplot(qlogis(res$pval) ~ Z)
abline(h = qlogis(max(res$pval[which(res$pBH<alpha)])), col=2)

# Storey-like FdR 
orderPval = order(res$pval); rankPval = rank(res$pval); 
p0 = 1 - crossProb[2, 2]
storey = p0*p*res$pval[orderPval] / (1:p)
res$storey = storey[rankPval]

# Kerfdr for the union pvalue
Ker = F_KerFdr(res$pval, p0=1-crossProb[2, 2], rescaleKernWidth=TRUE)
orderTau = order(1-Ker$tau); rankTau = rank(1-Ker$tau) 
fdr = cumsum(1-Ker$tau[orderTau]) / (1:p)
res$kerfdr =  fdr[rankTau]
plot(res$pBH, res$kerfdr); abline(0, 1)

###############################################################################
# Alternative union p-value
# FdR
F_Ker2FdR <- function(P, Tau){
   orderP = order(P); rankP = rank(P) 
   FdR = cumsum(Tau[orderP]) / (1:length(P))
   return(FdR[rankP])
}
Ker1$FdR = F_Ker2FdR(P[, 1], Ker1$tau); Ker2$FdR = F_Ker2FdR(P[, 2], Ker2$tau)
# FdR
F_Ker2Cdf <- function(P, Tau){
   orderP = order(P); rankP = rank(P) 
   Cdf = cumsum(Tau[orderP]) / sum(Tau)
   return(Cdf[rankP])
}
Ker1$Cdf = F_Ker2Cdf(P[, 1], Ker1$tau); Ker2$Cdf = F_Ker2Cdf(P[, 2], Ker2$tau)
# # Cheating empirical FdR under H1
# eCdf1 = ecdf(P[which((Z==3) | (Z==4)), 1])(P[, 1])
# eCdf2 = ecdf(P[which((Z==2) | (Z==4)), 2])(P[, 2])
# # Cheating exact FdR under H1
# eCdf1 = 1 - pnorm((qnorm(1 - P[, 1])-delta[1]))
# eCdf2 = 1 - pnorm((qnorm(1 - P[, 2])-delta[2]))

# Pvalue
# PvalUnion = (crossProb[1, 1] * pmax^2 + crossProb[1, 2] * pmax * Ker2$FdR + 
#    crossProb[2, 1] * pmax * Ker1$FdR) / (1 - crossProb[2, 2])
# PvalUnion = (crossProb[1, 1] * pmax^2 + crossProb[1, 2] * pmax * eCdf2 + 
#                 crossProb[2, 1] * pmax * eCdf2) / (1 - crossProb[2, 2])
PvalUnion = (crossProb[1, 1] * pmax^2 + crossProb[1, 2] * pmax * Ker2$Cdf + 
                crossProb[2, 1] * pmax * Ker1$Cdf) / (1 - crossProb[2, 2])
hist(PvalUnion, breaks=sqrt(p))
boxplot(PvalUnion ~ Z)   

KKer = F_KerFdr(PvalUnion, p0=1-crossProb[2, 2], rescaleKernWidth=TRUE)
orderTau = order(1-KKer$tau); rankTau = rank(1-KKer$tau) 
fdr = cumsum(1-KKer$tau[orderTau]) / (1:p)
KerUnion =  fdr[rankTau]
plot(pmax, res$kerfdr, xlim=c(0, alpha)); abline(0, 1); 
points(pmax, res$pBH, col=3); points(pmax, res$storey, col=4); 
points(pmax, KerUnion, col=2); points(pmax, p.adjust(PvalUnion, method='fdr'), col=5)

table((Z==K), (pmaxBH<alpha))
table((Z==K), (res$pBH<alpha))
table((Z==K), (res$storey<alpha))
table((Z==K), (res$kerfdr<alpha))
table((Z==K), (p.adjust(PvalUnion, method='fdr')<alpha))
table((Z==K), (KerUnion<alpha))

###############################################################################
# Estimate of the alternative (i.e. intersection of H1's) Cdf of Pmax
eCdfPmax = ecdf(pmax)(pmax)
eCdfPmaxH0 = PvalUnion
eCdfPmaxH1 = F_Ker2Cdf(pmax, Ker1$tau*Ker2$tau)
trueCdfPmaxH1 = ecdf(pmax[which(Z==4)])(pmax)
plot(pmax, eCdfPmax); points(pmax, eCdfPmaxH0, col=2)
# points(pmax, (eCdfPmax - (1 - crossProb[2, 2])*eCdfPmaxH0)/crossProb[2, 2], col=4)
points(pmax, trueCdfPmaxH1, col=8)
# points(pmax, (1-crossProb[2, 2])*eCdfPmaxH0 + crossProb[2, 2]*trueCdfPmaxH1, col=6)
points(pmax, eCdfPmaxH1, col=5)
points(pmax, (1-crossProb[2, 2])*eCdfPmaxH0 + crossProb[2, 2]*eCdfPmaxH1, col=6)

