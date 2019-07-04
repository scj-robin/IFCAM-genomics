# KerFdr

rm(list=ls())
library(KernSmooth); library(mclust); library(Matrix)

# # Parms
# n = 1e5; pi0 = .9; pi1 = 1-pi0

# # Data
# Z = rbinom(n, 1, pi1)
# X = rep(0, n)
# X[which(Z==0)] = rnorm(sum(Z==0))
# X[which(Z==1)] = rgamma(sum(Z==1), 2, 1)
# X = sort(X)
# f.true = dgamma(X, 2, 1)
# P = pnorm(X, lower.tail=F)

load('../../../Data-NotUpload/TTests-LogExp.Rdata'); 
P = PE; 
load('../../../Data-NotUpload/BetaRegTests-Methylation.Rdata'); 
PM = PM[-which(PM==0)]; P = PM; 

source('F_KerFdr.R')
F_KerFdr(P)

