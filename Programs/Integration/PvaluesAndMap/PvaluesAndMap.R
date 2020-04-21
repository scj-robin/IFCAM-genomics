# Analysis of pvalues & map data

rm(list=ls())
library(ks); library(mclust); library(spatstat)
source('../2020/Functions.R')

# Data
dataDir <- '../../../Results/PvaluesAndMap/'
load(paste0(dataDir, 'pvalue.map.I.II.IV.Rdata'))
head(mapPval); n <- nrow(mapPval)

# Configurations : at least comp with significant exp and meth
Q <- 6; expQ <- Q/2
mapConfig <- GetHinfo(Q, Q)$Hconfig

# Fitting marginals
mapPvalAlone <- as.matrix(mapPval[, -(1:5)])
mapPvalAlone[which(mapPvalAlone > 1-1e-14)] <- 1 - 1e-14 
# par(mfrow=c(Q/2, 4))
# sapply(1:Q, function(q){
#    hist(-qnorm(mapPvalAlone[, q]), breaks=sqrt(n), main=colnames(mapPvalAlone)[q])
#    try(FastKerFdr(mapPvalAlone[, q], plotting=TRUE))
#    })

# Fitting marginals and proportions
par(mfrow=c(Q/2, 2))
if(!file.exists(paste0(dataDir, 'pvalue.map.I.II.IV-marginals.Rdata'))){
   marginalFit <- MixtModProcedure(mapPvalAlone, mapConfig)
   save(marginalFit, file=paste0(dataDir, 'pvalue.map.I.II.IV-marginals.Rdata'))
}
load(paste0(dataDir, 'pvalue.map.I.II.IV-marginals.Rdata'))
postConfig <- apply(marginalFit$posterior, 1, which.max)
par(mfrow=c(1, 1)); image(1:2^Q, 1:2^Q, t(marginalFit$posterior)%*%marginalFit$posterior)
image(1:Q, 1:Q, t(marginalFit$tauKer)%*%marginalFit$tauKer)

# At least one significant pair exp-meth
composedName <- 'atLeastOnePair'
atLeastOnePairH1 <- which(unlist(lapply(mapConfig, function(config){
   sum(config[1:expQ] * config[expQ+(1:expQ)]) >=1
})))
atLeastOnePairH1
t(sapply(atLeastOnePairH1, function(c){mapConfig[[c]]}))

atLeastOnePairRes <- PerformMultipleTestingEM(posterior=marginalFit$posterior, Hconfig.H1=atLeastOnePairH1, Alpha=.05)
H1list <- which(atLeastOnePairRes$Rejection==1)
atLeastOnePairSel <- mapPval[H1list, ]
atLeastOnePairSel$lFdr <- atLeastOnePairRes$lFDR[H1list]
configH1 <- t(sapply(H1list, function(i){mapConfig[[postConfig[i]]]}))
atLeastOnePairSelconfig <- cbind(atLeastOnePairSel, postConfig[H1list], configH1)
names(atLeastOnePairSelconfig)[13] <- 'bestConfig'
save(atLeastOnePairSelconfig, file=paste0(dataDir, 'pvalue.map.I.II.IV-', composedName, '.Rdata'))
table(atLeastOnePairSel$gene)
table(atLeastOnePairSel$cpg)

# Exactly one significant pair exp-meth
composedName <- 'exactlyOnePair'
exactlyOnePairH1 <- which(unlist(lapply(mapConfig, function(config){
   (prod(config == c(1, 0, 0, 1, 0, 0))==1) | 
       (prod(config == c(0, 1, 0, 0, 1, 0))==1) | (prod(config == c(0, 0, 1, 0, 0, 1))==1)
})))
exactlyOnePairH1
t(sapply(exactlyOnePairH1, function(c){mapConfig[[c]]}))

exactlyOnePairRes <- PerformMultipleTestingEM(posterior=marginalFit$posterior, Hconfig.H1=exactlyOnePairH1, Alpha=.05)
H1list <- which(exactlyOnePairRes$Rejection==1)
exactlyOnePairSel <- mapPval[H1list, ]
exactlyOnePairSel$lFdr <- exactlyOnePairRes$lFDR[H1list]
configH1 <- t(sapply(H1list, function(i){mapConfig[[postConfig[i]]]}))
exactlyOnePairSelconfig <- cbind(exactlyOnePairSel, postConfig[H1list], configH1)
names(exactlyOnePairSelconfig)[13] <- 'bestConfig'
save(exactlyOnePairSelconfig, file=paste0(dataDir, 'pvalue.map.I.II.IV-', composedName, '.Rdata'))
table(exactlyOnePairSel$gene)
table(exactlyOnePairSel$cpg)
