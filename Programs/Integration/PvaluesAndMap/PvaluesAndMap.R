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
mapConfigH1 <- which(unlist(lapply(mapConfig, function(config){
   sum(config[1:expQ] * config[expQ+(1:expQ)]) >=1
   })))
mapConfigH1
t(sapply(mapConfigH1, function(c){mapConfig[[c]]}))

# Fitting marginals
mapPvalAlone <- as.matrix(mapPval[, -(1:5)])
mapPvalAlone[which(mapPvalAlone > 1-1e-14)] <- 1 - 1e-14 
# par(mfrow=c(Q/2, 4))
# sapply(1:Q, function(q){
#    hist(-qnorm(mapPvalAlone[, q]), breaks=sqrt(n), main=colnames(mapPvalAlone)[q])
#    try(FastKerFdr(mapPvalAlone[, q], plotting=TRUE))
#    })

# Whole procedure
ResMMP <- MixtModProcedure(mapPvalAlone, mapConfig)
names(ResMMP); dim(ResMMP$posterior); dim(mapPval)
finalRes <- PerformMultipleTestingEM(posterior=ResMMP$posterior, Hconfig.H1=mapConfigH1, Alpha=.05)
length(finalRes)
H1list <- which(finalRes==1)
mapH1 <- mapPval[which(finalRes==1), 1:3]
table(mapH1$gene)
table(mapH1$cpg)

           