# KerFdr

rm(list=ls())
library(KernSmooth); library(mclust); library(Matrix)
source('Functions/F_KerFdr.R')
DataDir = ("/home/robin/Bureau/RECHERCHE/EXPRESSION/IFCAM/Data-NotUpload/")

# load(paste0(DataDir, 'TTests-LogExp.Rdata')) 
# KFE = F_KerFdr(PE, plotting=F)
# save(PE, KFE, file=paste0(DataDir, 'TTests-LogExp-KerFdr.Rdata')) 
load(paste0(DataDir, 'TTests-LogExp-KerFdr.Rdata'))

# load(paste0(DataDir, 'BetaRegTests-Methylation.Rdata'))
# PM = PM[-which(PM==0)]; KFM = F_KerFdr(PM, plotting=F)
# save(PM, KFM, file=paste0(DataDir, 'BetaRegTests-Methylation-KerFdr.Rdata')) 
load(paste0(DataDir, 'BetaRegTests-Methylation-KerFdr.Rdata'))

# Observed stat products
LogExp = readRDS(paste0(DataDir, "LogExp.rds"))
MethInfo = readRDS(paste0(DataDir, "MethInfo.rds"))
XE = KFE$X; XM = KFM$X
CpGlist = ProdStat = ProdProb = list()
sapply(1:length(rownames(LogExp)), function(i){
   if(i%%100==0){cat(i, '')}
   CpGlist[[i]] <<- which(MethInfo$Gene==rownames(LogExp))
   CpGlist[[i]] <<- CpGlist[[i]][-which(CpGlist[[i]] > length(PM))]
   ProdStat[[i]] <<- XE[i]%*%XM[CpGlist[[i]]]
   ProdProb[[i]] <<- KFE$tau[i]%*%KFM$tau[CpGlist[[i]]]
})

# Null for the joint test stat
B = 1e4
XEsample = sample(1:length(PE), size=B, prob=1-KFE$tau, replace=T)
XMsample = sample(1:length(PM), size=B, prob=1-KFM$tau, replace=T)
hist(XEsample*XMsample)

