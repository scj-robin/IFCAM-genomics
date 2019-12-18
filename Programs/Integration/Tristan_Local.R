## ---------------------------
## Script name: 
## Purpose of script:
## Date of creation: 2019-12-10
## ---------------------------

rm(list=ls())
setwd("D:/IFCAM/IFCAM-genomics/Programs/Integration/")      
## For data generation
library(mvtnorm) 
library(gtools) 
## For analysis
library(ks)
library(mclust)
library(spatstat)
source('./Functions.R')


        #### Parameters


## For data generation
n = 1e5
Q = 6
DirichletProp = c(2,8)
MinPropForHconfig.H1 = 0.01
MinNbOfUnitPerHconfig.H1 = 5
Delta = sqrt(3*(1:Q))
Rho0= 0
Rho1 = 0

## For analysis
AtLeast = Q-2
Alpha = .1



        #### Data generation



## Get the list of H configuration, and the ones matching H1
Tmp <- GetHinfo(Q,AtLeast)
Hconfig <- Tmp$Hconfig
Hconfig.H1 <- Tmp$Hconfig.H1
rm(Tmp)

## Generate the pvalues
Tmp <- GeneratePvalues(n,Q,DirichletProp,Hconfig,Hconfig.H1,
                       MinPropForHconfig.H1,Delta,
                       Rho0,Rho1,MinNbOfUnitPerHconfig.H1)
Truth <- Tmp$Truth
PvalMat <- Tmp$P

## Compute the pvalue for the desired H1
P.H1 <- GetPH1(PvalMat,AtLeast)

## Have a look at the data
table(Truth)
boxplot(P.H1 ~ Truth)
hist(P.H1, breaks=sqrt(n))



        #### Naive fdr procedure applied to pmax



P.H1.BH = p.adjust(P.H1, method='fdr')



        #### Second naive procedure : apply fdr per series and cross lists



FDRnaive <- sapply(1:Q, function(q) p.adjust(PvalMat[,q],method = 'fdr'))
ClassifFDRnaive <- apply(FDRnaive,1, function(p) sum(p<Alpha)>=AtLeast)


        
        #### Compute pvalue for the global null hypothesis



## Compute the pvalues
PvalUnion <- ComputePValue(PvalMat,Hconfig,Hconfig.H1,P.H1)  
if(AtLeast == Q){
  PvalUnion_Qmax <- ComputePValue_Qmax(PvalMat,Hconfig,Hconfig.H1,P.H1)  
  plot(log(PvalUnion),log(PvalUnion_Qmax))
  abline(0,1,col=2)
}
if(Q==2){
  PvalUnion_2 <- ComputePValue_2(PvalMat,Hconfig,Hconfig.H1,P.H1)  
  plot(PvalUnion,PvalUnion_2)
}
boxplot(PvalUnion ~ Truth)   
hist(PvalUnion,sqrt(length(PvalUnion)))



        #### Multiple testing correction


N0=2*sum(PvalUnion>0.5)
Th=1e-5
FDR = Th*N0/sum(PvalUnion<Th);FDR
PvalUnion.fdr <- p.adjust(PvalUnion, method='fdr') 
table(Truth%in%Hconfig.H1, P.H1.BH<Alpha)
table(Truth%in%Hconfig.H1,ClassifFDRnaive)
table(Truth%in%Hconfig.H1, PvalUnion.fdr<Alpha)
if(AtLeast == Q){
  PvalUnion_Qmax.fdr <- p.adjust(PvalUnion_Qmax, method='fdr') 
  table(Truth%in%Hconfig.H1, PvalUnion_Qmax.fdr<Alpha)
}
if(Q==2){
  PvalUnion_2.fdr <- p.adjust(PvalUnion_2, method='fdr') 
  print(table(Truth%in%Hconfig.H1, PvalUnion_2.fdr<Alpha))
  plot(log10(PvalUnion_2),log10(PvalUnion_Qmax))
  abline(a=0,b=1,col=2)
}
max(abs(PvalUnion_2-PvalUnion_Qmax))



# ## nQ*Q^{2*(Q-K+1)}
# library(tidyverse)
# Qval = 6:10
# Kval = 1:10
# for(K in 1:length(Kval)){
#         for(Q in 1:length(Qval)){
#                 if(Qval[Q]>=Kval[K]){
#                         Cost[K,Q] = min((Qval[Q])**(2*(Kval[K])),(Qval[Q]**(2*(Qval[Q]-Kval[K]))))        
#                 } else {
#                         Cost[K,Q] = 0
#                 }
#                 
#         }
# }
# 
# 
# colnames(Cost)=6:10
# row.names(Cost)=1:10
# Cost
