## The goal here is to perform the joint anlysis of the Methylation and 
## Expression data for each gene the methylation. for each couple 
## (Gene exp)x(Methyl probe) the 2 pvalues are combined into a single test statistic.
## The true H0 distribution of this statistic is unknown, and has to be infered 
## using a kernel strategy (the one of kerfdr, see ref). Once this H0 distribution is 
## derived, one can compute pvalues associated to the test statistic by resampling.

rm(list=ls())
library(tidyverse)
library(KernSmooth)
library(mclust)
#library(Matrix)
# source('D:/IFCAM/IFCAM-genomics/Programs/Integration/Functions/F_KerFdr.R')
source('Functions/F_KerFdr.R')
source('Functions/F_EmpPvalues.R')


        #### Parameters


# DataRep <- 'D:/IFCAM/IFCAM-genomics/Data/PrivateIntegration/'
DataRep <- '../../Data/PrivateIntegration/'
TypeOfTest <-  'Welch'  #  'Wilcoxon' " 'Student' # 


        #### Read the data


Tests <- readRDS(paste0(DataRep,'ExpMeth_WelchAndBeta.rds'))
dim(Tests)
Tests

## Get the pvalues
Tests <- Tests %>%
  mutate(Pval.E = map_dbl(Test, ~.x$PExp),
         Pval.M = purrr::map(Test, ~.x$PMeth),
         NbProbePerGene = map_int(Meth,nrow),
         Qstat.E = map_dbl(Pval.E, ~ -qnorm(.x)),
         Qstat.M = purrr::map(Pval.M, ~ -qnorm(.x)),
         MinStat = map2(Qstat.E,Qstat.M, ~ sapply(.y, function(uu) min(.x,uu))),
         ProdStat = map2(Qstat.E,Qstat.M, ~ .y*.x)
         ) 


        #### Get the H1 posterior distributions


## Get the posteriors
GetH1.E <- Tests$Pval.E %>% F_KerFdr
GetH1.M <- Tests %>%
  select(Pval.M) %>% 
  unnest %>% 
  pull %>% 
  F_KerFdr

## Get the posterior products
DataForSampling <- Tests %>% 
  mutate(Tau.E = GetH1.E$tau) %>% 
  select(Gene,Tau.E,MinStat,Pval.E,Pval.M) %>% 
  unnest %>% 
  mutate(Tau.M = GetH1.M$tau,
         Tau.H0.joint = 1-(Tau.E*Tau.M))
Stephane <- DataForSampling %>% as.data.frame


        #### Perform sampling and plot the results


NbDraws = 1e8
H0Sample <- sample(1:nrow(DataForSampling), replace = T, prob = DataForSampling$Tau.H0.joint,size = NbDraws)

## Plot the actual and H0 distribution
H0.MinStat <- DataForSampling$MinStat[H0Sample]
Dens <- density(H0.MinStat)
Area <- mean(DataForSampling$Tau.H0.joint) 
Dens$y <- Dens$y*Area

par(mfrow=c(1,1))
hist(DataForSampling$MinStat,freq = FALSE,100)
lines(Dens, col=2, lwd=2)

## In case one does want to add the N(0,1)
Absc <- seq(-5,10,length.out = 1000)
Normal <- dnorm(Absc)*Area
lines(Absc, Normal,col=4,new=FALSE, lwd=2)


      #### Computing the p-values
# Use the function Functions/F_EmpPvalues.R only for B <= 1e8
# For B > 1e8, use the follogin loop

NbDraws = 1e6; BlockSize = 1e5; # BlockSize <= 1e7
NbBlock = round(NbDraws/BlockSize)
ObservedStat = Stephane$MinStat; NbStat = length(ObservedStat)
rankObservedStat = length(ObservedStat) - rank(ObservedStat) + 1
sortObservedStat = sort(ObservedStat, decreasing=T)
Pval = rep(0, NbStat)
sapply(1:NbBlock, function(b){
   cat(b, ': ')
   S0 = sample(1:nrow(DataForSampling), replace=T, 
                prob = DataForSampling$Tau.H0.joint, size=BlockSize)
   X0 = DataForSampling$MinStat[S0]
   Pval <<- Pval + F_EmpPvalues(sortObservedStat, X0)
   # H = hist(sortObservedStat, breaks=sqrt(NbStat))
   # lines(sortObservedStat, max(H$counts)/2*Pval/b, type='l', lwd=2, col=2)
   plot(sortObservedStat, Pval/b, type='l', lwd=2, col=2, log='y')
})
Pval = Pval[rankObservedStat] / NbBlock

H = hist(ObservedStat, breaks=sqrt(NbStat))
points(ObservedStat, max(H$counts)/2*Pval/b, pch=20, lwd=2, col=2)
