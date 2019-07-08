rm(list=ls())
library(KernSmooth)
library(mclust)
#library(Matrix)
# source('D:/IFCAM/IFCAM-genomics/Programs/Integration/Functions/F_KerFdr.R')
source('Functions/F_KerFdr.R')


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
  select(Gene,Tau.E,Qstat,Pval.E,Pval.M) %>% 
  unnest %>% 
  mutate(Tau.M = GetH1.M$tau,
         Tau.H0.joint = 1-(Tau.E*Tau.M))
Stephane <- DataForSampling %>% as.data.frame


        #### Perform sampling and plot the results


NbDraws = 1e8
H0Sample <- sample(1:nrow(DataForSampling), replace = T, prob = DataForSampling$Tau.H0.joint,size = NbDraws)

## Plot the actual and H0 distribution
H0.Qstat <- DataForSampling$Qstat[H0Sample]
Dens <- density(H0.Qstat)
Area <- mean(DataForSampling$Tau.H0.joint) 
Dens$y <- Dens$y*Area

par(mfrow=c(1,1))
hist(DataForSampling$Qstat,freq = FALSE,100)
lines(Dens,col=2)

## In case one does want to add the N(0,1)
Normal <- curve(dnorm,from=-5,to=10,new=FALSE)
lines(Normal$x, Normal$y*Area,col=4)


      #### Write the loops to perfom parallel computation according to Test/EmpiricalPvalues.R
