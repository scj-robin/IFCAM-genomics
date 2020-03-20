## ---------------------------
## Script name: 
## Purpose of script:
## Date of creation: 2019-12-10
## ---------------------------

rm(list=ls())
setwd("D:/IFCAM/IFCAM-genomics/Programs/Integration/2020")      

## For data generation
library(mvtnorm) 
library(gtools) 

## For analysis
library(ks)
library(mclust)
library(spatstat)
source('./Functions.R')

## For hype
library(tidyverse)
library(data.table)


        #### Parameters


## For data generation
NbObsList = c(1e4,1e5,1e6)
QList = c(2,4,8)
DirichletProp = c(8,2)
MinPropForHconfig.H1 = 0.03
MinNbOfUnitPerHconfig.H1 = 5
Rho0= 0
Rho1 = 0

## For analysis
AlphaList = c(0.01,0.05,0.1)

## Others
NbSim=20
ResRep <- 'D:/IFCAM/IFCAM-genomics/Results/Pmax/'



        #### Build the grid



ConfigList <- expand.grid(NbObsList,QList,1:NbSim)
NbConfig <- nrow(ConfigList)



        #### Run the loop that produces all required quantities


# ## Stoppé à 63
# sapply(1:NbConfig, function(ii){
# 
#   print(ii)
# 
# 
#           #### Local parameters
# 
# 
#   n = ConfigList[ii,1]
#   Q = ConfigList[ii,2]
#   Sim = ConfigList[ii,3]
#   Delta = sqrt(2*(1:Q))
#   AtLeast = Q
#   Name = paste0('PmaxStudy_NbObs',n,'_Q',Q,'_Sim',Sim)
# 
# 
#           #### Data generation
# 
# 
#   ## Get the list of H configuration, and the ones matching H1
#   Tmp <- GetHinfo(Q,AtLeast)
#   Hconfig <- Tmp$Hconfig
#   Hconfig.H1 <- Tmp$Hconfig.H1
#   rm(Tmp)
# 
#   ## Generate the pvalues
#   Tmp <- GeneratePvalues(n,Q,DirichletProp,Hconfig,Hconfig.H1,
#                          MinPropForHconfig.H1,Delta,
#                          Rho0,Rho1,MinNbOfUnitPerHconfig.H1)
#   Truth <- Tmp$Truth
#   PvalMat <- Tmp$P
# 
#   ## Save all pvalues
#   fwrite(x = data.frame(Truth=Truth),file = paste0(ResRep,Name,'_Truth.txt'))
#   fwrite(x = as.data.frame(PvalMat),file = paste0(ResRep,Name,'_PvalMat.txt'))
# 
#   ## Compute the pvalue for the desired H1
#   P.H1 <- GetPH1(PvalMat,AtLeast)
# 
# 
#           #### Compute marginal and configuration posteriors
# 
# 
#   ## Infer mixture marginally on each of the Q test series,
#   ## and get the posterior probability to be H1.
#   Tau.MixtMarg <- sapply(1:Q, function(q){
#     Res <- FastKerFdr(Pval=PvalMat[, q])
#     print(Res$p0)
#     return(Res$tau)
#   })
#   LogTau.MixtMarg <- log(Tau.MixtMarg)
#   LogOneMinusTau.MixtMarg <- log(1-Tau.MixtMarg)
# 
#   ## Get estimates of the Hconfig posteriors:
#   ## Post.Hconfig[i] = Prod_{q in H0} (1-Post.MixtMarg[i,q]) * Prod_{q in H1} Post.MixtMarg[i,q]
#   Tau.Hconfig <- sapply(Hconfig, function(h){
#     Tau <- rep(0,nrow(Tau.MixtMarg))
#     if (length(which(h==1))>0){
#       Tau <- Tau+rowSums(LogTau.MixtMarg[,which(h==1),drop=FALSE])
#     }
#     if (length(which(h==0))>0){
#       Tau <- Tau+rowSums(LogOneMinusTau.MixtMarg[,which(h==0),drop=FALSE])
#     }
#     return(exp(Tau))
#   })
# 
#   ## Infer priors of each Hconfig
#   Prior.Hconfig <- colMeans(Tau.Hconfig)
# 
# 
#           #### Compute pvalue for the global null hypothesis
# 
# 
#   ## Compute the pvalues
#   PvalUnion.H1 <- ComputePValue(PvalMat,LogTau.MixtMarg,Hconfig,Hconfig.H1,Prior.Hconfig,P.H1,Method='H1')
#   PvalUnion.H0 <- ComputePValue(PvalMat,LogTau.MixtMarg,Hconfig,Hconfig.H1,Prior.Hconfig,P.H1,Method='H0')
# 
#   ## Save  pvalues
#   fwrite(x = as.data.frame(PvalUnion.H1),file = paste0(ResRep,Name,'_PvalUnionH1.txt'))
#   fwrite(x = as.data.frame(PvalUnion.H0),file = paste0(ResRep,Name,'_PvalUnionH0.txt'))
# 
# })



        #### Run the loop for the analyzes


alpha = 0.05
Results <- lapply(1:NbConfig, function(ii){

  ## Don't get bored
  print(ii)
  
  ## Get local parameters        
  n = ConfigList[ii,1]
  Q = ConfigList[ii,2]
  Sim = ConfigList[ii,3]
  Delta = sqrt(2*(1:Q))
  AtLeast = Q
  
  ## Load the data
  Name = paste0('PmaxStudy_NbObs',n,'_Q',Q,'_Sim',Sim)
  Truth <- fread(paste0(ResRep,Name,'_Truth.txt')) %>% as.matrix
  PvalMat <- fread(paste0(ResRep,Name,'_PvalMat.txt')) %>% as.matrix
  PvalUnion.H0 <- fread(file = paste0(ResRep,Name,'_PvalUnionH0.txt')) %>% as.matrix
  PvalUnion.H1 <- fread(file = paste0(ResRep,Name,'_PvalUnionH1.txt')) %>% as.matrix

  ## Compute the H1 pvalue (eg pmax if AtLeast = Q)
  P.H1 <- GetPH1(PvalMat,AtLeast)

  ## Get the list of H configuration, and the ones matching H1
  Tmp <- GetHinfo(Q,AtLeast)
  Hconfig <- Tmp$Hconfig
  Hconfig.H1 <- Tmp$Hconfig.H1
  rm(Tmp)


          #### Multiple testing correction procedures

  
  MTMethods <- list() 
  TrueH1 <- as.numeric(Truth %in% Hconfig.H1)

  ## Version 1: naive Bonferroni on pmax
  MTMethods$P.H1.Bonf = p.adjust(P.H1, method='bonferroni')

  ## Version 2: naive Bonferroni on pmax
  MTMethods$P.H1.BH = p.adjust(P.H1, method='BH')
  
  ## Version 3: union-intersection, via H1
  MTMethods$PvalUnion.H1.fdr <- p.adjust(PvalUnion.H1, method='BH')

  ## Version 4: union-intersection, via H0
  MTMethods$PvalUnion.H0.fdr <- p.adjust(PvalUnion.H0, method='BH')

  ## Version 5: naive FDR per series + cross of lists
  IntersectFDR <- sapply(1:Q, function(q) p.adjust(PvalMat[,q],method = 'fdr'))
  
  
          #### Evaluate their performance
  
  EvalPerf <- function(Ref,Meth,alpha){
    head(Ref);head(Meth)
    NbReject <- length(which(Meth<=alpha))
    NbFP <- length(which((Meth<=alpha) & !Ref))
    NbTrueH1 <- sum(Ref)
    return(list(NbReject=NbReject,NbFP=NbFP,NbTrueH1=NbTrueH1))
  }
  EvalPerfIntersect <- function(Ref,Meth,alpha){
    #Meth = IntersectFDR
    Rejected <- Meth %>% `<`(alpha) %>% rowSums %>% `==`(2)
    NbReject <- length(which(Rejected))
    NbFP <- length(which(Rejected & !Ref))  
    NbTrueH1 <- sum(Ref)
    return(list(NbReject=NbReject,NbFP=NbFP,NbTrueH1=NbTrueH1))
  }
  Res <- map(MTMethods,~EvalPerf(TrueH1,.x,alpha))
  Res[['IntersectFDR']] <- EvalPerfIntersect(TrueH1,IntersectFDR,alpha)
  
  return(Res)
})
saveRDS(Results,paste0('SummaryPmaxStudy_',alpha,'.rds'))


        #### Synthesis

Division <- function(num,den){
  if (den>0){res=num/den} else {res = 0}
}
FDR <- function(elem){
  Division(elem$NbFP,elem$NbReject)
}
Power <- function(elem){
  (elem$NbReject-elem$NbFP)/elem$NbTrueH1
}


Summary <- readRDS(paste0('SummaryPmaxStudy_',alpha,'.rds')) %>% 
  tibble(Summary=.) %>% 
  mutate(NbObs=ConfigList[,1], Q=ConfigList[,2], Sim=ConfigList[,3]) %>% 
  mutate(Pmax.Bonf = map(Summary,~.x$P.H1.Bonf),
         Pmax.BH = map(Summary,~.x$P.H1.BH),
         PvalUnion.H1.fdr = map(Summary,~.x$PvalUnion.H1.fdr),
         PvalUnion.H0.fdr = map(Summary,~.x$PvalUnion.H1.fdr),
         IntersectFDR = map(Summary,~.x$IntersectFDR)) %>% 
 mutate(Pmax.Bonf_FDR = map_dbl(Pmax.Bonf,FDR),
        Pmax.BH_FDR = map_dbl(Pmax.BH,FDR),
        PvalUnion.H1.fdr_FDR = map_dbl(PvalUnion.H1.fdr,FDR),
        PvalUnion.H0.fdr_FDR = map_dbl(PvalUnion.H0.fdr,FDR),
        IntersectFDR_FDR = map_dbl(IntersectFDR,FDR)) %>% 
  mutate(Pmax.Bonf_Power = map_dbl(Pmax.Bonf,Power),
       Pmax.BH_Power = map_dbl(Pmax.BH,Power),
       PvalUnion.H1.fdr_Power = map_dbl(PvalUnion.H1.fdr,Power),
       PvalUnion.H0.fdr_Power = map_dbl(PvalUnion.H0.fdr,Power),
       IntersectFDR_Power = map_dbl(IntersectFDR,Power))  %>% 
  select(NbObs, Q, Sim,contains('_FDR'),contains('_Power')) %>% 
  select(-starts_with('Pmax.Bonf')) %>% 
  group_by(NbObs,Q) %>%
  summarise_at(vars(starts_with('Pmax'),starts_with('Pval'),starts_with('Intersect')),mean)
Summary  
write.table(Summary,file = 'PremiersResultats.txt', row.names = FALSE)

Results <- read.table(file = 'PremiersResultats.txt', header=TRUE)
round(Results,4) %>% select(-PvalUnion.H0)


