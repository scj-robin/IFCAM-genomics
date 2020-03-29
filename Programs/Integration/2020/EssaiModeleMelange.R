rm(list=ls());
library(data.table); 
library(ks); 
library(mclust); 
library(spatstat); 
library(ROCR)
library(tidyverse); 
library(data.table)
setwd('D:/IFCAM/IFCAM-genomics/Programs/Integration/2020/')
source('./Functions.R')


          #### Load the data


simDir <- '../../../Results/Pmax/';
simName <- 'PmaxStudy_NbObs10000_Q8_Sim1';
pValMat <- as.matrix(fread(paste0(simDir, simName, '_PvalMat.txt'), header=TRUE))
config <- as.matrix(fread(paste0(simDir, simName, '_Truth.txt'), header=TRUE))
Q <- ncol(pValMat); 
n <- nrow(pValMat)


        #### Configurations 


atLeast = Q; 
Tmp <- GetHinfo(Q, atLeast); 
Hconfig <- Tmp$Hconfig; 
Hconfig.H1 <- Tmp$Hconfig.H1
NbConfig <- length(Hconfig)


## True H1 & pMax
H1 <- 1*(config==Hconfig.H1); 
pMax <- pValH1 <- GetPH1(pValMat, atLeast)


        #### Marginal density estimation


## Fit kerFdr for each test
## to get the marginal densities
par(mfrow=c(3, 3), mex=.6)
tauMat <- f1Mat <- matrix(0, n, Q); p0 <- rep(0, Q)
for(q in 1:Q){
  ker <- FastKerFdr(pValMat[, q], plotting=TRUE)
  tauMat[, q] <- ker$tau; 
  f1Mat[,q] <- ker$f1
  p0[q] <- ker$p0
}
f0Mat <- matrix(dnorm(-qnorm(pValMat)),ncol=Q)

##Just to check
plot(-qnorm(pValMat[,1]),f0Mat[,1]*p0[1]+ f1Mat[,1]*(1-p0[1]),cex=0.7)
points(-qnorm(pValMat[,1]),f0Mat[,1]*p0[1],col=4,cex=0.7)
points(-qnorm(pValMat[,1]),f1Mat[,1]*(1-p0[1]),col=2,cex=0.7)
par(mfrow=c(1,1))

        #### From marginal to config densities


## Go from marginal to config densities
Logf0Mat <- log(f0Mat); 
Logf1Mat <- log(f1Mat); 
f.Hconfig <- sapply(Hconfig, function(h){
  f <- rep(0,nrow(Logf0Mat))
  if (length(which(h==1)) > 0){f <- f + rowSums(Logf1Mat[, which(h==1), drop=FALSE])}
  if (length(which(h==0)) > 0){f <- f + rowSums(Logf0Mat[, which(h==0), drop=FALSE])}
  return(exp(f))
})
dim(f.Hconfig)


        #### First step: cheat a beat


## Build the posteriors
priorHconfigTrue <- sapply(1:length(Hconfig), function(c){sum(config==c)})/n
posteriors <- f.Hconfig*(tcrossprod(rep(1:n),priorHconfigTrue))
posteriors <- posteriors/rowSums(posteriors)

## Perform multiple testing
localFDR <- 1-posteriors[,length(Hconfig)]
Order <- order(localFDR)
FDR <- cumsum(localFDR[Order])/(1:n)
plot(FDR)
NbReject <- max(which(FDR<=0.05))
Rejection <- rep(0,n)
Rejection[Order[1:NbReject]] <- 1
table(config==length(Hconfig),Rejection)

##So far so good !


        #### Second step: naive version with marginal priors


## Simple product of marginal priors estimator
priorHconfigFit <- sapply(1:length(Hconfig), function(c){
  prod(p0[which(Hconfig[[c]]==0)]) * prod(1-p0[which(Hconfig[[c]]==1)])
})

## How good is it ?
plot(priorHconfigTrue,priorHconfigFit)
points(priorHconfigTrue[NbConfig],priorHconfigFit[NbConfig],col=2,pch=16)
abline(0,1,col=4,lwd=2)
#Not that good

posteriors <- f.Hconfig*(tcrossprod(rep(1:n),priorHconfigFit))
posteriors <- posteriors/rowSums(posteriors)

## Perform multiple testing
localFDR <- 1-posteriors[,length(Hconfig)]
Order <- order(localFDR)
FDR <- cumsum(localFDR[Order])/(1:n)
plot(FDR)
NbReject <- max(which(FDR<=0.05))
NbReject
#We lose everything...


        #### Third step: EM for prior estimation


## Initialization
NotOK <- TRUE
Precision <- 1e-6
NoLowerThan <- 1e-7
NewPrior <- priorHconfigFit
while(NotOK){
  
  ## E step
  Tau <- f.Hconfig*(tcrossprod(rep(1:n),NewPrior))
  Tau <- Tau/rowSums(Tau)
  
  ## M step 
  OldPrior <- NewPrior
  NewPrior <- colMeans(Tau)
  NewPrior[NewPrior<NoLowerThan] <- NoLowerThan
  NewPrior <- NewPrior/sum(NewPrior)
  NotOK <- max((OldPrior-NewPrior)^2) > Precision
  
}

## How good is it ?
priorHconfigEM <- NewPrior
plot(priorHconfigTrue,priorHconfigEM)
points(priorHconfigTrue[NbConfig],priorHconfigEM[NbConfig],col=2,pch=16)
abline(0,1,col=4,lwd=2)

posteriors <- f.Hconfig*(tcrossprod(rep(1:n),priorHconfigEM))
posteriors <- posteriors/rowSums(posteriors)

## Perform multiple testing
localFDR <- 1-posteriors[,length(Hconfig)]
Order <- order(localFDR)
FDR <- cumsum(localFDR[Order])/(1:n)
plot(FDR)
lines(cumsum(as.numeric(config[Order]!=NbConfig))/(1:n),col=2,lwd=2)
NbReject <- max(which(FDR<=0.05))
NbReject
Rejection <- rep(0,n)
Rejection[Order[1:NbReject]] <- 1
table(config==length(Hconfig),Rejection)
#We get some guys.





