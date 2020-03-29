## Alternative try for the estimation of FDR for combined tests

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
library(data.table); library(ks); library(mclust); library(spatstat); library(ROCR); 
library(gtools); library(mvtnorm)
library(tidyverse); library(data.table)
source('./Functions.R')
setwd('/home/robin/RECHERCHE/EXPRESSION/IFCAM/IFCAM-genomics/Programs/Integration/2020/')
simDir <- '../../../Results/Pmax/'; 

## Simul results
simName <- 'PmaxStudy_NbObs10000_Q8_Sim1'; 
pValMat <- as.matrix(fread(paste0(simDir, simName, '_PvalMat.txt'), header=TRUE))
config <- as.matrix(fread(paste0(simDir, simName, '_Truth.txt'), header=TRUE))
Q <- ncol(pValMat); n <- nrow(pValMat)

## Configurations : test 'pMax' = H1 for all tests
atLeast = Q; Tmp <- GetHinfo(Q, atLeast); Hconfig <- Tmp$Hconfig; Hconfig.H1 <- Tmp$Hconfig.H1

## New simulations
n <- 1e4; Q <- 8; seed <- 1; set.seed(seed); simName <- paste0('NewSimSR_NbObs', n, '_Q', Q, '_Sim', seed);
DirichletProp <- c(8, 2); MinPropForHconfig.H1 <- 0.03; MinNbOfUnitPerHconfig.H1 <- 5; Rho0 <- 0; Rho1 <- 0;
Delta <- sqrt(2*(1:Q))
simul <- GeneratePvalues(n, Q, DirichletProp, Hconfig, Hconfig.H1, MinPropForHconfig.H1, Delta,
                         Rho0, Rho1, MinNbOfUnitPerHconfig.H1)
pValMat <- simul$P; config <- simul$Truth

## True H1 & pMax
H1 <- 1*(config==Hconfig.H1); pMax <- pValH1 <- GetPH1(pValMat, atLeast)

## Fit kerFdr for each test
if(!file.exists(paste0(simDir, simName,'_TauMat.txt'))){
   par(mfrow=c(3, 3), mex=.6)
   tauMat <- matrix(0, n, Q); p0 <- rep(0, Q)
   for(q in 1:Q){
      ker <- FastKerFdr(pValMat[, q], plotting=TRUE)
      tauMat[, q] <- ker$tau; p0[q] <- ker$p0
   }
   fwrite(x=data.frame(tauMat=tauMat), file=paste0(simDir, simName,'_TauMat.txt'))
   fwrite(x=data.frame(p0=p0), file=paste0(simDir, simName,'_p0KerFdR.txt'))
}
tauMat <- as.matrix(fread(paste0(simDir, simName, '_TauMat.txt'), header=TRUE))
p0fit <- as.matrix(fread(paste0(simDir, simName, '_p0KerFdR.txt'), header=TRUE))
tauProd <- apply(tauMat, 1, prod)
## -> Doubts about the sim parms: increasingly easy with q

## Get estimates of the Hconfig posteriors:
## Post.Hconfig[i] = Prod_{q in H0} (1-Post.MixtMarg[i,q]) * Prod_{q in H1} Post.MixtMarg[i,q]
logTauMat <- log(tauMat); log1_TauMat <- log(1-tauMat)
tauHconfig <- sapply(Hconfig, function(h){
   tau <- rep(0,nrow(tauMat))
   if (length(which(h==1)) > 0){tau <- tau + rowSums(logTauMat[, which(h==1), drop=FALSE])}
   if (length(which(h==0)) > 0){tau <- tau + rowSums(log1_TauMat[, which(h==0), drop=FALSE])}
   return(exp(tau))
})
priorHconfigOld <- colMeans(tauHconfig)
tauH1 <- tauHconfig[, Hconfig.H1]

## Combined p-values
pValUnionH1 <- ComputePValue(pValMat, logTauMat, Hconfig, Hconfig.H1, priorHconfigOld, pValH1, Method='H1')
pValUnionH0 <- ComputePValue(pValMat, logTauMat, Hconfig, Hconfig.H1, priorHconfigOld, pValH1, Method='H0')
par(mfrow=c(2, 2))
hist(pValH1, breaks=sqrt(n)); hist(log(tauH1), breaks=sqrt(n))
hist(pValUnionH1, breaks=sqrt(n)); hist(pValUnionH0, breaks=sqrt(n))

## Compare criteria 
crit <- as.data.frame(cbind(tauProd, tauH1, pMax, pValUnionH0, pValUnionH1))
colnames(crit) <- c('tauProd', 'tauH1', 'pMax', 'pValUnionH0', 'pValUnionH1')
# plot(crit, log='xy')
## -> tauProd = tauH1, same ordering for pMax, pValUnionH1, pValUnionH0

## ROC curves
par(mfrow=c(1, 1))
plot(performance(prediction(-pMax, H1), 'tpr', 'fpr'), xlim=c(0, .3), ylim=c(.5, 1), lwd=2)
plot(performance(prediction(-pValUnionH1, H1), 'tpr', 'fpr'), col=2, add=TRUE, lwd=2)
plot(performance(prediction(-pValUnionH0, H1), 'tpr', 'fpr'), col=4, add=TRUE, lwd=2)
plot(performance(prediction(tauProd, H1), 'tpr', 'fpr'), col=7, add=TRUE, lwd=2, lty=2)
plot(performance(prediction(tauH1, H1), 'tpr', 'fpr'), col=6, add=TRUE, lwd=2, lty=2)
## -> much better discrimination with tau's than with combined p-values

###############################################################################
## Changing the estimation of the priors proportions
###############################################################################
## True proportions for each configuration : CHEATING !!!!
priorHconfigTrue <- sapply(1:length(Hconfig), function(c){sum(config==c)})/n
###############################################################################
## Alternative estimates of the prior
priorHconfigNew <- sapply(1:length(Hconfig), function(c){
   prod(p0fit[which(Hconfig[[c]]==0)]) * prod(1-p0fit[which(Hconfig[[c]]==1)])
})
plot(1e-4+priorHconfigTrue, 1e-4+priorHconfigOld, log='xy'); abline(0, 1)
points(1e-4+priorHconfigTrue, 1e-4+priorHconfigNew, col=2)

###############################################################################
## Same analysis as before choosing the estimates
## Combining p-values
priorTest <- priorHconfigNew # priorHconfigTrue, priorHconfigOld, priorHconfigNew
pValUnionH1Test <- ComputePValue(pValMat, logTauMat, Hconfig, Hconfig.H1, priorTest, pValH1, Method='H1')
pValUnionH0Test <- ComputePValue(pValMat, logTauMat, Hconfig, Hconfig.H1, priorTest, pValH1, Method='H0')
# par(mfrow=c(2, 2))
# hist(pValH1, breaks=sqrt(n)); hist(log(tauH1), breaks=sqrt(n))
# hist(pValUnionH1Test, breaks=sqrt(n)); hist(pValUnionH0Test, breaks=sqrt(n))

## FDR estimation
# pdf(paste0(simDir, simName, 'FDRpriorTest.pdf'))
par(mfrow=c(1, 1))
orderPMax <- order(pMax); rankPMax <- rank(pMax)
trueFDR <- (cumsum(1-H1[orderPMax])/(1:n))[rankPMax]; plot(pMax, trueFDR, ylim=c(0, 1))
tauH1FDR <- (cumsum(1-tauH1[orderPMax])/(1:n))[rankPMax]; points(pMax, tauH1FDR, col=6)
pi0 <- sum(pValUnionH0Test>.5)/(n/2)
pValUnionH1TestFDR <- pi0*p.adjust(pValUnionH1Test, method='BH'); points(pMax, pValUnionH1TestFDR, col=2)
pValUnionH0TestFDR <- pi0*p.adjust(pValUnionH0Test, method='BH'); points(pMax, pValUnionH0TestFDR, col=4)
# dev.off()
