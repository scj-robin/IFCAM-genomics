# Sampling from the joint H0

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
source('../Functions/F_KerFdr.R')
source('../Functions/F_JointPvalueFromObserved.R')

# Parms
nA = 1e3; nB = 1e3; p0A = .8; p0B = .8; mu1A = 1; mu1B = 2; P = 1e4

# Test statistics
HA = rbinom(nA, 1, 1-p0A); statA = rnorm(nA) + mu1A*HA 
HB = rbinom(nB, 1, 1-p0B); statB = rnorm(nB) + mu1B*HB 

# Pairs to be tested
jointTests = as.data.frame(cbind(sample(1:nA, P, replace=TRUE), sample(1:nB, P, replace=TRUE)))
jointTests = jointTests[-which(duplicated(jointTests)), ]; 
names(jointTests) = c('A', 'B')

# Pvalues
pvalA = pnorm((statA), lower.tail=FALSE); hist(pvalA, breaks=sqrt(nA))
pvalB = pnorm((statB), lower.tail=FALSE); hist(pvalB, breaks=sqrt(nB))

# Fitting KerFdr 
KerFDRA = F_KerFdr(pvalA, rescaleKernWidth=TRUE); tauA = KerFDRA$tau; KerFDRA$p0
KerFDRB = F_KerFdr(pvalB, rescaleKernWidth=TRUE); tauB = KerFDRB$tau; KerFDRB$p0

# Make joint table
jointTests$tau = tauA[jointTests$A]*tauB[jointTests$B]
jointTests$stat = apply(cbind(statA[jointTests$A], statB[jointTests$B]), 1, min)
jointTests$H = HA[jointTests$A]*HB[jointTests$B]

# Pvalues based on observed pairs (not consistant with the assumption of conditional independence)
jointTests = F_JointPvalueFromObserved(jointTests)

# Check results
head(jointTests); hist(jointTests$pval, breaks=sqrt(P))
table((jointTests$pval < .05), jointTests$H)
KerFDRjoint = F_KerFdr(jointTests$pval, rescaleKernWidth=TRUE, plotting=TRUE)
KerFDRjoint$p0
