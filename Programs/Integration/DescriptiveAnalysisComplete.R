# Descriptive statistics
# Product of 2 gaussian: https://en.wikipedia.org/wiki/Product_distribution#Independent_central-normal_distributions
library(MASS)

rm(list=ls())
setwd('/home/robin/Bureau/RECHERCHE/EXPRESSION/IFCAM/IFCAM-genomics/Programs/Integration')
DataDir = ("../../../Data-NotUpload/")
library(dplyr)

# Function
F_StatMANOVA <- function(x, y, s){
   y = y-mean(y); x = x-mean(x); n = length(y)
   MANOVA = manova(cbind(x, y) ~ -1 + s, x=T)
   Parm = MANOVA$coefficients[1, ]
   VarParm = (n-1) / MANOVA$df * cov(MANOVA$residuals) / (n/2)
   # X = as.matrix(MANOVA$x)
   # VarParm = sum(((solve(crossprod(X))%*%t(X))[1, ])^2) * cov(MANOVA$residuals)
   Prod = prod(Parm)
   VarProd = (Parm[2:1]%*%VarParm%*%Parm[2:1])[1, 1]
   return((Prod/sqrt(VarProd)))
}

# Importation
Expression = readRDS(paste0(DataDir, "Expression.rds"))
Methylation = readRDS(paste0(DataDir, "Methylation.rds"))
MethInfo = readRDS(paste0(DataDir, "MethInfo.rds"))

# Log-transform for expression
eps = 0; LogExp = log(eps+Expression)
LogExp = t(scale(t(LogExp)))
# Methylation = t(scale(t(Methylation)))

# Remove NA
LogExp = LogExp[-which(rowSums(is.na(LogExp))>= ncol(Expression)/2), ]

# Status
TumorIndex = colnames(LogExp) %>% substr(., nchar(.), nchar(.)) %>% `==`('T') %>% which
NormalIndex = colnames(LogExp) %>% substr(., nchar(.), nchar(.)) %>% `==`('N') %>% which
Status = rep(0, ncol(LogExp)); Status[TumorIndex] = 1; Status = as.factor(Status)

# test CpGassoc
library(CpGassoc)
X = as.vector(as.matrix(Methylation[103, ]))
CPG = cpg.perm(Methylation[1:3, ], 1*(Status==1), nperm=1e1)

# Association methylation-expression: product of 2 tests
Pmanova = Ttest = Freg = matrix(NA, nrow(LogExp), max(table(MethInfo$Gene)))
sapply(1:1000, function(i){
   if(i%%10==0)(cat(i, ''))
   yDiff = LogExp[i, TumorIndex] - LogExp[i, NormalIndex]
   Texp <<- t.test(LogExp[i, TumorIndex], LogExp[i, NormalIndex])$stat
   yDiff = LogExp[i, TumorIndex] - LogExp[i, NormalIndex]
   MethTmp = Methylation[which(MethInfo$Gene==row.names(LogExp)[i]), ]
   invisible(sapply(1:nrow(MethTmp), function(j){
      Ttest[i, j] <<- Texp * 
         t.test(MethTmp[j, TumorIndex], MethTmp[j, NormalIndex])$stat
      Pmanova[i, j] <<- F_StatMANOVA(LogExp[i, ], as.vector(as.matrix(MethTmp[j, ])), Status)
      xDiff = as.vector(as.matrix(MethTmp[j, TumorIndex] - MethTmp[j, NormalIndex]))
      LM = lm(yDiff ~ xDiff)
      Freg[i, j] <<- anova(LM)[1, 4]
   }))
})
row.names(Freg) = row.names(Pmanova) = row.names(Ttest) = row.names(LogExp)
save(Freg, Pmanova, Ttest, file=paste0(DataDir, "Tests-LogExp-Meth.Rdata"))
par(mfcol=c(3, 2), mex=.6, pch=20)
H = hist(Ttest, breaks=sqrt(prod(dim(Ttest))), main='')
H = hist(Pmanova, breaks=sqrt(prod(dim(Pmanova))), main='')
H = hist(Freg, breaks=sqrt(prod(dim(Freg))), main='')
plot(abs(Ttest), abs(Pmanova), log='xy')
plot(abs(Ttest), Freg, log='xy')
plot(abs(Pmanova), Freg, log='xy')
# H = hist(log(abs(TtestVec)), breaks=sqrt(length(TtestVec)))
w = mean(diff(H$breaks)); B = length(TtestVec)
x.list = seq(-50, 50, length.out=1001)
# lines(x.list, B/2*w*exp(-abs(sqrt(2)*x.list)), col=2, lwd=2)
lines(x.list, B*w*besselK(abs(x.list), 0)/pi, col=4, lwd=2)

# Nul prod distribution
B = 1e4; df = ncol(Expression)-2
T00 = rt(B, df=df)*rt(B, df=df)
H = hist(T00, breaks=sqrt(B))
w = mean(diff(H$breaks))
x.list = seq(-5, 5, length.out=101)
lines(x.list, B*w*besselK(abs(x.list), 0)/pi, col=4, lwd=2)

# Nice example : PARP12 / 11
GeneName = 'PARP12'; ProbeNb = 11
y = as.vector(as.matrix(Expression[which(rownames(Expression)==GeneName),]))
y = log(y)
yy = as.vector(as.matrix(LogExp[which(rownames(LogExp)==GeneName),]))
xTmp = Methylation[which(MethInfo$Gene==GeneName), ]; 
x = as.vector(as.matrix(xTmp[ProbeNb, ]))
plot(x, y, col=Status)

# Test stat Fmanova
B= 1e3; Fs = rep(0, B); s = Status
for(b in 1:B){Fs[b] = F_StatMANOVA(rnorm(80), rnorm(80), s)}
sd(Fs)
