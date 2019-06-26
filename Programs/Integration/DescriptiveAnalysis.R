# Descriptive statistics
library(MASS)

rm(list=ls())
setwd('/home/robin/Bureau/RECHERCHE/EXPRESSION/IFCAM/IFCAM-genomics/Programs/Integration')
library(dplyr)

# Importation
Expression = readRDS("../../Data/Integration/expression.Rds")
Methylation = readRDS("../../Data/Integration/methylation.Rds")
Mutation = readRDS("../../Data/Integration/mutation.Rds")

# Functions
F_PlotSelectMethExp <- function(MethTab, ExpTab, location, diff=FALSE){
   # Plots expression, vs mean methylation for a given probe location type
   # MethTab=MethDiff; ExpTab=LogExpDiff; location = 'Body'
   if(location=='All'){
      MethTabSelect = MethTab %>% group_by(GeneName) %>% 
         summarise_at(colnames((MethTab)[-(1:3)]), mean)
   }else{
      MethTabSelect = MethTab %>% filter(ProbeLoc==location) %>% group_by(GeneName) %>% 
         summarise_at(colnames((MethTab)[-(1:3)]), mean)
   }
   rownames(MethTabSelect) = MethTabSelect$GeneName
   MethTabSelect$GeneName = c()
   MethTabSelect = as.matrix(MethTabSelect)
   ExpTabSelect = as.matrix(ExpTab[rownames(ExpTab)%in%rownames(MethTabSelect), ])
   GeneColor = (1:nrow(ExpTabSelect))%o%rep(1, ncol(ExpTabSelect))
   PatientColor = rep(1, nrow(ExpTabSelect))%o%(1:ncol(ExpTabSelect))
   if(diff){
      StatusColor = rep(1, nrow(ExpTabSelect))%o%(rep(1, ncol(ExpTabSelect)))
      }else{
         StatusColor = rep(1, nrow(ExpTabSelect))%o%(rep(c(1, 2), ncol(ExpTabSelect)/2))
         }
   par(mfrow=c(2, 2))
   plot(MethTabSelect, ExpTabSelect, col=GeneColor, pch=GeneColor, main=location)
   plot(MethTabSelect, ExpTabSelect, col=PatientColor, pch=PatientColor)
   plot(MethTabSelect, ExpTabSelect, col=StatusColor, pch=StatusColor)
   reg.coef = sapply(1:nrow(MethTabSelect), function(g){
      try(lm(ExpTabSelect[g, ] ~ MethTabSelect[g, ])$coef[2])
      })
   reg.coef = as.numeric(reg.coef[-which(substr(reg.coef, 1, 3)=='Err')])
   hist(reg.coef, breaks=sqrt(length(reg.coef)), main=length(reg.coef))
}

# Log-transform for expression
eps = 0; LogExp = log(eps+Expression)

# Normalisation : scaling (set mean=0, sd=1) across patients
LogExp = t(scale(t(LogExp)))
Methylation[, -(1:3)] = t(scale(t(Methylation[, -(1:3)])))
                              
# Methylation : a structure can be seen without scaling
MethHP = heatmap(as.matrix(Methylation[, -c(1:3)]), keep.dendro=T)
ClusThres = 8
plot(MethHP$Rowv); abline(h=ClusThres)
MethGroups = cutree(as.hclust(MethHP$Rowv), h=ClusThres)
heatmap(as.matrix(Methylation[, -c(1:3)]), RowSideColors=as.character(MethGroups))
# heatmap(as.matrix(Methylation[, -c(1:3)]), RowSideColors=as.character(as.vector(as.numeric(as.factor(Methylation$ProbeLoc)))))
ChosenGroup = 5
plot(MethGroups[MethHP$rowInd]); abline(h=ChosenGroup)
# Group of interest : 5
table(MethGroups, Methylation$ProbeLoc)[ChosenGroup, ]

# Genes from group 5
MethylationChosenGrp = Methylation[which(MethGroups==ChosenGroup), ]
GeneChosenGrp = unique(MethylationChosenGrp$GeneName)

# Mean methylation by probe location
LocationList = c(unique(Methylation$ProbeLoc), 'All')
sapply(LocationList, function(location){F_PlotSelectMethExp(Methylation, LogExp, location)})

# Computes difference between normal and tumor for LogExp and methylation 
TumorIndex = colnames(LogExp) %>% substr(., nchar(.), nchar(.)) %>% `==`('T') %>% which
NormalIndex = colnames(LogExp) %>% substr(., nchar(.), nchar(.)) %>% `==`('N') %>% which
LogExpDiff = as.matrix(LogExp[, TumorIndex] - LogExp[, NormalIndex])
MethDiff = as.matrix(Methylation[, 3+TumorIndex] - Methylation[, 3+NormalIndex])
MethDiff = cbind(Methylation[, 1:3], MethDiff)
invisible(sapply(LocationList, function(location){F_PlotSelectMethExp(MethDiff, LogExpDiff, location, diff=TRUE)}))

# Analysis gene by gene
GeneNum = 1; ChosenGene = rownames(Expression)[GeneNum]
MethylGene = Methylation[which(Methylation$GeneName==ChosenGene), ]
y = as.numeric(LogExp[GeneNum, ])
s = as.factor(substr(colnames(Expression), 6, 6))
par(mfrow=c(ceiling(sqrt(nrow(MethylGene))), round(sqrt(nrow(MethylGene)))))
boxplot(y ~ s)
sapply(1:nrow(MethylGene), function(i){
   x = as.numeric(MethylGene[i, -(1:3)])
   # boxplot(x ~ s, main=MethylGene$ProbeLoc[i])
   # Symbol = rep(Mutation[which(row.names(Mutation)==ChosenGene), ], each=2)
   plot(x,  y, col = s) # pch = 1+Symbol)
   LM = anova(lm(y ~ x*s))[3, 5]
})

# Association methylation-expression
s = as.factor(substr(colnames(Expression), 6, 6))
PvalJoint = PvalX = PvalY = Prod = list()
lapply(1:nrow(Expression), function(g){
   ChosenGene = rownames(Expression)[g]
   cat(g, ':')
   MethylGene = Methylation[which(Methylation$GeneName==ChosenGene), ]
   y = as.numeric(Expression[g, ])
   y = (y - mean(y))/sd(y)
   PvalJoint[[g]] <<- PvalX[[g]] <<- PvalY[[g]] <<- Prod[[g]] <<- rep(0, nrow(MethylGene))
   sapply(1:nrow(MethylGene), function(i){
      cat(i, '')
      x = as.numeric(MethylGene[i, -(1:3)])
      x = (x - mean(x))/sd(x)
      PvalX[[g]][i] <<- anova(lm(x ~ s))[1, 5]
      PvalY[[g]][i] <<- anova(lm(y ~ s))[1, 5]
      PvalJoint[[g]][i] <<- anova(manova(cbind(x, y) ~ s))[2, 6]
      Prod[[g]][i] <<- (t.test(x ~ s, var.equal=T)$statistic*t.test(y ~ s, var.equal=T)$statistic)^2
   })
   cat('\n')
})
df = length(y)-2; ProdThresh = qt(.975, df)^2
hist(log(unlist(Prod)), sqrt(length(unlist(Prod)))); abline(v=ProdThresh, col=2, lwd=2)
Color = 1 + 1*((unlist(Prod)^2)>ProdThresh)
hist(unlist(PvalJoint), sqrt(length(unlist(PvalJoint))))
hist(unlist(PvalX), sqrt(length(unlist(PvalX))))
hist(unlist(PvalY), sqrt(length(unlist(PvalY))))
plot(unlist(PvalX), unlist(PvalJoint), log='xy', col=Color)
plot(unlist(PvalY), unlist(PvalJoint), log='xy', col=Color)
plot(unlist(PvalX), unlist(PvalY), log='xy', col=Color)

