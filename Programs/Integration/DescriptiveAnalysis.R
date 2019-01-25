# Descriptive statistics

rm(list=ls())
setwd('/home/robin/Bureau/RECHERCHE/EXPRESSION/IFCAM/IFCAM-genomics/Programs/Integration')
library(dplyr)

# Importation
Expression = readRDS("../../Data/Integration/expression.Rds")
Methylation = readRDS("../../Data/Integration/methylation.Rds")
Mutation = readRDS("../../Data/Integration/mutation.Rds")

# Functions
F_PlotSelectMethExp <- function(MethTab, ExpTab, location){
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
   StatusColor = rep(1, nrow(ExpTabSelect))%o%(rep(c(1, 2), ncol(ExpTabSelect)/2))
   par(mfrow=c(2, 2))
   plot(MethTabSelect, ExpTabSelect, col=GeneColor, pch=GeneColor, main=location)
   plot(MethTabSelect, ExpTabSelect, col=PatientColor, pch=PatientColor)
   plot(MethTabSelect, ExpTabSelect, col=StatusColor, pch=StatusColor)
}

# Log-transform for expression
eps = 0; LogExp = log(eps+Expression)

# Normalisation : scaling (set mean=0, sd=1) across patients
LogExp = t(scale(t(LogExp)))
Methylation[, -(1:3)] = t(scale(t(Methylation[, -(1:3)])))
                              
# Methylation : a stucture can be seen without scaling
MethHP = heatmap(as.matrix(Methylation[, -c(1:3)]), keep.dendro=T)
plot(MethHP$Rowv); abline(h=1.2)
MethGroups = cutree(as.hclust(MethHP$Rowv), h=1.2)
heatmap(as.matrix(Methylation[, -c(1:3)]), RowSideColors=as.character(MethGroups))
plot(MethGroups[MethHP$rowInd]); abline(h=5)
# Group of interest : 5
table(MethGroups, Methylation$ProbeLoc)[5, ]

# Genes from group 5
MethylationGrp5 = Methylation[which(MethGroups==5), ]
GeneGrp5 = unique(MethylationGrp5$GeneName)

# Mean methylation by probe location
LocationList = c(unique(Methylation$ProbeLoc), 'All')
sapply(LocationList, function(location){F_PlotSelectMethExp(Methylation, LogExp, location)})

# Computes difference between normal and tumor for LogExp and methylation 
TumorIndex = colnames(LogExp) %>% substr(., nchar(.), nchar(.)) %>% `==`('T') %>% which
NormalIndex = colnames(LogExp) %>% substr(., nchar(.), nchar(.)) %>% `==`('N') %>% which
LogExpDiff = as.matrix(LogExp[, TumorIndex] - LogExp[, NormalIndex])
MethDiff = as.matrix(Methylation[, 3+TumorIndex] - Methylation[, 3+NormalIndex])
MethDiff = cbind(Methylation[, 1:3], MethDiff)
invisible(sapply(LocationList, function(location){F_PlotSelectMethExp(MethDiff, LogExpDiff, location)}))

