# Reading expression, methylation and mutation data
# Creating an Rds file for each

rm(list=ls())
setwd('/home/robin/Bureau/RECHERCHE/EXPRESSION/IFCAM/IFCAM-genomics/Programs/Integration')
library(dplyr)

# Expression data
DataExp <- read.csv("../../Data/Integration/expression.txt", sep='\t', header=T, as.is=T)
PatientIndexExp <- DataExp %>% colnames %>% substr(.,1,1) %>% `==`('X') %>% which
PatientNb = length(PatientIndexExp)/2
Expression = DataExp[, PatientIndexExp]; 
row.names(Expression) = DataExp$gene_short_name
colnames(Expression) = gsub('_0', '', colnames(Expression)) %>% gsub('X', 'Pt', .)
TumorIndex = colnames(Expression) %>% substr(., nchar(.), nchar(.)) %>% `==`('T') %>% which
# NormalIndex = colnames(Expression) %>% substr(., nchar(.), nchar(.)) %>% `==`('N') %>% which

# Methylation data
DataMeth <- read.csv("../../Data/Integration/methylation.txt", sep='\t', header=T, as.is=T)
PatientIndexMeth <- DataMeth %>% colnames %>% substr(.,1,2) %>% `==`('pt') %>% which
Methylation = DataMeth[, c(1, 2, ncol(DataMeth), (3:(2+length(PatientIndexMeth))))]
colnames(Methylation) = gsub('pt_', 'Pt', colnames(Methylation))
colnames(Methylation)[1] = 'GeneName'
colnames(Methylation)[2] = 'ProbeId'
colnames(Methylation)[3] = 'ProbeLoc'

# Mutation data
DataMut <- read.csv("../../Data/Integration/mutation.txt", sep='\t', header=T, as.is=T)
dim(DataMut)
MutationList = DataMut %>% group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>% summarise(N=n())
MutationList$Tumor_Sample_Barcode = gsub('Patient', 'Pt', MutationList$Tumor_Sample_Barcode)
Mutation = matrix(0, nrow(Expression), length(TumorIndex)); 
row.names(Mutation) = row.names(Expression)
colnames(Mutation) = colnames(Expression)[TumorIndex] %>% gsub('T', '', .)
invisible(sapply(1:nrow(MutationList), function(i){
   Mutation[which(row.names(Mutation)==MutationList$Hugo_Symbol[i]), which(colnames(Mutation)==MutationList$Tumor_Sample_Barcode[i])] <<- 
      Mutation[which(row.names(Mutation)==MutationList$Hugo_Symbol[i]), which(colnames(Mutation)==MutationList$Tumor_Sample_Barcode[i])] + MutationList$N[i]
   }))

# Exportation
saveRDS(Expression, "../../Data/Integration/expression.Rds")
saveRDS(Methylation, "../../Data/Integration/methylation.Rds")
saveRDS(Mutation, "../../Data/Integration/mutation.Rds")

