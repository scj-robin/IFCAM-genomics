# Reading and mergin of pvalues & map data

rm(list=ls())
library(data.table)

# Data
dataDir <- '../../../Results/PvaluesAndMap/'
expPval <- fread(paste0(dataDir, 'pvalue.expression.I.II.IV.txt'), data.table=FALSE)
methPval <- fread(paste0(dataDir, 'pvalue.methylation.I.II.IV.txt'), data.table=FALSE)
map <- fread(paste0(dataDir, 'map.gene.cpg.txt'), data.table=FALSE)
head(expPval); head(methPval); head(map)

# Renaming pvalues
names(expPval)[-1] <- paste0(names(expPval)[-1], '.exp')
names(methPval)[-1] <- paste0(names(methPval)[-1], '.meth')
rownames(expPval) <- expPval$gene; expPval <- as.matrix(expPval[, -1])
rownames(methPval) <- methPval$cpg; methPval <- as.matrix(methPval[, -1])
head(expPval); head(methPval); head(map)

# Histograms p-values
pdf(paste0(dataDir, 'pvaluesHist.pdf'))
par(mfcol=c(ncol(expPval), 2), mex=.5)
sapply(1:ncol(expPval), function(j){hist(expPval[, j], breaks=sqrt(nrow(expPval)), main=colnames(expPval)[j])})
sapply(1:ncol(methPval), function(j){hist(methPval[, j], breaks=sqrt(nrow(methPval)), main=colnames(methPval)[j])})
dev.off()

# Merging p-values
map$geneNum <- rep(NA, nrow(map))
map$cpgNum <- rep(NA, nrow(map))
for(i in 1:nrow(map)){
   if(i%%round(sqrt(nrow(map)))==0){cat(i, '')}
   geneNum <-  which(rownames(expPval)==map$gene[i])
   cpgNum <-  which(rownames(methPval)==map$cpg[i])
   if(length(geneNum)==1 & length(cpgNum)==1){
      map$geneNum[i] <- geneNum; map$cpgNum[i] <- cpgNum
   }
}
head(map)
mapPval <- map[-which(rowSums(is.na(map))>0), ]
head(mapPval)
mapPval <- cbind(mapPval, expPval[map$geneNum, ], methPval[map$cpgNum, ])
save(mapPval, file=paste0(dataDir, 'pvalue.map.I.II.IV.Rdata'))
