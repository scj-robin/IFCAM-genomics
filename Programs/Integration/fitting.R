
####################################################################
## Fitting distribution to gene expression pvalues
####################################################################

library(tidyverse)
library(KernSmooth)
library(mclust)

source('F_KerFdr.R')

a.GeneExp<-read.table("result.gene.expr.logcount.txt",header=T)
GetH1.E <- a.GeneExp$t.pvalue %>% F_KerFdr

a.Meth<-read.table("result.methylation.txt",header=T,sep="\t")
GetH1.M <- a.Meth$t.pvalue %>% F_KerFdr


