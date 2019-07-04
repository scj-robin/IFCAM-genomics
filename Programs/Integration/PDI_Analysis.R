rm(list=ls())
library(betareg)


        #### Parameters


DataRep <- 'D:/IFCAM/IFCAM-genomics/Data/PrivateIntegration/'
TypeOfTest <-  'Welch'  #  'Wilcoxon' " 'Student' # 


        #### Functions


TtestFunction <- function(x,y){
  
  ## Check x for -Inf
  Labs <- Labels
  ToBeRemoved <- which(x == -Inf)
  if(length(ToBeRemoved)>0){
    x <- x[-ToBeRemoved]
    Labs <- Labels[-ToBeRemoved]
    y <- y[,-ToBeRemoved]
  }
  ## Version wilcoxon
  if (TypeOfTest == 'Wilcoxon'){
    Res <- list(TExp = wilcox.test(x~Labs)$p.value %>% qnorm, TMeth = sapply(1:nrow(y), function(probe){ wilcox.test(as.numeric(y[probe,])~Labs)$p.value %>% qnorm } ))
  }
  ## Version student
  if (TypeOfTest == 'Student'){
    Res <- list(TExp = t.test(x~Labs,var.equal=TRUE)$p.value, TMeth = sapply(1:nrow(y), function(probe){ t.test(as.numeric(y[probe,])~Labs,var.equal=TRUE)$p.value } ))
  }
  ## Version welch
  if (TypeOfTest == 'Welch'){
    Res <- list(TExp = t.test(x~Labs,var.equal=FALSE)$p.value, TMeth = sapply(1:nrow(y), function(probe){ t.test(as.numeric(y[probe,])~Labs,var.equal=FALSE)$p.value } ))
  }
  return(Res)
}


TtestWelchAndBeta <- function(x,y){
  
  ## Check x for -Inf
  Labs <- Labels
  ToBeRemoved <- which(x == -Inf)
  if(length(ToBeRemoved)>0){
    x <- x[-ToBeRemoved]
    Labs <- Labels[-ToBeRemoved]
    y <- y[,-ToBeRemoved]
  }

  PExp = t.test(x~Labs,var.equal=FALSE)$p.value
  PMeth = purrr::map(1:nrow(y), function(probe){ betareg(as.numeric(y[probe,])~Labs) }) %>% 
    purrr::map(summary) %>%
    purrr::map(~ .x$coef$mean[2,4]) %>% 
    reduce(c)
  Res <- list(PExp=PExp,PMeth=PMeth)  
  return(Res)
}

        #### Load the data


Expression <- readRDS(paste0(DataRep,'Expression.rds'))
Methylation <- readRDS(paste0(DataRep,'Methylation.rds'))
MethylationInfo <- readRDS(paste0(DataRep,'MethInfo.rds'))
GeneInfo <- readRDS(paste0(DataRep,'GeneInfo.rds'))

## Get labels 
Labels <- purrr::map(colnames(Expression), ~ substr(.x,nchar(.x),nchar(.x))) %>% 
  reduce(c)


        #### Filter expression


MaxForT <- apply(Expression[,Labels=='T'],1,max) 
table(MaxForT==0)
NonExpT <- names(MaxForT)[MaxForT==0]
MaxForN <- apply(Expression[,Labels=='N'],1,max) 
table(MaxForN==0)
NonExpN <- names(MaxForN)[MaxForN==0]
NonExp <- union(NonExpN,NonExpT)


        #### Tidy the data


## Get gene info
GeneDF <- GeneInfo %>% 
  data.frame %>% 
  setNames("Pos") %>% 
  rownames_to_column(var = 'Gene') %>% 
  mutate(Chr = str_replace(str_extract(Pos,pattern = ".+(?=\\:)"),"chr",""),
         GeneStart = str_extract(Pos,pattern = "(?<=\\:).+(?=\\-)"),
         GeneEnd = str_extract(Pos,pattern = "(?<=\\-).+")
         ) %>% 
  select(-Pos)

## Get a tidy Methylation dataset
TidyMeth <- Methylation %>% 
  as.tibble %>% 
  mutate(Gene = MethylationInfo$Gene) %>% 
  nest(-Gene) %>% 
  rename(Meth=data) %>% 
  mutate(Meth = purrr::map(Meth, as.data.frame))

## Get a tidy Expression dataset
TidyExp <- Expression %>% 
  as.tibble %>% 
  rownames_to_column('Gene') %>% 
  nest(-Gene) %>% 
  rename(Exp=data) %>% 
  mutate(Exp = purrr::map(Exp, ~ as.data.frame(.x) %>% as.numeric),
         LogExp = purrr::map(Exp,log))

## Get a global tidy dataset
ExpMeth <- TidyMeth %>% 
  left_join(y=TidyExp) %>% 
  left_join(y=GeneDF) %>% 
  mutate(NbInfExp = map_dbl(LogExp, ~ sum(.x==-Inf))) %>% 
  filter(!(Gene %in%NonExp), NbInfExp < 30) 
ExpMeth


        #### Perform analysis


# ExpMethTest <- ExpMeth %>% 
#   mutate(Test = map2(Exp,Meth,TtestWelchAndBeta))
# saveRDS(ExpMethTest, paste0(DataRep,'ExpMeth_',TypeOfTest,'.rds'))

ExpMeth3.5 <- ExpMeth %>% 
  mutate(NbProbes = purrr::map(Meth,nrow)) %>% 
  filter(NbProbes >=3, NbProbes <=5) %>% 
  mutate(Test = map2(Exp,Meth,TtestWelchAndBeta))
saveRDS(ExpMeth3.5, paste0(DataRep,'ExpMeth_WelchAndBeta.rds'))


        #### Garbage


# ## Have a look at the number of probes per gene
# TidyMeth %>% 
#   mutate(NbProbesPerGene = map_int(Meth,nrow)) %>% 
#   summarise_at(.vars = c('NbProbesPerGene'), .funs = c('min','mean','median','max'))
# 
# 
# ## Focus on a gene
# NbProbe = 9
# Example <- TidyMeth %>% 
#   mutate(NbProbesPerGene = map_int(Meth,nrow)) %>% 
#   filter(NbProbesPerGene == NbProbe) %>% 
#   slice(2:2)
# 
# Example <- TidyMeth %>% 
#   filter(Gene == 'CASP8') 
# 
# 
# DF <- t(Example$Meth[[1]])[,1:20]
# Y <- log(1+Expression[Example$Gene,])
# x11()
# par(mfrow=c(ceiling(sqrt(ncol(DF))),round(sqrt(ncol(DF)))))
# for (ii in 1:ncol(DF)){
#   plot(DF[,ii],Y, main = paste("Probe",ii),col=as.factor(Labels))
# }


# Summary <- ExpMeth %>% 
#   select(Gene,Test) %>% 
#   unnest %>% 
#   arrange(Test) %>% 
#   filter(abs(Test) > 100)
# 
# Example <- TidyMeth %>% 
#   filter(Gene == 'PARP12') # LRP3
# 
# DF <- t(Example$Meth[[1]])[,1:min(20,nrow(Example$Meth[[1]]))]
# Y <- log(Expression[Example$Gene,])
# x11()
# par(mfrow=c(ceiling(sqrt(ncol(DF))),round(sqrt(ncol(DF)))))
# for (ii in 1:ncol(DF)){
#   plot(DF[,ii],Y, main = paste("Probe",ii),col=as.factor(Labels))
# }
# 
# Ydiff <- Y[seq(1,80,2)]-Y[seq(2,80,2)]
# DFdiff <- DF[seq(1,80,2),]-DF[seq(2,80,2),]
# plot(DFdiff[,ii],Ydiff, main = paste("Probe",ii),pch=16)
