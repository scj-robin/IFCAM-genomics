## The goal here is to build a tibble combining for each gene the methylation 
## and expression data, along with the pvalues of the tests that are performed 
## on each type of data. The resulting tibble is ExpMeth (or ExpMeth3.5 if only
## a subsample of the complete data set is wanted) that is store in a .rds format.

rm(list=ls())
library(betareg); 
library(tidyverse); 


        #### Parameters


DataRep <- 'D:/IFCAM/IFCAM-genomics/Data/PrivateIntegration/'
#DataRep <- '../../Data/PrivateIntegration/'
TypeOfTest <-  'Welch'  #  'Wilcoxon' " 'Student' # 


        #### Functions


## Corresponds mostly to initial attemps, should be obsolete now
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

## Corresponds to the function that performs a Welch unpaired T-test for expression,
## and a beta-regression test for methylation.
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


## Check for each class (T and N) whether there is some expression or not.
## If not, the gene will be removed later.
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
  merge(MethylationInfo,.,by="row.names") %>% 
  rename(ProbeName = Row.names) %>% 
  as.tibble %>% 
  nest(-Gene) %>% 
  mutate(Meth = purrr::map(data, ~ .x[,3:ncol(.x)] %>% as.data.frame),
         ProbeInfo = purrr::map(data, ~ .x[,1:2] %>% as.data.frame)
         ) %>% 
  select(-data)

## Get a tidy Expression dataset
TidyExp <- Expression %>% 
   rownames_to_column('Gene') %>% 
   as.tibble %>% 
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


## Complete datasets
# ExpMethTest <- ExpMeth %>% 
#   mutate(Test = map2(Exp,Meth,TtestWelchAndBeta))
# saveRDS(ExpMethTest, paste0(DataRep,'ExpMeth_',TypeOfTest,'.rds'))

## Only genes with 3 to 5 associated methylation probes
ExpMeth3.5 <- ExpMeth %>% 
  mutate(NbProbes = purrr::map(Meth,nrow)) %>% 
  filter(NbProbes >=3, NbProbes <=5) %>% 
  mutate(Test = map2(Exp,Meth,TtestWelchAndBeta))
saveRDS(ExpMeth3.5, paste0(DataRep,'ExpMeth_WelchAndBeta.rds'))


 