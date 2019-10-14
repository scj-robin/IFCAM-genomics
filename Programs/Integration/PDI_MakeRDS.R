## The goal of the present program is to shape the different datasets
## to be used in later programs. For now only the Methylation and Expression 
## data are processed.
## Four datasets are produced:
## Expression, GeneInfo, Methylation, MethylationInfo
## where the two "Info" files correspond to the meta-data describing 
## the genes and methylation probes. All these datasets are stored as .rds files

library(data.table)


        #### Parameters


# DataRep <- 'D:/IFCAM/IFCAM-genomics/Data/PrivateIntegration/'
DataRep <- '../../Data/PrivateIntegration/'
MutationFileName <- 'Mutation.txt'
MethylationFileName <- "Methylation.txt"
ExpressionFileName <- "fpkm_values_all_genes_from_40_patients_for_ISI_single_entry_JUNE.txt"


        #### Load the data


Methylation <-  fread(paste0(DataRep,MethylationFileName),data.table = FALSE)
Expression <-  fread(paste0(DataRep,ExpressionFileName),data.table = FALSE)
dim(Expression)
dim(Methylation)

        
        #### Expression


## Remove micro Rnas
GeneNames <- Expression$gene_short_name
WhichMiRNA <- which(substr(GeneNames,1,3) == 'MIR')
GeneNames <- GeneNames[-WhichMiRNA]
Expression <- Expression[-WhichMiRNA,]

## Remove gene names and information
GeneInfo <- Expression$locus
Expression <- Expression %>% 
  select(starts_with('RADS'), starts_with('pt_'))

## Add gene names as row names
row.names(Expression) <- GeneNames
names(GeneInfo) <- GeneNames

## Unify tissue status designation
## Patients were collected in two separated rounds. In the first round the tissues 
## were named according to the sample ID + 'N' for normal or 'T' for tumor. A similar 
## coding was used in the second round, except that tumors are coded with a 'D'. The following 
## chunk of code unifies the tissue designation by changing the 'D' into 'T'.
Unify <- function(Names){
  tibble(Names = Names) %>% 
    mutate(Last = substr(Names,nchar(Names),nchar(Names)),
         NewNames = ifelse(test = Last=='D',yes = paste0(substr(Names,1,nchar(Names)-1),'T'),no = Names)) %>% 
    pull(NewNames) 
     
}
colnames(Expression) = Unify(colnames(Expression))
Expression[1:10,1:10]
dim(Expression)


        #### Methylation


## First target methylation sites that correspond to genes whose expression
## is available, and vice-versa
ExpButNotMeth <- setdiff(GeneNames,Methylation$Gene)
MethButNotExp <- setdiff(Methylation$Gene,GeneNames)
Expression <- Expression[-which(GeneNames %in% ExpButNotMeth),]
GeneInfo <- GeneInfo[-which(GeneNames %in% ExpButNotMeth)]
Methylation <- Methylation[-which(Methylation$Gene %in% MethButNotExp),]

## Get Methylation probe names
ProbeNames <- Methylation$V1

## Get Methylation info
MethylationInfo <- Methylation %>% 
  select(Gene,Group)
row.names(MethylationInfo) <- ProbeNames

## Unify tissue designation
Methylation <- Methylation %>% 
  select(starts_with('RADS'),starts_with('pt_'))
colnames(Methylation) <- Unify(colnames((Methylation)))
colnames(Methylation)
row.names(Methylation) <- ProbeNames
dim(Methylation)

## CAUTION: reorder columns in Expression so that ordering is the one in Methylation
Expression[1:5,1:10]
Expression <- Expression[,order(match(colnames(Expression),colnames(Methylation)))]
Expression[1:5,1:10]


        #### Save all files


saveRDS(Expression,paste(DataRep,'Expression.rds'))
saveRDS(GeneInfo,paste(DataRep,'GeneInfo.rds'))
saveRDS(Methylation,paste(DataRep,'Methylation.rds'))
saveRDS(MethylationInfo,paste(DataRep,'MethInfo.rds'))
