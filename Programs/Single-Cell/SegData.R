dataDir <- '../../Data/Single-Cell/'

tab <- read.table(paste0(dataDir, 'T17225-counts.tsv'), sep='\t', header=TRUE, row.names=1)
dim(tab)
tab[1:5, 1:5]
hist(log10(colSums(tab)), breaks=sqrt(ncol(tab)))
hist(log10(rowSums(tab)), breaks=sqrt(nrow(tab)))
summary(colSums(tab))
summary(rowSums(tab))

library(Segmentor3IsBack); library(gfpop)
n = 1e5; Y = rpois(n, 5); Y <- tab[, 2]
myGraph <- graph(penalty = 10*log(n), type = "std")
seg <- gfpop(data=Y, mygraph=myGraph, type="poisson")
seg$parameters; 
length(seg$parameters); sum(abs(diff(seg$parameters))!=0)


# Pour des donnees multivairees ponderees, cf EM
plot(seg$changepoints, log(1+seg$parameters), type='s')
n <- 10; Y <- as.matrix(tab[, 1:n]); tau <- runif(n)
seg <- gfpop(data=as.vector(round(Y%*%tau)), mygraph=myGraph, type="poisson")
seg$parameters <- seg$parameters / sum(tau)
