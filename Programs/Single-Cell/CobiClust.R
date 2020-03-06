rm(list=ls())
library(data.table); library(cobiclust)

# Data
dataDir <- '../../Data/Single-Cell/'
tab <- as.data.frame(fread(paste0(dataDir, 'T17225-counts.tsv')))
row.names(tab) <- tab[, 1]; tab <- tab[, -1]
dim(tab); tab[1:5, 1:5]; 
par(mfrow=c(2, 1))

# Filtering
thres <- 1
data <- as.matrix(tab[which(rowMeans(tab)>=thres), which(colMeans(tab)>=thres)])
dim(data); n <- nrow(data); p <- ncol(data)
prod(dim(data)) / prod(dim(tab))
heatmap(log(1+data), Colv=NA)
heatmap(log(diag(1/rowMeans(1+data))%*%(1+data)), Colv=NA)

# Clustering
Kmax <- 5; Gmax <- 5; par(mfrow=c(Kmax, Gmax), mex=.5)
ICL <- LB <- BIC <- a <- matrix(0, Kmax, Gmax)
for(k in 1:Kmax){
   for(g in 1:Gmax){
      res <- cobiclust(data, K=k, G=g, nu_j=colMeans(data))
      expected <- res$info$s_ik %*% log(res$parameters$alpha) %*% t(res$info$t_jg)
      image(1:p, 1:n, t(expected[order(res$info$s_ik%*%(1:k)), ]), main=paste(k, g), xlab='',)
      crit <- selection_criteria(res, k, g)
      ICL[k, g] <- crit[1]
      BIC[k, g] <- crit[2]
      LB[k, g] <- crit[4]
      a[k, g] <- res$parameters$a
      print(c(k, g))
      cat('pi =', res$parameters$pi, '\n')
      cat('rho =', res$parameters$rho, '\n')
      print(log(res$parameters$alpha))
   }
}
BIC
ICL
a


