F_KerFdr <- function(P){
   library(KernSmooth); library(mclust); 

   # Data and sorting
   n = length(P); X = -qnorm(P); Xorder = order(X); Xrank = rank(X); 
   par(mfrow=c(2, 2), pch=20, lwd=1); 

   # Estimate of pi0 and phi0
   p0 = 2*sum((P > .5))/n; phi0 = dnorm(X)

   # Choice of w + Gram matrix : Gaussian kernel (small n)
   # w = dpik(X); G = exp(-as.matrix(dist(X, diag=T)^2)/2/w^2)/w/sqrt(2*pi)
   
   # Choice of w + Gram matrix : histogramm  (large n)
   h = dpih(X); h.inv = 1/h
   Xneighbors = list()
   invisible(sapply(1:n, function(i){Xneighbors[[i]] <<- which(abs(X-X[i])<=h/2)}))
   
   # Initialize
   GM = Mclust(X, G=2, modelNames='E')
   tau = GM$z[, which.max(GM$parameters$mean)]; 

   # Algo
   par(mfrow=c(4, 4), mex=.6, pch=20, cex=.5, lwd=2)
   hist(P, breaks=sqrt(n), main='')
   tol = 1e-5; diff = 2*tol; iter = 0
   while(diff > tol){
      iter = iter+1; # cat(iter, '')
      # # f : Gaussian kernel
      # f = tau%*%G / sum(tau)
      # f : rectangular kernel
      f = rep(0, n)
      sapply(1:n, function(i){f[Xneighbors[[i]]] <<- f[Xneighbors[[i]]] + tau[i]*h.inv})
      f = f / sum(tau)
      if(iter%%10==0){
         # plot(X[Xorder], (1-p0)*f[Xorder], type='l', xlab='', ylab='', main=iter); 
         # lines(X[Xorder], p0*phi0[Xorder], col=4)
         cat(i, '(', diff, ') ')
      }
      # Tau
      g = p0*phi0 + (1 - p0) * f
      tauNew = (1 - p0)*f / g
      # Test
      diff = max(abs(tau - tauNew))
      tau = tauNew
   }

   # Post-processing tau
   tauSort = tau[Xorder]
   tauMinIndex = max(which(tauSort <= 1e-2))
   tauSort[1:tauMinIndex] = 1e-2
   tau = tauSort[Xrank]
   
   return(list(p0=p0, f1=f, tau=tau))
}

