F_KerFdr <- function(P, p0=NULL, plotting=F, rescaleKernWidth=FALSE){
   # Fits a non parametric mixture for p-values
   # First estimate pi0
   # Then use probit transform X = -qnorm(P)
   # Then fit a mixture pi0 N(0, 1) + (1- pi0) F, 
   # to get a kernel estimate of the density f of distribution F
   library(KernSmooth); library(mclust); 
   
   if(plotting){par(mfrow=c(4, 4), mex=.6, pch=20)}

   # Data and sorting
   n = length(P); X = -qnorm(P); Xorder = order(X); Xrank = rank(X); 
   if(plotting){hist(P, breaks=sqrt(n)); plot(sort(P), (1:n)/n); hist(X, breaks=sqrt(n))}
   
   # Estimate of pi0 and phi0
   if(is.null(p0)){p0 = 2*sum((P > .5))/n}
   phi0 = dnorm(X)

   # Choice of w + Gram matrix : Gaussian kernel (small n)
   # w = dpik(X); G = exp(-as.matrix(dist(X, diag=T)^2)/2/w^2)/w/sqrt(2*pi)
   
   # Choice of w + Gram matrix : rectangular kernel  (large n)
   h = dpih(X); 
   if(rescaleKernWidth){h = h/sqrt(1-p0)}
   h.inv = 1/h
   Xneighbors = list()
   invisible(sapply(1:n, function(i){Xneighbors[[i]] <<- which(abs(X-X[i])<=h/2)}))
   
   # Initialize
   GM = Mclust(X, G=2, modelNames='E')
   tau = GM$z[, which.max(GM$parameters$mean)]; 

   # Algo
   # hist(P, breaks=sqrt(n), main='')
   tol = 1e-5; diff = 2*tol; iter = 0
   while(diff > tol){
      iter = iter+1; # cat(iter, '')
      # # f : Gaussian kernel
      # f = tau%*%G / sum(tau)
      # f : rectangular kernel
      f = rep(0, n)
      sapply(1:n, function(i){f[Xneighbors[[i]]] <<- f[Xneighbors[[i]]] + tau[i]*h.inv})
      f[which(X < -1)] = 0
      f = f / sum(tau)
      if(iter%%10==0){
         if(plotting){
            plot(X[Xorder], p0*phi0[Xorder], col=4, type='l', xlab='', ylab='', main=iter, xlim=c(min(X), max(X)));
            lines(X[Xorder], (1-p0)*f[Xorder], col=2)
         }
         cat(iter, '(', diff, ') ')
      }
      # Tau
      g = p0*phi0 + (1 - p0) * f
      tauNew = (1 - p0)*f / g
      # Test
      diff = max(abs(tau - tauNew))
      tau = tauNew
   }

   # # Post-processing tau
   # tauSort = tau[Xorder]
   # tauMinIndex = max(which(tauSort <= 1e-2))
   # tauSort[1:tauMinIndex] = 1e-2
   # tau = tauSort[Xrank]
   if(plotting){plot(X, tau)}
   
   return(list(X=X, p0=p0, f1=f, tau=tau))
}

