# KerFdr
library(KernSmooth); library(mclust); library(Matrix)

# Parms
n = 1e5; pi0 = .9; pi1 = 1-pi0

# Data
Z = rbinom(n, 1, pi1)
X = rep(0, n)
X[which(Z==0)] = rnorm(sum(Z==0))
X[which(Z==1)] = rgamma(sum(Z==1), 2, 1)
X = sort(X)
f.true = dgamma(X, 2, 1)
P = pnorm(X, lower.tail=F)

# Estimate of pi0 and phi0
p0 = 2*sum((P > .5))/n; phi0 = dnorm(X)

# Choice of w + Gram matrix : Gaussian kernel
# w = dpik(X); G = exp(-as.matrix(dist(X, diag=T)^2)/2/w^2)/w/sqrt(2*pi)

# Choice of w + Gram matrix : histogramm
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
   iter = iter+1; cat(iter, '')
   # # f : Gaussian kernel
   # f = tau%*%G / sum(tau)
   # f : rectangular kernel
   f = rep(0, n)
   sapply(1:n, function(i){f[Xneighbors[[i]]] <<- f[Xneighbors[[i]]] + tau[i]*h.inv})
   f = f / sum(tau)
   if(iter%%10==0){
      plot(X, (1-p0)*f, type='l', xlab='', ylab='', main=iter); 
      lines(X, (1-p0)*f.true, col=2)
      lines(X, p0*phi0, col=4)
      cat('(', diff, ') ')
   }
   # Tau
   g = p0*phi0 + (1 - p0) * f
   tauNew = (1 - p0)*f / g
   # Test
   diff = max(abs(tau - tauNew))
   tau = tauNew
}
plot(X, tau, type='l')

.5*diff(X)%*%(f.true[-n] + f.true[-1])
.5*diff(X)%*%(f[-n] + f[-1])
