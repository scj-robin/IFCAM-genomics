F_EmpPvalues <- function(Xsort, X0, W0=rep(1, length(X0))){
   # Computes empirical p-values P(X0 >= X) for a (large) set of observed statistics X 
   # based on (large) sample X0 from H0
   # Observed X MUST BE RANKED IN DECREASING ORDER
   # n = 100; X = sort(rgamma(n, 1, 1), decreasing=T); B0 = 1e4; X0 = rnorm(B0); W0 = rep(1, B0)
   W0 = W0/sum(W0); 
   W0 = W0[order(X0, decreasing=T)]; X0 = sort(X0, decreasing=T)
   n = length(Xsort); B0 = length(X0); Pval = rep(0, n); b = 1
   sapply(1:n, function(i){
      if(i%%round(sqrt(n))==0){cat(i, '')}
      if(i > 1){Pval[i] <<- Pval[i-1]}
      while((X0[b] > Xsort[i]) & (b < B0)){
         Pval[i] <<- Pval[i]+W0[b]; 
         b <<- b + 1
      }
   })
   cat('\n')
   # plot(Xsort, Pval, pch=20); lines(Xsort, pnorm(Xsort, lower.tail=FALSE), col=2, lwd=2)
   return(Pval)
}

