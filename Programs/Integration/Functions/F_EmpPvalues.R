F_EmpPvalues <- function(Xsort, X0){
   # Computes empirical p-values P(X0 >= X) for a (large) set of observed statistics X 
   # based on (large) sample X0 from H0
   # Observed X MUST BE RANKED IN DECREASING ORDER
   # n = 100; X = sort(rgamma(n, 1, 1), decreasing=T); B0 = 1e6; X0 = rnorm(B0)
   X0 = sort(X0, decreasing=T)
   n = length(X); B0 = length(X0); Pval = rep(0, n); b = 1
   sapply(1:n, function(i){
      # cat(i, b)
      while((X0[b] > X[i]) & (b < B0)){b = b + 1}
      Pval[i] <<- b-1
      # cat(b, '\n')
   })
   Pval = Pval / B0
   # plot(X, Pval, pch=20); lines(X, pnorm(X, lower.tail=FALSE), col=2, lwd=2)
   return(Pval)
}

