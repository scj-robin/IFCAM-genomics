# X0 = sample from H0
# X = vector of observed statistics
# Aim : compute P(X > x | H0) for each x
# Pb : X0 and X are huge

rm(list=ls()); par(pch=20)

# Dim
B = 1e8; n = 1e5
X0 = rnorm(B); X = rgamma(n, 1, 1)

# # Naive : too long
# Pval = sapply(1:n, function(i){if(i%%round(sqrt(n))==0){cat(i, '')}; Pval[i] <<- sum(X0 >= X[i])})/B

# # Empirical cdf : OK for B = 1e8, n = 1e5
# F0 = ecdf(X0); Pval = 1 - F0(X)

# Sorting : OK for B = 1e8, n = 1e6 but R cannot create a vector of length 1e9
# Indranil: split the 1e9 or 1e10 nulls into blocks of 1e8 and combine afterward
X0 = sort(X0, decreasing=TRUE)
Xorder = order(X, decreasing=TRUE); Xrank = n-rank(X)+1;
X = sort(X, decreasing=TRUE)
Pval = rep(0, n); i = 1; b = 1
sapply(1:n, function(i){
   while(X0[b] > X[i]){b <<- b+1}; Pval[i] <<- b-1
   if(i %% round(sqrt(n))==0){cat(i, '')}
})
Pval = Pval / B

Ptheo = 1-pnorm(X); Stheo = sqrt(Ptheo*(1-Ptheo)/B)
# plot(X, Pval); points(X, Ptheo, col=2)
plot(Ptheo, 2*Stheo, col=2, ylim=c(-2*max(Stheo), 2*max(Stheo)), type='l'); 
lines(Ptheo, -2*Stheo, col=2); abline(h=0, col=2)
lines(Ptheo, Pval-Ptheo)
