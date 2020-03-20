

########################   GeneratePvalues   ####################

## Function to generate the pvalues according to the model

GeneratePvalues <- function(n,Q,DirichletProp,Hconfig,Hconfig.H1,
                            MinPropForHconfig.H1,Delta,
                            Rho0,Rho1,MinNbOfUnitPerHconfig.H1){

  ## Build H configurations
  Hnumber = length(Hconfig)

  ## Build priors
  # Generate H0 and H1 proportion per class of test
  piq <- rdirichlet(Q, DirichletProp)
  # Generate proportions per Hconfig
  pi <- sapply(Hconfig, function(h){
    prod((piq[,1]**(1-h))*(piq[,2]**h))
  })
  # Make sure configs of interest are represented enough
  if (min(pi[Hconfig.H1])<MinPropForHconfig.H1){
    pi[Hconfig.H1] <- apply(cbind(pi[Hconfig.H1],MinPropForHconfig.H1),1,max)
  }
  # And shake a bit
  pi <- rdirichlet(1, n*pi)

  ## Build mean matrix per Hconfig
  mu = t(sapply(Hconfig, function(h) h*Delta))
  
  ## Build covariance matrix per Hconfig
  Sigma0 = matrix(Rho0, Q, Q); diag(Sigma0) = 1
  Sigma1 = matrix(Rho1, Q, Q); diag(Sigma1) = 1
  Sigma = lapply(Hconfig, function(h){
    Res <- Sigma0
    Res[which(h==1), which(h==1)] <- Sigma1[which(h==1), which(h==1)]
    return(Res)
  })

  ## Generate test statistics
  NotEnoughObsInHconfig.H1 <- TRUE
  while(NotEnoughObsInHconfig.H1){
    Z = sample(1:Hnumber,size = n,replace = TRUE,prob = pi)
    Tbl <- table(Z)
    NbObsMin <- min(Tbl[names(Tbl) %in% Hconfig.H1])
    NotEnoughObsInHconfig.H1 <- NbObsMin < MinNbOfUnitPerHconfig.H1 
  }
  Y = matrix(0, n, Q)
  sapply(1:Hnumber, function(h){
    if(length(which(Z==h))>0){
      Y[which(Z==h), ] <<- rmvnorm(sum(Z==h), mean=mu[h, ], sigma=Sigma[[h]])
    }
  })
  P = pnorm(Y, lower.tail=FALSE)

  ## Return the result
  return(list(P=P, Truth=Z,pi=pi,mu=mu,Sigma=Sigma,Y=Y))
}



########################   GetHinfo   ####################

GetHinfo <- function(Q,AtLeast){

  ## Build H configurations
  Hconfig <- as.matrix(expand.grid(lapply(1:Q, function(q) 0:1)))
  Hconfig <- split(Hconfig, seq(2^Q))
  
  ## Find the ones that match H1
  MatchingH1 <- sapply(Hconfig, function(h){sum(h)>=AtLeast})
  Hconfig.H1 <- which(MatchingH1)
  names(Hconfig.H1) <- NULL
  
  ## Collect results
  return(list(Hconfig=Hconfig,Hconfig.H1=Hconfig.H1))
  
}  



########################   GetPH1   ####################

GetPH1 <- function(PvalMat,AtLeast){
  
  Q = ncol(PvalMat)
  if(AtLeast==Q){
    Res <- apply(PvalMat, 1, max)
  } else {
    Res <- sapply(1:nrow(PvalMat), function(x) sort(PvalMat[x,])[AtLeast])
  }
  return(Res)
}



########################   FastKerFdr   ####################

FastKerFdr <- function(Pval,p0=NULL,plotting=FALSE,NbKnot=1e5,tol = 1e-5){

  ## Transform pvalues into N(0,1) quantiles
  n = length(Pval)
  X = -qnorm(Pval)

  ## Get a p0 estimate
  if(is.null(p0)){
    p0 = 2*sum(X < 0)/n; p1 = 1 - p0    
  }
  
  ## Knots, counts and initialization (using Mclust)
  if(length(X)>NbKnot){
    Hist = hist(X, breaks=NbKnot, plot=FALSE) 
    Knots = Hist$mids; ActualNbKnot = length(Knots); Counts = Hist$counts;
    Xsample = sample(X, NbKnot)
    GM = Mclust(Xsample, G=3, modelNames='E'); 
    mu = max(GM$parameters$mean) 
  } else {
    Knots = X; ActualNbKnot = length(X); Counts = rep(1,ActualNbKnot); 
    GM = Mclust(X, G=3, modelNames='E'); 
    mu = max(GM$parameters$mean)
  }
  if (plotting){
    Order <- order(Knots)
    Knots <- Knots[Order]
    Counts <- Counts[Order]
  }
  
  ## Initialize the taus using GM
  phi = dnorm(Knots); f1 = dnorm(Knots, mean=mu, sd=1)
  tau = p1*f1/(p0*phi + p1*f1)

  ## Get the weighted kernel density estimate
  diff = 2*tol; iter = 0
  while(diff > tol){
    iter = iter + 1
    weights = tau*Counts; weights = ActualNbKnot * weights / sum(weights) 
    f1 = kde(x=Knots, w=weights, eval.points=Knots)$estimate
    tauNew = p1*f1/(p0*phi + p1*f1)
    diff = max(abs(tau - tauNew))
    tau = tauNew
  }
  if(plotting){
    Max=max(p0*phi+p1*f1)
    plot(Knots, p0*phi, type='l', main='', col=4,lwd=2,ylim=c(0,Max),
         xlab="Q-transformed pvalues",ylab="Densities"); 
    lines(Knots, p1*f1, col=2,lwd=2); 
    lines(Knots, p0*phi+p1*f1 + 0.01*Max,lwd=2,lty=2)
    legend("topright", legend=c("H0 dist", "H1 dist","Mixture Dist"),
           col=c("blue","red", "black"), lty=c(1,1,2), cex=0.8)  
  }
  
  ## Get the Tau
  KDE = kde(x=Knots, w=weights, eval.points=X)
  f1 = KDE$estimate
  tau = p1*f1 / (p1*f1 + p0*dnorm(X))
  
  return(list(p0=p0,tau=tau))
}  



########################   GetTrueMargTau   ####################

GetTrueMargTau <- function(Tmp,Hconfig){
  pi.config <-  Tmp$pi
  Y <- Tmp$Y
  mu <- Tmp$mu
  Res <- sapply(1:ncol(Y), function(k){
    
    ## Get the marginal prior of being H1
    pi1.marg <- sapply(Hconfig, function(x) x[k]==1) %>% 
      which %>% 
      pi.config[.] %>% 
      sum
    print(pi1.marg)
    
    ## Compute H0 and H1 densities
    Densities <- sapply(unique(mu[,k]), function(mu.loc){
      dnorm(Y[,k],mean=mu.loc,sd=1)
    })
    
    ## Post. prob 
    Tau <- pi1.marg*Densities[,2] / ((1-pi1.marg)*Densities[,1] + pi1.marg*Densities[,2])
    return(Tau)
  })
  return(Res)
}



########################   Tau2Cdf   ####################

Tau2Cdf <- function(P, Tau, Quant){
  ewcdf(P, weights=Tau)(Quant)
}



########################   ComputePValue   ####################


ComputePValue <- function(PvalMat,LogTau.MixtMarg,Hconfig,Hconfig.H1,Prior.Hconfig,P.H1,Method='H0'){

  ## Basic stuff
  Q <- ncol(PvalMat)
  n <- nrow(PvalMat)
  Prior.Hconfig.H1 <- Prior.Hconfig[Hconfig.H1]

  ######   First case: full H1   ######
  if(length(Hconfig.H1)==1){

    if(Method == 'H1'){
  
      ## Infer cdf of each Hconfig.H1
      # Cdf.Hconfig.H1 <- sapply(Hconfig[Hconfig.H1], function(h){
      #   P.H1.h = apply(PvalMat[, which(h==1), drop=FALSE], 1, max)
      #   tau.h <- rowSums(LogTau.MixtMarg[, which(h==1), drop=FALSE]) %>% exp
      #   H1comp <- Tau2Cdf(P.H1.h, tau.h, P.H1)
      #   Res<- H1comp * (P.H1^sum(h==0))
      # })
      Cdf.Hconfig.H1 <- LogTau.MixtMarg %>% 
        rowSums %>% 
        exp %>% 
        Tau2Cdf(P.H1, ., P.H1) 
      Num <- Tau2Cdf(P.H1,rep(1/n,n),P.H1) - Cdf.Hconfig.H1*Prior.Hconfig.H1
      Den <-1-sum(Prior.Hconfig.H1)
      PvalUnion =  Num/Den 
    }
    
    if(Method == 'H0'){
      
      Cdf.Hconfig <- sapply(Hconfig[-Hconfig.H1], function(h){
        if(sum(h) > 0){
          P.H1.h = apply(PvalMat[, which(h==1), drop=FALSE], 1, max)
          tau.h <- rowSums(LogTau.MixtMarg[, which(h==1), drop=FALSE]) %>% exp
          H1comp <- Tau2Cdf(P.H1.h, tau.h, P.H1)
          Res<- H1comp * (P.H1^sum(h==0))
        } else {
          Res <- P.H1^Q
        }
      })
      PvalUnion = Cdf.Hconfig%*%Prior.Hconfig[-Hconfig.H1] / sum(Prior.Hconfig[-Hconfig.H1])
      
    }
  
  ######   Second case: composite H1   ######
  } else {
    
    if(Method == 'H1'){
      
      ## Get the marginal H1 cumulative distribution functions
      ## estimated at p = Pval_q for each serie q
      # Cdf.MixtMarg <- sapply(1:Q, function(q){
      #   Tau2Cdf(PvalMat[, q], exp(LogTau.MixtMarg[,q]), P.H1)
      # })
      LogCdf.MixtMarg <- sapply(1:Q, function(q){
        Tau2Cdf(PvalMat[, q], exp(LogTau.MixtMarg[,q]), P.H1) %>% log
      })
      
      ## Infer cdf of each Hconfig.H1
      Cdf.Hconfig.H1 <- sapply(Hconfig[Hconfig.H1], function(h.cfg){
        Tmp <- sapply(Hconfig[Hconfig.H1], function(h.lower){
          
          #Deal with the H0 (lower and higher)
          NbH0Lower <- sum((h.cfg==0) & (h.lower==1))
          NbH0Higher <- sum((h.cfg==0) & (h.lower==0))
          LogCdf <- NbH0Lower*log(P.H1) + NbH0Higher*log(1-P.H1)
          #Deal with H1 and lower
          WhichH1Lower <- which((h.cfg==1) & (h.lower==1))
          if(length(WhichH1Lower)>0){
            LogCdf <- LogCdf + rowSums(LogCdf.MixtMarg[,WhichH1Lower,drop=FALSE])
          }
          #Deal with H1 and higher
          WhichH1Higher <- which((h.cfg==1) & (h.lower==0))
          if(length(WhichH1Higher)>0){
            LogCdf <- LogCdf + rowSums(log(1-exp(LogCdf.MixtMarg[,WhichH1Higher,drop=FALSE])))
          }
          return(exp(LogCdf))
        })
        return(rowSums(Tmp))
      })
      # Cdf.Hconfig.H1.bis <- sapply(Hconfig[Hconfig.H1], function(h.cfg){
      #   Tmp <- sapply(Hconfig[Hconfig.H1], function(h.lower){
      #     
      #     #Deal with the H0 (lower and higher)
      #     NbH0Lower <- sum((h.cfg==0) & (h.lower==1))
      #     NbH0Higher <- sum((h.cfg==0) & (h.lower==0))
      #     Cdf <- (P.H1^NbH0Lower)*((1-P.H1)^NbH0Higher)
      #     #Deal with H1 and lower
      #     WhichH1Lower <- which((h.cfg==1) & (h.lower==1))
      #     if(length(WhichH1Lower)>0){
      #       Cdf <- Cdf*apply(Cdf.MixtMarg[,WhichH1Lower,drop=FALSE],1,prod)
      #     }
      #     #Deal with H1 and higher
      #     WhichH1Higher <- which((h.cfg==1) & (h.lower==0))
      #     if(length(WhichH1Higher)>0){
      #       Cdf <- Cdf*apply(1-Cdf.MixtMarg[,WhichH1Higher,drop=FALSE],1,prod)
      #     }
      #     return(Cdf)
      #   })
      #   return(rowSums(Tmp))
      # })
      
      ## Next two lines are equivalent
      Num <- rank(P.H1)/length(P.H1) - Cdf.Hconfig.H1%*%Prior.Hconfig.H1
      #Num <- Tau2Cdf(P.H1,1/length(P.H1),P.H1) - Cdf.Hconfig.H1%*%Prior.Hconfig.H1
      Den <-1-sum(Prior.Hconfig.H1)
      PvalUnion =  Num/Den 
      
    }
      
  }
  
  return(PvalUnion)
}



########################   ComputePValue_2   ####################

ComputePValue_2 <- function(PvalMat,Tau.MixtMarg,Hconfig,Hconfig.H1,Prior.Hconfig,P.H1){
  
  ## Basic stuff
  Q <- ncol(PvalMat)
  if(Q !=2){
    stop("ComputePvalue_2 only works when Q=2 (Q=Nb of series of tests")
  }
  
  ## Infer cdf of each Hconfig
  Cdf.Hconfig <- sapply(Hconfig, function(h){
    if(sum(h) > 0){
      P.H1.h = apply(PvalMat[, which(h==1), drop=FALSE], 1, max)
      tau.h = rep(1, n)
      sapply(which(h==1), function(q){
        tau.h <<- tau.h * Tau.MixtMarg[,q]
      })
      H1comp <- Tau2Cdf(P.H1.h, tau.h, P.H1)
      Res<- H1comp * (P.H1^sum(h==0))
    } else {
      Res <- P.H1^Q
    }
  })
  PvalUnion = Cdf.Hconfig[,-Hconfig.H1]%*%Prior.Hconfig[-Hconfig.H1] / sum(Prior.Hconfig[-Hconfig.H1])
  return(PvalUnion)
}

