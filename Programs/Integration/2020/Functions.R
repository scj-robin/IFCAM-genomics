

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

GetHinfo <- function(Q,AtLeast,Consecutive=FALSE){

  ## Build H configurations
  Hconfig <- as.matrix(expand.grid(lapply(1:Q, function(q) 0:1)))
  Hconfig <- split(Hconfig, seq(2^Q))
  
  ## Find the ones that match H1
  if (!Consecutive){
    MatchingH1 <- sapply(Hconfig, function(h){sum(h)>=AtLeast})
    Hconfig.H1 <- which(MatchingH1)
    names(Hconfig.H1) <- NULL
  } else {
    Consec <- paste(rep(1,AtLeast),collapse='')
    Hconcat <- sapply(Hconfig, function(hh){paste(hh,collapse='')})
    Hconfig.H1 <- grep(pattern = Consec,x = Hconcat)
  }
  
  ## Collect results
  return(list(Hconfig=Hconfig,Hconfig.H1=Hconfig.H1))
  
}  

GetHinfoEqual <- function(Q,Equal,Consecutive=FALSE){
  
  ## Build H configurations
  Hconfig <- as.matrix(expand.grid(lapply(1:Q, function(q) 0:1)))
  Hconfig <- split(Hconfig, seq(2^Q))
  
  ## Find the ones that match H1
  if (!Consecutive){
    MatchingH1 <- sapply(Hconfig, function(h){sum(h)==Equal})
    Hconfig.H1 <- which(MatchingH1)
    names(Hconfig.H1) <- NULL
  } else {
    Consec <- paste(rep(1,Equal),collapse='')
    ConsecP1 <- paste(rep(1,Equal+1),collapse='')
    Hconcat <- sapply(Hconfig, function(hh){paste(hh,collapse='')})
    Hconfig.H1 <- intersect(grep(pattern = Consec,x = Hconcat),which(sapply(Hconfig, function(h){sum(h)==Equal})))
  }
  
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
   ## Pval <- pValMat[, 1]; NbKnot=1e5; tol = 1e-5

  ## Transform pvalues into N(0,1) quantiles
  n = length(Pval)
  X = -qnorm(Pval)

  ## Get a p0 estimate
  if(is.null(p0)){
    p0 = min(2*sum(X < 0)/n,1-1/n);     
  }
  p1 = 1 - p0
  
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
    ## Dirty job 1: get rid of the f1 mass on the left
    tauNew[Knots< -3] <- 0
    diff = max(abs(tau - tauNew))
    tau = tauNew
  }
  if(plotting){
    Hist.fig <- hist(X, freq=TRUE, breaks=sqrt(n), main='', border=8, 
                     xlab="Q-transformed pvalues", ylab="Densities")
    bin.width <- mean(diff(Hist.fig$breaks))
    lines(Knots, n*bin.width*p0*phi, type='l', col=4, lwd=2); 
    lines(Knots, n*bin.width*p1*f1, col=2,lwd=2); 
    lines(Knots, n*bin.width*(p0*phi+p1*f1), lwd=2)
    legend("topright", legend=c("H0 dist", "H1 dist","Mixture Dist"),
           col=c("blue","red", "black"), lty=c(1,1,2), cex=0.8)  
  }
  
  
  ## Dirty job 1: get rid of the f1 mass on the left
  #plot(Knots, f1, col=2,lwd=2); 
  # MassToReaffect <- sum(weights[which(Knots<= -3)])
  # TotalMass <- sum(weights)
  # weights[which(Knots<= -3)] <- 0
  # weights[which(Knots> -3)] <- weights[which(Knots> -3)]*(TotalMass/(TotalMass-MassToReaffect))

  ## Now get the f1 estimate
  KDE = kde(x=Knots, w=weights, eval.points=X)
  f1 = KDE$estimate
  
  ## Dirty job 2: get rid of numeric problems
  f1[f1<0] <- 1e-30
  tau = p1*f1 / (p1*f1 + p0*dnorm(X))
  
  return(list(p0=p0,tau=tau,f1=f1))
}  

PlotPvalHistDist <- function(Pval, p0, f1){
   
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


        #### Full mixture model procedure


MixtModProcedure <- function(pValMat,Hconfig){
  
  n <- nrow(pValMat)
  Q <- ncol(pValMat)
  
  #### Step 1: Marginal density estimation
  
  ## Get p0 estimates
  p0 <- rep(0, Q)
  for (q in 1:Q){
    p0[q] = min(2*sum(pValMat[,q] > 0.5)/n,1-1/n);     
  }
  SomeH1 <-  which(p0<1)
  NoH1 <- which(p0==1-1/n); 
  if(length(NoH1)==1){
    message(paste("Pvalue serie",NoH1, "may have very few H1 (or a weird distribution)"))
  }
  if(length(NoH1)>1){
    message(paste("Pvalue series",paste(NoH1,collapse=' '), "may have very few H1 (or a weird distribution)"))
  }
  
  ## Fit a 2-component mixture to each test serie using kerFdr
  f1Mat <- tauKerMat <- matrix(1, n, Q); 
  for(q in SomeH1){
    ker <- FastKerFdr(pValMat[, q], p0=p0[q], plotting=TRUE)
    f1Mat[, q] <- ker$f1
    tauKerMat[, q] <- ker$tau
  }
  f0Mat <- matrix(dnorm(-qnorm(pValMat)),ncol=Q)
  
  #### Step 2: transform marginal densities into config densities
  
  Logf0Mat <- log(f0Mat); 
  Logf1Mat <- log(f1Mat); 
  f.Hconfig <- sapply(Hconfig, function(h){
    f <- rep(0,nrow(Logf0Mat))
    if (length(which(h==1)) > 0){f <- f + rowSums(Logf1Mat[, which(h==1), drop=FALSE])}
    if (length(which(h==0)) > 0){f <- f + rowSums(Logf0Mat[, which(h==0), drop=FALSE])}
    return(exp(f))
  })
  
  #### Step 3: Infer prior estimation using an EM procedure
  
  ## Initialization: simple product of marginal priors estimator
  NewPrior <- sapply(1:length(Hconfig), function(c){
    prod(p0[which(Hconfig[[c]]==0)]) * prod(1-p0[which(Hconfig[[c]]==1)])
  })
  PriorsAt0 <- which(NewPrior==0)
  
  ## EM calibration
  NotOK <- TRUE
  Precision <- 1e-6
  NoLowerThan <- 1e-7
  while(NotOK){
    
    ## E step
    Tau <- f.Hconfig*(tcrossprod(rep(1:n),NewPrior))
    Tau <- Tau/rowSums(Tau)
    
    ## M step 
    OldPrior <- NewPrior
    NewPrior <- colMeans(Tau)
    if(length(PriorsAt0)==0){
      NewPrior[NewPrior<NoLowerThan] <- NoLowerThan
    } else {
      NoLowerCoord <- setdiff(which(NewPrior<NoLowerThan),PriorsAt0)
      if(length(NoLowerCoord)>0){
        NewPrior[NoLowerCoord] <- NoLowerThan  
      }
    }
    NewPrior <- NewPrior/sum(NewPrior)
    NotOK <- max((OldPrior-NewPrior)^2) > Precision
    
  }
  priorHconfigEM <- NewPrior

  #### Step 4: Posterior computation
  
  posterior <- f.Hconfig*(tcrossprod(rep(1:n),priorHconfigEM))
  posterior <- posterior/rowSums(posterior)
  
  #### Last but not least: output results
  Res <- list(prior=priorHconfigEM, posterior=posterior, Tau=Tau, tauKer=tauKerMat)
  return(Res)
}


PerformMultipleTestingEM <- function(posterior,Hconfig.H1,Alpha=0.05){
  
  n <- nrow(posterior)
  localFDR <- 1-rowSums(posterior[,Hconfig.H1,drop=FALSE])
  Order <- order(localFDR)
  FDR <- cumsum(localFDR[Order])/(1:n)
  NbReject <- max(which(FDR<=Alpha))
  Rejection <- rep(0,n)
  if (NbReject>0){
    Rejection[Order[1:NbReject]] <- 1
  }
  return(list(Rejection=Rejection,lFDR=localFDR))
  
}

