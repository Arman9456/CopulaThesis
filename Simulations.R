
#SIMUL1 IS FINISHED, STILL NEED TO FINISH SIMUL2!
# WHY CANT I PUT NEW STUFF IN SIMUL1, ITS LIKE THE CODE ONLY RUNS WITH THE OLD VALUES
# ->>> PRESS SOURCE, THAT WILL FIX THE ABOVE

Simul1 <- function(T,N,rho){
  
  #intceptvecOLS <- rep(0,N) # intercept coefficients to fill up
  alphavecOLS <- rep(0,N) # endogenous variable coefficients to fill up
  betavecOLS <- rep(0,N) # exogenous variable coefficients to fill up
  
  #intceptvecGC <- rep(0,N) # intercept coefficients to fill up
  alphavecGC <- rep(0,N) # endogenous variable coefficients to fill up
  betavecGC <- rep(0,N) # exogenous variable coefficients to fill up
  
  #intceptvecIV <- rep(0,N) # intercept coefficients to fill up
  alphavecIV<- rep(0,N) # endogenous variable coefficients to fill up
  betavecIV <- rep(0,N) # exogenous variable coefficients to fill up
  
  for(i in 1:N){
    # Generate bivariate normally generated error terms between error term of original equation and error term of equation for endogenous regressor P
    # S.t. the error terms are correlated and so P and the error term of the original equaton are correlated
    bivariate_errors = mvrnorm(n=T,mu=c(0, 0),Sigma=matrix(c(1, rho,rho, 1), ncol=2))
    # IF YOU PUT HIGH VARIANCE (SIGMA = C(3,rho,rho,2), THEN GC PERFORMS WELL)
    
    
    # Construct endogenous variable that has instrument Z which is perfectly correlated with P but uncorrelated to error term of original equation
    # Generate endogenous regressor P with instrument Z~gamma(1,1) s.t. P is nonnormal just like in GC assumption
    Z = rgamma(n=T, shape = 1, rate = 0.5) # could also use uniform distribution just like in Park & Gupta (or Gamma like in other thesis)
    P = Z + bivariate_errors[,2] # should I add intercept to this equation?
    
    #cov(bivariate_errors[,1], P) #check whether covariance between P and error term of original equation equals rho (TRUE)
    # Generate X_t's 
    X = rnorm(n=T,mean=0,sd=1)
    
    # Coefficient values of original equations
    #intcpt = rep(0.5,T)
    alpha = -1
    beta = 2
    
    Y = alpha * P + beta*X + bivariate_errors[,1]
    
    
    # OLS estimates
    #intceptvecOLS[i] = OLS1(Y,X,P,T)[1]
    alphavecOLS[i]= OLS1(Y,X,P,T)[1]
    betavecOLS[i] = OLS1(Y,X,P,T)[2]

    
    
    #GC estimates
    
    #intceptvecGC[i] = GC1(Y,X,P)[[1]]
    fit = GC1(Y,X,P)
    alphavecGC[i] = fit[[2]]
    betavecGC[i] = fit[[1]]

    
    #IV estimates
    #intceptvecIV[i] = IV1(Y,X,P,Z)[[1]]
    alphavecIV[i] = IV1(Y,X,P,Z)[[2]]
    betavecIV[i] = IV1(Y,X,P,Z)[[1]]
    
    
   
  }
  BiasOLSEstimates = c(mean(alphavecOLS)-alpha,mean(betavecOLS)-beta)
  BiasGCstimates = c(mean(alphavecGC)-alpha,mean(betavecGC)-beta)
  BiasIVEstimates = c(mean(alphavecIV)-alpha,mean(betavecIV)-beta)
  
  hist(alphavecOLS-alpha,main = expression(paste("Distribution of Bias",(delta))),xlab="Bias OLS")
  hist(alphavecGC-alpha,main = expression(paste("Distribution of Bias",(delta))),xlab="Bias GC")
  hist(alphavecIV-alpha,main = expression(paste("Distribution of Bias",(delta))),xlab="Bias IV")
  
  hist(betavecOLS-beta,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias OLS")
  hist(betavecGC-beta,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias GC")
  hist(betavecIV-beta,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias IV")
  
  
  MatrixOfMeanEstimates = matrix(c(BiasOLSEstimates,BiasGCstimates,BiasIVEstimates),ncol=3)
  rownames(MatrixOfMeanEstimates) <- c("alpha","beta")
  colnames(MatrixOfMeanEstimates) <- c("Bias OLS", "Bias GC","Bias IV")
  return(MatrixOfMeanEstimates)
}




Simul2 <- function(T,N,rho_12,rho_13,rho_23){
  
  #intceptvecOLS <- rep(0,N) # intercept coefficients to fill up
  alpha1vecOLS <- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecOLS <- rep(0,N)
  betavecOLS <- rep(0,N) # exogenous variable coefficients to fill up
  
  #intceptvecGC <- rep(0,N) # intercept coefficients to fill up
  alpha1vecGC <- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecGC <- rep(0,N)
  betavecGC <- rep(0,N) # exogenous variable coefficients to fill up
  
  #intceptvecIV <- rep(0,N) # intercept coefficients to fill up
  alpha1vecIV<- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecIV<- rep(0,N)
  betavecIV <- rep(0,N) 
  
  #intceptvecIVCOP <- rep(0,N) # intercept coefficients to fill up
  alpha1vecIVCOP <- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecIVCOP <- rep(0,N)
  betavecIVCOP <- rep(0,N) # exogenous variable coefficients to fill up
  
  alpha1vecIVCOPREV <- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecIVCOPREV <- rep(0,N)
  betavecIVCOPREV <- rep(0,N) # exogenous variable coefficients to fill up
  
  
  
  for(i in 1:N){
  
  #Create errors that are correlated across all equations (the two endogenous regressor equations and the original equation) just like in Park & Gupta p.32
  # Because if both endogeneous regressors are correlated with the error term from the original equation, they must be correlated to eachother
  
  
  trivariate_errors = mvrnorm(n=T,mu=c(0, 0,0),Sigma=matrix(c(1, rho_12,rho_13,
                                                              rho_12,1,rho_23,
                                                              rho_13,rho_23,1), ncol=3,byrow=TRUE))
  
  X2 = rnorm(n=T,mean=0,sd=1) # create an exogenous variable, since only 2 endogenous variables gave weird dresults
  beta2= 0.5
  
  alpha1 = -1
  alpha2 = 1
  
  Z1 = rgamma(n=T, shape = 1, rate = 1) # instrument for endogenous regressor P1 s.t. P1 is not normally distributed
  Z2 = rgamma(n=T,shape=1,0.5)  # could also use uniform distribution just like in Park & Gupta
  
  P2 = Z2 + trivariate_errors[,3] 
  P1 = Z1 + trivariate_errors[,2]
  
  #intcpt2 = rep(0.3,T)
  # Construct sample from DGP
  Y2= (alpha1 * P1) + (alpha2*P2) + beta2*X2 + trivariate_errors[,1]
  

  # Iv for one of them (P1) and GC for other one (P2)
  
  # first estimate P1 equation by OLS as first stage in 2SLS, I don't think I should include P2 as a regressor in that equation, because P2 
  # is correlated with the original error term. (If P2 were to be an exogenous regressor, I would have to include it as a first stage instrument as well)
  #https://timeseriesreasoning.com/contents/two-stage-least-squares-estimation/
  
  
  #IVCOPMIX1(Y2,X2,P1,P2,Z1,Z2) # Damn, pretty decent estimates
  
  #intceptvecOLS[i] = OLS2(Y2,X2,P1,P2,T)[1]
  alpha1vecOLS[i]= OLS2(Y2,X2,P1,P2,T)[1]
  alpha2vecOLS[i]= OLS2(Y2,X2,P1,P2,T)[2]
  betavecOLS[i] = OLS2(Y2,X2,P1,P2,T)[3]
  
  
  
  #GC estimates
  #intceptvecGC[i] = GC2(Y2,X2,P1,P2)[[1]]
  fitGC = GC2(Y2,X2,P1,P2)
  alpha1vecGC[i] = fitGC[[1]]
  alpha2vecGC[i] = fitGC[[2]]
  betavecGC[i] = fitGC[[3]]
  

  
  #IV estimates
  #intceptvecIV[i] =   IV2(Y2,X2,P1,P2,Z1,Z2)[[1]]
  alpha1vecIV[i] =   IV2(Y2,X2,P1,P2,Z1,Z2)[[1]]
  alpha2vecIV[i] =   IV2(Y2,X2,P1,P2,Z1,Z2)[[2]]
  betavecIV[i] =   IV2(Y2,X2,P1,P2,Z1,Z2)[[3]]
  
 
  
  #intceptvecIVCOP[i] =   IVCOPMIX1(Y2,X2,P1,P2,Z1,Z2)[[1]]
  fitIVCOP = IVCOPMIX1(Y2,X2,P1,P2,Z1,Z2)
  alpha1vecIVCOP[i] =   fitIVCOP[[1]]
  alpha2vecIVCOP[i] =   fitIVCOP[[2]]
  betavecIVCOP[i] =   fitIVCOP[[3]]
  
  fitIVCOPREV = IVCOPMIX1REVERSE(Y2,X2,P1,P2,Z1,Z2)
  alpha1vecIVCOPREV[i] =   fitIVCOPREV[[1]]
  alpha2vecIVCOPREV[i] =   fitIVCOPREV[[2]]
  betavecIVCOPREV[i] =   fitIVCOPREV[[3]]
  

  }
  
  MeanOLSEstimates = c(mean(alpha1vecOLS)-alpha1,mean(alpha2vecOLS)-alpha2,mean(betavecOLS)-beta2)
  MeanGCstimates = c(mean(alpha1vecGC)-alpha1,mean(alpha2vecGC)-alpha2,mean(betavecGC)-beta2)
  MeanIVEstimates = c(mean(alpha1vecIV)-alpha1,mean(alpha2vecIV)-alpha2,mean(betavecIV)-beta2)
  MeanIVCOPEstimates = c(mean(alpha1vecIVCOP)-alpha1,mean(alpha2vecIVCOP)-alpha2,mean(betavecIVCOP)-beta2)
  MeanIVCOPREVEstimates = c(mean(alpha1vecIVCOPREV)-alpha1,mean(alpha2vecIVCOPREV)-alpha2,mean(betavecIVCOPREV)-beta2)
  
  
  hist(alpha1vecOLS-alpha1,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias OLS")
  hist(alpha1vecGC-alpha1,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias GC")
  hist(alpha1vecIV-alpha1,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias IV")
  hist(alpha1vecIVCOP-alpha1,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias IVCOP")
  hist(alpha1vecIVCOPREV-alpha1,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias IVCOP Reverse")
  
  hist(alpha2vecOLS-alpha2,main = expression(paste("Distribution of Bias",(delta[2]))),xlab="Bias OLS")
  hist(alpha2vecGC-alpha2,main = expression(paste("Distribution of Bias",(delta[2]))),xlab="Bias GC")
  hist(alpha2vecIV-alpha2,main = expression(paste("Distribution of Bias",(delta[2]))),xlab="Bias IV")
  hist(alpha2vecIVCOP-alpha2,main = expression(paste("Distribution of Bias",(delta[2]))),xlab="Bias IVCOP")
  hist(alpha2vecIVCOPREV-alpha2,main = expression(paste("Distribution of Bias",(delta[2]))),xlab="Bias IVCOP Reverse")
  
  hist(betavecOLS-beta2,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias OLS")
  hist(betavecGC-beta2,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias GC")
  hist(betavecIV-beta2,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias IV")
  hist(betavecIVCOP-beta2,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias IVCOP")
  hist(betavecIVCOPREV-beta2,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias IVCOP Reverse")
  
  
MatrixOfMeanEstimates = matrix(c(MeanOLSEstimates,MeanGCstimates,MeanIVEstimates,MeanIVCOPEstimates,MeanIVCOPREVEstimates),ncol=5)
rownames(MatrixOfMeanEstimates) <- c("alpha1","alpha2","beta")
colnames(MatrixOfMeanEstimates) <- c("Bias OLS", "Bias GC","Bias IV","Bias IVCOP", "Bias IVCOP Reverse")
return(MatrixOfMeanEstimates)
  
  
}

# IV on 2 endo, GC on one endo (1exo and 3 endo's in total)
Simul3 <- function(T,N,rho_12,rho_13,rho_23,rho_14,rho_24,rho_34){
  
  #intceptvecOLS <- rep(0,N) # intercept coefficients to fill up
  alpha1vecOLS <- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecOLS <- rep(0,N)
  alpha3vecOLS <- rep(0,N)
  betavecOLS <- rep(0,N) # exogenous variable coefficients to fill up
  
  #intceptvecGC <- rep(0,N) # intercept coefficients to fill up
  alpha1vecGC <- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecGC <- rep(0,N)
  alpha3vecGC <- rep(0,N)
  betavecGC <- rep(0,N) # exogenous variable coefficients to fill up
  
  #intceptvecIV <- rep(0,N) # intercept coefficients to fill up
  alpha1vecIV<- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecIV<- rep(0,N)
  alpha3vecIV <- rep(0,N)
  betavecIV <- rep(0,N) 
  
  #intceptvecIVCOP <- rep(0,N) # intercept coefficients to fill up
  alpha1vecIVCOP <- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecIVCOP <- rep(0,N)
  alpha3vecIVCOP <- rep(0,N)
  betavecIVCOP <- rep(0,N) # exogenous variable coefficients to fill up
  
  
  for(i in 1:N){
    
    #Create errors that are correlated across all equations (the two endogenous regressor equations and the original equation) just like in Park & Gupta p.32
    # Because if both endogeneous regressors are correlated with the error term from the original equation, they must be correlated to eachother
    
    
    multivariate_errors = mvrnorm(n=T,mu=c(0,0,0,0),Sigma=matrix(c(1, rho_12,rho_13, rho_14,
                                                                rho_12,1,rho_23,rho_24,
                                                                rho_13,rho_23,1,rho_34, 
                                                                rho_14,rho_24,rho_34,1),ncol=4,byrow=TRUE))
    
    X2 = rnorm(n=T,mean=0,sd=1) # create an exogenous variable, since only 2 endogenous variables gave weird dresults
    beta2= 0.5
    
    alpha1 = -1
    alpha2 = 1
    alpha3 = 0.7
    
    Z1 = rgamma(n=T, shape = 1, rate = 1) # instrument for endogenous regressor P1 s.t. P1 is not normally distributed
    Z2 = rgamma(n=T, shape = 1, rate = 0.5)  # If I put the rate to 2 it gives bad solutions
    Z3 = rgamma(n=T,shape=0.5,rate = 0.5)
    
    P3 = Z3 + multivariate_errors[,4]
    P2 = Z2 + multivariate_errors[,3] 
    P1 = Z1 + multivariate_errors[,2]
    
    #intcpt2 = rep(0.3,T)
    # Construct sample from DGP
    Y2= (alpha1 * P1) + (alpha2*P2) + alpha3*P3+ beta2*X2 + multivariate_errors[,1]
    
    
    # IV for two of them (P1,P2) and GC for other one (P3)
    
    # first estimate P1,P2 equation by OLS as first stage in 2SLS , because P2 
    # is correlated with the original error term. (If P2 were to be an exogenous regressor, I would have to include it as a first stage instrument as well)
    #https://timeseriesreasoning.com/contents/two-stage-least-squares-estimation/
    
    
    #IVCOPMIX1(Y2,X2,P1,P2,Z1,Z2) # Damn, pretty decent estimates
    
    #intceptvecOLS[i] = OLS2(Y2,X2,P1,P2,T)[1]
    alpha1vecOLS[i]= OLS3(Y2,X2,P1,P2,P3,T)[1]
    alpha2vecOLS[i]= OLS3(Y2,X2,P1,P2,P3,T)[2]
    alpha3vecOLS[i]= OLS3(Y2,X2,P1,P2,P3,T)[3]
    betavecOLS[i] = OLS3(Y2,X2,P1,P2,P3,T)[4]
    
    
    
    #GC estimates
    #intceptvecGC[i] = GC2(Y2,X2,P1,P2)[[1]]
    fitGC = GC3(Y2,X2,P1,P2,P3)
    alpha1vecGC[i] = fitGC[[1]]
    alpha2vecGC[i] = fitGC[[2]]
    alpha3vecGC[i] = fitGC[[3]]
    betavecGC[i] = fitGC[[4]]
    
    
    
    #IV estimates
    #intceptvecIV[i] =   IV2(Y2,X2,P1,P2,Z1,Z2)[[1]]
    alpha1vecIV[i] =   IV3(Y2,X2,P1,P2,P3,Z1,Z2,Z3)[[1]]
    alpha2vecIV[i] =   IV3(Y2,X2,P1,P2,P3,Z1,Z2,Z3)[[2]]
    alpha3vecIV[i] = IV3(Y2,X2,P1,P2,P3,Z1,Z2,Z3)[[3]]
    betavecIV[i] =   IV3(Y2,X2,P1,P2,P3,Z1,Z2,Z3)[[4]]
    
   
    
    #intceptvecIVCOP[i] =   IVCOPMIX1(Y2,X2,P1,P2,Z1,Z2)[[1]]
    fitCOPIV = IVCOPMIX2(Y2,X2,P1,P2,P3,Z1,Z2,Z3)
    alpha1vecIVCOP[i] =   fitCOPIV[[1]]
    alpha2vecIVCOP[i] =   fitCOPIV[[2]]
    alpha3vecIVCOP[i] = fitCOPIV[[3]]
    betavecIVCOP[i] =   fitCOPIV[[4]]
    
   
  }
  
  MeanOLSEstimates = c(mean(alpha1vecOLS)-alpha1,mean(alpha2vecOLS)-alpha2,mean(alpha3vecOLS)-alpha3,mean(betavecOLS)-beta2)
  MeanGCstimates = c(mean(alpha1vecGC)-alpha1,mean(alpha2vecGC)-alpha2,mean(alpha3vecGC)-alpha3,mean(betavecGC)-beta2)
  MeanIVEstimates = c(mean(alpha1vecIV)-alpha1,mean(alpha2vecIV)-alpha2,mean(alpha3vecIV)-alpha3,mean(betavecIV)-beta2)
  MeanIVCOPEstimates = c(mean(alpha1vecIVCOP)-alpha1,mean(alpha2vecIVCOP)-alpha2,mean(alpha3vecIVCOP)-alpha3,mean(betavecIVCOP)-beta2)
  
  hist(alpha1vecOLS-alpha1,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias OLS")
  hist(alpha1vecGC-alpha1,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias GC")
  hist(alpha1vecIV-alpha1,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias IV")
  hist(alpha1vecIVCOP-alpha1,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias IVCOP")
  
  hist(alpha2vecOLS-alpha2,main = expression(paste("Distribution of Bias",(delta[2]))),xlab="Bias OLS")
  hist(alpha2vecGC-alpha2,main = expression(paste("Distribution of Bias",(delta[2]))),xlab="Bias GC")
  hist(alpha2vecIV-alpha2,main = expression(paste("Distribution of Bias",(delta[2]))),xlab="Bias IV")
  hist(alpha2vecIVCOP-alpha2,main = expression(paste("Distribution of Bias",(delta[2]))),xlab="Bias IVCOP")
  
  
  hist(alpha3vecOLS-alpha3,main = expression(paste("Distribution of Bias",(delta[3]))),xlab="Bias OLS")
  hist(alpha3vecGC-alpha3,main = expression(paste("Distribution of Bias",(delta[3]))),xlab="Bias GC")
  hist(alpha3vecIV-alpha3,main = expression(paste("Distribution of Bias",(delta[3]))),xlab="Bias IV")
  hist(alpha3vecIVCOP-alpha3,main = expression(paste("Distribution of Bias",(delta[3]))),xlab="Bias IVCOP")
  
  
  hist(betavecOLS-beta2,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias OLS")
  hist(betavecGC-beta2,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias GC")
  hist(betavecIV-beta2,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias IV")
  hist(betavecIVCOP-beta2,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias IVCOP")
  
  MatrixOfMeanEstimates = matrix(c(MeanOLSEstimates,MeanGCstimates,MeanIVEstimates,MeanIVCOPEstimates),ncol=4)
  rownames(MatrixOfMeanEstimates) <- c("alpha1","alpha2","alpha3","beta")
  colnames(MatrixOfMeanEstimates) <- c("Bias OLS", "Bias GC","Bias IV","Bias IVCOP")
  return(MatrixOfMeanEstimates)
  
  
}

Simul4 <- function(T,N,rho_12,rho_13,rho_23,rho_14,rho_24,rho_34){
  
  #intceptvecOLS <- rep(0,N) # intercept coefficients to fill up
  alpha1vecOLS <- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecOLS <- rep(0,N)
  alpha3vecOLS <- rep(0,N)
  betavecOLS <- rep(0,N) # exogenous variable coefficients to fill up
  
  #intceptvecGC <- rep(0,N) # intercept coefficients to fill up
  alpha1vecGC <- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecGC <- rep(0,N)
  alpha3vecGC <- rep(0,N)
  betavecGC <- rep(0,N) # exogenous variable coefficients to fill up
  
  #intceptvecIV <- rep(0,N) # intercept coefficients to fill up
  alpha1vecIV<- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecIV<- rep(0,N)
  alpha3vecIV <- rep(0,N)
  betavecIV <- rep(0,N) 
  
  #intceptvecIVCOP <- rep(0,N) # intercept coefficients to fill up
  alpha1vecIVCOP <- rep(0,N) # endogenous variable coefficients to fill up
  alpha2vecIVCOP <- rep(0,N)
  alpha3vecIVCOP <- rep(0,N)
  betavecIVCOP <- rep(0,N) # exogenous variable coefficients to fill up
  
  
  for(i in 1:N){
    
    #Create errors that are correlated across all equations (the two endogenous regressor equations and the original equation) just like in Park & Gupta p.32
    # Because if both endogeneous regressors are correlated with the error term from the original equation, they must be correlated to eachother
    
    
    multivariate_errors = mvrnorm(n=T,mu=c(0,0,0,0),Sigma=matrix(c(1, rho_12,rho_13, rho_14,
                                                                   rho_12,1,rho_23,rho_24,
                                                                   rho_13,rho_23,1,rho_34, 
                                                                   rho_14,rho_24,rho_34,1),ncol=4,byrow=TRUE))
    
    X2 = rnorm(n=T,mean=0,sd=1) # create an exogenous variable, since only 2 endogenous variables gave weird dresults
    beta2= 0.5
    
    alpha1 = -1
    alpha2 = 1
    alpha3 = 0.7
    
    Z1 = rgamma(n=T, shape = 1, rate = 1) # instrument for endogenous regressor P1 s.t. P1 is not normally distributed
    Z2 = rgamma(n=T, shape = 1, rate = 0.5)  # If I put the rate to 2 it gives bad solutions
    Z3 = rgamma(n=T,shape=0.5,rate = 0.5)
    
    P3 = Z3 + multivariate_errors[,4]
    P2 = Z2 + multivariate_errors[,3] 
    P1 = Z1 + multivariate_errors[,2]
    
    #intcpt2 = rep(0.3,T)
    # Construct sample from DGP
    Y2= (alpha1 * P1) + (alpha2*P2) + alpha3*P3+ beta2*X2 + multivariate_errors[,1]
    
    
    # IV for one of them P1 and GC for other two P2,P3
    
    # first estimate P1, equation by OLS as first stage in 2SLS , 
    #https://timeseriesreasoning.com/contents/two-stage-least-squares-estimation/
    
    
  
    
    #intceptvecOLS[i] = OLS2(Y2,X2,P1,P2,T)[1]
    alpha1vecOLS[i]= OLS3(Y2,X2,P1,P2,P3,T)[1]
    alpha2vecOLS[i]= OLS3(Y2,X2,P1,P2,P3,T)[2]
    alpha3vecOLS[i]= OLS3(Y2,X2,P1,P2,P3,T)[3]
    betavecOLS[i] = OLS3(Y2,X2,P1,P2,P3,T)[4]
    
    
    
    #GC estimates
    #intceptvecGC[i] = GC2(Y2,X2,P1,P2)[[1]]
    fitGC = GC3(Y2,X2,P1,P2,P3)
    alpha1vecGC[i] = fitGC[[1]]
    alpha2vecGC[i] = fitGC[[2]]
    alpha3vecGC[i] = fitGC[[3]]
    betavecGC[i] = fitGC[[4]]
    
    
    
    #IV estimates
    #intceptvecIV[i] =   IV2(Y2,X2,P1,P2,Z1,Z2)[[1]]
    alpha1vecIV[i] =   IV3(Y2,X2,P1,P2,P3,Z1,Z2,Z3)[[1]]
    alpha2vecIV[i] =   IV3(Y2,X2,P1,P2,P3,Z1,Z2,Z3)[[2]]
    alpha3vecIV[i] = IV3(Y2,X2,P1,P2,P3,Z1,Z2,Z3)[[3]]
    betavecIV[i] =   IV3(Y2,X2,P1,P2,P3,Z1,Z2,Z3)[[4]]
    
    
    
    #intceptvecIVCOP[i] =   IVCOPMIX1(Y2,X2,P1,P2,Z1,Z2)[[1]]
    fitCOPIV = IVCOPMIX3(Y2,X2,P1,P2,P3,Z1,Z2,Z3)
    alpha1vecIVCOP[i] =   fitCOPIV[[1]]
    alpha2vecIVCOP[i] =   fitCOPIV[[2]]
    alpha3vecIVCOP[i] = fitCOPIV[[3]]
    betavecIVCOP[i] =   fitCOPIV[[4]]
    
    
  }
  
  MeanOLSEstimates = c(mean(alpha1vecOLS)-alpha1,mean(alpha2vecOLS)-alpha2,mean(alpha3vecOLS)-alpha3,mean(betavecOLS)-beta2)
  MeanGCstimates = c(mean(alpha1vecGC)-alpha1,mean(alpha2vecGC)-alpha2,mean(alpha3vecGC)-alpha3,mean(betavecGC)-beta2)
  MeanIVEstimates = c(mean(alpha1vecIV)-alpha1,mean(alpha2vecIV)-alpha2,mean(alpha3vecIV)-alpha3,mean(betavecIV)-beta2)
  MeanIVCOPEstimates = c(mean(alpha1vecIVCOP)-alpha1,mean(alpha2vecIVCOP)-alpha2,mean(alpha3vecIVCOP)-alpha3,mean(betavecIVCOP)-beta2)
  

  hist(alpha1vecIVCOP-alpha1,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias IVCOP")
  hist(alpha2vecIVCOP-alpha2,main = expression(paste("Distribution of Bias",(delta[2]))),xlab="Bias IVCOP")
  hist(alpha3vecIVCOP-alpha3,main = expression(paste("Distribution of Bias",(delta[3]))),xlab="Bias IVCOP")
  hist(betavecIVCOP-beta2,main = expression(paste("Distribution of Bias",(beta))),xlab="Bias IVCOP")
  
  MatrixOfMeanEstimates = matrix(c(MeanOLSEstimates,MeanGCstimates,MeanIVEstimates,MeanIVCOPEstimates),ncol=4)
  rownames(MatrixOfMeanEstimates) <- c("alpha1","alpha2","alpha3","beta")
  colnames(MatrixOfMeanEstimates) <- c("Bias OLS", "Bias GC","Bias IV","Bias IVCOP")
  return(MatrixOfMeanEstimates)
  
  
}


