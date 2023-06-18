# Author: Arman Behbood
# Code for all the functions that will be used to estimate OLS,GC and IVCOP without intercept

OLS1 <- function(Y,X,P,T){
  
 #VectorOf1 = rep(1,T)
 Xmatrix = matrix(c(P,X),ncol=2)  # merge intcpt, P and X into one matrix of regressors PX, because OLS doesn't distinguish between exogenous and endogenous regressors
 theta_OLS = inv(t(Xmatrix)%*%Xmatrix) %*% t(Xmatrix) %*% Y # alternative: lm(Y~P+X)
 
 return(theta_OLS)
}

GC1 <- function(Y,X,P){
  dataMatrix = data.frame(Y,X,P)
  theta_cop = copulaCorrection(Y~X+P-1|continuous(P),data= dataMatrix,num.boots=2)
  return(theta_cop$coefficients)
}

IV1 <- function(Y,X,P,Z){
  dataMatrix = data.frame(Y,X,P,Z)
  theta_IV = ivreg(Y ~ X + P-1| X + Z,data = dataMatrix) # X and Z are used as instruments (X is also exogenous regressor which should be included as instrument)
  return(theta_IV$coefficients)
}


OLS2 <- function(Y,X,P1,P2,T){
  #VectorOf1 = rep(1,T)
  Xmatrix = matrix(c(P1,P2,X),ncol=3)
  theta_ols = inv(t(Xmatrix)%*%(Xmatrix)) %*% (t(Xmatrix) %*% Y)
  return(theta_ols)
}

GC2 <- function(Y,X,P1,P2){
  dataMatrix = data.frame(Y,P1,P2,X)
  theta_cop = copulaCorrection(Y~P1+P2+X-1|continuous(P1,P2), data= dataMatrix,num.boots=2)
  return(theta_cop$coefficients)
}

IV2 <- function(Y,X,P1,P2,Z1,Z2){
  dataMatrix = data.frame(Y,P1,P2,X)
  theta_IV = ivreg(Y~P1+P2+X-1|Z1+Z2+X,data=dataMatrix)# Z1 and Z2 are used as instruments and X2 as exogenous variable instrument, since you regress endogenous regressors on exogenous variables as well
  return(theta_IV$coefficients)
}

IVCOPMIX1 <- function(Y,X,P1,P2,Z1,Z2){
  
firstStageRegressors = matrix(c(Z1,Z2,X),ncol=3) # regress P1 on all instruments as well as regressor X2 which is also used as instrument
# sO P1 = beta1_1sls * Z1 + beta2_1sls * Z2 + beta3_1sls * X2

beta_2SLS = inv(t(firstStageRegressors)%*%(firstStageRegressors)) %*% (t(firstStageRegressors) %*% P1)  
P1_hat = beta_2SLS[1]*Z1+beta_2SLS[2]*Z2+beta_2SLS[3]*X # form estimated endogenous variable P1

#estimate new original equation with P1_hat plugged into it instead of P1 by GC for endogenous regressor P2 since P1_hat is now considered exogenous

dataMatrix = data.frame(Y,P1_hat,P2,X)
theta_cop3 = copulaCorrection(Y~P1_hat+P2+X-1|continuous(P2),data= dataMatrix,num.boots=2)
return(theta_cop3$coefficients) 


}

IVCOPMIX1REVERSE <- function(Y,X,P1,P2,Z1,Z2){
  
  firstStageRegressors = matrix(c(Z1,Z2,X),ncol=3) # regress P2 on all instruments as well as regressor X2 which is also used as instrument
  # sO P2 = beta1_1sls * Z1 + beta2_1sls * Z2 + beta3_1sls * X2
  
  beta_2SLS = inv(t(firstStageRegressors)%*%(firstStageRegressors)) %*% (t(firstStageRegressors) %*% P2)  
  P2_hat = beta_2SLS[1]*Z1+beta_2SLS[2]*Z2+beta_2SLS[3]*X # form estimated endogenous variable P1
  
  #estimate new original equation with P1_hat plugged into it instead of P1 by GC for endogenous regressor P2 since P1_hat is now considered exogenous
  
  dataMatrix = data.frame(Y,P1,P2_hat,X)
  theta_cop3 = copulaCorrection(Y~P1+P2_hat+X-1|continuous(P1),data= dataMatrix,num.boots=2)
  return(theta_cop3$coefficients) 
  
  
}


OLS3 <- function(Y,X,P1,P2,P3,T){
  
  #VectorOf1 = rep(1,T)
  Xmatrix = matrix(c(P1,P2,P3,X),ncol=4)
  theta_ols = inv(t(Xmatrix)%*%(Xmatrix)) %*% (t(Xmatrix) %*% Y)
  return(theta_ols)
}

GC3 <- function(Y,X,P1,P2,P3){
  dataMatrix = data.frame(Y,P1,P2,P3,X)
  theta_cop = copulaCorrection(Y~P1+P2+P3+X-1|continuous(P1,P2,P3), data= dataMatrix,num.boots=2)
  return(theta_cop$coefficients)
}

IV3 <- function(Y,X,P1,P2,P3,Z1,Z2,Z3){
  dataMatrix = data.frame(Y,P1,P2,P3,X)
  theta_IV = ivreg(Y~P1+P2+P3+X-1|Z1+Z2+Z3+X,data=dataMatrix)# Z1 and Z2,Z3 are used as instruments and X2 as exogenous variable instrument, since you regress endogenous regressors on exogenous variables as well
  return(theta_IV$coefficients)
}

# IV for two of them (P1,P2) and GC for other one (P3)
IVCOPMIX2 <- function(Y,X,P1,P2,P3,Z1,Z2,Z3){
  
  firstStageRegressors = matrix(c(Z1,Z2,Z3,X),ncol=4) # regress P1 on all instruments as well as regressor X2 which is also used as instrument
  # sO P1 = beta1_1sls * Z1 + beta2_1sls * Z2 + beta3_1sls * Z3+ beta3_1sls * X2
  
  beta_2SLS = inv(t(firstStageRegressors)%*%(firstStageRegressors)) %*% (t(firstStageRegressors) %*% P1)  
  P1_hat = beta_2SLS[1]*Z1+beta_2SLS[2]*Z2+beta_2SLS[3]*Z3+beta_2SLS[4]*X # form estimated endogenous variable P1
  beta_2SLS_P2 = inv(t(firstStageRegressors)%*%(firstStageRegressors)) %*% (t(firstStageRegressors) %*% P2)  
  P2_hat =  beta_2SLS_P2[1]*Z1+beta_2SLS_P2[2]*Z2+beta_2SLS_P2[3]*Z3+beta_2SLS_P2[4]*X 
  #estimate new original equation with P1_hat plugged into it instead of P1 by GC for endogenous regressor P2 since P1_hat is now considered exogenous
  
  dataMatrix = data.frame(Y,P1_hat,P2_hat,P3,X)
  theta_cop3 = copulaCorrection(Y~P1_hat+P2_hat+P3+X-1|continuous(P3),data= dataMatrix,num.boots=2)
  return(theta_cop3$coefficients) 
  
  
}


# IV for one of them (P1) and GC for other two (P2,P3)
IVCOPMIX3 <- function(Y,X,P1,P2,P3,Z1,Z2,Z3){
  
  firstStageRegressors = matrix(c(Z1,Z2,Z3,X),ncol=4) # regress P1 on all instruments as well as regressor X2 which is also used as instrument
  # sO P1 = beta1_1sls * Z1 + beta2_1sls * Z2 + beta3_1sls * Z3+ beta3_1sls * X2
  
  beta_2SLS = inv(t(firstStageRegressors)%*%(firstStageRegressors)) %*% (t(firstStageRegressors) %*% P1)  
  P1_hat = beta_2SLS[1]*Z1+beta_2SLS[2]*Z2+beta_2SLS[3]*Z3+beta_2SLS[4]*X # form estimated endogenous variable P1

  #estimate new original equation with P1_hat plugged into it instead of P1 by GC for endogenous regressors P2,P3 since P1_hat is now considered exogenous
  
  dataMatrix = data.frame(Y,P1_hat,P2,P3,X)
  theta_cop3 = copulaCorrection(Y~P1_hat+P2+P3+X-1|continuous(P2,P3),data= dataMatrix,num.boots=2)
  return(theta_cop3$coefficients) 
  
  
}

