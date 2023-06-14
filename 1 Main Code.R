#install.packages("REndo")
library("REndo")
library("MASS")
library("matlib")
library("AER")
library("boot")
rm(list = ls())

source("Functions.R")
source("Simulations.R")
source("IntceptFunctions.R")
source("IntCeptSimulations.R")
set.seed(6187807) # make results reproducable
# Case 1 endogenous regressor and 1 exogenous regressor Y_t = alpha*P_t + beta*X_t + eps_t (or with intercept? but park & gupta dont look at intercept)
# But copula function in R does include intercept
# X_t is N(0,1)



N = 100 #Number of Simulations 
T = 100 # size of sample
corr = 0.9 # correlation between error terms
cov = corr* (sqrt(1))

bivariate_errors = mvrnorm(n=T,mu=c(0, 0),Sigma=matrix(c(1, cov,cov, 1), ncol=2))
# IF YOU PUT HIGH VARIANCE (SIGMA = C(3,rho,rho,2), THEN GC PERFORMS WELL)


# Construct endogenous variable that has instrument Z which is perfectly correlated with P but uncorrelated to error term of original equation
# Generate endogenous regressor P with instrument Z~gamma(1,1) s.t. P is nonnormal just like in GC assumption
Z = rgamma(T,shape=1) # could also use uniform distribution just like in Park & Gupta (or Gamma like in other thesis)
P = Z + bivariate_errors[,2] # should I add intercept to this equation?

#cov(bivariate_errors[,1], P) #check whether covariance between P and error term of original equation equals rho (TRUE)
# Generate X_t's 
X1 = rnorm(n=T,mean=0,sd=1)
#X2 = rgamma(T,shape=1)

# Coefficient values of original equations
#intcpt = rep(0.5,T)
alpha = -1
beta1 = 2



Y =  alpha * P + beta1*X1 + bivariate_errors[,1]

dataMatrix = data.frame(Y,X1,P)
theta_cop = copulaCorrection(Y~X1+P-1|continuous(P),data= dataMatrix,num.boots=2) #-1 makes sure there is no intercept in the model
theta_cop$coefficients


################# BOOTSTRAP CONFIDENCE INTERVAL ATTEMPT (NOT REALLY WORKING) #####################
# GC1boot <- function(data){
#   
#   Y = dataMatrix[,1]
#   X = dataMatrix[,2]
#   P = dataMatrix[,3]
#   dataMatrix = data.frame(Y,X,P)
#   theta_cop = copulaCorrection(Y~X+P-1|continuous(P),data= dataMatrix,num.boots=2)
#   return(theta_cop$coefficients)
# }
# GC1boot(dataMatrix)
# 
# bootsamples = boot(data=dataMatrix,statistic=GC1boot,R=10,sim="parametric")
# boot.ci(bootsamples)


testdiff = rnorm(100)
hist(testdiff,breaks=20,main = expression(paste("Distribution of Bias",(delta[1]))),xlab="Bias GC")
#Mean OLS/GC/IV estimates of case 1 endo 1 exo
results1 = Simul1(T,N,cov)
results1

results1Intcept = Simul1INT(T,N,cov) # changing NR of bootstraps to 1000 gives even worse bias
results1Intcept


GC1(Y,X1,P)[[2]]
OLS1(Y,X1,P,T)[1]
IV1(Y,X1,P,Z)[[2]]
# GC combined with IV -> Not possible yet I think with only 1 endogenous regressor since you want to apply IV to one endogenous regressor and GC to other endogenous regressor
# For Combined GC with IV, I think we regress the IV endogenous variable P on its strong instrument and plug the estimated P the original equation
# Then we treat this estimated P as exogenous and apply GC to the original equation with this exogenous variable P_hat and the other endogenous variable


################# Case Two endogenous regressors ###################
corr12 = 0.9
corr13 = 0.5

rho_12 = corr12* (sqrt(1)*sqrt(1))
rho_13 = corr13 *sqrt(1)*sqrt(1)
rho_23 = 0.2


#Create errors that are correlated across all equations (the two endogenous regressor equations and the original equation) just like in Park & Gupta p.32
# Because if both endogeneous regressors are correlated with the error term from the original equation, they must be correlated to eachother
trivariate_errors = mvrnorm(n=T,mu=c(0, 0,0),Sigma=matrix(c(1, rho_12,rho_13,
                                                          rho_12,1,rho_23,
                                                          rho_13,rho_23,1), ncol=3,byrow=TRUE))

X2 = rnorm(n=T,mean=0,sd=1) # create an exogenous variable, since only 2 endogenous variables gave weird dresults
beta2= 0.5

alpha1 = -1
alpha2 = 1

Z1 = rgamma(n=T, shape = 2, rate = 0.5) # instrument for endogenous regressor P1 s.t. P1 is not normally distributed
Z2 = rgamma(n=T,shape=2,1)  # could also use uniform distribution just like in Park & Gupta

P2 = Z2 + trivariate_errors[,3] 
P1 = Z1 + trivariate_errors[,2]

intcpt2 = rep(0.3,T)
# Construct sample from DGP
Y2=  (alpha1 * P1) + (alpha2*P2) + beta2*X2 + trivariate_errors[,1]

plot(density(Z1))

# OLS
OLS2(Y2,X2,P1,P2,T)

# GC estimates
GC2(Y2,X2,P1,P2)

GC2INT(Y2,X2,P1,P) # with intcept

#IV estimates
IV2(Y2,X2,P1,P2,Z1,Z2)

IV2INT(Y2,X2,P1,P2,Z1,Z2)
# Iv for one of them (P1) and GC for other one (P2)

# first estimate P1 equation by OLS as first stage in 2SLS, I don't think I should include P2 as a regressor in that equation, because P2 
# is correlated with the original error term. (If P2 were to be an exogenous regressor, I would have to include it as a first stage instrument as well)
#https://timeseriesreasoning.com/contents/two-stage-least-squares-estimation/


IVCOPMIX1(Y2,X2,P1,P2,Z1,Z2) #  Pretty decent estimates
IVCOPMIX1INT(Y2,X2,P1,P2,Z1,Z2) # With intcpt

IVCOPMIX1REVERSE(Y2,X2,P1,P2,Z1,Z2) 
IVCOPMIX1REVERSEINT(Y2,X2,P1,P2,Z1,Z2) 

results2 = Simul2(T,N,rho_12,rho_13,rho_23)
results2

results2Intcpt = Simul2INT(T,N,rho_12,rho_13,rho_23)
results2Intcpt


# CASE 3 ENDOGENOUS REGRESSORS, 1 EXOGENOUS REGRESSOR
rho_14 = 0.4
rho_24 = 0.2
rho_34 = 0.1

# just to check if sigma is positive definite
multivariate_errors = mvrnorm(n=T,mu=c(0, 0,0,0),Sigma=matrix(c(1, rho_12,rho_13,rho_14,
                                                                              rho_12,1,rho_23,rho_24,
                                                                              rho_13,rho_23,1,rho_34
                                                            ,rho_14,rho_24,rho_34,1), ncol=4,byrow=TRUE))


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

IVCOPMIX3(Y2,X2,P1,P2,P3,Z1,Z2,Z3)

# Case IV for P1,P2 and GC for P3
results3 = Simul3(T,N,rho_12,rho_13,rho_23,rho_14,rho_24,rho_34)
results3

results3INT = Simul3INT(T,N,rho_12,rho_13,rho_23,rho_14,rho_24,rho_34)
results3INT

# Case IV for P1 and GC for P2,P3
results4 = Simul4(T,N,rho_12,rho_13,rho_23,rho_14,rho_24,rho_34)
results4

results4INT = Simul4INT(T,N,rho_12,rho_13,rho_23,rho_14,rho_24,rho_34)
results4INT

