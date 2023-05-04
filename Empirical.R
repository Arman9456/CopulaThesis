# Empirical application
# Retail price is endogenous, Bonus promotion as well
# Prices of different stores are correlated and often used as IV for eachother
# Direct price reduction variables is considered exogenous
# following Park & Gupta
# Use retail price of the other store as strong instrument for price
# Use cop for Bonus promotion since no strong instrument
library("REndo")
library("MASS")
library("matlib")
library("AER")
library("boot")
library("ivDiag")

load("C:/Users/arman/Desktop/Thesis MSC Copula/Empirical application R/store_level_data.RData")


Storedata = store_dat

Store1data = subset(Storedata,STORE==115)
Store2data = subset(Storedata,STORE==126)

# Delete empty 0 rows (also cant take log of 0's)
Store1data = subset(Store1data,PRICE!=0)
Store2data = subset(Store2data,PRICE!=0)

# Make sure that after deleting empty rows, the weeks of the two stores are matching
Store1data = Store1data[Store1data$WEEK %in% Store2data$WEEK,]
Store2data = Store2data[Store2data$WEEK %in% Store1data$WEEK,]

# some of the variables are in logs, others are not like in Yang paper
Y1 = log(Store1data$SALES)
PriceRed1 = Store1data$PRICEREDU # exogenous
RetailPrice1 = log(Store1data$PRICE) # endogenous with IV
BonusProm1 = Store1data$BONUS # endogenous without IV
Instrument1 = log(Store2data$PRICE) # instrument
Instrument2 = Store2data$BONUS # extra Instrument for 2SLS

plot(RetailPrice1,type="l")
plot(Instrument1,type="l")


Timepoints = 1:length(RetailPrice1)
DetrendModel1 = lm(RetailPrice1~Timepoints)
predicted1 = predict(DetrendModel1)

DetrendModel2 = lm(Instrument1~Timepoints)
predicted2 = predict(DetrendModel2)

#Detrended version of log(retailprice1) and log(instrument1) because we saw that they had a trend in data
RetailPrice1 = RetailPrice1 - predicted1
Instrument1 = Instrument1 - predicted2

plot(RetailPrice1,type="l")
plot(Instrument1,type="l")

DataFrame = data.frame(Y1,PriceRed1,RetailPrice1,BonusProm1,Instrument1,Instrument2)

# Check strength of instrument
ivDiag(data=DataFrame,Y="Y1",D="RetailPrice1",Z="Instrument1",controls="PriceRed1") # Anderson-Rubin test has H0: instrument is weak. We see that with an F-statistic of 32 and P-value of 0 we reject the null hypothesis of a weak instrument
ivDiag(data=DataFrame,Y="Y1",D="BonusProm1",Z="Instrument2",controls="PriceRed1") # F-stat (0.7399) and p-value=0.3902 of AR test for the instrument for bonusprom shows that this is a weak instrument




hist(RetailPrice1)

############# IVCOP below##################

firstStageRegressors = matrix(c(Instrument1,PriceRed1),ncol=2) # regress RetailPrice on all instruments as well as exogenous regressors
# sO RetailPrice1 = beta1_1sls * Instrument1 + beta2_1sls * PriceRed1

beta_2SLS = inv(t(firstStageRegressors)%*%(firstStageRegressors)) %*% (t(firstStageRegressors) %*% RetailPrice1)  
RetailPrice1_hat = beta_2SLS[1]*Instrument1+beta_2SLS[2]*PriceRed1 # form estimated endogenous variable RetailPrice

#estimate new original equation with P1_hat plugged into it instead of P1 by GC for endogenous regressor P2 since P1_hat is now considered exogenous

dataMatrix = data.frame(Y1,RetailPrice1_hat,BonusProm1,PriceRed1)
theta_cop3 = copulaCorrection(Y1~RetailPrice1_hat+BonusProm1+PriceRed1-1|continuous(BonusProm1),data= dataMatrix,num.boots=2)
theta_cop3_const = copulaCorrection(Y1~RetailPrice1_hat+BonusProm1+PriceRed1|continuous(BonusProm1),data= dataMatrix,num.boots=2) # with constant


alpha1IVCOP =   theta_cop3$coefficients[[1]] # coeff for retailprice1_hat
alpha2IVCOP =   theta_cop3$coefficients[[2]] # coeff for BonusProm1
betaIVCOP =   theta_cop3$coefficients[[3]] # coeff for PriceRed1

IVCOPconst = theta_cop3_const$coefficients[[1]] # constant
alpha1IVCOPconst =   theta_cop3_const$coefficients[[2]] # coeff for retailprice1_hat
alpha2IVCOPconst =   theta_cop3_const$coefficients[[3]] # coeff for BonusProm1
betaIVCOPconst =   theta_cop3_const$coefficients[[4]] # coeff for PriceRed1


########### IV(2SLS) below ###########

# BonusProm1 will have a negative coefficient with 2SLS with constant which seems unlikely, confirming that we had a weak instrument for bonusprom1
# Thus, IVCOP definitely does better since it has a positive coefficient for BonusProm1

theta_IV = ivreg(Y1~RetailPrice1+BonusProm1+PriceRed1-1|Instrument1+Instrument2+PriceRed1,data=DataFrame)
theta_IV_Const = ivreg(Y1~RetailPrice1+BonusProm1+PriceRed1|Instrument1+Instrument2+PriceRed1,data=DataFrame)

alpha1vecIV =   theta_IV$coefficients[[1]]
alpha2vecIV =   theta_IV$coefficients[[2]]
betavecIV =  theta_IV$coefficients[[3]]

IVconst = theta_IV_Const$coefficients[[1]]
alpha1IVconst =   theta_IV_Const$coefficients[[2]]
alpha2IVconst =   theta_IV_Const$coefficients[[3]]
betaIVconst =  theta_IV_Const$coefficients[[4]]


######### OLS below ################

Xmatrix = matrix(c(RetailPrice1,BonusProm1,PriceRed1),ncol=3)
theta_ols = inv(t(Xmatrix)%*%(Xmatrix)) %*% (t(Xmatrix) %*% Y1)
theta_ols

c = rep(1,length(RetailPrice1))

XmatrixConst = matrix(c(c,RetailPrice1,BonusProm1,PriceRed1),ncol=4)
theta_ols_const = inv(t(XmatrixConst)%*%(XmatrixConst)) %*% (t(XmatrixConst) %*% Y1)
theta_ols_const


############ Copula Below ##########

theta_cop = copulaCorrection(Y1~RetailPrice1+BonusProm1+PriceRed1-1|continuous(RetailPrice1,BonusProm1), data= DataFrame,num.boots=2)
theta_cop$coefficients

theta_cop_const = copulaCorrection(Y1~RetailPrice1+BonusProm1+PriceRed1|continuous(RetailPrice1,BonusProm1), data= DataFrame,num.boots=2)
theta_cop_const$coefficients



