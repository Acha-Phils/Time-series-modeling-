################################################################################
#      Time series modeling 
#       April 2024
################################################################################

library(rugarch)
library(forecast)
library(tseries)
library(fGarch)

################### Data Extraction ############################################
# Extract the data using get.hist.quote function
ftmc_dat <- get.hist.quote(instrument = "^FTMC", start="2006-01-01", end ="2023-03-30", 
                            quote = c("Open", "High", "Low", "Close", "Volume"),
                            provider="yahoo", retclass = "zoo", compression="d")
str(ftmc_dat)
ftmc<-na.omit(ftmc_dat)
str(ftmc)

# randomly select 4000 samples
set.seed(123)
n<-nrow(ftmc)
ftmc_4000<- ftmc[1:4000, ]
str(ftmc_4000)

# Extract the returns from the close column
ftmc_return<-ftmc_4000[,4]
str(ftmc_return)

###########  1 #########################################################

# plotting the time series of the data over time
par(mfrow=c(1,3))
plot(ftmc_return, xlab= "Year", ylab= "Daily Returns",
     main="Plot of daily returns over time", col="skyblue")
# Does not look to be stationary 

# Testing for Stationarity using the adf test
# Using augmented Dickey–Fuller test (ADF)test Ho: data is not stationary
adf.test(ftmc_return) 
# p-value=0.2584 implies not stationary data
log_return<-log(ftmc_return)

# plotting the log returns of the time series
plot(log_return, main= "Plot of log returns", ylab="Log returns",
     col="skyblue")
#Indicates non stationary

# ADF test of stationarity
adf.test(log_return)

# Differencing the log returns of the time series to achieve stationarity
df_log_return<-diff(log_return)
# Plotting the first difference of the series
plot(df_log_return,xlab="Year",ylab="Daily Returns",
     main="Plot of log difference of the time series", col="skyblue")
par(mfrow=c(1,1))

# Testing for stationarity using the (ADF) test
adf.test(df_log_return)

####################  2 #################################################
# Selecting a suitable ARMA process with auto.arima function
auto.arima(df_log_return, max.p=5, max.q=5, max.P=5, ic='bic')

#AR(2) looks to be the best model  
return.AR2<-arima(df_log_return,c(2,0,0))
residAR2<-na.remove(return.AR2$residuals)
residAR2_sq<-residAR2^2

# Residual plots
par(mfrow=c(1,2))
#ACF plot
acf(residAR2,cex.axis=1.5,cex.lab=1.5,main="") 
title("ACF plot of AR(2) residuals")
#PACF plot
pacf(residAR2,cex.axis=1.5,cex.lab=1.5,main="") 
title("PACF plot of AR(2) squared residuals",cex.main=1.5)
par(mfrow=c(1,1))

# Performing Box test on the ARMA(2,0) model
Ljung.test.ar2<-Box.test(residAR2, lag=10, type="Ljung-Box")
Pierce.test.ar2<-Box.test(residAR2, lag=10, type="Box-Pierce")
Ljung.test.ar2
Pierce.test.ar2

###################################################
#Fitting an AR(2)-GARCH(1,1) process with Gaussian noise
garch.spec.norm=ugarchspec(variance.model=list(garchOrder=c(1,1)),
                             mean.model=list(armaOrder=c(2,0)))
return.garch.norm = ugarchfit(spec=garch.spec.norm, data=df_log_return)

coef(return.garch.norm)
infocriteria(return.garch.norm)

## Obtain residuals and Model diagnostics#######
residuals.n <- return.garch.norm@fit$residuals
sq_residuals.n<-residuals.n^2

# Plot ACF of residuals
par(mfrow=c(2,2))
acf(residuals.n, main="ACF of ARMA(2,0)-GARCH(1,1)")

# Plot PACF of residuals
pacf(residuals.n, main="PACF of ARMA(2,0)-GARCH(1,1)")

# Normality of Residuals
#QQ-plot
qqnorm(residuals.n, main="ARMA(2,0)-GARCH(1,1)")
qqline(residuals.n, col='red', )

# hist of residuals
h<-hist(residuals.n, col="skyblue", main="ARMA(2,0)-GARCH(1,1) Residuals")
xfit<-seq(min(residuals.n),max(residuals.n),length=length(residuals.n))
yfit<-dnorm(xfit,mean=mean(residuals.n),sd=sd(residuals.n))
yfit <- yfit*diff(h$mids[1:2])*length(residuals.n)
lines(xfit, yfit, col='red', lwd=2)
par(mfrow=c(1,1))

#Test  for normality
shapiro.test(residuals.n)

# Plot ACF and PACF of squared residuals
acf(sq_residuals.n, main="ACF of ARMA(2,0)-GARCH(1,1) squared residuals")
pacf(sq_residuals.n, main="ACF of ARMA(2,0)-GARCH(1,1) squared residuals")

# Plot standardized residuals
residuals.std <- residuals.n/sigma(return.garch.norm)
plot(residuals.std, type = "l", main = "Standardized Residuals")
abline(h = 0, col = "red")

# Plot squared standardized residuals
plot(residuals.std^2, type = "l", main = "Squared Standardized Residuals")
abline(h = 1, col = "red")


plot(sigma(return.garch.norm),
     main="AR(2)-GARCH(1,1) Conditional Volatility")

#########################################################################################

#Fitting an ARMA(2,0)-GARCH(1,1) process with t-distributed noise
garch.spec.t = ugarchspec(variance.model=list(garchOrder=c(1,1)), 
                          mean.model=list(armaOrder=c(2,0)),distribution.model="std")
return.garch.t = ugarchfit(spec=garch.spec.t, data=df_log_return)
coef(return.garch.t)
infocriteria(return.garch.t)

## Obtain residuals and Model diagnostics#######
residuals.t <- return.garch.t@fit$residuals
sq_residuals.t<-residuals.t^2
residuals.std <- residuals.t/sigma(return.garch.t)

# Plot ACF of residuals
par(mfrow=c(2,2))
acf(residuals.t, main="ARMA(2,0)-GARCH(1,1)") 

# Plot PACF of residuals
pacf(residuals.t, main="ARMA(2,0)-GARCH(1,1)")
length(residuals.t)

#ACF of squared residuals
acf(sq_residuals.t, main="ACF of ARMA(2,0)-GARCH(1,1) squared residuals")

#QQplot of residuals
df.est = coef(return.garch.t)["shape"]
qqplot(qt(ppoints(3999), df.est), residuals.t, col='skyblue', 
       main="QQplot ofARMA(2,0)-GARCH(1,1)")
par(mfrow=c(1,1))

#note the increase in volatility in the series
plot(residuals.std, type = "l", ylab = "Standardized Residuals",
     main = "ARMA(2,0)-GARCH(1,1) Model")
abline(h = 0, col = "red")

# Plot ACF and PACF of squared residuals
acf(sq_residuals.t) #title("ACF of ARMA(2,0)-GARCH(1,1) squared residuals with t-dist")
pacf(sq_residuals.t) #title("PACF of ARMA(2,0)-GARCH(1,1) squared residuals with t-dist")

# Perform Ljung-Box test
ljung_box_test.t <- Box.test(residuals.t, lag = 10, type = "Ljung-Box")
ljung_box_test.t
# Perform Box-Pierce test
box_pierce_test.t <- Box.test(residuals.t, lag = 10, type = "Box-Pierce")
box_pierce_test.t

plot(sigma(return.garch.t),main="ARMA(2,0)-GARCH(1,1) Conditional Volatility")
#note the increase in volatility in the series

########################################################################################
#           To capture asymmetric volatility
#######################################################################################
# Fit ARMA(2,0)-iGARCH(1,1) Model 
igarch_spec <- ugarchspec(variance.model=list(model = "iGARCH", garchOrder=c(1,1)), 
                          mean.model = list(armaOrder=c(2,0)), distribution.model="std")
fit.igarch<- ugarchfit(igarch_spec, data = df_log_return)

#Fit ARMA(2,0)-eGARCH model
spec <- ugarchspec(mean.model = list(armaOrder = c(2, 0)), 
                   distribution.model = "std", 
                   variance.model = list(model = "eGARCH", garchOrder = c(1, 1)))
fit.egarch <- ugarchfit(spec, df_log_return)


####### select the better fitting model ###########

infocriteria(fit.igarch)
infocriteria(fit.egarch)
# based on BIC eGARCH model is better

### Best model is ARMA(2,0)-eGARCH ################################


## Obtain residuals and Model diagnostics#######
residuals.egarch <- fit.egarch@fit$residuals
sq_residuals.egarch<-residuals.egarch^2
residuals.estd <- residuals.egarch/sigma(fit.egarch)

# Plot ACF of residuals
par(mfrow=c(2,2))
acf(residuals.egarch, main="ARMA(2,0)-eGARCH(1,1)")

# Plot PACF of residuals
pacf(residuals.egarch, main="ARMA(2,0)-eGARCH(1,1)")

#QQplot of residuals
df = coef(fit.egarch)["shape"]
qqplot(qt(ppoints(3999), df), residuals.egarch, col='skyblue', 
       main="QQplot ofARMA(2,0)-eGARCH(1,1)")
#note the increase in volatility in the series
plot(residuals.estd, type = "l", ylab = "Standardized Residuals",
     main = "ARMA(2,0)-GARCH(1,1) Model")
abline(h = 0, col = "red")
par(mfrow=c(1,1))




shapiro.test(residuals.egarch)


# Plot ACF and PACF of squared residuals
acf(sq_residuals.egarch); title("ACF of ARMA(2,0)-GARCH(1,1) squared residuals")
pacf(sq_residuals.egarch); title("PACF of ARMA(2,0)-GARCH(1,1) squared residuals")
par(mfrow=c(1,1))


# Perform Ljung-Box test
ljung_box_test.e <- Box.test(residuals.egarch, lag = 10, type = "Ljung-Box")
ljung_box_test.e
# Perform Box-Pierce test
box_pierce_test.e <- Box.test(residuals.egarch, lag = 10, type = "Box-Pierce")
box_pierce_test.e

# Conditional Volatility plot
plot(sigma(fit.egarch),main="ARMA(2,0)-eGARCH(1,1) Conditional Volatility",
     col='skyblue')

# Compare the volatility of the two final models
vole_garch<-sigma(return.garch.t)
vole_egarch<-sigma(fit.egarch)

# Plot volatility of both models on the same plot
plot(vole_garch, col="blue", main="Volatility Comparison")
lines(vole_egarch, col = "red")
legend("top", legend=c("GARCH", "eGARCH"), col = c("blue", "red"), lty =1)

# Model criteria selection for all models considered
infocriteria(return.garch.norm)
infocriteria(return.garch.t)
infocriteria(fit.egarch)
infocriteria(fit.igarch)


cor(residuals.t,residuals.egarch)




# Check how close the two processes are
cor(vole_garch,vole_egarch)









