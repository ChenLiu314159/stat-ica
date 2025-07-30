library(VineCopula)
library(goftest)
library(KScorrect)
library(fGarch)

## 1.Data Overview 
# Read-in csv data of stock prices
data <- read.csv("CAC40_Nk225.csv")
CAC40 <- data$CAC40
Nk225 <- data$Nikkei225

# Calculate log-returns 
rt_CAC40 <- diff(log(CAC40))
rt_Nk225 <- diff(log(Nk225))

# Plot time series graphs
# We decide to observe only first 10 lags.
plot_ts <- function(x){
  par(mfrow=c(2,2))
  acf(x,lag.max= 10, col="green")
  pacf(x,lag.max= 10, col="blue")
  acf(x^2,lag.max= 10, col="red")
  par(mfrow=c(1,1))
}
plot_ts(rt_CAC40)
plot_ts(rt_Nk225)
# It is observed that, for CAC 40 indices, all lags in the ACF plot lie within 
# the significance level apart from the first lag. The tail off pattern is 
# significant. Its pacf graph represents a cut off at lag 1. Therefore, an
# AR(1) model is suggested.
# In terms of Nikkei 225, the exponential decay in the acf graph is significant.
# However, we couldn't detect any cut off in the pacf graph, and the autocorrelation
# test above indicates an independence among data. None AR model is suggested.
# The acf plots of squared returns for both indices reveal some serial correlation
# among variances with time, which suggest GARCH models for their conditional 
# variances.

#Test for autocorrelation/independency
Box1 <- Box.test(rt_CAC40, lag = 10, type = c("Ljung-Box"), fitdf = 1) 
Box2 <- Box.test(rt_Nk225, lag = 10, type = c("Ljung-Box"), fitdf = 1) 
Autocorrelation=rbind(c(Box1$p.value,Box2$p.value))
rownames(Autocorrelation) <- c("p-value")
colnames(Autocorrelation) <- c("CAC40 Log Return","Nikkei225 Log Return")
Autocorrelation
# The p-values of log-returns for CAC 40 and Nikkei 225 are 0.01 and
# 0.87 respectively. The p-value of CAC 40 lies below the significance level,
# which indicates an autocorrelation for CAC 40 data. On the other hand, 
# Nikkei 225 stock prices are considered as independently distributed.


## 2.Model Selection
# Model for CAC 40 first.

# Initialise an empty data frame to store the results
results <- data.frame(order = character(), AIC = numeric(), Distribution = character(), stringsAsFactors = FALSE)
distributions <- c("norm","snorm", "ged", "sged", "std", "sstd")

# Nested for loops to iterate over GARCH orders and distributions for CAC40
for (i in 1:3) {
  for (j in 1:3) {
    for (dist in distributions) {
      # Construct the formula for the current GARCH model order
      formula <- as.formula(paste("~ arma(1,0) + garch(", i, ",", j, ")", sep = ""))
      
      # Fit the GARCH model with the current order and distribution
      model <- garchFit(formula = formula, data = rt_CAC40, trace = F, cond.dist = dist)
      
      # Extract the AIC value
      aic_value <- model@fit$ics["AIC"]
      
      # Add the results to the dataframe
      results <- rbind(results, data.frame(order = paste("garch(", i, ",", j, ")", sep = ""),
                                           AIC = aic_value, Distribution = dist))

      
    }
  }
}

# Order the results by AIC
ordered_results_CAC40 <- results[order(results$AIC), ]

# Output the best model with the smallest value of AIC.
best_model_CAC40 <- head(ordered_results_CAC40, 1)
print(best_model_CAC40)
rm(formula, model, aic_value, results)# Remove repeated data
# The best fit model has been worked out. It has the lowest AIC value of  -4.510153
# among all models tested,
# and this GARCH(1,1) model follows the standardized skew Student-t distribution.

# Now it turns to Nikkei 225.
# Nested for loops to iterate over GARCH orders and distributions for Nikkei225
results <- data.frame(order = character(), AIC = numeric(), Distribution = character(), stringsAsFactors = FALSE)
distributions <- c("norm","snorm", "ged", "sged", "std", "sstd")

for (i in 1:3) {
  for (j in 1:3) {
    for (dist in distributions) {
      # Construct the formula for the current GARCH model order
      formula <- as.formula(paste("~garch(", i, ",", j, ")", sep = ""))
      
      # Fit the GARCH model with the current order and distribution
      model <- garchFit(formula = formula, data = rt_Nk225, trace = F, cond.dist = dist)
      
      # Extract the AIC value
      aic_value <- model@fit$ics["AIC"]
      
      # Add the results to the dataframe
      results <- rbind(results, data.frame(order = paste("garch(", i, ",", j, ")", sep = ""),
                                           AIC = aic_value, Distribution = dist))
      
      
    }
  }
}

# Order the results by AIC and Output the best model with the smallest AIC.
ordered_results_Nk225 <- results[order(results$AIC), ]
best_model_Nk225 <- head(ordered_results_Nk225, 1)
print(best_model_Nk225)
rm(formula, model, aic_value, results) # Remove repeated data
# The best fit model for Nikkei225 GARCH(1,1) with a lowest AIC of  -4.510153 and
# conditional variances follow the standardized skew Student-t distribution.

# Rename best fit models for CAC40 and Nikkei225.
model_CAC40 = garchFit(formula = ~ arma(1,0) +garch(1,1), data = rt_CAC40, trace = F, cond.dist="sstd")
model_Nk225 = garchFit(formula = ~ garch(1,1), data = rt_Nk225, trace = F, cond.dist="sstd")

# Check residuals.
res1 <- residuals(model_CAC40, standardize=TRUE)
res2 <- residuals(model_Nk225, standardize=TRUE)

plot_res <- function(x){
  par(mfrow=c(2,1))
  acf(x, main= "Residuals", lag.max= 10, col="green", lwd=2)
  acf(x^2, main= "Residuals squared", lag.max= 10, col="red", lwd=2)
  par(mfrow=c(1,1))
}

# Then check models fit.
plot_res(res1)
plot_res(res2)
# The graphs demonstrate ACF plots of residuals and residuals squared for two
# models. It is obvious that all plots appear a white noise pattern with
# lags lie within the significance level. This indicates that residuals appear 
# randomly, and model_CAC40 and model_Nk225 are of good fit.

Box.test(res1, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res2^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
# The two Ljung-Box tests conduct p-values of 0.2449 and 0.1533 for model_CAC40,
# and conduct p-values of 0.9343 and 0.5077 for model_Nk225.
# All p-values are higher than 0.05.
# Hence, there is no evidence to reject the null hypothesis that residuals
# are independently distributed. model_CAC40 and model_Nk225 pass Ljung-Box
# tests.

## 3.Constructing Copula 

# CAC40
# Check histogram
u1<-psstd(res1, mean=0, sd=1, nu=coef(model_CAC40)["shape"], xi=coef(model_CAC40)["skew"])[2:length(rt_CAC40)]
hist(u1,  main = "Histogram of CAC40 Model Residuals after PIT")
# Fairly uniform.

# Kolmogorov-Smirnov Test
KStest1<-LcKS(u1, cdf = "punif")
KStest1$p.value
# Pass

# Anderson-Darling test
ADtest1<-ad.test(u1, null="punif")
ADtest1$p.value
# Pass

# Nikkei225
# Check histogram
u2<-psstd(res2,mean=0, sd=1,nu=coef(model_Nk225)["shape"], xi=coef(model_Nk225)["skew"])
hist(u2,  main = "Histogram of Nikkei225 Model Residuals after PIT")

# Kolmogorov-Smirnov Test
KStest2<-LcKS(u2, cdf = "punif")
KStest2$p.value
# Pass the test.

# Anderson-Darling test
ADtest2<-ad.test(u2, null="punif")
ADtest2$p.value
# Pass the test

# A SYMMARY TABLE OF GOODNIESS OF FIT TEST FOR TWO MODELS
fittest=rbind(c(KStest1$p.value,ADtest1$p.value),c(KStest2$p.value,ADtest2$p.value))
rownames(fittest) <- c("CAC40","Nikkei225")
colnames(fittest) <- c("KS","AD")
fittest

# We pass the test for uniformity for both transformed log-returns, so we can proceed 
# to copula modelling.

# Fit Copula 
# After removing N/A values, we could find a difference in observations of 2 between
# two data. After conducting standard skew student-t distributions, there were 
# still one difference in observation numbers between CAC40 and Neikkei225.
# Therefore, we have to remove one observation from CAC40 in order to match lengths
# of u1 and u2.
length(u1)
length(u2)
u2 <-u2[-1]

# Using BiCopSelect function fit various copulas to the dataset and select
# the copula that provides the best fit based on the AIC criterion.
model=BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05)
model
# Bivariate copula: 
# Rotated Tawn type 2 180 degrees (par = 1.15, par2 = 0.27, tau = 0.06) 

# Value-at-Risk uisng Monte Carlo simulation
N=10000
set.seed(0123)
model=BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05)

par(mfrow=c(1,2))
plot(u1, u2 , ylab="CAC40",xlab="Neikki225", pch=20,
     col = rgb(128/255, 0/255, 128/255, 0.4),main="Observed Copula")
plot(u_sim[,1],u_sim[,2], ylab="CAC40",xlab="Neikki225", pch=20,
     col = rgb(128/255, 0/255, 128/255, 0.2),main="Simulated Copula")
par(mfrow=c(1,1))

# IPIT
rt_CAC40_sim<-qsstd(u_sim[,1], nu=coef(model_CAC40)["shape"], xi=coef(model_CAC40)["skew"])
rt_Nk225_sim<-qsstd(u_sim[,2], nu=coef(model_Nk225)["shape"], xi=coef(model_Nk225)["skew"])

## 4.Reintroducing AR-GARCH effects 

# Re-introduce AR-GARCH for CAC40
alpha1<-coef(model_CAC40)["ar1"]
mu<-coef(model_CAC40)["mu"]
alpha<-coef(model_CAC40)["alpha1"]
beta<-coef(model_CAC40)["beta1"]
omega <- coef(model_CAC40)["omega"]
res <-tail(model_CAC40@residuals,1)
epsilon <- rt_CAC40_sim
sigma<-tail(model_CAC40@sigma.t,1)

# Conditional variance at time t.
cond_var<- omega + alpha * res^2 + beta * sigma^2

CAC40_last<- mu + alpha1*tail(rt_CAC40,1) + sqrt(cond_var)*epsilon

# Remove
rm(alpha1,mu,alpha,beta,omega,res,epsilon,sigma,cond_var)


# Re-introduce GARCH for Nk225
mu<-coef(model_Nk225)["mu"]
alpha<-coef(model_Nk225)["alpha1"]
beta<-coef(model_Nk225)["beta1"]
omega <- coef(model_Nk225)["omega"]
res<-tail(model_Nk225@residuals,1)
epsilon <- rt_Nk225_sim
sigma<-tail(model_Nk225@sigma.t,1)


cond_var<- omega + alpha * res^2 + beta * sigma^2
Nk225_last<- mu + sqrt(cond_var)*epsilon

rm(mu,alpha,beta,omega,res,epsilon,sigma,cond_var)
par(mfrow=c(1,1))
plot(rt_CAC40, rt_Nk225, ylab="Neikki 225 Log Return",xlab="CAC 40 Log Return",
     main="Original and Simulated log-returns", pch=20,col = rgb(1, 165/255,0, 0.4))
points(CAC40_last, Nk225_last, pch = 20, col = rgb(128/225, 0, 128/255, 0.06))
legend(-0.02,-0.2, pch = 20, col = rgb(c(1, 128/255), c(165/255,0), c(0/255,128/255), 0.7), 
       legend = c("Original Log Return","Simulated Log Rerurn"), cex = 1)

# Compute VaR
port_sim <- matrix(0, nrow = N, ncol = 1)
var_sim <- matrix(0, nrow = 1, ncol = 2)
port_sim=log(1+((exp(CAC40_last)-1)+(exp(Nk225_last)-1))*(1/2))
var_sim=quantile(port_sim,c(0.01,0.05))
var_sim
port_sim <- matrix(0, nro= N, ncol = 1)
var_sim <- matrix(0, nrow = 1, ncol = 2)

port_sim=log(1+((exp(CAC40_last)-1)+(exp(Nk225_last)-1))*(1/2))

# The estimated 99% and 95% Value-at-Risk estimates are:
var_sim=quantile(port_sim,c(0.01,0.05))
var_sim



