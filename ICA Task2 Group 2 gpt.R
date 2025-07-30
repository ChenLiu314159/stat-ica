# Load necessary libraries
library(rugarch)
library(copula)
library(zoo)
library(tseries)

# CAC 40 Data Preparation
cac40 <- read.csv("C:\\Users\\吴睿仪\\Desktop\\UCL\\2024 term 2\\stat0011\\project\\gpt\\^FCHI.csv")
cac40 <- na.omit(cac40)
cac40_returns <- ts(cac40[,-1], start = c(2000, 1), frequency = 52)
CAC40_returns <- diff(log(cac40_returns[,5]))

# Nikkei 225 Data Preparation
nikkei <- read.csv("C:\\Users\\吴睿仪\\Desktop\\UCL\\2024 term 2\\stat0011\\project\\gpt\\^N225.csv") 
nikkei <- na.omit(nikkei)
nikkei_returns <- ts(nikkei[,-1], start = c(2000, 1), frequency = 52)
NIKKEI_returns <- diff(log(nikkei_returns[,5]))

# Fit GARCH(1,1) models to each index's returns
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                   mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
                   distribution.model = "sstd")

garch_cac40 <- ugarchfit(spec, CAC40_returns)
garch_nikkei <- ugarchfit(spec, NIKKEI_returns)

# Extract standardized residuals and conditional volatilities
std_resid_CAC40 <- residuals(garch_cac40, standardize = TRUE)
std_resid_N225 <- residuals(garch_nikkei, standardize = TRUE)

# Simulate future standardized residuals from the fitted GARCH models
n_sim <- 10000
sim_cac40 <- ugarchsim(garch_cac40, n.sim = n_sim, m.sim = 1)
sim_nikkei <- ugarchsim(garch_nikkei, n.sim = n_sim, m.sim = 1)

sim_cac40_resid <- as.vector(sigma(sim_cac40) * rnorm(n_sim))
sim_nikkei_resid <- as.vector(sigma(sim_nikkei) * rnorm(n_sim))

# Combine the residuals and transform them using empirical CDF (pobs)
u_cac40 <- pobs(sim_cac40_resid)
u_nikkei <- pobs(sim_nikkei_resid)
u_combined <- cbind(u_cac40, u_nikkei)

# Fit a copula to the transformed residuals
cop <- normalCopula(dim = 2)
fit_copula <- fitCopula(cop, u_combined, method = "ml")

# Generate samples from the copula
copula_samples <- rCopula(n_sim, copula = fit_copula@copula)

# Assuming equal weights for CAC 40 and Nikkei 225 in the portfolio
portfolio_returns_sim <- rowMeans(copula_samples) # Simplified portfolio returns from copula samples

# Conceptual adjustment for calculating losses
portfolio_losses = 1 - (1 + portfolio_returns_sim) # Reflecting losses from an initial investment of 1

# Calculate Value-at-Risk (VaR) based on losses
VaR_99 <- quantile(portfolio_losses, 0.99) # Calculating as negative losses
VaR_95 <- quantile(portfolio_losses, 0.95) # Calculating as negative losses

# Print VaR results
print(paste("99% VaR:", VaR_99))
print(paste("95% VaR:", VaR_95))

plot_ts <- function(x){
  par(mfrow=c(2,2))
  acf(x,lag.max= 10, col="green")
  pacf(x,lag.max= 10, col="blue")
  acf(x^2,lag.max= 10, col="red")
  par(mfrow=c(1,1))
}
plot_ts(NIKKEI_returns)
plot_ts(CAC40_returns)

library(VineCopula) # Add a pachage that chatgpt doesn't use
u1<-pnorm(std_resid_CAC40,mean=0, sd=1)
u2<-pnorm(std_resid_N225,mean=0, sd=1)
model=BiCopSelect(ts(u1), ts(u2), familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05)
model # student t copula is the best
u_sim=BiCopSim(N=10000, family= model$family, model$par,  model$par2)
par(mfrow=c(1,2))
plot(ts(u1), ts(u2) , ylab="CAC40",xlab="Neikki225", pch=20,
     col = rgb(128/255, 0/255, 128/255, 0.4),main="Observed Copula")
plot(u_sim[,1],u_sim[,2], ylab="CAC40",xlab="Neikki225", pch=20,
     col = rgb(128/255, 0/255, 128/255, 0.2),main="Simulated Copula")
par(mfrow=c(1,1))