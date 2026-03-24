#-------------------------------------------------------------------------------------------
# 0. Loading libraries
#-------------------------------------------------------------------------------------------

library(readxl)
library(quantmod)
library(PerformanceAnalytics)
library(robustbase)
library(rugarch)
library(Cairo)
library(RJDemetra)
library(tidyverse)
library(BEKKs)
library(rmgarch)
library(DescTools)
library(FinTS)
library(forecast)
library(MCS)

#-------------------------------------------------------------------------------------------
# 1. Exploratory analysis over 2021-2025
#-------------------------------------------------------------------------------------------

# Downloading the SSE Composite Index data from Yahoo
getSymbols(
    Symbols = "000001.SS", src = "yahoo",
    from = as.Date("2021-01-01"), to = as.Date("2025-12-31")
)
class(`000001.SS`)
head(`000001.SS`)
dim(`000001.SS`)
price <- `000001.SS`[, 4]

# The four price series together
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(`000001.SS`[, 1:4], legend.loc = "bottomleft", main = "SSE Composite Index Prices", col = rainbow(4))

#-------------------------------------------------------------------------------------------

# Computing the returns and squared returns of the SSE Composite Index
return <- dailyReturn(price)
creturn <- Return.clean(return, method = "boudt")

# Identification of corrected outliers
outliers_idx <- which(return != creturn)
outliers_df <- data.frame(
    Date = index(return)[outliers_idx],
    Original = coredata(return)[outliers_idx],
    Cleaned = coredata(creturn)[outliers_idx]
)
outliers_df$Difference <- outliers_df$Cleaned - outliers_df$Original
outliers_df <- outliers_df[order(abs(outliers_df$Difference), decreasing = TRUE), ]
outliers_df

creturn2 <- creturn^2

# Plot of the closing price and returns of the SSE Composite Index
par(mfrow = c(2, 1))
plot.xts(`000001.SS`[, 4], , main = "SSE Composite Index Prices", col = rainbow(4))
plot.xts(return, main = "SSE Composite Index Returns", col = "blue")

# Plot of the returns and cleaned returns of the SSE Composite Index
par(mfrow = c(1, 1))
data <- cbind(creturn, return)
colnames(data) <- c("cleaned returns", "returns")
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(data, legend.loc = "topleft", main = NULL, col = rainbow(4))

# Plot of the returns and squared returns of the SSE Composite Index
par(mfrow = c(2, 1))
plot.xts(creturn, main = "Cleaned Daily Returns", col = rainbow(4))
plot.xts(creturn2, main = "Cleaned Squared Daily Returns", col = "blue")


#-------------------------------------------------------------------------------------------
# 	Correlogram + histogram
#-------------------------------------------------------------------------------------------
# Correlograms of returns
par(mfrow = c(2, 2))
acf(creturn, main = "ACF of Returns")
pacf(creturn, main = "PACF of Returns")
acf(creturn2, main = "ACF of Squared Returns")
pacf(creturn2, main = "PACF of Squared Returns")

# Histogram of returns
par(mfrow = c(1, 1))
# chart.Histogram(creturn, xlab = "Cleaned Daily Returns", ylab = "Frequencies", methods = c("add.density", "add.normal"))

# 1. Margins: c(bottom, left, top, right) = c(5, 4, 4, 2) by default
par(mar = c(4, 4, 1, 1))
# 2. Histogram without automatic titles
chart.Histogram(creturn,
    methods = c("add.density", "add.normal"),
    ann = FALSE
)
title(xlab = "Cleaned Daily Returns", ylab = "Density")

#-------------------------------------------------------------------------------------------
# Descriptive statistics
#-------------------------------------------------------------------------------------------
# descriptive stat on returns with PerformanceAnalytics package (warning: give Ex.Kurtosis)
stats <- table.Stats(creturn * 100)
stats
table.Distributions(creturn)

# Skewness
skew <- stats["Skewness", "daily.returns"] # -0.2449

# Excess Kurtosis ("Kurtosis" dans le tableau)
excess_kur <- stats["Kurtosis", "daily.returns"] # 3.2566

# Sample size (T)
T_obs <- stats["Observations", "daily.returns"] # 5035

# Skewness test statistic
nu_1 <- skew / sqrt(6 / T_obs)
nu_1

# Kurtosis test statistic
nu_2 <- excess_kur / sqrt(24 / T_obs)
nu_2

# Jarque-Bera normality test with DescTools package
JarqueBeraTest(creturn, robust = FALSE, method = "chisq")

# autocorrelation Box-Pierce and Ljung-Box tests with 10 lags
Box.test(creturn, lag = 10, type = c("Box-Pierce", "Ljung-Box"), fitdf = 0)

# ARCH test with 5 and 10 lags
ArchTest(creturn, lag = 5)
ArchTest(creturn, lag = 10)

#-------------------------------------------------------------------------------------------
# 2. Volatility model estimation over 2021-2024 (on cleaned data)
#-------------------------------------------------------------------------------------------


# 2.1. Normal distribution
#-------------------------------------------------------------------------------------------

#----------------------------------------------------------------
# GARCH modelling
spec <- ugarchspec(variance.model = list(model = "sGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))
mod <- ugarchfit(data = creturn, spec = spec)
mod

# Compute persistence
pers <- persistence(mod)
show(pers)

# Compute half-life
hl <- halflife(mod)
show(hl)

#----------------------------------------------------------------
# IGARCH modelling
spec <- ugarchspec(variance.model = list(model = "iGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))
mod <- ugarchfit(data = creturn, spec = spec)
mod

#----------------------------------------------------------------
# Risk modelling
spec <- ugarchspec(variance.model = list(model = "iGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), fixed.pars = list(omega = 0, alpha1 = 0.06, beta1 = 0.94))
mod <- ugarchfit(data = creturn, spec = spec)
mod

#----------------------------------------------------------------
# GJR-GARCH modelling
spec <- ugarchspec(variance.model = list(model = "gjrGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))
mod <- ugarchfit(data = creturn, spec = spec)
mod


# 2.2. Student's t-distribution
#-------------------------------------------------------------------------------------------

#----------------------------------------------------------------
# GARCH modelling
spec <- ugarchspec(variance.model = list(model = "sGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std")
mod <- ugarchfit(data = creturn, spec = spec)
mod

# Compute persistence
pers <- persistence(mod)
show(pers)

# Compute half-life
hl <- halflife(mod)
show(hl)

#----------------------------------------------------------------
# IGARCH modelling
spec <- ugarchspec(variance.model = list(model = "iGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std")
mod <- ugarchfit(data = creturn, spec = spec)
mod

#----------------------------------------------------------------
# Risk modelling
spec <- ugarchspec(variance.model = list(model = "iGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), fixed.pars = list(omega = 0, alpha1 = 0.06, beta1 = 0.94), distribution.model = "std")
mod <- ugarchfit(data = creturn, spec = spec)
mod

#----------------------------------------------------------------
# GJR-GARCH modelling
spec <- ugarchspec(variance.model = list(model = "gjrGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std")
mod <- ugarchfit(data = creturn, spec = spec)
mod

#-------------------------------------------------------------------------------------------
# Volatility and Value at Risk (VaR) forecasting over 2025
#-------------------------------------------------------------------------------------------

y <- creturn

show(return[969, ]) # 2024-12-31	-0.01630692
b <- nrow(return) # b=1211
estim <- 969 # number of return observations from 2021 to 2024
h <- b - 970 + 1 # number of return observations in 2025
original <- return[970:1211, 1] # original returns in 2025

#----------------------------------------------------------
# Matrix initialization
foremat <- matrix(nrow = h, ncol = 1) # matrix containing the variance forecasts
varmat <- matrix(nrow = h, ncol = 1) # matrix containing the VaR forecasts
esmat <- matrix(nrow = h, ncol = 1)

#----------------------------------------------------------
# Estimating GARCH(1,1) with only an intercept in the mean
spec <- ugarchspec(variance.model = list(model = "sGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))

for (i in 1:h)
{
    yy <- y[i:(estim - 1 + i), 1]
    fit <- ugarchfit(data = yy, spec = spec)
    forc <- ugarchforecast(fit, n.ahead = 1)
    foremat[i, 1] <- sigma(forc)^2
    varmat[i, 1] <- qnorm(0.05) * sigma(forc)
    esmat[i, 1] <- -dnorm(qnorm(0.05)) / 0.05 * sigma(forc)
}

error_sGARCH <- original^2 - foremat
mse_sGARCH <- mean(error_sGARCH^2)
mse_sGARCH * 100000

mean(varmat) # Mean VaR
mean(esmat) # Expected Shortfall

# VaR figure
data <- cbind(original, varmat)
colnames(data) <- c("returns", "VaR")
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(data, legend.loc = "bottomleft", main = "", col = rainbow(4))

# Kupiec and Engle-Manganelli VaR tests
print(VaRTest(0.05, as.numeric(original), as.numeric(varmat)))

#----------------------------------------------------------
# Estimating GJR(1,1) with only an intercept in the mean
spec <- ugarchspec(variance.model = list(model = "gjrGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))

for (i in 1:h)
{
    yy <- y[1:(estim - 1 + i), 1]
    fit <- ugarchfit(data = yy, spec = spec)
    forc <- ugarchforecast(fit, n.ahead = 1)
    foremat[i, 1] <- sigma(forc)^2
    varmat[i, 1] <- qnorm(0.05) * sigma(forc)
    esmat[i, 1] <- -dnorm(qnorm(0.05)) / 0.05 * sigma(forc)
}

error_gjrGARCH <- original^2 - foremat
mse_gjrGARCH <- mean(error_gjrGARCH^2)
mse_gjrGARCH * 100000
dm.test(error_sGARCH, error_gjrGARCH, h = 1)

mean(varmat) # Mean VaR
mean(esmat) # Expected Shortfall

# VaR figure
data <- cbind(original, varmat)
colnames(data) <- c("returns", "VaR")
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(data, legend.loc = "bottomleft", main = "", col = rainbow(4))

# Kupiec and Engle-Manganelli VaR tests
print(VaRTest(0.05, as.numeric(original), as.numeric(varmat)))

#----------------------------------------------------------
# Estimating GARCH(1,1) Student with only an intercept in the mean
spec <- ugarchspec(variance.model = list(model = "sGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std")

for (i in 1:h)
{
    yy <- y[1:(estim - 1 + i), 1]
    fit <- ugarchfit(data = yy, spec = spec)
    forc <- ugarchforecast(fit, n.ahead = 1)
    foremat[i, 1] <- sigma(forc)^2
    varmat[i, 1] <- qnorm(0.05) * sigma(forc)
    esmat[i, 1] <- -dnorm(qnorm(0.05)) / 0.05 * sigma(forc)
}

error_sGARCH_std <- original^2 - foremat
mse_sGARCH_std <- mean(error_sGARCH_std^2)
mse_sGARCH_std * 100000
dm.test(error_sGARCH, error_sGARCH_std, h = 1)

mean(varmat) # Mean VaR
mean(esmat) # Expected Shortfall

# VaR figure
data <- cbind(original, varmat)
colnames(data) <- c("returns", "VaR")
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(data, legend.loc = "bottomleft", main = "", col = rainbow(4))

# Kupiec and Engle-Manganelli VaR tests
print(VaRTest(0.05, as.numeric(original), as.numeric(varmat)))

#----------------------------------------------------------
# Estimating GJR(1,1) Student with only an intercept in the mean
spec <- ugarchspec(variance.model = list(model = "gjrGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std")

for (i in 1:h)
{
    yy <- y[1:(estim - 1 + i), 1]
    fit <- ugarchfit(data = yy, spec = spec)
    forc <- ugarchforecast(fit, n.ahead = 1)
    foremat[i, 1] <- sigma(forc)^2
    varmat[i, 1] <- qnorm(0.05) * sigma(forc)
    esmat[i, 1] <- -dnorm(qnorm(0.05)) / 0.05 * sigma(forc)
}

error_gjrGARCH_std <- original^2 - foremat
mse_gjrGARCH_std <- mean(error_gjrGARCH_std^2)
mse_gjrGARCH_std * 100000
dm.test(error_sGARCH, error_gjrGARCH_std, h = 1)

mean(varmat) # Mean VaR
mean(esmat) # Expected Shortfall

# VaR figure
data <- cbind(original, varmat)
colnames(data) <- c("returns", "VaR")
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(data, legend.loc = "bottomleft", main = "", col = rainbow(4))

# Kupiec and Engle-Manganelli VaR tests
print(VaRTest(0.05, as.numeric(original), as.numeric(varmat)))

#-------------------------------------------------------------------------------------------
# Global forecast evaluation via the Model Confidence Set (MCS)

# 1. Construction of the loss matrix
# The loss function used here is the quadratic loss (squared errors)
Loss_matrix <- cbind(
    error_sGARCH^2,
    error_gjrGARCH^2,
    error_sGARCH_std^2,
    error_gjrGARCH_std^2
)
colnames(Loss_matrix) <- c("sGARCH", "gjrGARCH", "sGARCH_std", "gjrGARCH_std")

# 2. Execution of the MCS procedure
# The confidence level alpha (typically 0.1 or 0.05) represents the risk level of rejecting a well-performing model.
# B = 5000: number of bootstrap replications
MCS_result <- MCSprocedure(Loss = Loss_matrix, alpha = 0.10, B = 5000, statistic = "Tmax", cl = NULL)
