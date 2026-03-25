# ==============================================================================
# 0. Loading libraries
# ==============================================================================

library(quantmod) # Section 1.1
library(PerformanceAnalytics) # Section 1.1
library(robustbase) # Section 1.1
library(DescTools) # Section 1.3
library(FinTS) # Section 1.3
library(rugarch) # Section 2
library(forecast) # Section 3
library(MCS) # Section 4

# ==============================================================================
# 1. Exploratory and descriptive analysis (2021-2025)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1.1. Data and return series
# ------------------------------------------------------------------------------

# Downloading the SSE Composite Index data from Yahoo
getSymbols(
    Symbols = "000001.SS", src = "yahoo",
    from = as.Date("2021-01-01"), to = as.Date("2025-12-31")
)
price <- `000001.SS`[, 4]

# The four price series together
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(`000001.SS`[, 1:4], legend.loc = "bottomleft", main = "SSE Composite Index Prices", col = rainbow(4))

# Computing the returns and squared returns of the SSE Composite Index
returns <- dailyReturn(price)
clean_returns <- Return.clean(returns, method = "boudt")

# Identification of corrected outliers
outliers_idx <- which(returns != clean_returns)
outliers_df <- data.frame(
    Date = index(returns)[outliers_idx],
    Original = coredata(returns)[outliers_idx],
    Cleaned = coredata(clean_returns)[outliers_idx]
)
outliers_df$Difference <- outliers_df$Cleaned - outliers_df$Original
outliers_df <- outliers_df[order(abs(outliers_df$Difference), decreasing = TRUE), ]
outliers_df

clean_returns_sq <- clean_returns^2

# Plot of the closing price and returns of the SSE Composite Index
par(mfrow = c(2, 1))
plot.xts(`000001.SS`[, 4], main = "SSE Composite Index Prices", col = rainbow(4))
plot.xts(returns, main = "SSE Composite Index Returns", col = "blue")

# Plot of the returns and cleaned returns of the SSE Composite Index
par(mfrow = c(1, 1))
plot_data <- cbind(clean_returns, returns)
colnames(plot_data) <- c("Cleaned Returns", "Returns")
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(plot_data, legend.loc = "topleft", main = NULL, col = rainbow(4))

# Plot of the returns and squared returns of the SSE Composite Index
par(mfrow = c(2, 1))
plot.xts(clean_returns, main = "Cleaned Daily Returns", col = rainbow(4))
plot.xts(clean_returns_sq, main = "Cleaned Squared Daily Returns", col = "blue")


# ------------------------------------------------------------------------------
# 1.2. Correlogram + histogram
# ------------------------------------------------------------------------------
# Correlograms of returns
par(mfrow = c(2, 2))
acf(clean_returns, main = "ACF of Returns")
pacf(clean_returns, main = "PACF of Returns")
acf(clean_returns_sq, main = "ACF of Squared Returns")
pacf(clean_returns_sq, main = "PACF of Squared Returns")

# Histogram of returns
par(mfrow = c(1, 1))

# 1. Margins: c(bottom, left, top, right) = c(5, 4, 4, 2) by default
par(mar = c(4, 4, 1, 1))
# 2. Histogram without automatic titles
chart.Histogram(clean_returns,
    methods = c("add.density", "add.normal"),
    ann = FALSE
)
title(xlab = "Cleaned Daily Returns", ylab = "Density")

# ------------------------------------------------------------------------------
# 1.3. Descriptive statistics
# ------------------------------------------------------------------------------
# descriptive stat on returns with PerformanceAnalytics package (warning: give Ex.Kurtosis)
stats <- table.Stats(clean_returns * 100)
stats
table.Distributions(clean_returns)

# Skewness
skew <- stats["Skewness", "daily.returns"] # -0.2449

# Excess Kurtosis ("Kurtosis" in the table)
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
JarqueBeraTest(clean_returns, robust = FALSE, method = "chisq")

# autocorrelation Box-Pierce and Ljung-Box tests with 10 lags
Box.test(clean_returns, lag = 10, type = c("Box-Pierce", "Ljung-Box"), fitdf = 0)

# ARCH test with 5 and 10 lags
ArchTest(clean_returns, lag = 5)
ArchTest(clean_returns, lag = 10)

# ==============================================================================
# 2. Volatility model estimation (2021-2024)
# ==============================================================================


# ------------------------------------------------------------------------------
# 2.1. Normal distribution
# ------------------------------------------------------------------------------

#----------------------------------------------------------------
# GARCH modelling
spec <- ugarchspec(variance.model = list(model = "sGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))
mod <- ugarchfit(data = clean_returns, spec = spec)
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
mod <- ugarchfit(data = clean_returns, spec = spec)
mod

#----------------------------------------------------------------
# Risk modelling
spec <- ugarchspec(variance.model = list(model = "iGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), fixed.pars = list(omega = 0, alpha1 = 0.06, beta1 = 0.94))
mod <- ugarchfit(data = clean_returns, spec = spec)
mod

#----------------------------------------------------------------
# GJR-GARCH modelling
spec <- ugarchspec(variance.model = list(model = "gjrGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))
mod <- ugarchfit(data = clean_returns, spec = spec)
mod


# ------------------------------------------------------------------------------
# 2.2. Student's t-distribution
# ------------------------------------------------------------------------------

#----------------------------------------------------------------
# GARCH modelling
spec <- ugarchspec(variance.model = list(model = "sGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std")
mod <- ugarchfit(data = clean_returns, spec = spec)
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
mod <- ugarchfit(data = clean_returns, spec = spec)
mod

#----------------------------------------------------------------
# Risk modelling
spec <- ugarchspec(variance.model = list(model = "iGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), fixed.pars = list(omega = 0, alpha1 = 0.06, beta1 = 0.94), distribution.model = "std")
mod <- ugarchfit(data = clean_returns, spec = spec)
mod

#----------------------------------------------------------------
# GJR-GARCH modelling
spec <- ugarchspec(variance.model = list(model = "gjrGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std")
mod <- ugarchfit(data = clean_returns, spec = spec)
mod

# ==============================================================================
# 3. Volatility and VaR forecasting (2025)
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1. Out-of-sample forecasting and backtesting
# ------------------------------------------------------------------------------

y <- clean_returns

# Defining the in-sample and out-of-sample boundaries
estim <- which(index(returns) == as.Date("2024-12-31")) # number of return observations from 2021 to 2024
b <- nrow(returns)
h <- b - estim # number of return observations in 2025
original <- returns[(estim + 1):b, 1] # original returns in 2025

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
plot_data <- cbind(original, varmat)
colnames(plot_data) <- c("Returns", "VaR")
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(plot_data, legend.loc = "bottomleft", main = "", col = rainbow(4))

# Kupiec and Engle-Manganelli VaR tests
print(VaRTest(0.05, as.numeric(original), as.numeric(varmat)))

#----------------------------------------------------------
# Estimating GJR(1,1) with only an intercept in the mean
spec <- ugarchspec(variance.model = list(model = "gjrGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))

for (i in 1:h)
{
    yy <- y[i:(estim - 1 + i), 1]
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
plot_data <- cbind(original, varmat)
colnames(plot_data) <- c("Returns", "VaR")
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(plot_data, legend.loc = "bottomleft", main = "", col = rainbow(4))

# Kupiec and Engle-Manganelli VaR tests
print(VaRTest(0.05, as.numeric(original), as.numeric(varmat)))

#----------------------------------------------------------
# Estimating GARCH(1,1) Student with only an intercept in the mean
spec <- ugarchspec(variance.model = list(model = "sGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std")

for (i in 1:h)
{
    yy <- y[i:(estim - 1 + i), 1]
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
plot_data <- cbind(original, varmat)
colnames(plot_data) <- c("Returns", "VaR")
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(plot_data, legend.loc = "bottomleft", main = "", col = rainbow(4))

# Kupiec and Engle-Manganelli VaR tests
print(VaRTest(0.05, as.numeric(original), as.numeric(varmat)))

#----------------------------------------------------------
# Estimating GJR(1,1) Student with only an intercept in the mean
spec <- ugarchspec(variance.model = list(model = "gjrGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), distribution.model = "std")

for (i in 1:h)
{
    yy <- y[i:(estim - 1 + i), 1]
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
plot_data <- cbind(original, varmat)
colnames(plot_data) <- c("Returns", "VaR")
options(repr.plot.res = 300, repr.plot.height = 4.4)
plot.xts(plot_data, legend.loc = "bottomleft", main = "", col = rainbow(4))

# Kupiec and Engle-Manganelli VaR tests
print(VaRTest(0.05, as.numeric(original), as.numeric(varmat)))

# ==============================================================================
# 4. Forecast evaluation via the Model Confidence Set (MCS)
# ==============================================================================

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
