# Conditional Volatility Modeling and Risk Measurement (VaR): Application to the SSE Composite Index
*This project models conditional volatility and forecasts market risk on the Shanghai Stock Exchange Composite Index.*

---

## 🎯 Overview

This project applies GARCH-family econometric models to the daily returns of the SSE Composite Index (2021–2025), with the dual objective of estimating conditional volatility in-sample and forecasting Value-at-Risk (VaR) and Expected Shortfall (ES) out-of-sample.

**Objectives**
- Characterize the distributional and autocorrelation properties of SSE Composite Index daily returns
- Estimate and compare eight GARCH-type specifications under normal and Student-*t* distributional assumptions
- Forecast one-step-ahead conditional VaR (5%) and Expected Shortfall over the 2025 out-of-sample horizon
- Evaluate forecast accuracy via the Model Confidence Set (MCS) and validate capital coverage through Kupiec backtesting

---

## 🗄️ Data

- **Source:** Yahoo Finance (`quantmod::getSymbols`, ticker `000001.SS`)
- **Time Period / Size:** January 1, 2021 – December 31, 2025, daily frequency; 1,211 adjusted observations (in-sample: 2021–2024)
- **Target Variable:** Daily log-returns of the SSE Composite Index closing price
- **Key Predictors / Features:** Lagged squared returns (ARCH component), lagged conditional variance (GARCH component), asymmetry indicator for negative shocks (GJR leverage parameter γ)
- **Preprocessing:** Outlier detection and winsorization via `Return.clean` (Boudt method); three extreme observations corrected (2024-09-30: −0.0514, 2024-10-09: +0.0370, 2025-04-07: +0.0443)
- **Data Availability:** Publicly available via Yahoo Finance API; retrieved programmatically through `quantmod`

---

## 🧠 Methodology

- **Theoretical Approach:** Conditional heteroskedasticity modeling via GARCH(1,1), IGARCH(1,1), RiskMetrics (EWMA with fixed parameters α = 0.06, β = 0.94), and GJR-GARCH(1,1), estimated under both Gaussian and Student-*t* distributional assumptions
- **Mathematical Framework:** Maximum likelihood estimation of the conditional variance equation σ²<sub>t</sub> = ω + αε²<sub>t−1</sub> + γε²<sub>t−1</sub>𝟙<sub>{ε<sub>t−1</sub><0}</sub> + βσ²<sub>t−1</sub>; model selection via log-likelihood, AIC, and Hannan-Quinn (HQ) criteria; residual validation through Ljung-Box Q(5), Q²(5), and LM-ARCH(5) tests
- **Evaluation Strategy:** Rolling-window out-of-sample forecasting (h = 1 step ahead) over 2025; MSE and R²<sub>OOS</sub> for variance accuracy; Diebold-Mariano (DM) test and Model Confidence Set (MCS, α = 0.10, B = 5,000 bootstrap replications, T<sub>max</sub> statistic) for equal predictive ability; Kupiec unconditional coverage (POF) test for VaR backtesting at the 5% level

---

## ⚙️ Features

- **Detect and Neutralize Outliers:** Apply the Boudt winsorization method to remove extreme return observations prior to estimation
- **Estimate GARCH-Family Models:** Fit eight distinct volatility specifications under two distributional assumptions via maximum likelihood using `rugarch`
- **Compare In-Sample Specifications:** Rank valid models by log-likelihood, AIC, and HQ to identify the optimal conditional variance architecture
- **Forecast Conditional VaR and ES:** Generate rolling one-step-ahead 5% VaR and Expected Shortfall forecasts over the 2025 out-of-sample period
- **Backtest Risk Projections:** Validate regulatory capital coverage through the Kupiec unconditional coverage test and count empirical VaR exceptions
- **Evaluate Forecast Accuracy:** Apply the Diebold-Mariano test and the MCS procedure to formally test equal predictive ability across model candidates

---

## 🧰 Tech Stack

- **Language:** R 4.5.3
- **Econometrics & Statistical Inference:** DescTools, FinTS, robustbase
- **Time Series Analysis:** forecast, rugarch, MCS
- **Quantitative Finance:** quantmod, PerformanceAnalytics
- **Reporting & Documentation:** LaTeX (scrreprt)
- **MLOps & Infrastructure:** renv 1.2.0

---

## 📦 Installation

```bash
git clone https://github.com/floriancrochet/master-year2-advanced-financial-econometrics.git
cd master-year2-advanced-financial-econometrics
Rscript -e 'install.packages(c("quantmod", "PerformanceAnalytics", "robustbase", "DescTools", "FinTS", "rugarch", "forecast", "MCS"))'
```

---

## 💻 Usage Example

### Reproducing the Analysis / Execution Pipeline

```r
library(quantmod)
library(rugarch)
library(MCS)

# 1. Download SSE Composite Index data
getSymbols("000001.SS", src = "yahoo", from = as.Date("2021-01-01"), to = as.Date("2025-12-31"))
price        <- `000001.SS`[, 4]
returns      <- dailyReturn(price)
clean_returns <- PerformanceAnalytics::Return.clean(returns, method = "boudt")

# 2. Estimate GJR-GARCH(1,1) Student-t (optimal in-sample model)
spec <- ugarchspec(
  variance.model  = list(model = "gjrGARCH"),
  mean.model      = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
mod <- ugarchfit(data = clean_returns, spec = spec)
show(mod)

# 3. Run full pipeline (OOS forecasting + MCS) via project.R
source("project.R")
```

---

## 📂 Project Structure

```text
master-year2-advanced-financial-econometrics/
│
├── report/
│   └── report.pdf
├── .Rprofile
├── .gitignore
├── LICENSE
├── README.md
├── master-year2-advanced-financial-econometrics.Rproj
├── project.R                                            # GARCH/GJR-GARCH estimation and VaR forecasting pipeline
├── renv.lock
```

---

## 📈 Results

### Performance Metrics

| Model | Distribution | Persistence | MSE (×10⁻⁵) | Exceptions (OOS) | Kupiec p-value |
|---|---|---|---|---|---|
| GARCH(1,1) | Normal | 0.9503 | 0.0127 | 7 | 0.1040 |
| GJR-GARCH(1,1) | Normal | 0.9447 | 0.0127 | 8 | 0.1990 |
| GARCH(1,1) | Student-*t* | 0.9539 | 0.0127 | 7 | 0.1040 |
| **GJR-GARCH(1,1)** | **Student-*t*** | **0.9451** | **0.0127** | **8** | **0.1990** |

### Key Findings

- **Optimal In-Sample Model:** GJR-GARCH(1,1) under the Student-*t* distribution achieved the highest log-likelihood and strictly minimized both AIC and HQ criteria, driven by the significance of the leverage parameter γ (negative shocks generate higher variance) and the shape parameter (leptokurtic tail absorption).
- **Leverage Effect Confirmed:** The statistically significant γ coefficient in GJR specifications confirms that negative return shocks on the SSE Index generate greater conditional variance than positive shocks of equal magnitude, consistent with empirical leverage effect theory.
- **High Volatility Persistence:** All valid GARCH specifications display persistence estimates between 0.9447 and 0.9539, indicating that volatility shocks decay slowly and that market risk is highly autocorrelated over time.
- **OOS Forecast Homogeneity:** All four validated models produce a rigorously identical MSE of 0.0127 × 10⁻⁵ over 2025; neither the DM test (p-values: 0.0554–0.8931) nor the MCS procedure (p-values: 0.7486–1.0000) were able to statistically discriminate among them.
- **Regulatory Capital Adequacy:** With 7–8 empirical exceptions against a theoretical expectation of 12 (at 5%), all models systematically overestimate volatility out-of-sample; Kupiec p-values (0.1040–0.1990) confirm unconditional coverage cannot be rejected at the 5% level.
- **Exogenous Shock Limitation:** The April 7, 2025 tariff-induced crash (daily return below −7%) breached the 5% VaR threshold across all models simultaneously, exposing the structural inability of autoregressive GARCH smoothing to anticipate sudden geopolitical information shocks.
- **Expected Shortfall Sensitivity:** The mean ES (approximately −1.6600%, ranging from −1.6500% to −1.6700% across specifications) is strongly inflated in absolute value by the isolated April 7, 2025 event, underscoring the high sensitivity of empirical ES to exogenous tail observations.

---

## 🚧 Limitations & Future Work

- **Parametric Shock Blindness:** GARCH-family models cannot anticipate or adequately quantify the magnitude of unexpected geopolitical shocks (e.g., sudden tariff escalations); coupling with historical or Monte Carlo stress tests is required.
- **Post-Shock Adjustment Inertia:** The VaR threshold adjusts only moderately following extreme events (dropping to approximately −2% after the April 2025 crash), reflecting GARCH smoothing limitations under sudden regime changes.
- **RiskMetrics Misspecification:** Fixed EWMA parameters (α = 0.06, β = 0.94) fail to filter residual ARCH effects specific to the SSE Index, as evidenced by Q²(5) p-values of 0.0430–0.0450.
- **OOS Distributional Indifference:** The Student-t assumption, which proves superior in-sample, does not translate into measurable predictive gains in variance forecasting at the 5% horizon, limiting its out-of-sample justification beyond regulatory compliance.
- **Integrate Stress Testing:** Complement GARCH forecasts with historical or Monte Carlo stress scenarios to guard against unpredictable extreme shocks
- **Extend to Conditional Coverage:** Apply the Christoffersen conditional independence test to detect clustering of VaR exceptions over time
- **Explore Semi-Parametric Models:** Evaluate Extreme Value Theory (EVT) or filtered historical simulation approaches to improve tail risk estimation beyond parametric GARCH assumptions

---

## 📜 License

This project is released under the MIT License.  
© 2026 Florian Crochet

---

## 👤 Author

**Florian Crochet**  
[GitHub Profile](https://github.com/floriancrochet)

*Master 2 – Econometrics and Statistics, Applied Econometrics Track*  
*Nantes Université – IAE Nantes*

---

## 🤝 Acknowledgments

This work was conducted as part of the Advanced Financial Econometrics course, supervised by Mr. Olivier Darné.