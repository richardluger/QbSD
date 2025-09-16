# Illustrative example: QbSD Value-at-Risk (VaR) and Expected Shortfall (ES) for the S&P 500.

# Clean the environment
rm(list = ls())


# Check if required packages are installed; install them if missing, then load
packages <- c("quantreg", "Rcpp", "RcppDist")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}


# Load R and C++ functions
source("../Utilities/Functions.R")       # R wrapper for QbSD
sourceCpp("../Utilities/Functions.cpp")  # C++ core computation


########################################################################################
## Load the data and assign series names
## 'Data.Rdata' contains returns for 8 international stock indices and 6 US bank stocks:
## Stock indices: S&P 500, DJIA, NASDAQ, EURO STOXX 50, FTSE 100, DAX, CAC 40, TSX
## Banks: JPM, BAC, CITI, WFC, GS, MS
########################################################################################

load("../Data/Data.Rdata")
series_names <- c("S&P 500", "DJIA", "NASDAQ", "EURO STOXX 50", "FTSE 100", "DAX", "CAC 40", "TSX", "JPM", "BAC", "CITI", "WFC", "GS", "MS")
names(Data) <- series_names


#####################################################################################
## Extract one series and apply the QbSD method over the full sample for illustration
#####################################################################################

series <- "S&P 500"
rets <- as.numeric(Data[[series]])
dates <- as.Date(names(Data[[series]]))


# Set risk level (alpha) for VaR and ES forecasts
alpha <- 0.01


# Choose QbSD model configuration
location_type <- "QAR"	# or "zero"
CAViaR_type   <- "AS"	# or "SAV"


# Compute in-sample estimates and OOS forecasts
results <- QbSD(rets, location=location_type, CAViaR=CAViaR_type, alpha=alpha)


# Extract in-sample VaR and ES estimates
VaR_in <- results$VaR_in
ES_in  <- results$ES_in


# Extract OOS forecasts (1-step-ahead forecasts at the end of the sample)
VaR_oos <- results$VaR_oos
ES_oos  <- results$ES_oos


# Plot returns with in-sample VaR and ES
ylim_range <- range(c(rets, VaR_in, ES_in), na.rm = TRUE)

plot(dates, rets, type = "l", col = "black", lwd = 1.0, main = paste(series, "returns with in-sample VaR and ES at", alpha * 100, "%"), xlab = "", ylab = "", ylim = ylim_range, axes = TRUE)

lines(dates, VaR_in, col = "#D55E00", lwd = 1.2, lty = 2)  
lines(dates, ES_in,  col = "#0072B2", lwd = 1.2, lty = 3)  
abline(h = 0, col = "gray40")

legend("bottomleft", inset = 0.01, legend = c("Returns", paste0("VaR (", alpha * 100, "%)"), paste0("ES (", alpha * 100, "%)")), col = c("black", "#D55E00", "#0072B2"), lty = c(1, 2, 3), lwd = c(1.0, 1.2, 1.2), bg = "white", cex = 0.85, bty = "n")




