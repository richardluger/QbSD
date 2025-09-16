# QbSD
**Q**uantile-**b**ased **S**cale **D**ynamics

This repository contains the code accompanying the paper  
*“Quantile-based modeling of scale dynamics in financial returns for Value-at-Risk and Expected Shortfall forecasting”*  
by **Xiaochun Liu** and **Richard Luger**.

Assembled on: 15 September 2025

---

## Repository Structure

### `Data/`
* **Data.Rdata** – Publicly available daily returns for 8 international stock indices and 6 U.S. bank stocks, obtained from Yahoo Finance.  
  *Stock indices*: S&P 500, DJIA, NASDAQ, EURO STOXX 50, FTSE 100, DAX, CAC 40, TSX  
  *Banks*: JPM, BAC, CITI, WFC, GS, MS  
  **Note**: Only the index returns are used in the paper.  


### `Example_QbSD/`
* **Run_example.R** – Illustrative example: QbSD Value-at-Risk (VaR) and Expected Shortfall (ES) for the S&P 500.

### `Tables_1-4/`
* **Run_simulations.R** – Reproduces the simulation results reported in Tables 1–4.

### `Tables_5-10/`
Contains the code and outputs used to generate the forecasting results reported in Tables 5–10.  
Key subfolders:
* **Full_forecasting/** – Driver scripts (`01_run_SP500.R` … `08_run_TSX.R`) to run the full forecasting experiment for each index.  
* **Full_forecasting_output/** – **Intermediate generated datasets**: saved forecast objects produced by the scripts in `Full_forecasting/`. These objects are included so users can reproduce the tables without re-running the full forecasting experiment.  
* **Gen_tables/** – Script (`Run_Tables_5-10.R`) to process outputs and generate the summary tables.

### `Utilities/`
Shared helper code used across the project.  
* **Functions.R** – R helper functions and wrappers (e.g., the `QbSD` interface).  
* **Functions.cpp** – C++ implementations of core routines (called from R via Rcpp).


## Requirements

### Software
- **R** (version 4.0 or higher recommended; **tested with 4.5.1**)

### R packages  
CRAN packages required (**tested versions in parentheses**):  
`quantreg` (**6.1**), `Rcpp` (**1.1.0**), `RcppDist` (**0.1.1.1**), `MCS` (**0.1.3**)



Install with: 
```r
install.packages(c("quantreg", "Rcpp", "RcppDist", "MCS")) 
```

### C++ toolchain
Compiled code via `Rcpp::sourceCpp()` requires:
- **Windows**: Rtools  
- **macOS**: Xcode Command Line Tools (`xcode-select --install`)  
- **Linux**: GNU build tools (e.g., `sudo apt install build-essential`)


## Runtime

**Computer type**: MacBook Pro   
**CPU**: Apple M3 Max, 14-core CPU  
**Memory**: 96 GB unified RAM

**Measured runtimes**
- Illustrative example (`Example_QbSD/Run_example.R`): **< 1 minute**
- Monte Carlo simulations (for a single DGP configuration: **T = 2500, v = 20, λ = 0, α = 0.01; 1000 reps**): **~1.6 hours**
- Full forecasting results for a **single index**: **~67.5 hours** (**~2.8 days**)
- Generate tables (`Tables_5-10/Gen_tables/Run_Tables_5-10.R`): **~1 hour**




## Contact
Richard Luger — <richard.luger@fsa.ulaval.ca>
