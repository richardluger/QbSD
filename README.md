# QbSD
**Q**uantile-**b**ased **S**cale **D**ynamics

This repository contains the code accompanying the paper  
*“Quantile-based modeling of scale dynamics in financial returns for Value-at-Risk and Expected Shortfall forecasting”*  
by **Xiaochun Liu** and **Richard Luger**.

Updated on: December 8, 2025

---

## Repository Structure


### `Data/`
* **Data.Rdata** – Daily return series for 8 international stock indices, publicly available from Yahoo Finance.  
  *Indices:* S&P 500, DJIA, NASDAQ, EURO STOXX 50, FTSE 100, DAX, CAC 40, TSX



### `Example_QbSD/`
* **run_example.R** – Illustrative example: QbSD Value-at-Risk (VaR) and Expected Shortfall (ES) for the S&P 500.


### `Tables_1-4/`
* **run_Tables_1-2.R** and **run_Tables_3-4.R** – Scripts to reproduce the simulation results for Tables 1–4.

### `Tables_5-10/`
Contains the code and outputs used to generate the forecasting results reported in Tables 5–10 and in the Supplementary Tables C1–C12.  
Key subfolders:
* **Full_forecasting/** – Driver scripts (`01_run_SP500.R` … `08_run_TSX.R`) to run the full forecasting experiment for each index.  
* **Full_forecasting_output/** – **Intermediate generated datasets**: saved forecast objects produced by the scripts in `Full_forecasting/`. These objects are included so users can reproduce the tables without re-running the full forecasting experiment.  
* **Gen_tables/** – Scripts (`run_Tables_5-10.R`, `run_Tables_C1-C6.R`, `run_Tables_C7-C12.R`) that process outputs and generate the summary tables.

### `Utilities/`
Shared helper code used across the project.  
* **functions.R** – R helper functions and wrappers (e.g., the `QbSD` interface).  
* **functions.cpp** – C++ implementations of core routines (called from R via Rcpp).


> **Note on file paths**  
> All scripts use **relative paths** that assume the same folder hierarchy as this repository.  To run the code without editing the `load()`, `source()`, or `sourceCpp()` commands,  clone or download the repository and keep the directory structure exactly as shown above. Alternatively, adjust the paths in the scripts if you choose a different working directory.


## Requirements

### Software
- **R** (version 4.0 or higher recommended; **tested with 4.5.0**)

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

All experiments were run on a High-Performance Computing (HPC) system of the Digital Research Alliance of Canada (Rorqual cluster).

**Node type**: Dell Inc. PowerEdge R6625  
**CPUs**: 2 × AMD EPYC 9654 (“Genoa”), 96 cores each (192 cores per node)  
**Memory**: ~755 GB RAM per node




**Measured runtimes**
- Illustrative example (`Example_QbSD/run_example.R`): **< 1 min**
- Monte Carlo simulations (for a single DGP configuration: **T = 1250, v = 20, λ = 0, α = 0.01; 1000 reps**): **~2 h**
- Full forecasting results for a single index: **~5 days 16 h**
- Generate tables (`Tables_5-10/Gen_tables/run_Tables_5-10.R`): **~1 h 45 min**


## Contact
Richard Luger — <richard.luger@fsa.ulaval.ca>
