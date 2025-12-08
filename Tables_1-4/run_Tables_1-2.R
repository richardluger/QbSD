# Clean the environment
rm(list=ls())  

# Check if required packages are installed; install them if missing, then load
packages <- c("quantreg", "Rcpp", "RcppDist")

for (p in packages){
  if (!requireNamespace(p, quietly = TRUE)){
    install.packages(p)
  }
  library(p, character.only = TRUE)
}


# Load R and C++ functions
source("../Utilities/functions.R")
sourceCpp("../Utilities/functions.cpp")	

	
filename <- 'output1-2.txt'             # plain-text output file

totrep <- 1000                         	# Monte Carlo replications
Tsize <- 1250                           # sample length 

v_set <- c(20, 5)						# df of Hansen (1994) skew-t: 20 (light), 5 (heavy)
lambda_set <- c(0, -0.5)				# skewness: 0 (sym), -0.5 (left-skew)
alpha_set <- c(0.01, 0.025, 0.05)		# risk levels for VaR/ES: 0.01, 0.025, 0.05


for (v in v_set) for (lambda in lambda_set) for (alpha in alpha_set){

	
true_VaR <- numeric(totrep)

QbSD_gSAV_VaR_oos <- numeric(totrep)  
QbSD_gAS_VaR_oos  <- numeric(totrep) 

AL_mult_SAV_VaR_oos <- numeric(totrep) 
AL_ar_SAV_VaR_oos   <- numeric(totrep) 
AL_mult_AS_VaR_oos  <- numeric(totrep) 
AL_ar_AS_VaR_oos    <- numeric(totrep)

GAS_VaR_oos <- numeric(totrep)

GARCH_norm_VaR_oos   <- numeric(totrep)
GARCH_t_VaR_oos      <- numeric(totrep)
GARCH_skew_t_VaR_oos <- numeric(totrep)

GJR_norm_VaR_oos   <- numeric(totrep)
GJR_t_VaR_oos      <- numeric(totrep)
GJR_skew_t_VaR_oos <- numeric(totrep)

EGARCH_VaR_oos <- numeric(totrep)


true_ES <- numeric(totrep)

QbSD_gSAV_ES_oos <- numeric(totrep) 
QbSD_gAS_ES_oos  <- numeric(totrep) 

AL_mult_SAV_ES_oos <- numeric(totrep) 
AL_ar_SAV_ES_oos   <- numeric(totrep) 
AL_mult_AS_ES_oos  <- numeric(totrep)
AL_ar_AS_ES_oos    <- numeric(totrep)

GAS_ES_oos <- numeric(totrep)

GARCH_norm_ES_oos   <- numeric(totrep)
GARCH_t_ES_oos      <- numeric(totrep)
GARCH_skew_t_ES_oos <- numeric(totrep)

GJR_norm_ES_oos   <- numeric(totrep)
GJR_t_ES_oos      <- numeric(totrep)
GJR_skew_t_ES_oos <- numeric(totrep)

EGARCH_ES_oos <- numeric(totrep)



set.seed(123456)
 
for (irep in 1:totrep){ 
					
	print(irep)
  	
	omega <- 0.05
	beta  <- 0.85
	gamma <- 0.10
	delta <- 1.5
	theta <- 0

	burnin <- 500
	TT <- Tsize + burnin

	sigma <- numeric(TT)
	r     <- numeric(TT)

	# crude, stable-ish init (burn-in will wash this out)
	sigma[1] <- ( omega / (1 - min(0.99, gamma + beta)) )^(1/delta)

	epsilon <- vapply(runif(TT), Hansen_qt, numeric(1), v = v, lambda = lambda)

	r[1] <- sigma[1] * epsilon[1]

	for (t in 2:TT){
  		sigma[t] <- ( omega + beta * sigma[t-1]^delta + gamma * ( abs(r[t-1]) - theta * r[t-1] )^delta  )^(1/delta)
		r[t] <- sigma[t] * epsilon[t]
	}

	sigma <- sigma[(burnin + 1):TT]
	r <- r[(burnin + 1):TT]

	true_VaR[irep] <- sigma[Tsize] * Hansen_qt(alpha, v, lambda)
	true_ES[irep]  <- sigma[Tsize] * Hansen_ES(alpha, v, lambda)$ES

	returns <- r[1:(Tsize - 1)]

 
	## QbSD

	results <- QbSD(returns, location="zero", CAViaR="gSAV", alpha)                        
	QbSD_gSAV_VaR_oos[irep] <- results$VaR_oos	
	QbSD_gSAV_ES_oos[irep]  <- results$ES_oos    

	results <- QbSD(returns, location="zero", CAViaR="gAS", alpha)        
	QbSD_gAS_VaR_oos[irep] <- results$VaR_oos
	QbSD_gAS_ES_oos[irep]  <- results$ES_oos    


	## AL 
			
	results <- AL(returns, location="zero", CAViaR="SAV", ES="mult", alpha)              
	AL_mult_SAV_VaR_oos[irep] <- results$VaR_oos
	AL_mult_SAV_ES_oos[irep]  <- results$ES_oos

	results <- AL(returns, location="zero", CAViaR="SAV", ES="ar", alpha)                                         
	AL_ar_SAV_VaR_oos[irep] <- results$VaR_oos
	AL_ar_SAV_ES_oos[irep]  <- results$ES_oos

	results <- AL(returns, location="zero", CAViaR="AS", ES="mult", alpha)                            
	AL_mult_AS_VaR_oos[irep] <- results$VaR_oos
	AL_mult_AS_ES_oos[irep]  <- results$ES_oos
		
	results <- AL(returns, location="zero", CAViaR="AS", ES="ar", alpha)                            
	AL_ar_AS_VaR_oos[irep] <- results$VaR_oos
	AL_ar_AS_ES_oos[irep]  <- results$ES_oos
		    

	## GAS 
				
	results <- GAS(returns, location="zero", alpha)
	GAS_VaR_oos[irep] <- results$VaR_oos
	GAS_ES_oos[irep]  <- results$ES_oos

				
	## GARCH family 
        
	results <- GARCH_norm(returns, location="zero", alpha)  	# GARCH + normal
	GARCH_norm_VaR_oos[irep] <- results$VaR_oos
	GARCH_norm_ES_oos[irep]  <- results$ES_oos

	results <- GARCH_t(returns, location="zero", alpha)     	# GARCH + t
	GARCH_t_VaR_oos[irep] <- results$VaR_oos
	GARCH_t_ES_oos[irep]  <- results$ES_oos

	results <- GARCH_skew_t(returns, location="zero", alpha)  	# GARCH + skew-t
	GARCH_skew_t_VaR_oos[irep] <- results$VaR_oos
	GARCH_skew_t_ES_oos[irep]  <- results$ES_oos

	results <- GJR_norm(returns, location="zero", alpha)     	# GJR + normal
	GJR_norm_VaR_oos[irep] <- results$VaR_oos
	GJR_norm_ES_oos[irep]  <- results$ES_oos

	results <- GJR_t(returns, location="zero", alpha)        	# GJR + t
	GJR_t_VaR_oos[irep] <- results$VaR_oos
	GJR_t_ES_oos[irep]  <- results$ES_oos

	results <- GJR_skew_t(returns, location="zero", alpha)    	# GJR + skew-t
	GJR_skew_t_VaR_oos[irep] <- results$VaR_oos
	GJR_skew_t_ES_oos[irep]  <- results$ES_oos				

	results <- EGARCH(returns, location="zero", alpha)  		# EGARCH + normal
	EGARCH_VaR_oos[irep] <- results$VaR_oos
	EGARCH_ES_oos[irep]  <- results$ES_oos
											
}	


h1   <- function(title) cat(title, "\n", file = filename, append = TRUE)
rule <- function(txt)   cat(paste0("\n", txt, "\n"), file = filename, append = TRUE)

.pad_right <- function(x, width){
  w <- nchar(x, type = "width", allowNA = FALSE, keepNA = FALSE)
  paste0(x, strrep(" ", pmax.int(0L, width - w)))
}

w_align_block <- function(labels, values, header = NULL){
  if (!is.null(header)) rule(header)
  colw  <- max(nchar(labels, type = "width"))
  lab   <- .pad_right(labels, colw)
  lines <- sprintf("%s  %12.6f", lab, values)
  cat(paste0(lines, collapse = "\n"), "\n", file = filename, append = TRUE)
}

w_block_mae <- function(labels, series, truth, header){
  mae  <- function(x) mean(abs(x - truth))
  vals <- vapply(series, mae, numeric(1), USE.NAMES = FALSE)
  w_align_block(labels, vals, header)
}

w_block_rmse <- function(labels, series, truth, header){
  rmse <- function(x) sqrt(mean((x - truth)^2))
  vals <- vapply(series, rmse, numeric(1), USE.NAMES = FALSE)
  w_align_block(labels, vals, header)
}

h1("==== DGP =========================================================")
cat(sprintf("Tsize = %d   v = %.4f   lambda = %.4f   alpha = %.4f   theta = %.4f\n\n", Tsize, v, lambda, alpha, theta), file = filename, append = TRUE)

labels <- c(
  "QbSD — gSAV",
  "QbSD — gAS",
  "AL — SAV (mult)",
  "AL — SAV (AR)",
  "AL — AS (mult)",
  "AL — AS (AR)",
  "GAS",
  "GARCH (Normal)",
  "GARCH (t)",
  "GARCH (skew-t)",
  "GJR (Normal)",
  "GJR (t)",
  "GJR (skew-t)",
  "EGARCH (Normal)"
)

series_VaR <- list(
  QbSD_gSAV_VaR_oos,
  QbSD_gAS_VaR_oos,
  AL_mult_SAV_VaR_oos,
  AL_ar_SAV_VaR_oos,
  AL_mult_AS_VaR_oos,
  AL_ar_AS_VaR_oos,
  GAS_VaR_oos,
  GARCH_norm_VaR_oos,
  GARCH_t_VaR_oos,
  GARCH_skew_t_VaR_oos,
  GJR_norm_VaR_oos,
  GJR_t_VaR_oos,
  GJR_skew_t_VaR_oos,
  EGARCH_VaR_oos
)

series_ES <- list(
  QbSD_gSAV_ES_oos,
  QbSD_gAS_ES_oos,
  AL_mult_SAV_ES_oos,
  AL_ar_SAV_ES_oos,
  AL_mult_AS_ES_oos,
  AL_ar_AS_ES_oos,
  GAS_ES_oos,
  GARCH_norm_ES_oos,
  GARCH_t_ES_oos,
  GARCH_skew_t_ES_oos,
  GJR_norm_ES_oos,
  GJR_t_ES_oos,
  GJR_skew_t_ES_oos,
  EGARCH_ES_oos
)


w_block_mae (labels, series_VaR, true_VaR, "VaR ************* MAE ************************************** VaR")
w_block_rmse(labels, series_VaR, true_VaR, "VaR ************* RMSE ************************************* VaR")

w_block_mae (labels, series_ES, true_ES, "ES ************* MAE *************************************** ES")
w_block_rmse(labels, series_ES, true_ES, "ES ************* RMSE ************************************** ES")
					 
cat(c("\n"), file=filename, append=TRUE)
cat(c("\n"), file=filename, append=TRUE)
cat(c("\n"), file=filename, append=TRUE)
			 


}
	
	
			 



