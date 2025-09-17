rm(list=ls())

# Check if required packages are installed; install them if missing, then load
packages <- c("quantreg", "Rcpp", "RcppDist")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}


# Load R and C++ functions
source("../../Utilities/Functions.R")       # R wrapper for QbSD_oos
sourceCpp("../../Utilities/Functions.cpp")  # C++ core computation


####################################################################################################################################
## Load the data and assign series names
## 'Data.Rdata' contains returns for 8 international stock indices: S&P 500, DJIA, NASDAQ, EURO STOXX 50, FTSE 100, DAX, CAC 40, TSX
####################################################################################################################################

load("../../Data/Data.Rdata")

bnam <- names(Data)
nb <- length(bnam) 

# Define start date for each series (in year + month/day format)
bdate <- matrix(c(2002,250*10/12+18, # 10 months and 18 days into 2002
                  2002,250*10/12+18,
                  2002,250*10/12+18,
                  2002,250*10/12+18,
                  2002,250*10/12+18,
                  2002,250*10/12+18,
                  2002,250*10/12+18,
                  2002,250*10/12+18),
                nb, 2, byrow = TRUE)
dimnames(bdate)[[1]] <- bnam       
datt <- Data

################################################################################################
## Set global parameters for the forecasting experiment
################################################################################################

rw <- c(250, 1250, 2500)            	# rolling window sizes
ALP <- c(0.01, 0.025, 0.05)         	# VaR and ES confidence levels

# Settings for different methods
RLOC <- c("zero", "QAR")           # location for robust QbSD_oos: constant vs. quantile autoregression
CAV  <- c("SAV", "AS")              # CAViaR types: symmetric absolute value, asymmetric slope
RFO  <- c("VaRf", "ESf")           # forecasted VaR and ES
LOC  <- c("zero", "AR")            # location for other models: constant vs. AR(1)
ESD  <- c("basic", "AR")           # ES dynamics in Taylor_oos model
gmo  <- c("GARCH_oos-n", "GARCH_oos-t", "GARCH_oos-skewt", "GJR-n", "GJR-t", "GJR-skewt", "EGARCH_oos-n")  # GARCH_oos variants

# Initialize containers for forecast outputs 
QbSD_oos <- list()
Taylor_oos <- list()
Patton_oos <- list()
GARCH_oos  <- list()
 
 
################################################################################################
## Forecasting loop for one series (DAX)
################################################################################################

n <- 6  # index for DAX 
 
dat <- as.numeric(datt[[n]])
TT <- length(dat)

# Initialize nested lists for each model 
QbSD_oos[[n]] <- list()
Taylor_oos[[n]] <- list()
Patton_oos[[n]] <- list()
GARCH_oos[[n]]  <- list()
 
for (w in 3) {
	  
	win <- rw[w]          # current rolling window size
  fp <- TT - win        # number of out-of-sample forecasts
  
	QbSD_oos[[n]][[w]] <- list()  
	Taylor_oos[[n]][[w]] <- array(0,c(fp,2,2,2,2,length(ALP) ),dimnames=list(c(),RFO,ESD,CAV,LOC,ALP ))
	Patton_oos[[n]][[w]] <- array(0,c(fp,2,length(ALP) ),dimnames=list( c(), RFO,ALP ))
	GARCH_oos[[n]][[w]]  <- array(0,c(fp,length(gmo),2,2,length(ALP) ),dimnames=list(c(),gmo,RFO,LOC,ALP ))

 
 	for (a in 1:length(ALP)) {
		alpha <- ALP[a]  # VaR/ES probability level

		# Forecasts from QbSD method

		QbSD_oos[[n]][[w]][[a]] <- array(0,c(fp, 2, 2, 2), dimnames=list(NULL, RFO, CAV, RLOC))

 		for (i in 1:fp) {

 			print(paste(bnam[n],"-",rw[w],"-",alpha,"-",i))
  
       		rets <- dat[i:(i + win - 1)]  # current rolling window

	        ##########################################################################################
      		## Quantile-based Scale Dynamics (QbSD_oos) forecasts
      		##########################################################################################

			for (ir in 1:2) {
    			for (jr in 1:2) {
    				     	
      				results <- QbSD(rets, location=RLOC[ir], CAViaR=CAV[jr], alpha=alpha )

					QbSD_oos[[n]][[w]][[a]][i,1,jr,ir] <- results$VaR_oos
					QbSD_oos[[n]][[w]][[a]][i,2,jr,ir] <- results$ES_oos  			        			            
			    }
			}

 
       		##########################################################################################
      		## Taylor_oos (2019) joint VaR and ES models
      		##########################################################################################

			for (itay in 1:2) {
				for (jtay in 1:2) {
			        for (stay in 1:2)  {
          	
	  	              	results <- Taylor(rets, location=LOC[itay], CAViaR=CAV[jtay], ES=ESD[stay], alpha)	  	              			                
		                Taylor_oos[[n]][[w]][i,,stay,jtay,itay,a] <- c(results$VaR_oos,results$ES_oos)
		           }
			   }
			}


	        ##########################################################################################
      		## Patton_oos et al. (2019) approach
      		##########################################################################################

			results <- Patton(rets, alpha)
			Patton_oos[[n]][[w]][i,,a] <- c(results$VaR_oos, results$ES_oos)

     		##########################################################################################
      		## GARCH family with various innovations and location dynamics
      		##########################################################################################

			for (ig in 1:2) {
  
				results <- GARCHnorm(rets, location=LOC[ig], alpha)
				NormalVaRf <- results$VaR_oos
				NormalESf  <- results$ES_oos - results$mu_oos 

				results <- GARCHt(rets, location=LOC[ig], alpha)
				StudentVaRf <- results$VaR_oos
				StudentESf  <- results$ES_oos - results$mu_oos 

				results <- Hansen(rets, location=LOC[ig], alpha)    
			  	SkewtVaRf <- results$VaR_oos
			  	SkewtESf  <- results$ES_oos - results$mu_oos 
		
				results <- GJRnorm(rets, location=LOC[ig], alpha)
				GJRNormalVaRf <- results$VaR_oos
				GJRNormalESf  <- results$ES_oos - results$mu_oos 

				results <- GJRt(rets, location=LOC[ig], alpha)
				GJRStudentVaRf <- results$VaR_oos
				GJRStudentESf  <- results$ES_oos - results$mu_oos 

				results <- GJRskewt(rets, location=LOC[ig], alpha)
				GJRskewtVaRf <- results$VaR_oos
				GJRskewtESf  <- results$ES_oos - results$mu_oos 

				results <- EGARCHnorm(rets, location=LOC[ig], alpha)
				EGARCHNormalVaRf <- results$VaR_oos
				EGARCHNormalESf  <- results$ES_oos - results$mu_oos 

				GARCH_oos[[n]][[w]][i,1,,ig,a] <- c(as.numeric(NormalVaRf),as.numeric(NormalESf))	
				GARCH_oos[[n]][[w]][i,2,,ig,a] <- c(as.numeric(StudentVaRf),as.numeric(StudentESf) )
				GARCH_oos[[n]][[w]][i,3,,ig,a] <- c(as.numeric(SkewtVaRf),as.numeric(SkewtESf) )
				GARCH_oos[[n]][[w]][i,4,,ig,a] <- c(as.numeric(GJRNormalVaRf),as.numeric(GJRNormalESf))
				GARCH_oos[[n]][[w]][i,5,,ig,a] <- c(as.numeric(GJRStudentVaRf),as.numeric(GJRStudentESf))
				GARCH_oos[[n]][[w]][i,6,,ig,a] <- c(as.numeric(GJRskewtVaRf),as.numeric(GJRskewtESf))
				GARCH_oos[[n]][[w]][i,7,,ig,a] <- c(as.numeric(EGARCHNormalVaRf),as.numeric(EGARCHNormalESf))
				
			}			
		}
	}
}


################################################################################################
## Save forecast objects to RData files
################################################################################################

save(QbSD_oos,file="../Full_forecasting_output/QbSD.DAX.Rdata")
save(Taylor_oos,file="../Full_forecasting_output/Taylor.DAX.Rdata")
save(Patton_oos,file="../Full_forecasting_output/Patton.DAX.Rdata")
save(GARCH_oos,file="../Full_forecasting_output/GARCH.DAX.Rdata")
 
 
