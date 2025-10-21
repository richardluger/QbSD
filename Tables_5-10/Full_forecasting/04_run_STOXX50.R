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
source("../../Utilities/functions.R")       # R wrapper for QbSD_oos
sourceCpp("../../Utilities/functions.cpp")  # C++ core computation


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
gCAV  <- c("gSAV", "gAS")          # global CAViaR types: global symmetric absolute value, global asymmetric slope
CAV  <- c("SAV", "AS")             # CAViaR types: symmetric absolute value, asymmetric slope
RFO  <- c("VaRf", "ESf")           # forecasted VaR and ES
LOC  <- c("zero", "AR")            # location for other models: constant vs. AR(1)
ESD  <- c("mult", "ar")            # ES dynamics in AL_oos model
gmo  <- c("GARCH_oos-n", "GARCH_oos-t", "GARCH_oos-skewt", "GJR-n", "GJR-t", "GJR-skewt", "EGARCH_oos-n")  # GARCH_oos variants

# Initialize containers for forecast outputs 
QbSD_oos <- list()
AL_oos <- list()
GAS_oos <- list()
GARCH_oos  <- list()
 
 
################################################################################################
## Forecasting loop for one series (STOXX50)
################################################################################################

n <- 4  # index for STOXX50 
 
dat <- as.numeric(datt[[n]])
TT <- length(dat)

# Initialize nested lists for each model 
QbSD_oos[[n]] <- list()
AL_oos[[n]] <- list()
GAS_oos[[n]] <- list()
GARCH_oos[[n]]  <- list()

for (w in 1:length(rw)) {
  
  win <- rw[w]          # current rolling window size
  fp <- TT - win        # number of out-of-sample forecasts
  
  QbSD_oos[[n]][[w]]  <- array(0,c(fp, 2, 2, 2, length(ALP)),dimnames = list(NULL, RFO, gCAV, RLOC, ALP))
  AL_oos[[n]][[w]]    <- array(0,c(fp,2,2,2,2,length(ALP) ),dimnames=list(NULL,RFO,ESD,CAV,LOC,ALP ))
  GAS_oos[[n]][[w]] <- array(0, c(fp, 2, 2, length(ALP)), dimnames = list(NULL, RFO, LOC, ALP))
  GARCH_oos[[n]][[w]] <- array(0,c(fp,length(gmo),2,2,length(ALP) ),dimnames=list(NULL,gmo,RFO,LOC,ALP ))
  
  for (a in 1:length(ALP)) {
    alpha <- ALP[a]  # VaR/ES probability level
    
    # Forecasts from QbSD method
    
    for (i in 1:fp) {
      
      print(paste(bnam[n],"-",rw[w],"-",alpha,"-",i))
      
      rets <- dat[i:(i + win - 1)]  # current rolling window
      
      ##########################################################################################
      ## Quantile-based Scale Dynamics (QbSD_oos) forecasts
      ##########################################################################################
      
      for (ir in 1:2) {
        for (jr in 1:2) {
          
          results <- QbSD(rets, location=RLOC[ir], CAViaR=gCAV[jr], alpha=alpha )
          
          QbSD_oos[[n]][[w]][i,1,jr,ir,a] <- results$VaR_oos
          QbSD_oos[[n]][[w]][i,2,jr,ir,a] <- results$ES_oos - results$mu_oos 			        			            
        }
      }
      
      ##########################################################################################
      ## AL joint VaR and ES forecasts
      ##########################################################################################
      
      for (itay in 1:2) {
        for (jtay in 1:2) {
          for (stay in 1:2)  {
            
            results <- AL(rets, location=LOC[itay], CAViaR=CAV[jtay], ES=ESD[stay], alpha=alpha)	  	              			                
            AL_oos[[n]][[w]][i,,stay,jtay,itay,a] <- c(results$VaR_oos, results$ES_oos - results$mu_oos)
          }
        }
      }
      
      ##########################################################################################
      ## GAS forecasts
      ##########################################################################################
      
      for (itay in 1:2) {
        results <- GAS(rets, location = LOC[itay], alpha = alpha)
        GAS_oos[[n]][[w]][i, , itay, a] <- c(results$VaR_oos, results$ES_oos - results$mu_oos)
      }
      
      ##########################################################################################
      ## GARCH family forecasts
      ##########################################################################################
      
      for (ig in 1:2) {
        
        results <- GARCH_norm(rets, location=LOC[ig], alpha=alpha)
        NormalVaRf <- results$VaR_oos
        NormalESf  <- results$ES_oos - results$mu_oos 
        
        results <- GARCH_t(rets, location=LOC[ig], alpha=alpha)
        StudentVaRf <- results$VaR_oos
        StudentESf  <- results$ES_oos - results$mu_oos 
        
        results <- GARCH_skew_t(rets, location=LOC[ig], alpha=alpha)    
        SkewtVaRf <- results$VaR_oos
        SkewtESf  <- results$ES_oos - results$mu_oos 
        
        results <- GJR_norm(rets, location=LOC[ig], alpha=alpha)
        GJRNormalVaRf <- results$VaR_oos
        GJRNormalESf  <- results$ES_oos - results$mu_oos 
        
        results <- GJR_t(rets, location=LOC[ig], alpha=alpha)
        GJRStudentVaRf <- results$VaR_oos
        GJRStudentESf  <- results$ES_oos - results$mu_oos 
        
        results <- GJR_skew_t(rets, location=LOC[ig], alpha=alpha)
        GJRskewtVaRf <- results$VaR_oos
        GJRskewtESf  <- results$ES_oos - results$mu_oos 
        
        results <- EGARCH(rets, location=LOC[ig], alpha=alpha)
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

save(QbSD_oos,file="../Full_forecasting_output/QbSD.STOXX50.Rdata")
save(AL_oos,file="../Full_forecasting_output/AL.STOXX50.Rdata")
save(GAS_oos,file="../Full_forecasting_output/GAS.STOXX50.Rdata")
save(GARCH_oos,file="../Full_forecasting_output/GARCH.STOXX50.Rdata")
 
 
