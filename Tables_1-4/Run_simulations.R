# Clean the environment
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
source("../Utilities/Functions.R")
sourceCpp("../Utilities/Functions.cpp")	

	
filename <- 'output.txt'                # plain-text output file

totrep <- 1000                          # Monte Carlo replications
TT <- 2500                              # sample length 

v <- 20									# df of Hansen (1994) skew-t: 20 (light), 5 (heavy)
lambda <- 0								# skewness: 0 (sym), -0.5 (left-skew)
alpha <- 0.01							# risk levels for VaR/ES: 0.01, 0.025, 0.05

	
TrueVaR <- matrix(0, totrep,1) 
TrueES <- matrix(0, totrep,1) 

SAVRobustVaRfmean <- matrix(0, totrep,1) 
SAVRobustESfmean <- matrix(0, totrep,1) 

ASRobustVaRfmean <- matrix(0, totrep,1) 
ASRobustESfmean <- matrix(0, totrep,1) 

SAVTaylorVaRf.basic <- matrix(0, totrep,1) 
SAVTaylorESf.basic <- matrix(0, totrep,1) 

SAVTaylorVaRf.AR <- matrix(0, totrep,1) 
SAVTaylorESf.AR <- matrix(0, totrep,1) 

ASTaylorVaRf.basic <- matrix(0, totrep,1) 
ASTaylorESf.basic <- matrix(0, totrep,1) 

ASTaylorVaRf.AR <- matrix(0, totrep,1) 
ASTaylorESf.AR <- matrix(0, totrep,1) 

PattonVaRf <- matrix(0, totrep,1) 
PattonESf  <- matrix(0, totrep,1) 

NormalVaRf <- matrix(0, totrep,1) 
NormalESf <- matrix(0, totrep,1) 

StudentVaRf <- matrix(0, totrep,1) 
StudentESf <- matrix(0, totrep,1) 

HansenVaRf <- matrix(0, totrep,1) 
HansenESf <- matrix(0, totrep,1) 

EGARCHNormalVaRf <- matrix(0, totrep,1) 
EGARCHNormalESf <- matrix(0, totrep,1) 

GJRNormalVaRf <- matrix(0, totrep,1) 
GJRNormalESf <- matrix(0, totrep,1) 

GJRStudentVaRf <- matrix(0, totrep,1) 
GJRStudentESf <- matrix(0, totrep,1) 

GJRskewtVaRf <- matrix(0, totrep,1) 
GJRskewtESf <- matrix(0, totrep,1) 

set.seed(123456)
 
for (irep in 1:totrep){ 
					
	print(irep)
  	
	y <- h <- epsilon <-  matrix(0, TT, 1)

	omegahat <- 0.05
	alphahat <- 0.10
	betahat  <- 0.85
	deltahat <- 1.5	
	gammahat <- 0

	h[1] <- omegahat/(1-alphahat-betahat)
	
	epsilon[1] <- Hansen_qt(runif(1), v, lambda)	
	y[1] <- h[1]* epsilon[1]

	for (t in 2:TT){	
				
		h[t] <- ( omegahat + alphahat*(abs(y[t-1]) - gammahat*y[t-1])^deltahat + betahat*h[t-1]^deltahat )^(1/deltahat)

		epsilon[t] <- Hansen_qt(runif(1), v, lambda)
		
		y[t] <- h[t]*epsilon[t]
	}
		       
	TrueVaR[irep] <- h[TT]*Hansen_qt(alpha, v, lambda)       
	TrueES[irep]  <- h[TT]*Hansen_ES(alpha, v, lambda)$ES 
 
	rets <- y[1:(TT-1)]
 	

	## QbSD

	results <- QbSD(rets, location="zero", CAViaR="SAV", alpha)                        
	SAVRobustVaRfmean[irep] <- results$VaR_oos
	SAVRobustESfmean[irep]  <- results$ES_oos    

	results <- QbSD(rets, location="zero", CAViaR="AS", alpha)        
	ASRobustVaRfmean[irep] <- results$VaR_oos
	ASRobustESfmean[irep]  <- results$ES_oos    


	## Taylor (2019) 
				
	results <- Taylor(rets, location="zero", CAViaR="SAV", ES="basic", alpha)              
	SAVTaylorVaRf.basic[irep] <- results$VaR_oos
	SAVTaylorESf.basic[irep]  <- results$ES_oos

	results <- Taylor(rets, location="zero", CAViaR="SAV", ES="AR", alpha)                                         
	SAVTaylorVaRf.AR[irep] <- results$VaR_oos
	SAVTaylorESf.AR[irep]  <- results$ES_oos

	results <- Taylor(rets, location="zero", CAViaR="AS", ES="basic", alpha)                            
	ASTaylorVaRf.basic[irep] <- results$VaR_oos
	ASTaylorESf.basic[irep]  <- results$ES_oos

	results <- Taylor(rets, location="zero", CAViaR="AS", ES="AR", alpha)                            
	ASTaylorVaRf.AR[irep] <- results$VaR_oos
	ASTaylorESf.AR[irep]  <- results$ES_oos
				    

	## Patton et al. (2019) 
				
	results <- Patton(rets, alpha)
	PattonVaRf[irep] <- results$VaR_oos
	PattonESf[irep]  <- results$ES_oos

				
	## GARCH family 
        
    results <- Hansen(rets, location="zero", alpha)  # GARCH + skew-t
	HansenVaRf[irep] <- results$VaR_oos
	HansenESf[irep]  <- results$ES_oos

	results <- GARCHnorm(rets, location="zero", alpha)  # GARCH + normal
	NormalVaRf[irep] <- results$VaR_oos
	NormalESf[irep]  <- results$ES_oos

	results <- GARCHt(rets, location="zero", alpha)     # GARCH + t
	StudentVaRf[irep] <- results$VaR_oos
	StudentESf[irep]  <- results$ES_oos

	results <- EGARCHnorm(rets, location="zero", alpha)  # EGARCH + normal
	EGARCHNormalVaRf[irep] <- results$VaR_oos
	EGARCHNormalESf[irep]  <- results$ES_oos

	results <- GJRnorm(rets, location="zero", alpha)     # GJR + normal
	GJRNormalVaRf[irep] <- results$VaR_oos
	GJRNormalESf[irep]  <- results$ES_oos

	results <- GJRt(rets, location="zero", alpha)        # GJR + t
	GJRStudentVaRf[irep] <- results$VaR_oos
	GJRStudentESf[irep]  <- results$ES_oos

	results <- GJRskewt(rets, location="zero", alpha)    # GJR + skew-t
	GJRskewtVaRf[irep] <- results$VaR_oos
	GJRskewtESf[irep]  <- results$ES_oos				
											
}	

			 
cat(c('==== DGP =========================================================', "\n"), file=filename, append=TRUE)
cat(c('Tsize =', TT, 'v = ', v, 'lambda = ', lambda, 'alpha = ', alpha, 'gammahat = ', gammahat, "\n"), file=filename, append=TRUE)
cat(c("\n"), file=filename, append=TRUE)
cat(c('VaR ************* MAE ************************************** VaR', "\n"), file=filename, append=TRUE)


cat(c('mean( abs(SAVRobustVaRfmean - TrueVaR))',   mean( abs(SAVRobustVaRfmean - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(ASRobustVaRfmean - TrueVaR))',    mean( abs(ASRobustVaRfmean - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(SAVTaylorVaRf.basic - TrueVaR))', mean( abs(SAVTaylorVaRf.basic - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(SAVTaylorVaRf.AR - TrueVaR))',    mean( abs(SAVTaylorVaRf.AR - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(ASTaylorVaRf.basic - TrueVaR))',  mean( abs(ASTaylorVaRf.basic - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(ASTaylorVaRf.AR - TrueVaR))',     mean( abs(ASTaylorVaRf.AR - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(PattonVaRf - TrueVaR))',          mean( abs(PattonVaRf - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(NormalVaRf - TrueVaR))',          mean( abs(NormalVaRf - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(StudentVaRf - TrueVaR))',         mean( abs(StudentVaRf - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(HansenVaRf - TrueVaR))',          mean( abs(HansenVaRf - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(GJRNormalVaRf - TrueVaR))',       mean( abs(GJRNormalVaRf - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(GJRStudentVaRf - TrueVaR))',      mean( abs(GJRStudentVaRf - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(GJRskewtVaRf - TrueVaR))',        mean( abs(GJRskewtVaRf - TrueVaR)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(EGARCHNormalVaRf - TrueVaR))',    mean( abs(EGARCHNormalVaRf - TrueVaR)),  "\n" ), file=filename, append=TRUE)

cat(c("\n"), file=filename, append=TRUE)
cat(c('VaR ************* RMSE ************************************** VaR', "\n"), file=filename, append=TRUE)

cat(c('sqrt(mean((SAVRobustVaRfmean - TrueVaR)^2))',   sqrt(mean((SAVRobustVaRfmean - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((ASRobustVaRfmean - TrueVaR)^2))',    sqrt(mean((ASRobustVaRfmean - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((SAVTaylorVaRf.basic - TrueVaR)^2))', sqrt(mean((SAVTaylorVaRf.basic - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((SAVTaylorVaRf.AR - TrueVaR)^2))',    sqrt(mean((SAVTaylorVaRf.AR - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((ASTaylorVaRf.basic - TrueVaR)^2))',  sqrt(mean((ASTaylorVaRf.basic - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((ASTaylorVaRf.AR - TrueVaR)^2))',     sqrt(mean((ASTaylorVaRf.AR - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((PattonVaRf - TrueVaR)^2))',          sqrt(mean((PattonVaRf - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((NormalVaRf - TrueVaR)^2))',          sqrt(mean((NormalVaRf - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((StudentVaRf - TrueVaR)^2))',         sqrt(mean((StudentVaRf - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((HansenVaRf - TrueVaR)^2))',          sqrt(mean((HansenVaRf - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((GJRNormalVaRf - TrueVaR)^2))',       sqrt(mean((GJRNormalVaRf - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((GJRStudentVaRf - TrueVaR)^2))',      sqrt(mean((GJRStudentVaRf - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((GJRskewtVaRf - TrueVaR)^2))',        sqrt(mean((GJRskewtVaRf - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((EGARCHNormalVaRf - TrueVaR)^2))',    sqrt(mean((EGARCHNormalVaRf - TrueVaR)^2)),  "\n" ), file=filename, append=TRUE)

cat(c("\n"), file=filename, append=TRUE)
cat(c("\n"), file=filename, append=TRUE)
cat(c('ES ************ MAE *************************************** ES', "\n"), file=filename, append=TRUE)

cat(c('mean( abs(SAVRobustESfmean - TrueES))',   mean( abs(SAVRobustESfmean - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(ASRobustESfmean - TrueES))',    mean( abs(ASRobustESfmean - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(SAVTaylorESf.basic - TrueES))', mean( abs(SAVTaylorESf.basic - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(SAVTaylorESf.AR - TrueES))',    mean( abs(SAVTaylorESf.AR - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(ASTaylorESf.basic - TrueES))',  mean( abs(ASTaylorESf.basic - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(ASTaylorESf.AR - TrueES))',     mean( abs(ASTaylorESf.AR - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(PattonESf - TrueES))',          mean( abs(PattonESf - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(NormalESf - TrueES))',          mean( abs(NormalESf - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(StudentESf - TrueES))',         mean( abs(StudentESf - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(HansenESf - TrueES))',          mean( abs(HansenESf - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(GJRNormalESf - TrueES))',       mean( abs(GJRNormalESf - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(GJRStudentESf - TrueES))',      mean( abs(GJRStudentESf - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(GJRskewtESf - TrueES))',        mean( abs(GJRskewtESf - TrueES)),  "\n" ), file=filename, append=TRUE)
cat(c('mean( abs(EGARCHNormalESf - TrueES))',    mean( abs(EGARCHNormalESf - TrueES)),  "\n" ), file=filename, append=TRUE)

cat(c("\n"), file=filename, append=TRUE)
cat(c('ES ************ RMSE *************************************** ES', "\n"), file=filename, append=TRUE)

cat(c('sqrt(mean((SAVRobustESfmean - TrueES)^2))',   sqrt(mean((SAVRobustESfmean - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((ASRobustESfmean - TrueES)^2))',    sqrt(mean((ASRobustESfmean - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((SAVTaylorESf.basic - TrueES)^2))', sqrt(mean((SAVTaylorESf.basic - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((SAVTaylorESf.AR - TrueES)^2))',    sqrt(mean((SAVTaylorESf.AR - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((ASTaylorESf.basic - TrueES)^2))',  sqrt(mean((ASTaylorESf.basic - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((ASTaylorESf.AR - TrueES)^2))',     sqrt(mean((ASTaylorESf.AR - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((PattonESf - TrueES)^2))',          sqrt(mean((PattonESf - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((NormalESf - TrueES)^2))',          sqrt(mean((NormalESf - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((StudentESf - TrueES)^2))',         sqrt(mean((StudentESf - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((HansenESf - TrueES)^2))',          sqrt(mean((HansenESf - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((GJRNormalESf - TrueES)^2))',       sqrt(mean((GJRNormalESf - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((GJRStudentESf - TrueES)^2))',      sqrt(mean((GJRStudentESf - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((GJRskewtESf - TrueES)^2))',        sqrt(mean((GJRskewtESf - TrueES)^2)),  "\n" ), file=filename, append=TRUE)
cat(c('sqrt(mean((EGARCHNormalESf - TrueES)^2))',    sqrt(mean((EGARCHNormalESf - TrueES)^2)),  "\n" ), file=filename, append=TRUE)

cat(c("\n"), file=filename, append=TRUE)
cat(c("\n"), file=filename, append=TRUE)
cat(c("\n"), file=filename, append=TRUE)		



