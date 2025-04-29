# QAR1: Estimates a Quantile Autoregression model of order 1 at the median (tau = 0.5)
QAR1 <- function(rets){	
	Tsize <- length(rets) 
	temp <- rq(rets[2:Tsize]~rets[1:(Tsize-1)], tau=0.5)
	muhat <- as.numeric(temp$coef[1])
	phihat <- as.numeric(temp$coef[2]) 
	medianhat <- matrix(0, Tsize, 1)
	medianhat[2:Tsize] <- muhat + phihat*rets[1:(Tsize-1)]	
	medianf <- muhat + phihat*rets[Tsize]		
	return(list(medianhat=medianhat, medianf=medianf))		
}



#--------------------------------------------------------------------------------------------------
# Robust: Quantile-based Scale Dynamics (QbSD) model for VaR and ES forecasting
#
# Args:
#   rets     : A numeric vector of returns
#   location : Location specification ("zero" for constant, "QAR" for QAR(1) median)
#   CAViaR   : Scale dynamics specification ("SAV" or "AS")
#   alpha    : VaR/ES confidence level (e.g., 0.01, 0.025, 0.05)
#
# Returns:
#   A list of in-sample and out-of-sample (OOS) estimates for VaR and ES using the QbSD approach
#--------------------------------------------------------------------------------------------------

QbSD <- function(rets, location, CAViaR, alpha){

	#-----------------------------------------------------------------------------------------------
	# Scale(): Estimate the conditional scale using restricted quantile regressions
	#-----------------------------------------------------------------------------------------------  
	Scale <- function(centeredrets, p){

		probs <- c(p, 1-p) 
		Q0 <- quantile(centeredrets, probs)

		# Ensure scale is defined by checking that lower quantile is negative and upper is positive    			
		while ((Q0[1] >= 0) | (Q0[2] <= 0)){				
			p <- p - 0.05
			probs <- c(p, 1-p) 
			Q0 <- quantile(centeredrets, probs)				
		}
												
		if (CAViaR == "SAV"){
								
			# Grid search over initial values				
			beta0 <- seq(0.1, 0.9, length.out=3)
			gamma0 <- seq(0.1, 0.9, length.out=3)
								
			results <- matrix(0, length(beta0)*length(gamma0), 5)
				
			kk <- 1
				
			for (i in beta0){
				for (j in gamma0){ 
										
					pars0 <- c(Q0, log(i), log(j))			
						
					temp <- optim(par=pars0, fn=SAVscale_Cpp, centeredrets=centeredrets, tau=p, q0=Q0)
					pars1 <- temp$par			
					pars1[3:4] <- exp(pars1[3:4])
						
					results[kk,1:4] <- pars1
					results[kk, 5] <- temp$value
						
					kk <- kk + 1						
				}
			}		
				
   		# Select best result and compute fitted/predicted quantiles
			idx <- which.min(results[,5])				
			pars1 <- results[idx, 1:4]
																		
			SAVtemp <- SAVscalehatf_Cpp(centeredrets=centeredrets, q0=Q0, pars=pars1)

			Q1hat   <- SAVtemp[1:Tsize,1]				
			Q1f     <- SAVtemp[Tsize+1,1]
																
			Q2hat   <- SAVtemp[1:Tsize,2]				
			Q2f     <- SAVtemp[Tsize+1,2]				
		}	

		if (CAViaR == "AS"){	
								
			beta0 <- seq(0.1, 0.9, length.out=3)
			gamma10 <- seq(0.1, 0.9, length.out=3)
			gamma20 <- seq(0.1, 0.9, length.out=3)
								
			results <- matrix(0, length(beta0)*length(gamma10)*length(gamma20), 6)
				
			kk <- 1
				
			for (i in beta0){
				for (j in gamma10){ 
					for (k in gamma20 ){
								
						pars0 <- c(Q0, log(i), log(j), log(k))		
																
						temp <- optim(par=pars0, fn=ASscale_Cpp, centeredrets=centeredrets, tau=p, q0=Q0)
														
						pars1 <- temp$par			
						pars1[3:5] <- exp(pars1[3:5])
													
						results[kk,1:5] <- pars1
						results[kk, 6] <- temp$value
						
						kk <- kk + 1																					
					}							
				}					
			}

			idx <- which.min(results[,6])				
			pars1 <- results[idx, 1:5]

			AStemp <- ASscalehatf_Cpp(centeredrets=centeredrets, q0=Q0, pars=pars1)

			Q1hat   <- AStemp[1:Tsize,1]				
			Q1f     <- AStemp[Tsize+1,1]
																
			Q2hat   <- AStemp[1:Tsize,2]				
			Q2f     <- AStemp[Tsize+1,2]
		}	
			
		scalehat <- Q2hat - Q1hat					
		scalef <- Q2f - Q1f
		return(list(scalehat=scalehat, scalef=scalef))	
	}

	#-----------------------------------------------------------------------------------------------
  	# VaR(): Compute conditional quantile given standardized innovation
	#----------------------------------------------------------------------------------------------- 
	VaR <- function(){			
		QQ0 <- quantile(epsilon, alpha)		
		varhat <- muhat + scalehat* QQ0
		varf <- muf + scalef* QQ0
		return(list(varhat=varhat, varf=varf))
	}

	#-----------------------------------------------------------------------------------------------
	# ES(): Approximate expected shortfall as average of conditional quantiles below VaR
	#-----------------------------------------------------------------------------------------------  
	ES <- function(){										
		TargetTaus <- seq(from=alpha/N, to=alpha, by=alpha/N)							
		Nsize <- length(TargetTaus)
		QQ0 <- matrix(quantile(epsilon, TargetTaus), Nsize, 1)		
		QQhat <- matrix(0, Tsize, Nsize)		
		QQf   <- matrix(0, 1, Nsize)		

		for (i in 1:Nsize){
			QQhat[,i] <- muhat + scalehat* QQ0[i]
			QQf[,i]   <- muf + scalef* QQ0[i]				
		}	

		eshat <- apply(QQhat, 1, mean)	
		esf   <- mean(QQf)
		return(list(eshat=eshat, esf=esf))
	}
		
	#-----------------------------------------------------------------------------------------------
	# Main computation
	#-----------------------------------------------------------------------------------------------
  
	Tsize <- length(rets)
			
	P <- seq(from=0.05, to=0.25, by=0.05)
		
	# Compute conditional location
	if (location == "zero"){ 
		muhat <- matrix(0, Tsize, 1)		
		muf   <- 0
	}	

	if (location == "QAR"){
		temp <- QAR1(rets)
		muhat <- temp$medianhat
		muf   <- temp$medianf		
	}	
	
	centeredrets <- rets - muhat
		
	# Containers for multiple quantile levels
	ScaleHat <- matrix(0, Tsize, length(P))
	ScaleF   <- numeric(length(P))

	VaRp  <- matrix(0, Tsize, length(P))	
	VaRfp  <- matrix(0, 1, length(P))	
		
	for (j in 1:length(P)){		

		p <- P[j]

		temp <- Scale(centeredrets, p)
		ScaleHat[,j] <- scalehat <- temp$scalehat
		ScaleF[j]    <- scalef   <- temp$scalef
							
		epsilon <- (rets-muhat)/scalehat

		temp <- VaR()		
		VaRp[,j] <- temp$varhat		
		VaRfp[j] <- temp$varf					
	}

	VaR_in    <- apply(VaRp[,1:length(P)], 1, mean)
	VaR_oos   <- mean(VaRfp)

	#-----------------------------------------------------------------------------------------------
	# ES forecast loop: increase number of quantile grid points until convergence
	#-----------------------------------------------------------------------------------------------
  
	ESp   <- matrix(0, Tsize, length(P))	
	ESfp   <- matrix(0, 1, length(P))	

	tol <- 0.0001
				
	N_values <- 4:100
	
	results <- matrix(0, length(N_values), 2)
	
	for (idx in 1:length(N_values)){
	
		N <- N_values[idx]
	
		for (j in 1:length(P)){		
			p <- P[j]

			scalehat <- ScaleHat[,j]
			scalef <- ScaleF[j]
			
			epsilon <- (rets-muhat)/scalehat

			temp <- ES()					
			ESp[,j] <- temp$eshat		
			ESfp[j] <- temp$esf		
		}
	
		ES_in  <- apply(ESp[,1:length(P)], 1, mean)
		ES_oos  <- mean(ESfp)
	
		results[idx, 1] <- N
		results[idx, 2] <- ES_oos
				
		if ((idx >=2) && (abs(results[idx,2] - results[idx-1,2])< tol)) break	
	}
		
	#-----------------------------------------------------------------------------------------------
 	# Return VaR and ES in-sample estimates and OOS forecasts
	#-----------------------------------------------------------------------------------------------  													
	return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos))	
}	









