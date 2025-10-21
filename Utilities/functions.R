QAR1 <- function(returns){	
	Tsize <- length(returns) 
	temp <- rq(returns[2:Tsize]~returns[1:(Tsize-1)], tau=0.5)
	muhat <- as.numeric(temp$coef[1])
	phihat <- as.numeric(temp$coef[2]) 
	medianhat <- numeric(Tsize)
	medianhat[1] <- if (abs(phihat) < 0.999) muhat/(1 - phihat) else median(returns)
	medianhat[2:Tsize] <- muhat + phihat*returns[1:(Tsize-1)]	
	medianf <- muhat + phihat*returns[Tsize]		
	return(list(medianhat=medianhat, medianf=medianf))		
}



AR1 <- function(returns){	
  Tsize <- length(returns) 
  temp <- lm(returns[2:Tsize]~returns[1:(Tsize-1)])
  muhat <- as.numeric(temp$coef[1])
  phihat <- as.numeric(temp$coef[2]) 
  meanhat <- numeric(Tsize)  
  meanhat[1] <- if (abs(phihat) < 0.999) muhat/(1 - phihat) else mean(returns)  
  meanhat[2:Tsize] <- muhat + phihat*returns[1:(Tsize-1)]	
  meanf <- muhat + phihat*returns[Tsize]	
  return(list(meanhat=meanhat, meanf=meanf))		
}



QbSD <- function(returns, location, CAViaR, alpha){

	Scale <- function(centeredreturns, p){
		
		Tsize <- length(centeredreturns)
		
		probs <- c(p, 1-p) 
		Q0 <- quantile(centeredreturns, probs)

		p_min <- 0.010001

		while ((Q0[1] >= 0 || Q0[2] <= 0) && p > p_min){
		  p  <- max(p_min, p - 0.01)
		  Q0 <- quantile(centeredreturns, c(p, 1 - p), na.rm = TRUE, names = FALSE)
		}

		if (Q0[1] >= 0 || Q0[2] <= 0){
		  stop(sprintf("Endpoints did not straddle zero at p = %.4f", p))
		}

												
		if (CAViaR == "gSAV"){
								
			beta0   <- seq(0.85, 0.95, length.out=3)
			gamma0  <- seq(0.01, 0.1, length.out=3)								
			results <- matrix(0, length(beta0)*length(gamma0), 5)
				
			kk <- 1				
			for (i in beta0){
				for (j in gamma0){ 										
					pars0 <- c(Q0, log(i), log(j))									
					temp <- optim(par=pars0, fn=gSAV_Cpp, centeredreturns=centeredreturns, tau=p, q0=Q0)										
					pars1 <- temp$par									
					results[kk,1:4] <- pars1
					results[kk, 5] <- temp$value						
					kk <- kk + 1						
				}
			}		
				
			idx <- which.min(results[,5])				
			pars1 <- results[idx, 1:4]
																		
			SAVtemp <- gSAV_forecast_Cpp(centeredreturns=centeredreturns, q0=Q0, pars=pars1)

			Q1hat   <- SAVtemp[1:Tsize, 1]				
			Q1f     <- SAVtemp[Tsize+1, 1]
																
			Q2hat   <- SAVtemp[1:Tsize, 2]				
			Q2f     <- SAVtemp[Tsize+1, 2]				
		}	

		if (CAViaR == "gAS"){	
								
			beta0   <- seq(0.85, 0.95, length.out=3)
			gamma10 <- seq(0.01, 0.1, length.out=3)
			gamma20 <- seq(0.01, 0.1, length.out=3)								
			results <- matrix(0, length(beta0)*length(gamma10)*length(gamma20), 6)
				
			kk <- 1				
			for (i in beta0){
				for (j in gamma10){ 
					for (k in gamma20 ){								
						pars0 <- c(Q0, log(i), log(j), log(k))																		
						temp <- optim(par=pars0, fn=gAS_Cpp, centeredreturns=centeredreturns, tau=p, q0=Q0)														
						pars1 <- temp$par																
						results[kk,1:5] <- pars1
						results[kk, 6] <- temp$value						
						kk <- kk + 1																					
					}							
				}					
			}

			idx <- which.min(results[,6])				
			pars1 <- results[idx, 1:5]

			AStemp <- gAS_forecast_Cpp(centeredreturns=centeredreturns, q0=Q0, pars=pars1)

			Q1hat   <- AStemp[1:Tsize,1]				
			Q1f     <- AStemp[Tsize+1,1]
																
			Q2hat   <- AStemp[1:Tsize,2]				
			Q2f     <- AStemp[Tsize+1,2]
		}	
			
		scalehat <- pmax(Q2hat - Q1hat, 1e-8)
		scalef   <- max(Q2f - Q1f, 1e-8)

		return(list(scalehat=scalehat, scalef=scalef))	
	}

	VaR <- function(){			
		QQ0 <- quantile(epsilon, alpha)		
		varhat <- muhat + scalehat* QQ0
		varf <- mu_oos + scalef* QQ0
		return(list(varhat=varhat, varf=varf))
	}

	ES <- function(){										
		TargetTaus <- seq(from=alpha/N, to=alpha, by=alpha/N)							
		Nsize <- length(TargetTaus)
		QQ0 <- as.numeric(quantile(epsilon, TargetTaus, na.rm = TRUE, names = FALSE))		
		QQhat <- matrix(0, Tsize, Nsize)		
		QQf <- numeric(Nsize)

		for (i in 1:Nsize){
			QQhat[,i] <- muhat  + scalehat* QQ0[i]
			QQf[i]    <- mu_oos + scalef* QQ0[i]				
		}	

		eshat <- apply(QQhat, 1, mean)	
		esf   <- mean(QQf)
		return(list(eshat=eshat, esf=esf))
	}
		
  
	Tsize <- length(returns)
			
	P <- seq(from=0.05, to=0.25, by=0.05)
			
	if (location == "zero"){
	  muhat  <- rep(0, Tsize)
	  mu_oos <- 0
	} else if (location == "QAR"){
	  tmp    <- QAR1(returns)
	  muhat  <- as.numeric(tmp$medianhat)
	  mu_oos <- as.numeric(tmp$medianf)
	} else {
	  stop("location must be 'zero' or 'QAR'")
	}	
	
	centeredreturns <- returns - muhat
		
	ScaleHat <- matrix(0, Tsize, length(P))
	ScaleF   <- numeric(length(P))

	VaRp  <- matrix(0, Tsize, length(P))	
	VaRfp <- numeric(length(P))
		
	for (j in 1:length(P)){		

		p <- P[j]

		temp <- Scale(centeredreturns, p)
		ScaleHat[,j] <- scalehat <- temp$scalehat
		ScaleF[j]    <- scalef   <- temp$scalef
							
		epsilon <- (returns-muhat)/scalehat

		temp <- VaR()		
		VaRp[,j] <- temp$varhat		
		VaRfp[j] <- temp$varf					
	}

	VaR_in   <- apply(VaRp[,1:length(P)], 1, mean)
	VaR_oos  <- mean(VaRfp)
  
	ESp    <- matrix(0, Tsize, length(P))	
	ESfp  <- numeric(length(P))	

	tol <- 0.0001
				
	N_values <- 4:100
	
	results <- matrix(0, length(N_values), 2)
	
	for (idx in 1:length(N_values)){
	
		N <- N_values[idx]
	
		for (j in 1:length(P)){		
			p <- P[j]

			scalehat <- ScaleHat[,j]
			scalef <- ScaleF[j]
			
			epsilon <- (returns-muhat)/scalehat

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
		
	return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos, mu_oos=mu_oos))	
}	



AL <- function(returns, location, CAViaR, ES, alpha){	
  	Tsize <- length(returns)

	if (location == "zero"){
    	muhat  <- rep(0, Tsize)
    	mu_oos <- 0
	} else if (location == "AR"){
    	tmp    <- AR1(returns)
    	muhat  <- tmp$meanhat
    	mu_oos <- tmp$meanf
  	} else {
    	stop("location must be 'zero' or 'AR'")
  	}

  	if (!(alpha > 0 && alpha < 1)){
    	stop("alpha must be in (0,1). For 99% VaR use 0.01; for 95% use 0.05, etc.")
  	}
  	
	centeredreturns <- returns - muhat

  	prob <- alpha
  	Q0   <- as.numeric(quantile(centeredreturns, prob))
  	
  	
	if (identical(CAViaR, "SAV") && identical(ES, "mult")){

    	beta10 <- seq(0.85, 0.95, length.out = 3)
    	beta20 <- seq(-0.01, 0.01, length.out = 3)
    	grid   <- matrix(0, length(beta10) * length(beta20), 4)

    	kk <- 1
    	for (b1 in beta10) for (b2 in beta20){
      		pars0 <- c(Q0, log(b1), b2)
      		opt   <- optim(par = pars0, fn = SAV_Cpp, centeredreturns = centeredreturns, tau = prob, q0 = Q0)
      		grid[kk, 1:3] <- opt$par
      		grid[kk, 4]   <- opt$value
      		kk <- kk + 1
    	}
    	betahat0 <- grid[which.min(grid[, 4]), 1:3]

    	VaR_path_centered <- SAV_forecast_Cpp(centeredreturns = centeredreturns, tau=prob, q0=Q0, pars=betahat0)
    	VaRhat_centered   <- as.numeric(VaR_path_centered[1:Tsize])

    	opt1d     <- optimize(f = mult_Cpp, interval = c(-10, 10), centeredreturns = centeredreturns, alpha = prob, VaR = VaRhat_centered, maximum = TRUE)
    	gammahat0 <- opt1d$maximum

    	pars0 <- c(betahat0, gammahat0)
    	mle   <- optim(par = pars0, fn = AL_mult_SAV_Cpp, centeredreturns = centeredreturns, alpha = prob, q0 = Q0, control = list(fnscale = -1))
    	pars1 <- mle$par

    	fcst     <- AL_mult_SAV_forecast_Cpp(centeredreturns = centeredreturns, alpha=prob, q0=Q0, pars=pars1)
    	VaR_in   <- muhat + as.numeric(fcst[1:Tsize, 1])
    	VaR_oos  <- mu_oos + as.numeric(fcst[Tsize + 1, 1])
    	ES_in    <- muhat + as.numeric(fcst[1:Tsize, 2])
    	ES_oos   <- mu_oos + as.numeric(fcst[Tsize + 1, 2])

  	} else if (identical(CAViaR, "AS") && identical(ES, "mult")){

    	beta10 <- seq(0.85, 0.95, length.out = 3)
    	beta20 <- seq(-0.01, 0.01, length.out = 3)
    	beta30 <- seq(-0.01, 0.01, length.out = 3)
    	grid   <- matrix(0, length(beta10) * length(beta20) * length(beta30), 5)

    	kk <- 1
    	for (b1 in beta10) for (b2 in beta20) for (b3 in beta30){
      		pars0 <- c(Q0, log(b1), b2, b3)
      		opt   <- optim(par = pars0, fn = AS_Cpp, centeredreturns = centeredreturns, tau=prob, q0=Q0)
	      	grid[kk, 1:4] <- opt$par
      		grid[kk, 5]   <- opt$value
      		kk <- kk + 1
    	}
    	betahat0 <- grid[which.min(grid[, 5]), 1:4]

    	VaR_path_centered <- AS_forecast_Cpp(centeredreturns = centeredreturns, tau=prob, q0=Q0, pars=betahat0)
    	VaRhat_centered   <- as.numeric(VaR_path_centered[1:Tsize])

    	opt1d     <- optimize(f = mult_Cpp, interval = c(-10, 10), centeredreturns = centeredreturns, alpha = prob, VaR = VaRhat_centered, maximum = TRUE)
    	gammahat0 <- opt1d$maximum

    	pars0 <- c(betahat0, gammahat0)
    	mle   <- optim(par = pars0, fn = AL_mult_AS_Cpp, centeredreturns = centeredreturns, alpha = prob, q0 = Q0, control = list(fnscale = -1))
    	pars1 <- mle$par

    	fcst     <- AL_mult_AS_forecast_Cpp(centeredreturns = centeredreturns, alpha=prob, q0=Q0, pars=pars1)
    	VaR_in   <- muhat + as.numeric(fcst[1:Tsize, 1])
    	VaR_oos  <- mu_oos + as.numeric(fcst[Tsize + 1, 1])
    	ES_in    <- muhat + as.numeric(fcst[1:Tsize, 2])
    	ES_oos   <- mu_oos + as.numeric(fcst[Tsize + 1, 2])

  	} else if (identical(CAViaR, "SAV") && identical(ES, "ar")){
   
    	beta10 <- seq(0.85, 0.95, length.out = 3)
    	beta20 <- seq(-0.01, 0.01, length.out = 3)
    	grid   <- matrix(0, length(beta10) * length(beta20), 4)

    	kk <- 1
    	for (b1 in beta10) for (b2 in beta20){
      		pars0 <- c(Q0, log(b1), b2)
      		opt   <- optim(par = pars0, fn = SAV_Cpp, centeredreturns = centeredreturns, tau = prob, q0 = Q0)
      		grid[kk, 1:3] <- opt$par
      		grid[kk, 4]   <- opt$value
      		kk <- kk + 1
    	}
    	betahat0 <- grid[which.min(grid[, 4]), 1:3]
    	
    	VaR_path_centered <- SAV_forecast_Cpp(centeredreturns = centeredreturns, tau=prob, q0=Q0, pars=betahat0)
    	VaRhat_centered   <- as.numeric(VaR_path_centered[1:Tsize])

	    gamma00 <- seq(0.0001, 0.01, length.out = 3)
	    gamma10 <- seq(0.10,   0.20, length.out = 3)
	    gamma20 <- seq(0.60,   0.80, length.out = 3)

	    grid <- matrix(0, length(gamma00) * length(gamma10) * length(gamma20), 4)
	    kk <- 1
	    for (g0 in gamma00) for (g1 in gamma10) for (g2 in gamma20){
	    	pars0 <- c(log(g0), log(g1), log(g2))
      		opt   <- optim(par = pars0, fn = AR_Cpp, centeredreturns = centeredreturns, alpha = prob, VaR = VaRhat_centered, control = list(fnscale = -1)) 
      		grid[kk, 1:3] <- opt$par
      		grid[kk, 4]   <- opt$value
      		kk <- kk + 1
    	}
    	uhat0 <- (grid[which.max(grid[, 4]), 1:3])

    	pars0 <- c(betahat0, uhat0)
    	mle   <- optim(par = pars0, fn = AL_ar_SAV_Cpp, centeredreturns = centeredreturns, alpha = prob, q0 = Q0, control = list(fnscale = -1))
    	pars1 <- mle$par

    	fcst     <- AL_ar_SAV_forecast_Cpp(centeredreturns = centeredreturns, alpha=prob, q0=Q0, pars=pars1)
    	VaR_in   <- muhat + as.numeric(fcst[1:Tsize, 1])
    	VaR_oos  <- mu_oos + as.numeric(fcst[Tsize + 1, 1])
    	ES_in    <- muhat + as.numeric(fcst[1:Tsize, 2])
    	ES_oos   <- mu_oos + as.numeric(fcst[Tsize + 1, 2])

	  } else if (identical(CAViaR, "AS") && identical(ES, "ar")){

    	beta10 <- seq(0.85, 0.95, length.out = 3)
    	beta20 <- seq(-0.01, 0.01, length.out = 3)
    	beta30 <- seq(-0.01, 0.01, length.out = 3)	    	    
	    grid   <- matrix(0, length(beta10) * length(beta20) * length(beta30), 5)
    	
    	kk <- 1
    	for (b1 in beta10) for (b2 in beta20) for (b3 in beta30){
      		pars0 <- c(Q0, log(b1), b2, b3)
      		opt   <- optim(par = pars0, fn = AS_Cpp, centeredreturns = centeredreturns, tau=prob, q0=Q0)
	      	grid[kk, 1:4] <- opt$par
      		grid[kk, 5]   <- opt$value
      		kk <- kk + 1
    	}
    	betahat0 <- grid[which.min(grid[, 5]), 1:4]

	    VaR_path_centered <- AS_forecast_Cpp(centeredreturns = centeredreturns, tau=prob, q0=Q0, pars=betahat0)
	    VaRhat_centered   <- as.numeric(VaR_path_centered[1:Tsize])
	    
	    gamma00 <- seq(0.0001, 0.01, length.out = 3)
	    gamma10 <- seq(0.10,   0.20, length.out = 3)
	    gamma20 <- seq(0.60,   0.80, length.out = 3)
	    
		grid <- matrix(0, length(gamma00)*length(gamma10)*length(gamma20), 4)
		kk <- 1
		for (g0 in gamma00) for (g1 in gamma10) for (g2 in gamma20) {
	    	pars0 <- c(log(g0), log(g1), log(g2))
      		opt   <- optim(par = pars0, fn = AR_Cpp, centeredreturns = centeredreturns, alpha = prob, VaR = VaRhat_centered, control = list(fnscale = -1)) 
  			grid[kk, 1:3] <- opt$par
  			grid[kk, 4]   <- opt$value
  			kk <- kk + 1
		}
		uhat0 <- grid[which.max(grid[, 4]), 1:3]

    	pars0 <- c(betahat0, uhat0)
		
    	mle   <- optim(par = pars0, fn = AL_ar_AS_Cpp, centeredreturns = centeredreturns, alpha = prob, q0 = Q0, control = list(fnscale = -1))
    	pars1 <- mle$par			
    	    	
    	fcst     <- AL_ar_AS_forecast_Cpp(centeredreturns = centeredreturns, alpha=prob, q0=Q0, pars=pars1)
    	VaR_in   <- muhat + as.numeric(fcst[1:Tsize, 1])
    	VaR_oos  <- mu_oos + as.numeric(fcst[Tsize + 1, 1])
    	ES_in    <- muhat + as.numeric(fcst[1:Tsize, 2])
    	ES_oos   <- mu_oos + as.numeric(fcst[Tsize + 1, 2])

  	} else {
    	stop("Unsupported combination: CAViaR must be 'SAV' or 'AS' and ES must be 'mult' or 'ar'.")
  	}

  	list(VaR_in = VaR_in, ES_in = ES_in, VaR_oos = VaR_oos, ES_oos = ES_oos, mu_oos = mu_oos)
}



GAS <- function(returns, location, alpha){

  	Tsize <- length(returns)

	if (location == "zero"){
    	muhat  <- rep(0, Tsize)
    	mu_oos <- 0
	} else if (location == "AR"){
    	tmp    <- AR1(returns)
    	muhat  <- tmp$meanhat
    	mu_oos <- tmp$meanf
  	} else {
    	stop("location must be 'zero' or 'AR'")
  	}

  	if (!(alpha > 0 && alpha < 1)){
    	stop("alpha must be in (0,1). For 99% VaR use 0.01; for 95% use 0.05, etc.")
  	}
  	
	centeredreturns <- returns - muhat
  
  	a_grid_abs   <- c(0.5, 1.0, 2.0)        
	p1_grid      <- log(a_grid_abs)

  	delta_grid   <- c(0.10, 0.20, 0.50)     
  	p2_grid      <- log(delta_grid)

  	beta_grid    <- c(0.90, 0.95, 0.99)     
  	p3_grid      <- log(beta_grid)

  	gamma_grid   <- c(0.001, 0.01, 0.05)    
  	p4_grid      <- log(gamma_grid)

  	results <- matrix(NA, nrow = length(p1_grid)*length(p2_grid)*length(p3_grid)*length(p4_grid), ncol = 5)

  	kk <- 1
  	for (p1 in p1_grid) for (p2 in p2_grid) for (p3 in p3_grid) for (p4 in p4_grid) {
    	pars0 <- c(p1, p2, p3, p4)
    	opt   <- optim(par = pars0, fn = GAS_Cpp, centeredreturns = centeredreturns, alpha = alpha)
    	results[kk, 1:4] <- opt$par
    	results[kk, 5]   <- opt$value
    	kk <- kk + 1
  	}

  	best   <- which.min(results[, 5])
  	pars1  <- results[best, 1:4]

  	res <- GAS_forecast_Cpp(centeredreturns = centeredreturns, pars = pars1, alpha = alpha)

  	VaR_in  <- muhat + res[1:Tsize, 1]
  	VaR_oos <- mu_oos + res[Tsize + 1, 1]
  	ES_in   <- muhat + res[1:Tsize, 2]
  	ES_oos  <- mu_oos + res[Tsize + 1, 2]
  	
  	list(VaR_in = VaR_in, ES_in = ES_in, VaR_oos = VaR_oos, ES_oos = ES_oos, mu_oos = mu_oos)
}



GARCH_norm <- function(returns, location = "zero", alpha){
  	Tsize <- length(returns)

  	if (location == "zero"){
    	muhat  <- rep(0, Tsize)
    	mu_oos <- 0
  	} else if (location == "AR"){
    	temp   <- AR1(returns)
    	muhat  <- temp$meanhat
    	mu_oos <- temp$meanf
  	} else {
    	stop("location must be 'zero' or 'AR'")
  	}

  	if (!(alpha > 0 && alpha < 1)){
    	stop("alpha must be in (0, 1). For 99% VaR use alpha = 0.01, for 95% use 0.05, etc.")
  	}

  	centeredreturns <- returns - muhat

  	beta_to_theta <- function(beta){
    	beta0 <- beta[1]; b1 <- beta[2]; b2 <- beta[3]
    	c <- log(beta0)
    	den <- 1 - b1 - b2          
    	s1 <- b1 / den              
    	s2 <- b2 / den
    	a <- log(s1)                
    	b <- log(s2)
    	c(c = c, a = a, b = b)
  	}

  	pars0 <- beta_to_theta(c(0.01, 0.8, 0.1))
  	var0  <- var(centeredreturns)

  	temp <- optim(par = pars0, fn = GARCH_norm_Cpp,	centeredreturns = centeredreturns, var0 = var0,	control = list(fnscale = -1))
  	pars1 <- temp$par

  	sigma2_path     <- GARCH_norm_forecast_Cpp(centeredreturns= centeredreturns, pars=pars1, var0=var0)
  	sigma2_in       <- sigma2_path[1:Tsize] 
  	sigma2_forecast <- sigma2_path[Tsize + 1]  
  	sigma_in        <- sqrt(sigma2_in)       
  	sigma_forecast  <- sqrt(sigma2_forecast)  

  	z        <- qnorm(alpha)
  	VaR_in   <- as.numeric(muhat + sigma_in * z)
  	ES_in    <- as.numeric(muhat - sigma_in * dnorm(z) / alpha)
  	VaR_oos  <- as.numeric(mu_oos + sigma_forecast * z)
  	ES_oos   <- as.numeric(mu_oos - sigma_forecast * dnorm(z) / alpha)

  	list(VaR_in = VaR_in, ES_in = ES_in, VaR_oos = VaR_oos, ES_oos = ES_oos, mu_oos = mu_oos)
}



GARCH_t <- function(returns, location = "zero", alpha){		
  	Tsize <- length(returns)

  	if (location == "zero"){
    	muhat  <- rep(0, Tsize)
    	mu_oos <- 0
  	} else if (location == "AR"){
    	temp   <- AR1(returns)
    	muhat  <- temp$meanhat
    	mu_oos <- temp$meanf
  	} else {
    	stop("location must be 'zero' or 'AR'")
  	}

  	if (!(alpha > 0 && alpha < 1)){
    	stop("alpha must be in (0, 1). For 99% VaR use alpha = 0.01, for 95% use 0.05, etc.")
  	}

  	centeredreturns <- returns - muhat
  	var0 <- var(centeredreturns)
  
  	beta_to_theta <- function(beta){
    	beta0 <- beta[1]; b1 <- beta[2]; b2 <- beta[3]; v <- beta[4]
    	c  <- log(beta0)
	    den <- max(1e-12, 1 - b1 - b2)
	    s1 <- b1 / den
	    s2 <- b2 / den
	    a <- log(s1)
	    b <- log(s2)
	    g <- log(max(1e-12, v - 2))
	    c(c = c, a = a, b = b, g = g)
	  }

  	  beta_start <- c(beta0 = 0.01, beta1 = 0.8, beta2 = 0.1, v = 30)	       
	  pars0  <- beta_to_theta(beta_start)    

	  fit <- optim(par=pars0, fn=GARCH_t_Cpp, centeredreturns=centeredreturns, var0=var0, control=list(fnscale = -1))
	  pars1 <- fit$par
	  v <- 2 + exp(pars1[4])

	  sigma2_path     <- GARCH_t_forecast_Cpp(centeredreturns= centeredreturns, pars=pars1, var0=var0)
	  sigma2_in       <- sigma2_path[1:Tsize]
	  sigma2_forecast <- sigma2_path[Tsize + 1]
	  sigma_in        <- sqrt(sigma2_in)
	  sigma_forecast  <- sqrt(sigma2_forecast)

	  # VaR / ES under standardized t innovations with Var=1 scaling √((v-2)/v)
	  zt <- qt(alpha, df = v)
	  scale_t <- sqrt((v - 2) / v)

	  VaR_in  <- as.numeric(muhat  + sigma_in       * zt * scale_t)
	  ES_in   <- as.numeric(muhat  - sigma_in       * scale_t * ((v + zt^2) / (v - 1)) * dt(zt, df = v) / alpha)
	  VaR_oos <- as.numeric(mu_oos + sigma_forecast * zt * scale_t)
	  ES_oos  <- as.numeric(mu_oos - sigma_forecast * scale_t * ((v + zt^2) / (v - 1)) * dt(zt, df = v) / alpha)

	  list(VaR_in = VaR_in, ES_in = ES_in, VaR_oos = VaR_oos, ES_oos = ES_oos, mu_oos = mu_oos)
}



Hansen_qt <- function(u, v, lambda){			    
  	c <- gamma((v+1)/2)/( gamma(v/2)*sqrt(pi * (v-2) ) )    
  	a <- 4*lambda*c*( (v-2)/(v-1) )
  	b <- sqrt(1 + 3*lambda^2 - a^2)
  
  	if (u < (1-lambda)/2 ){        
    	answer <- ( (1-lambda) * sqrt( (v-2)/v ) * qt(u/(1-lambda), v) - a )/b    
	}else{        
    	answer <- ( (1+lambda) * sqrt( (v-2)/v ) * qt((u + lambda)/(1+lambda), v) - a )/b    
  	}
  	return(answer)
}



Hansen_pt <- function(z, v, lambda){	
  	c <- gamma((v+1)/2)/( gamma(v/2)*sqrt(pi * (v-2) ) )    
  	a <- 4*lambda*c*( (v-2)/(v-1) )
  	b <- sqrt(1 + 3*lambda^2 - a^2)
  
  	if (z < -a/b ){
    	x <- sqrt(v/(v-2) ) * ( b*z + a )/(1-lambda)     			        
    	answer <-  (1-lambda) * pt(x, v)    
  	} else {
    	x <- sqrt(v/(v-2) ) * ( b*z + a )/(1+lambda)         
    	answer <-  (1+lambda) * pt(x, v) - lambda   
  	}
  	return(answer)
}



Hansen_ES <- function(alpha, v, lambda){
	c <- gamma((v+1)/2)/( gamma(v/2)*sqrt(pi * (v-2) ) )    
  	a <- 4*lambda*c*( (v-2)/(v-1) )
  	b <- sqrt(1 + 3*lambda^2 - a^2)

    mode_cut <- -a/b
  	VaR_z <- Hansen_qt(alpha, v, lambda)     

  	if (VaR_z < mode_cut){    	  		  			
    	alphatil <- Hansen_pt( (b/(1-lambda)) * (VaR_z + a/b), v, 0 )    	  	
    	qq <- qt(alphatil, v)
    	ES_y <- -(1/alphatil)*dt(qq, v)*((v+qq^2)/(v-1))*sqrt((v-2)/v)    	    	
    	ES_z <- (alphatil/alpha) * (1-lambda) * ( -a/b + ((1-lambda)/b)*ES_y )    	
  	} else {    	
    	lambda2 <- -lambda
	    c2 <- c
	    a2 <- 4*lambda2*c2*( (v-2)/(v-1) )
	    b2 <- sqrt(1 + 3*lambda2^2 - a2^2)    	
	    VaR_z2 <- Hansen_qt(1-alpha, v, lambda2)     	
	    alphatil2 <- Hansen_pt( (b2/(1-lambda2)) * (VaR_z2 + a2/b2), v, 0 )		
	    qq2 <- qt(alphatil2, v)	    
	    ES_y2 <- -(1/alphatil2)*dt(qq2, v)*((v+qq2^2)/(v-1))*sqrt((v-2)/v)    				
		ES_z <- (alphatil2 / alpha) * (1 - lambda2) * ( -a2/b2 + ((1 - lambda2)/b2) * ES_y2 )	       	
	}
	return(list(VaR=VaR_z, ES=ES_z))    	
}



GARCH_skew_t <- function(returns, location = "zero", alpha){		
  	Tsize <- length(returns)

  	if (location == "zero"){
    	muhat  <- rep(0, Tsize)
    	mu_oos <- 0
  	} else if (location == "AR"){
    	tmp    <- AR1(returns)
    	muhat  <- tmp$meanhat
    	mu_oos <- tmp$meanf
  	} else {
    	stop("location must be 'zero' or 'AR'")
  	}

  	if (!(alpha > 0 && alpha < 1)){
    	stop("alpha must be in (0, 1). For 99% VaR use alpha = 0.01, for 95% use 0.05, etc.")
  	}

  	centeredreturns <- returns - muhat
  	var0 <- var(centeredreturns)

  	beta_to_theta <- function(beta){
    	beta0 <- beta[1]; b1 <- beta[2]; b2 <- beta[3]; v <- beta[4]; lambda <- beta[5]
    	den <- 1 - b1 - b2
    	c  <- log(beta0)
	    a  <- log(b1 / den)
	    b  <- log(b2 / den)
	    g  <- log(v - 2)
	    ell <- qlogis((lambda + 1)/2)
	    c(c = c, a = a, b = b, g = g, ell = ell)
	}

	beta_start <- c(beta0 = 0.01, beta1 = 0.8, beta2 = 0.1, v = 30, lambda = 0)
  	pars0 <- beta_to_theta(beta_start)

  	fit <- optim(par=pars0, fn=GARCH_skew_t_Cpp, centeredreturns = centeredreturns, var0=var0, control=list(fnscale = -1))
  	pars1 <- fit$par
    v <- 2 + exp(pars1[4])
    lambda <- 2 * ( exp(pars1[5])/(1 + exp(pars1[5])) ) - 1
  	  
  	sigma2_path <- GARCH_skew_t_forecast_Cpp(centeredreturns= centeredreturns, pars=pars1, var0=var0)
  	sigma2_in   <- sigma2_path[1:Tsize]
  	sigma2_f    <- sigma2_path[Tsize + 1]
  	sigma_in    <- sqrt(sigma2_in)
  	sigma_f     <- sqrt(sigma2_f)

  	hz <- Hansen_ES(alpha, v, lambda)
  	VaR_in  <- as.numeric(muhat + sigma_in * hz$VaR)
  	ES_in   <- as.numeric(muhat + sigma_in * hz$ES)
  	VaR_oos <- as.numeric(mu_oos + sigma_f * hz$VaR)
  	ES_oos  <- as.numeric(mu_oos + sigma_f * hz$ES)

  	list(VaR_in = VaR_in, ES_in = ES_in, VaR_oos = VaR_oos, ES_oos = ES_oos, mu_oos = mu_oos)
}



GJR_norm <- function(returns, location = "zero", alpha){
  	Tsize <- length(returns)

  	if (location == "zero"){
    	muhat  <- rep(0, Tsize)
    	mu_oos <- 0
  	} else if (location == "AR"){
    	tmp    <- AR1(returns)
    	muhat  <- tmp$meanhat
    	mu_oos <- tmp$meanf
  	} else {
    	stop("location must be 'zero' or 'AR'")
  	}

  	if (!(alpha > 0 && alpha < 1)){
    	stop("alpha must be in (0, 1). For 99% VaR use alpha = 0.01, for 95% use 0.05, etc.")
  	}

  	centeredreturns <- returns - muhat
  	var0 <- var(centeredreturns)

  	beta_to_theta <- function(beta){
	    beta0 <- beta[1]; b1 <- beta[2]; b2 <- beta[3]; b3 <- beta[4]
	    den <- 1 / max(1e-12, 1 - b1 - b2 - b3 / 2)
	    s1 <- b1*den
	    s2 <- b2*den
	    s3 <- (b3*den)/2
	    c <- log(beta0)
   		a <- log(s1)
    	b <- log(s2)
    	d <- log(s3)
    	c(c = c, a = a, b = b, d = d)
  	}
  
  	beta_start <- c(beta0 = 0.01, beta1 = 0.8, beta2 = 0.1, beta3 = 0.0001) 
  	pars0 <- beta_to_theta(beta_start)

  	fit <- optim(par = pars0, fn = GJR_norm_Cpp, centeredreturns = centeredreturns, var0 = var0, control = list(fnscale = -1))
  	pars1 <- fit$par

  	sigma2_path     <- GJR_norm_forecast_Cpp(centeredreturns= centeredreturns, pars=pars1, var0=var0)
  	sigma2_in       <- sigma2_path[1:Tsize]
  	sigma2_forecast <- sigma2_path[Tsize + 1]
	sigma_in        <- sqrt(sigma2_in)
  	sigma_forecast  <- sqrt(sigma2_forecast)

  	z   <- qnorm(alpha)
  	VaR_in  <- as.numeric(muhat  + sigma_in * z)
  	ES_in   <- as.numeric(muhat  - sigma_in * dnorm(z) / alpha)
  	VaR_oos <- as.numeric(mu_oos + sigma_forecast * z)
  	ES_oos  <- as.numeric(mu_oos - sigma_forecast * dnorm(z) / alpha)

  	list(VaR_in = VaR_in, ES_in = ES_in, VaR_oos = VaR_oos, ES_oos = ES_oos, mu_oos = mu_oos)
}



GJR_t <- function(returns, location = "zero", alpha){
  	Tsize <- length(returns)

  	if (location == "zero"){
    	muhat  <- rep(0, Tsize)
    	mu_oos <- 0
  	} else if (location == "AR"){
    	tmp    <- AR1(returns)
    	muhat  <- tmp$meanhat
    	mu_oos <- tmp$meanf
  	} else {
    	stop("location must be 'zero' or 'AR'")
  	}

  	if (!(alpha > 0 && alpha < 1)){
    	stop("alpha must be in (0, 1). For 99% VaR use alpha = 0.01, for 95% use 0.05, etc.")
  	}

  	centeredreturns <- returns - muhat
  	var0 <- var(centeredreturns)

	beta_to_theta <- function(beta){
	  beta0 <- beta[1]; b1 <- beta[2]; b2 <- beta[3]; b3 <- beta[4]; v <- beta[5]
	  rho <- b1 + b2 + 0.5 * b3
	  inv <- 1 / pmax(1e-12, 1 - rho)
	  s1  <- b1 * inv
	  s2  <- b2 * inv
	  s3  <- b3 * inv
	  eps <- 1e-8
	  c   <- log(beta0)
	  a   <- log(pmax(s1, eps))
	  b   <- log(pmax(s2, eps))
	  d   <- log(pmax(s3, eps))
	  g   <- log(pmax(v - 2, eps))
	  c(c = c, a = a, b = b, d = d, g = g)
	}

  	beta_start <- c(beta0 = 0.01, beta1 = 0.8, beta2 = 0.1, beta3 = 0, v = 30)
  	pars0 <- beta_to_theta(beta_start)

  	fit <- optim(par = pars0, fn = GJR_t_Cpp, centeredreturns = centeredreturns, var0 = var0, control = list(fnscale = -1))
  	pars1 <- fit$par
    v <- 2 + exp(pars1[5])
  	
	sigma2_path     <- GJR_t_forecast_Cpp(centeredreturns= centeredreturns, pars=pars1, var0=var0)
  	sigma2_in       <- sigma2_path[1:Tsize]
  	sigma2_forecast <- sigma2_path[Tsize + 1]
  	sigma_in        <- sqrt(sigma2_in)
  	sigma_forecast  <- sqrt(sigma2_forecast)

  	## VaR / ES for standardized t with unit variance (scale = sqrt((v-2)/v))
  	zt      <- qt(alpha, df = v)
  	scale_t <- sqrt((v - 2) / v)

  	VaR_in  <- as.numeric(muhat  + sigma_in       * zt * scale_t)
  	ES_in   <- as.numeric(muhat  - sigma_in       * scale_t * ((v + zt^2) / (v - 1)) * dt(zt, df = v) / alpha)
  	VaR_oos <- as.numeric(mu_oos + sigma_forecast * zt * scale_t)
  	ES_oos  <- as.numeric(mu_oos - sigma_forecast * scale_t * ((v + zt^2) / (v - 1)) * dt(zt, df = v) / alpha)

  	list(VaR_in = VaR_in, ES_in = ES_in, VaR_oos = VaR_oos, ES_oos = ES_oos, mu_oos = mu_oos)
}



GJR_skew_t <- function(returns, location = "zero", alpha){
 	Tsize <- length(returns)

  	if (location == "zero"){
    	muhat  <- rep(0, Tsize)
    	mu_oos <- 0
  	} else if (location == "AR"){
    	tmp    <- AR1(returns)
    	muhat  <- tmp$meanhat
    	mu_oos <- tmp$meanf
  	} else {
    	stop("location must be 'zero' or 'AR'")
  	}

  	if (!(alpha > 0 && alpha < 1)){
    	stop("alpha must be in (0,1)")
  	}

  	centeredreturns <- returns - muhat
  	var0 <- var(centeredreturns)

  	beta_to_theta <- function(beta){
    	beta0 <- beta[1]; b1 <- beta[2]; b2 <- beta[3]; b3 <- beta[4]; v <- beta[5]; lambda <- beta[6]
   		rho <- b1 + b2 + 0.5*b3
    	den <- 1 / (1 - rho)
    	s1 <- b1 * den
    	s2 <- b2 * den
    	s3 <- b3 * den
    	c  <- log(beta0)
    	a  <- log(pmax(1e-12, s1))
    	b  <- log(pmax(1e-12, s2))
    	d  <- log(pmax(1e-12, s3))
    	g  <- log(pmax(1e-12, v - 2))
    	ell <- qlogis((lambda + 1)/2)
    	c(c=c, a=a, b=b, d=d, g=g, ell=ell)
  	}

  	beta_start <- c(beta0=0.01, beta1=0.8, beta2=0.1, beta3=0, v=30, lambda=0)
  	pars0 <- beta_to_theta(beta_start)

  	fit <- optim(par = pars0, fn = GJR_skew_t_Cpp, centeredreturns = centeredreturns, var0 = var0, control = list(fnscale = -1))
  	pars1 <- fit$par
    v <- 2 + exp(pars1[5])
    lambda <- 2 * ( exp(pars1[6])/(1 + exp(pars1[6])) ) - 1

  	sigma2_path <- GJR_skew_t_forecast_Cpp(centeredreturns= centeredreturns, pars=pars1, var0=var0)

  	sigma2_in       <- sigma2_path[1:Tsize]
  	sigma2_forecast <- sigma2_path[Tsize + 1]
  	sigma_in        <- sqrt(sigma2_in)
  	sigma_forecast  <- sqrt(sigma2_forecast)

	# VaR / ES for Hansen skew-t (already standardized to Var=1)
  	hz <- Hansen_ES(alpha, v, lambda)  # returns list(VaR=z_alpha, ES=ES_alpha) on unit scale
  	VaR_in  <- as.numeric(muhat  + sigma_in       * hz$VaR)
  	ES_in   <- as.numeric(muhat  + sigma_in       * hz$ES)
  	VaR_oos <- as.numeric(mu_oos + sigma_forecast * hz$VaR)
  	ES_oos  <- as.numeric(mu_oos + sigma_forecast * hz$ES)

  	list(VaR_in = VaR_in, ES_in = ES_in, VaR_oos = VaR_oos, ES_oos = ES_oos, mu_oos = mu_oos)
}



EGARCH <- function(returns, location = "zero", alpha){			
  	Tsize <- length(returns)

  	if (location == "zero"){
    	muhat  <- rep(0, Tsize)
    	mu_oos <- 0
  	} else if (location == "AR"){
    	tmp    <- AR1(returns)
    	muhat  <- tmp$meanhat
    	mu_oos <- tmp$meanf
  	} else {
    	stop("location must be 'zero' or 'AR'")
  	}

  	if (!(alpha > 0 && alpha < 1)){
    	stop("alpha must be in (0, 1). For 99% VaR use alpha = 0.01, for 95% use 0.05, etc.")
  	}
  
  	centeredreturns <- returns - muhat
  	var0 <- var(centeredreturns)
  
	beta0_target <- 0.95
	b0_raw <- qlogis(beta0_target)   
	w0     <- (1 - beta0_target) * log(var0)  # so E[log σ²] ≈ log(var0)
	a0 <- 0.1
	g0 <- 0
  	pars0 <- c(w0, b0_raw, a0, g0)
  
  	fit <- optim(par = pars0, fn = EGARCH_Cpp, centeredreturns = centeredreturns, var0 = var0, control = list(fnscale = -1))
  	pars1 <- fit$par
  
  	sigma2_path     <- EGARCH_forecast_Cpp(centeredreturns= centeredreturns, pars=pars1, var0=var0)
  	sigma2_in       <- sigma2_path[1:Tsize]
  	sigma2_forecast <- sigma2_path[Tsize + 1]
  	sigma_in        <- sqrt(sigma2_in)
  	sigma_forecast  <- sqrt(sigma2_forecast)

  	z  <- qnorm(alpha)
  	VaR_in  <- as.numeric(muhat  + sigma_in * z)
  	ES_in   <- as.numeric(muhat  - sigma_in * dnorm(z) / alpha)
  	VaR_oos <- as.numeric(mu_oos + sigma_forecast * z)
  	ES_oos  <- as.numeric(mu_oos - sigma_forecast * dnorm(z) / alpha)

  	list(VaR_in = VaR_in, ES_in = ES_in, VaR_oos = VaR_oos, ES_oos = ES_oos, mu_oos = mu_oos)
}	







