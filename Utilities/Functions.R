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




AR1 <- function(rets){	
  Tsize <- length(rets) 
  temp <- lm(rets[2:Tsize]~rets[1:(Tsize-1)])
  muhat <- as.numeric(temp$coef[1])
  phihat <- as.numeric(temp$coef[2]) 
  meanhat <- matrix(0, Tsize, 1)
  meanhat[2:Tsize] <- muhat + phihat*rets[1:(Tsize-1)]	
  meanf <- muhat + phihat*rets[Tsize]	
  return(list(meanhat=meanhat, meanf=meanf))		
}




QbSD <- function(rets, location, CAViaR, alpha){

	Scale <- function(centeredrets, p){

		probs <- c(p, 1-p) 
		Q0 <- quantile(centeredrets, probs)

		while ((Q0[1] >= 0) | (Q0[2] <= 0)){				
			p <- p - 0.05
			probs <- c(p, 1-p) 
			Q0 <- quantile(centeredrets, probs)				
		}
												
		if (CAViaR == "SAV"){
								
			beta0   <- seq(0.1, 0.9, length.out=3)
			gamma0  <- seq(0.1, 0.9, length.out=3)								
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
				
			idx <- which.min(results[,5])				
			pars1 <- results[idx, 1:4]
																		
			SAVtemp <- SAVscalehatf_Cpp(centeredrets=centeredrets, q0=Q0, pars=pars1)

			Q1hat   <- SAVtemp[1:Tsize, 1]				
			Q1f     <- SAVtemp[Tsize+1, 1]
																
			Q2hat   <- SAVtemp[1:Tsize, 2]				
			Q2f     <- SAVtemp[Tsize+1, 2]				
		}	

		if (CAViaR == "AS"){	
								
			beta0   <- seq(0.1, 0.9, length.out=3)
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

	VaR <- function(){			
		QQ0 <- quantile(epsilon, alpha)		
		varhat <- muhat + scalehat* QQ0
		varf <- mu_oos + scalef* QQ0
		return(list(varhat=varhat, varf=varf))
	}

	ES <- function(){										
		TargetTaus <- seq(from=alpha/N, to=alpha, by=alpha/N)							
		Nsize <- length(TargetTaus)
		QQ0 <- matrix(quantile(epsilon, TargetTaus), Nsize, 1)		
		QQhat <- matrix(0, Tsize, Nsize)		
		QQf   <- matrix(0, 1, Nsize)		

		for (i in 1:Nsize){
			QQhat[,i] <- muhat + scalehat* QQ0[i]
			QQf[,i]   <- mu_oos + scalef* QQ0[i]				
		}	

		eshat <- apply(QQhat, 1, mean)	
		esf   <- mean(QQf)
		return(list(eshat=eshat, esf=esf))
	}
		
  
	Tsize <- length(rets)
			
	P <- seq(from=0.05, to=0.25, by=0.05)
		
	if (location == "zero"){ 
		muhat <- matrix(0, Tsize, 1)		
		mu_oos   <- 0
	}	

	if (location == "QAR"){
		temp <- QAR1(rets)
		muhat <- temp$medianhat
		mu_oos   <- temp$medianf		
	}	
	
	centeredrets <- rets - muhat
		
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

	VaR_in   <- apply(VaRp[,1:length(P)], 1, mean)
	VaR_oos  <- mean(VaRfp)
  
	ESp    <- matrix(0, Tsize, length(P))	
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
		
	return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos, mu_oos=mu_oos))	
}	




Taylor <- function(rets, location, CAViaR, ES, alpha){
  
  Tsize <- length(rets)	
  
  if (location == "zero"){ 
    muhat <- matrix(0, Tsize, 1)		
    mu_oos   <- 0
  }	
  
  if (location == "AR"){
    temp <- AR1(rets)
    muhat <- temp$meanhat
    mu_oos   <- temp$meanf		
  }	
  
  prob <- alpha
  
  if ((CAViaR == "SAV") & (ES == "basic")){	
    
    Q0 <- quantile(rets, prob)		
    beta10 <- seq(0.1, 0.9, length.out=3)
    beta20 <- seq(0.1, 0.9, length.out=3)    
    results <- matrix(0, length(beta10)*length(beta20), 4)
    
    kk <- 1    
    for (i in beta10){
      for (j in beta20){ 
        
        pars0 <- c(Q0, i, j)			        
        temp <- optim(pars0, SAV_Cpp, rets=rets, mu=muhat, tau=prob, q0=Q0)				
        pars1 <- temp$par		
        
        results[kk, 1:3] <- pars1
        results[kk, 4] <- temp$value        
        kk <- kk + 1						
      }
    }		
    
    idx <- which.min(results[,4])				
    betahat0 <- results[idx, 1:3]
    
    VaRhatf <- SAVhatf_Cpp(rets=rets, mu=muhat, tau=prob, q0=Q0, pars=betahat0)				
    VaRhat  <- VaRhatf[1:Tsize]				
    
    temp <- optimize(f=Basic_Cpp, interval=c(-10,10), rets=rets, mu=muhat, alpha= prob, VaR=VaRhat, maximum=TRUE)
    gammahat0 <- temp$maximum		
        
    pars0 <- c(betahat0, gammahat0)			
    
    temp <- optim(pars0, Taylor_SAV_Basic_Cpp, rets=rets, mu=muhat, alpha=prob, q0=Q0, control=list(fnscale=-1))				
    pars1 <- temp$par		
    
    temp <- Taylor_SAV_BasicHatf_Cpp(rets=rets, mu=muhat, alpha=prob, q0=Q0, pars=pars1 )				
    VaR_in  <- temp[1:Tsize, 1]				
    VaR_oos <- temp[Tsize+1, 1]	
    ES_in   <- temp[1:Tsize]		
    ES_oos  <- temp[Tsize+1]    
  }
  
  
  if ((CAViaR == "AS") & (ES == "basic")){	
    
    Q0 <- quantile(rets, prob)		
    beta10 <- seq(0.1, 0.9, length.out=3)
    beta20 <- seq(0.1, 0.9, length.out=3)
    beta30 <- seq(0.1, 0.9, length.out=3)    
    results <- matrix(0, length(beta10)*length(beta20)*length(beta30), 5)
    
    kk <- 1    
    for (i in beta10){
      for (j in beta20){ 
        for (k in beta30 ){
          
          pars0 <- c(Q0, i, j, k)		          
          temp <- optim(pars0, AS_Cpp, rets=rets, mu=muhat, tau=prob, q0=Q0)		          
          pars1 <- temp$par
          
          results[kk,1:4] <- pars1
          results[kk, 5] <- temp$value
          
          kk <- kk + 1																					
        }							
      }					
    }
    
    idx <- which.min(results[,5])				
    betahat0 <- results[idx, 1:4]
    
    VaRhatf <- AShatf_Cpp(rets=rets, mu=muhat, tau=prob, q0=Q0, pars= betahat0)			
    VaRhat  <- VaRhatf[1:Tsize]				
    
    temp <- optimize(f=Basic_Cpp, interval=c(-10,10), rets=rets, mu=muhat, alpha= prob, VaR=VaRhat, maximum=TRUE)
    gammahat0 <- temp$maximum		
    
    pars0 <- c(betahat0, gammahat0)			
    
    temp <- optim(pars0, Taylor_AS_Basic_Cpp, rets=rets, mu=muhat, alpha=prob, q0=Q0, control=list(fnscale=-1))				
    pars1 <- temp$par		
    
    temp <- Taylor_AS_BasicHatf_Cpp(rets=rets, mu=muhat, alpha=prob, q0=Q0, pars=pars1 )				
    VaR_in  <- temp[1:Tsize, 1]				
    VaR_oos <- temp[Tsize+1, 1]	
    ES_in   <- temp[1:Tsize]		
    ES_oos  <- temp[Tsize+1]    
  }
  
  
  
  if ((CAViaR == "SAV") & (ES == "AR")){	
    
    Q0 <- quantile(rets, prob)		
    beta10 <- seq(0.1, 0.9, length.out=3)
    beta20 <- seq(0.1, 0.9, length.out=3)
    
    results <- matrix(0, length(beta10)*length(beta20), 4)
    
    kk <- 1
    
    for (i in beta10){
      for (j in beta20){ 
        
        pars0 <- c(Q0, i, j)			        
        temp <- optim(pars0, SAV_Cpp, rets=rets, mu=muhat, tau=prob, q0=Q0)				
        pars1 <- temp$par		
        
        results[kk,1:3] <- pars1
        results[kk, 4] <- temp$value        
        kk <- kk + 1						
      }
    }		
    
    idx <- which.min(results[,4])				
    betahat0 <- results[idx, 1:3]
    
    VaRhatf <- SAVhatf_Cpp(rets=rets, mu=muhat, tau=prob, q0=Q0, pars=betahat0)				
    VaRhat  <- VaRhatf[1:Tsize]				
    
    gamma00 <- seq(0.001, 0.1, length.out=3)
    gamma10 <- seq(0.1, 0.9, length.out=3)
    gamma20 <- seq(0.1, 0.9, length.out=3)
    
    results <- matrix(0, length(gamma00)*length(gamma10)*length(gamma20), 4)
    
    kk <- 1
    
    for (i in gamma00){
      for (j in gamma10){ 
        for (k in gamma20 ){
          
          pars0 <- c(log(i), log(j), log(k))	          
          temp <- optim(pars0, AR_Cpp, rets=rets, mu=muhat, alpha=prob, VaR=VaRhat, control=list(fnscale=-1))
          pars1 <- temp$par																	
          
          results[kk,1:3] <- pars1
          results[kk, 4] <- temp$value          
          kk <- kk + 1																					
        }							
      }					
    }
    
    idx <- which.max(results[,4])				
    gammahat0 <- exp(results[idx, 1:3])
    
    tol <- 0.00001
    gammahat0[gammahat0 <= tol] <- tol
    
    pars0 <- c(betahat0, log(gammahat0))			
    
    temp <- optim(pars0, Taylor_SAV_AR_Cpp, rets=rets, mu=muhat, alpha=prob, q0=Q0, control=list(fnscale=-1))				
    pars1 <- temp$par
    pars1[4:6] <- exp(pars1[4:6])	
    
    temp <- Taylor_SAV_ARHatf_Cpp(rets=rets, mu=muhat, alpha=prob, q0=Q0, pars=pars1 )				
    VaR_in  <- temp[1:Tsize, 1]				
    VaR_oos <- temp[Tsize+1, 1]	
    ES_in   <- temp[1:Tsize]		
    ES_oos  <- temp[Tsize+1]    
  }
  
  
  if ((CAViaR == "AS") & (ES == "AR")){	
    
    
    Q0 <- quantile(rets, prob)		
    beta10 <- seq(0.1, 0.9, length.out=3)
    beta20 <- seq(0.1, 0.9, length.out=3)
    beta30 <- seq(0.1, 0.9, length.out=3)
    
    results <- matrix(0, length(beta10)*length(beta20)*length(beta30), 5)
    
    kk <- 1
    
    for (i in beta10){
      for (j in beta20){ 
        for (k in beta30 ){
          
          pars0 <- c(Q0, i, j, k)		
          
          temp <- optim(pars0, AS_Cpp, rets=rets, mu=muhat, tau=prob, q0=Q0)		
          
          pars1 <- temp$par
          
          results[kk,1:4] <- pars1
          results[kk, 5] <- temp$value
          
          kk <- kk + 1																					
        }							
      }					
    }
    
    idx <- which.min(results[,5])				
    betahat0 <- results[idx, 1:4]
    
    VaRhatf <- AShatf_Cpp(rets=rets, mu=muhat, tau=prob, q0=Q0, pars= betahat0)			
    VaRhat  <- VaRhatf[1:Tsize]				
    
    
    gamma00 <- seq(0.001, 0.1, length.out=3)
    gamma10 <- seq(0.1, 0.9, length.out=3)
    gamma20 <- seq(0.1, 0.9, length.out=3)
    
    results <- matrix(0, length(gamma00)*length(gamma10)*length(gamma20), 4)
    
    kk <- 1
    
    for (i in gamma00){
      for (j in gamma10){ 
        for (k in gamma20 ){
          
          pars0 <- c(log(i), log(j), log(k))	
          
          temp <- optim(pars0, AR_Cpp, rets=rets, mu=muhat, alpha=prob, VaR=VaRhat, control=list(fnscale=-1))
          pars1 <- temp$par																	
          
          results[kk,1:3] <- pars1
          results[kk, 4] <- temp$value
          
          kk <- kk + 1																					
        }							
      }					
    }
    
    idx <- which.max(results[,4])				
    gammahat0 <- exp(results[idx, 1:3])
    
    tol <- 0.00001
    gammahat0[gammahat0 <= tol] <- tol
    
    
    pars0 <- c(betahat0, log(gammahat0))			
    
    temp <- optim(pars0, Taylor_AS_AR_Cpp, rets=rets, mu=muhat, alpha=prob, q0=Q0, control=list(fnscale=-1))				
    pars1 <- temp$par
    pars1[5:7] <- exp(pars1[5:7])	
    
    temp <- Taylor_AS_ARHatf_Cpp(rets=rets, mu=muhat, alpha=prob, q0=Q0, pars=pars1 )				
    VaR_in  <- temp[1:Tsize, 1]				
    VaR_oos <- temp[Tsize+1, 1]	
    ES_in   <- temp[1:Tsize]		
    ES_oos  <- temp[Tsize+1]    
  }
  
  return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos, mu_oos=mu_oos))	    
}	




Patton <- function(rets, alpha){
  
  Tsize <- length(rets)		
  
  beta0 <- seq(0.9, 0.99, length.out=3)	
  gamma0 <- seq(0.001, 0.1, length.out=3)
  
  results <- matrix(0, length(beta0)*length(gamma0), 5)
  
  kk <- 1
  
  for (i in beta0){
    for (j in gamma0){
      
      pars0 <- c(0, 0, log(i/(1-i)), log(j))
      
      temp <- optim(par=pars0, fn=PattonGAS_Cpp, rets=rets, alpha=alpha)
      pars1 <- temp$par
      
      results[kk,1:4] <- pars1
      results[kk, 5] <- temp$value
      
      kk <- kk + 1																					
    }							
  }					
  
  idx <- which.min(results[,5])				
  pars1 <- results[idx, 1:4]
  
  results <- PattonGAShatf_Cpp(rets=rets, pars=pars1, alpha=alpha)   
  VaR_in  <- results[1:Tsize, 1]				
  VaR_oos <- results[Tsize+1, 1]		  
  ES_in   <- results[1:Tsize, 2]		
  ES_oos  <- results[Tsize+1, 2]
    
  return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos))	    
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
    answer =  (1-lambda) * pt(x, v)    
  }else{
    x <- sqrt(v/(v-2) ) * ( b*z + a )/(1+lambda)         
    answer =  (1+lambda) * pt(x, v) - lambda   
  }
  return(answer)
}




Hansen_ES <- function(alpha, v, lambda){
  c <- gamma((v+1)/2)/( gamma(v/2)*sqrt(pi * (v-2) ) )    
  a <- 4*lambda*c*( (v-2)/(v-1) )
  b <- sqrt(1 + 3*lambda^2 - a^2)
  VaR_z <- Hansen_qt(alpha, v, lambda)     
  if (VaR_z <= -a/b){    	
    alphatil <- Hansen_pt( (b/(1-lambda)) * (VaR_z + a/b), v, 0 )    	  	
    qq <- qt(alphatil, v)
    ES_y <- -(1/alphatil)*dt(qq, v)*((v+qq^2)/(v-1))*sqrt((v-2)/v)    	    	
    ES_z <- (alphatil/alpha) * (1-lambda) * ( -a/b + ((1-lambda)/b)*ES_y )    	
  }else{    	
    lambda2 <- -lambda
    v2 <- v
    c2 <- c
    a2 <- 4*lambda2*c2*( (v2-2)/(v2-1) )
    b2 <- sqrt(1 + 3*lambda2^2 - a2^2)    	
    VaR_z <- Hansen_qt(1-alpha, v2, lambda2)     	
    alphatil <- Hansen_pt( (b2/(1-lambda2)) * (VaR_z + a2/b2), v2, 0 )		
    qq <- qt(alphatil, v2)
    ES_y <- -(1/alphatil)*dt(qq, v)*((v2+qq^2)/(v2-1))*sqrt((v2-2)/v2)    	
    ES_z <- (alphatil/alpha) * (1-lambda2) * ( -a2/b2 + ((1-lambda2)/b2)*ES_y )    	
  }
  return(list(VaR=VaR_z, ES=ES_z))    
}




Hansen <- function(rets, location, alpha){
      
  Tsize <- length(rets)	
  
  if (location == "zero"){ 
    muhat <- matrix(0, Tsize, 1)		
    mu_oos   <- 0
  }	
  
  if (location == "AR"){
    temp <- AR1(rets)
    muhat <- temp$meanhat
    mu_oos   <- temp$meanf		
  }	
  
  pars0 <- c(log(0.01), log(0.1), log(0.85), log(30), 0)		
  
  temp <- optim(pars0, GARCHskewt_Cpp, rets=rets, mu = muhat, control=list(fnscale=-1))
  pars1 <- temp$par	
  
  sigma2hatf <- GARCHskewtHatf_Cpp(rets, muhat, pars1) 	
  sigma2hat <- sigma2hatf[1:Tsize]
  sigma2f <- sigma2hatf[Tsize]
  
  beta0 <- pars1[1]
  beta1 <- pars1[2]
  beta2 <- pars1[3]
  v <- pars1[4]
  lambda <- pars1[5]
    
  beta0 <- exp(beta0);
  beta1 <- exp(beta1);
  beta2 <- exp(beta2);
  v <- 2 + exp(v);
  lambda <- 2*( exp(lambda)/(1+exp(lambda)) ) -1;
  
  temp <- Hansen_ES(alpha, v, lambda)
  
  VaR_in <- muhat + sqrt(sigma2hat)*temp$VaR
  ES_in  <- muhat + sqrt(sigma2hat)*temp$ES
  
  VaR_oos <- mu_oos + sqrt(sigma2f)*temp$VaR
  ES_oos  <- mu_oos + sqrt(sigma2f)*temp$ES
    
  return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos, mu_oos=mu_oos))	        
}	




GARCHnorm <- function(rets, location, alpha){
  
  Tsize <- length(rets)
  if (location == "zero"){ 
    muhat <- matrix(0, Tsize, 1)		
    mu_oos   <- 0
  }	
  
  if (location == "AR"){
    temp <- AR1(rets)
    muhat <- temp$meanhat
    mu_oos   <- temp$meanf		
  }
  
  centeredrets <- rets - muhat
    
  pars0 <- c(  0.01, 0.7, 0.1)		
  
  temp <- optim(par=pars0, GARCHnorm_Cpp, centeredrets=centeredrets, control=list(fnscale=-1))
  pars1 <- temp$par	
  sigma2hatf <- GARCHnormhatf_Cpp(centeredrets, pars1) 
  sigma2hat <- sigma2hatf[1:Tsize]
  sigma2f <- sigma2hatf[Tsize]
  sigmahat <- sqrt(sigma2hat)
  sigmaf <- sqrt(sigma2f)
  
  VaRhat <- muhat + sigmahat*qnorm(alpha)
  VaR_in <- as.numeric(VaRhat)
  
  EShat <- muhat - sigmahat*dnorm(qnorm(alpha))/alpha
  ES_in <- as.numeric(EShat)
  
  VaR_oos <- mu_oos + sigmaf*qnorm(alpha)
  ES_oos  <- mu_oos - sigmaf*dnorm(qnorm(alpha))/alpha
    
  return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos, mu_oos=mu_oos))	    
}	




GARCHt <- function(rets, location, alpha){
  Tsize <- length(rets)	
  if (location == "zero"){ 
    muhat <- matrix(0, Tsize, 1)		
    mu_oos   <- 0
  }	
  
  if (location == "AR"){
    temp <- AR1(rets)
    muhat <- temp$meanhat
    mu_oos   <- temp$meanf		
  }
  
  centeredrets = rets - muhat;
  
  pars0 <- c(0.01 , 0.7, 0.1, log(10))		
  
  temp <- optim(par=pars0, GARCHt_Cpp, centeredrets=centeredrets, control=list(fnscale=-1))
  pars1 <- temp$par	
  sigma2hatf <- GARCHthatf_Cpp(centeredrets, pars1) 
  sigma2hat <- sigma2hatf[1:Tsize]
  sigma2f <- sigma2hatf[Tsize]
  sigmahat <- sqrt(sigma2hat)
  sigmaf <- sqrt(sigma2f)
  
  v <- 2+exp(pars1[length(pars1)])
  
  VaRhat <- muhat + sigmahat*qt(alpha, df=v) * sqrt((v-2)/v)
  VaR_in <- as.numeric(VaRhat)
  VaR_oos <- mu_oos + sigmaf*qt(alpha, df=v)* sqrt((v-2)/v)		
  
  qstar <- qt(alpha, df=v) #* sqrt((v-2)/v)		
  EShat <- muhat - sigmahat * sqrt((v-2)/v) *  ( (v + qstar^2) / (v-1) )*dt(qstar, df=v)/alpha
  ES_in <- as.numeric(EShat)  
  ES_oos <- mu_oos - sigmaf* sqrt((v-2)/v)*  ( (v + qstar^2) / (v-1) )*dt(qstar, df=v)/alpha
  
  return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos, mu_oos=mu_oos))	      
}	




EGARCHnorm <- function(rets, location, alpha){
  Tsize <- length(rets)
  if (location == "zero"){ 
    muhat <- matrix(0, Tsize, 1)		
    mu_oos   <- 0
  }	
  
  if (location == "AR"){
    temp <- AR1(rets)
    muhat <- temp$meanhat
    mu_oos   <- temp$meanf		
  }
  
  centeredrets = rets - muhat;
  
  
  pars0 <- c(log(var(centeredrets)), 0.5, 0.1, 0)		
  
  temp <- optim(par=pars0, EGARCHnorm_Cpp, centeredrets=centeredrets, control=list(fnscale=-1))
  pars1 <- temp$par	
  sigma2hatf <- EGARCHnormhatf_Cpp(centeredrets, pars1) 
  sigma2hat <- sigma2hatf[1:Tsize]
  sigma2f <- sigma2hatf[Tsize]
  sigmahat <- sqrt(sigma2hat)
  sigmaf <- sqrt(sigma2f)
  
  
  VaRhat <- muhat + sigmahat*qnorm(alpha)
  VaR_in <- as.numeric(VaRhat)
  
  EShat <- muhat - sigmahat*dnorm(qnorm(alpha))/alpha
  ES_in <- as.numeric(EShat)
  
  VaR_oos <- mu_oos + sigmaf*qnorm(alpha)
  ES_oos  <- mu_oos - sigmaf*dnorm(qnorm(alpha))/alpha
  
  return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos, mu_oos=mu_oos))	        
}	





GJRnorm <- function(rets, location, alpha){
  Tsize <- length(rets)
  if (location == "zero"){ 
    muhat <- matrix(0, Tsize, 1)		
    mu_oos   <- 0
  }	
  
  if (location == "AR"){
    temp <- AR1(rets)
    muhat <- temp$meanhat
    mu_oos   <- temp$meanf		
  }
  
  centeredrets = rets - muhat;
    
  pars0 <- c(0.01, 0.7, 0.1, 0)		
  
  temp <- optim(par=pars0, GJRnorm_Cpp, centeredrets=centeredrets, control=list(fnscale=-1))
  pars1 <- temp$par	
  sigma2hatf <- GJRnormhatf_Cpp(centeredrets, pars1) 	
  sigma2hat <- sigma2hatf[1:Tsize]
  sigma2f <- sigma2hatf[Tsize]
  sigmahat <- sqrt(sigma2hat)
  sigmaf <- sqrt(sigma2f)
  
  VaRhat <- muhat + sigmahat*qnorm(alpha)
  VaR_in <- as.numeric(VaRhat)
  
  EShat <- muhat - sigmahat*dnorm(qnorm(alpha))/alpha
  ES_in <- as.numeric(EShat)
  
  VaR_oos <- mu_oos + sigmaf*qnorm(alpha)
  ES_oos  <- mu_oos - sigmaf*dnorm(qnorm(alpha))/alpha
  
  return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos, mu_oos=mu_oos))	          	
}	




GJRt <- function(rets, location, alpha){
  Tsize <- length(rets)		
  if (location == "zero"){ 
    muhat <- matrix(0, Tsize, 1)		
    mu_oos   <- 0
  }	
  
  if (location == "AR"){
    temp <- AR1(rets)
    muhat <- temp$meanhat
    mu_oos   <- temp$meanf		
  }
  
  centeredrets = rets - muhat;
    
  pars0 <- c(0.01 , 0.7, 0.1, 0, log(10))		
  
  temp <- optim(par=pars0, GJRt_Cpp, centeredrets=centeredrets, control=list(fnscale=-1))
  pars1 <- temp$par	
  sigma2hatf <- GJRthatf_Cpp(centeredrets, pars1) 	
  sigma2hat <- sigma2hatf[1:Tsize]
  sigma2f <- sigma2hatf[Tsize]
  sigmahat <- sqrt(sigma2hat)
  sigmaf <- sqrt(sigma2f)
  
  v <- 2+exp(pars1[length(pars1)])
  
  VaRhat <- muhat + sigmahat*qt(alpha, df=v) * sqrt((v-2)/v)
  VaR_in <- as.numeric(VaRhat)
  VaR_oos <- mu_oos + sigmaf*qt(alpha, df=v)* sqrt((v-2)/v)		
  
  qstar <- qt(alpha, df=v) #* sqrt((v-2)/v)		
  EShat <- muhat - sigmahat * sqrt((v-2)/v) *  ( (v + qstar^2) / (v-1) )*dt(qstar, df=v)/alpha
  ES_in <- as.numeric(EShat)
  
  ES_oos <- mu_oos - sigmaf* sqrt((v-2)/v)*  ( (v + qstar^2) / (v-1) )*dt(qstar, df=v)/alpha
  
  return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos, mu_oos=mu_oos))	          	  
}	




GJRskewt <- function(rets, location, alpha){
    
  Tsize <- length(rets)	
  
  if (location == "zero"){ 
    muhat <- matrix(0, Tsize, 1)		
    mu_oos   <- 0
  }	
  
  if (location == "AR"){
    temp <- AR1(rets)
    muhat <- temp$meanhat
    mu_oos   <- temp$meanf		
  }	
    
  pars0 <- c(0.01, 0.7, 0.1, 0, log(10), 0)		
  
  temp <- optim(pars0, GJRskewt_Cpp, rets=rets, mu = muhat, control=list(fnscale=-1))
  pars1 <- temp$par	
  
  sigma2hatf <- GJRskewthatf_Cpp(rets, muhat, pars1) 	
  sigma2hat <- sigma2hatf[1:Tsize]
  sigma2f <- sigma2hatf[Tsize]
  
  v <- pars1[5]
  lambda <- pars1[6]
   
  v <- 2 + exp(v);
  lambda <- 2*( exp(lambda)/(1+exp(lambda)) ) -1;
  
  temp <- Hansen_ES(alpha, v, lambda)
  
  VaR_in <- muhat + sqrt(sigma2hat)*temp$VaR
  ES_in  <- muhat + sqrt(sigma2hat)*temp$ES
  
  VaR_oos <- mu_oos + sqrt(sigma2f)*temp$VaR
  ES_oos  <- mu_oos + sqrt(sigma2f)*temp$ES
  
  return(list(VaR_in=VaR_in, ES_in=ES_in, VaR_oos=VaR_oos, ES_oos=ES_oos, mu_oos=mu_oos))	          	      
}



