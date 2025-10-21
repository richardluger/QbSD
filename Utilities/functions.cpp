#include <RcppDist.h>
//[[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace arma;



//[[Rcpp::export]]
double check_Cpp(double u, double tau){
    int indic = 0;
    if (u < 0) indic = 1;
    return u*( tau - indic);
}



//[[Rcpp::export]]
double gSAV_Cpp(vec centeredreturns, double tau, vec q0, vec pars){
    double omega1 = pars(0), omega2 = pars(1), beta = exp(pars(2)), gamma = exp(pars(3));
    int Tsize = centeredreturns.n_rows;
    vec qlower(Tsize, fill::zeros), qupper(Tsize, fill::zeros);
    vec errorlower(Tsize, fill::zeros), errorupper(Tsize, fill::zeros);
    vec epsilon(Tsize, fill::zeros), scale(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    qlower(0) = q0(0), qupper(0) = q0(1);
    scale(0) = qupper(0) - qlower(0);
    epsilon(0) = centeredreturns(0)/scale(0);
    double answer = 1e8;
        
    if (omega1 < omega2){
        for (int t = 1; t <= (Tsize-1); ++t){
            qlower(t) = omega1 + beta * qlower(t-1) + gamma * qlower(t-1) * abs(epsilon(t-1));
            qupper(t) = omega2 + beta * qupper(t-1) + gamma * qupper(t-1) * abs(epsilon(t-1));
            scale(t) = qupper(t) - qlower(t);
            epsilon(t) = centeredreturns(t)/scale(t);
            errorlower(t) = centeredreturns(t)  - qlower(t);
            errorupper(t) = centeredreturns(t)  - qupper(t);
            objective(t) = check_Cpp( errorlower(t), tau ) + check_Cpp( errorupper(t), 1-tau );
        }
        answer =  sum( objective( span( 1,(Tsize-1) ) ) );
    }
    return answer;
}



//[[Rcpp::export]]
mat gSAV_forecast_Cpp(vec centeredreturns, vec q0, vec pars){
    double omega1 = pars(0), omega2 = pars(1), beta = exp(pars(2)), gamma = exp(pars(3));
    int Tsize = centeredreturns.n_rows;
    vec qlower(Tsize+1, fill::zeros), qupper(Tsize+1, fill::zeros);
    vec epsilon(Tsize+1, fill::zeros), scale(Tsize+1, fill::zeros);
    qlower(0) = q0(0), qupper(0) = q0(1);
    scale(0) = qupper(0) - qlower(0);
    epsilon(0) = centeredreturns(0)/scale(0);
                                
    for (int t = 1; t <= (Tsize-1); ++t){
        qlower(t) = omega1 + beta * qlower(t-1) + gamma * qlower(t-1) * abs(epsilon(t-1));
        qupper(t) = omega2 + beta * qupper(t-1) + gamma * qupper(t-1) * abs(epsilon(t-1));
        scale(t) = qupper(t) - qlower(t);
        epsilon(t) = centeredreturns(t)/scale(t);
    }
    
    int t = Tsize;
    qlower(t) = omega1 + beta * qlower(t-1) + gamma * qlower(t-1) * abs(epsilon(t-1));
    qupper(t) = omega2 + beta * qupper(t-1) + gamma * qupper(t-1) * abs(epsilon(t-1));
              
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = qlower;
    results(span(0, Tsize),1) = qupper;
    return results;
}



//[[Rcpp::export]]
double gAS_Cpp(vec centeredreturns, double tau, vec q0, vec pars){
    double omega1 = pars(0), omega2 = pars(1), beta = exp(pars(2)), gamma1 = exp(pars(3)), gamma2 = exp(pars(4));
    int Tsize = centeredreturns.n_rows;
    vec qlower(Tsize, fill::zeros), qupper(Tsize, fill::zeros);
    vec errorlower(Tsize, fill::zeros), errorupper(Tsize, fill::zeros);
    vec epsilon(Tsize, fill::zeros), scale(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    qlower(0) = q0(0), qupper(0) = q0(1);
    scale(0) = qupper(0) - qlower(0);
    epsilon(0) = centeredreturns(0)/scale(0);
    double answer = 1e8;
        
    if (omega1 < omega2){
        for (int t = 1; t <= (Tsize-1); ++t){
            int indic = 0;
            if ( centeredreturns(t-1) > 0 ) indic = 1;
            qlower(t) = omega1 + beta * qlower(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qlower(t-1) * abs(epsilon(t-1));
            qupper(t) = omega2 + beta * qupper(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qupper(t-1) * abs(epsilon(t-1));
            scale(t) = qupper(t) - qlower(t);
            epsilon(t) = centeredreturns(t)/scale(t);
            errorlower(t) = centeredreturns(t)  - qlower(t);
            errorupper(t) = centeredreturns(t)  - qupper(t);
            objective(t) = check_Cpp( errorlower(t), tau ) + check_Cpp( errorupper(t), 1-tau );
        }
        answer =  sum(   objective( span( 1,(Tsize-1) ) )   );
    }
    return answer;
}



//[[Rcpp::export]]
mat gAS_forecast_Cpp(vec centeredreturns, vec q0, vec pars){
    double omega1 = pars(0), omega2 = pars(1), beta = exp(pars(2)), gamma1 = exp(pars(3)), gamma2 = exp(pars(4));
    int Tsize = centeredreturns.n_rows;
    vec qlower(Tsize+1, fill::zeros), qupper(Tsize+1, fill::zeros);
    vec epsilon(Tsize+1, fill::zeros), scale(Tsize+1, fill::zeros);
    qlower(0) = q0(0), qupper(0) = q0(1);
    scale(0) = qupper(0) - qlower(0);
    epsilon(0) = centeredreturns(0)/scale(0);

    for (int t = 1; t <= (Tsize-1); ++t){
        int indic = 0;
        if ( centeredreturns(t-1) > 0 ) indic = 1;
        qlower(t) = omega1 + beta * qlower(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qlower(t-1) * abs(epsilon(t-1));
        qupper(t) = omega2 + beta * qupper(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qupper(t-1) * abs(epsilon(t-1));
        scale(t) = qupper(t) - qlower(t);
        epsilon(t) = centeredreturns(t)/scale(t);
    }
    
    int t = Tsize;
    int indic = 0;
    if ( centeredreturns(t-1) > 0 ) indic = 1;
     
    qlower(t) = omega1 + beta * qlower(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qlower(t-1) * abs(epsilon(t-1));
    qupper(t) = omega2 + beta * qupper(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qupper(t-1) * abs(epsilon(t-1));
              
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = qlower;
    results(span(0, Tsize),1) = qupper;
    return results;
}



//[[Rcpp::export]]
double SAV_Cpp(vec centeredreturns, double tau, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2);
    int Tsize = centeredreturns.n_rows;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    q(0) = q0;
                    
    for (int t = 1; t <= (Tsize-1); ++t){
        q(t) = beta0 + + beta1 * q(t-1) + beta2 * abs( centeredreturns(t-1) );
        e(t) = centeredreturns(t)  - q(t);
        objective(t) = check_Cpp( e(t), tau );
    }
        
    double answer =  sum(   objective( span( 1,(Tsize-1) ) )   );
    return answer;
}



//[[Rcpp::export]]
vec SAV_forecast_Cpp(vec centeredreturns, double tau, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2);
    int Tsize = centeredreturns.n_rows;
    vec q(Tsize+1, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    q(0) = q0;
                    
    for (int t = 1; t <= Tsize; ++t){
        q(t) = beta0 + beta1 * q(t-1) + beta2 * abs( centeredreturns(t-1) );
    }
    return q;
}



//[[Rcpp::export]]
double AS_Cpp(vec centeredreturns, double tau, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2), beta3= pars(3);
    int Tsize = centeredreturns.n_rows;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    q(0) = q0;
                    
    for (int t = 1; t <= (Tsize-1); ++t){
        int indic = 0;
        
        if ( centeredreturns(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 +  beta1 * q(t-1) + beta2 *indic*abs( centeredreturns(t-1) ) + beta3*(1-indic)*abs( centeredreturns(t-1) );
        e(t) = centeredreturns(t)  - q(t);
        objective(t) = check_Cpp( e(t), tau );
    }
        
    double answer =  sum( objective( span( 1,(Tsize-1) ) ) );
    return answer;
}



//[[Rcpp::export]]
vec AS_forecast_Cpp(vec centeredreturns, double tau, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2), beta3= pars(3);
    int Tsize = centeredreturns.n_rows;
    vec q(Tsize+1, fill::zeros);
    q(0) = q0;
                    
    for (int t = 1; t <= Tsize; ++t){
        int indic = 0;
        
        if ( centeredreturns(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 +  beta1 * q(t-1) + beta2 *indic*abs( centeredreturns(t-1) ) + beta3*(1-indic)*abs( centeredreturns(t-1) );
    }
    return q;
}



//[[Rcpp::export]]
double mult_Cpp(vec centeredreturns, double alpha, vec VaR, double pars){
    double gamma0 = pars;
    int Tsize = centeredreturns.n_rows;
    vec e(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    ES(0)= (1+exp(gamma0))*VaR(0);
    
    double penalty = -1e7;
    
    for (int t = 1; t <= (Tsize-1); ++t){
        ES(t)= (1+exp(gamma0))*VaR(t);
        e(t) = centeredreturns(t) - VaR(t);
        objective(t) = penalty;
        
        if ( ES(t) < 0 ){
            objective(t) =  log( (1-alpha)/( -ES(t) )  )  -  (  check_Cpp( e(t), alpha ) / ( alpha * ( -ES(t) ) )   );
        }
    }
        
    double answer =  sum(   objective( span( 1,(Tsize-1) ) )   );
    return answer;
}



//[[Rcpp::export]]
double AL_mult_SAV_Cpp(vec centeredreturns, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2);
    double gamma0 = pars(3);

    int Tsize = centeredreturns.n_rows;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    
    q(0) = q0;
    ES(0)= (1+exp(gamma0))*q(0);
    
    double penalty = -1e7;
        
    for (int t = 1; t <= (Tsize-1); ++t){
        q(t) = beta0 + beta1*q(t-1) + beta2 * abs( centeredreturns(t-1) );
        ES(t)= (1+exp(gamma0))*q(t);
        e(t) = centeredreturns(t) - q(t);
        
        objective(t) = penalty;
        
        if ( ES(t) < 0 ){
            objective(t) =  log( (1-alpha)/( -ES(t) )  )  -  ( check_Cpp( e(t), alpha ) / ( alpha * ( -ES(t) ) )  );
        }
    }
        
    double answer =  sum(   objective( span( 1,(Tsize-1) ) )   );
    return answer;
}



//[[Rcpp::export]]
mat AL_mult_SAV_forecast_Cpp(vec centeredreturns, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2);
    double gamma0 = pars(3);

    int Tsize = centeredreturns.n_rows;
    vec q(Tsize+1, fill::zeros);
    vec ES(Tsize+1, fill::zeros);
    
    q(0) = q0;
    ES(0)= (1+exp(gamma0))*q(0);
        
    for (int t = 1; t <= Tsize; ++t){
        q(t) = beta0 + beta1 * q(t-1) + beta2 * abs( centeredreturns(t-1) );
        ES(t)= (1+exp(gamma0))*q(t);
    }
        
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = q;
    results(span(0, Tsize),1) = ES;
    return results;
}



//[[Rcpp::export]]
double AR_Cpp(vec centeredreturns, double alpha, vec VaR, vec pars){
    double gamma0 = exp(pars(0));
    double a = exp(pars(1));
    double b = exp(pars(2));
    double den = 1.0 + a + b;
    double gamma1 = a/den, gamma2 = b/den;
    
    int Tsize = centeredreturns.n_rows;
    vec e(Tsize, fill::zeros);
    vec x(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);

    x(0) = gamma0;
    ES(0)= VaR(0) - x(0);
    
    double penalty = -1e7;
    
    for (int t = 1; t <= (Tsize-1); ++t){
        x(t) = x(t-1);
        
        if ( VaR(t-1) - centeredreturns(t-1) >= 0 ){
            x(t) = gamma0 + gamma1*( VaR(t-1) - centeredreturns(t-1) ) + gamma2*x(t-1);
        }
                    
        ES(t)= VaR(t) - x(t);
        e(t) = centeredreturns(t) - VaR(t);
        
        objective(t) = penalty;
        
        if ( ES(t) < 0 ){
            objective(t) =  log( (1-alpha)/( -ES(t) )  )  -  (  check_Cpp( e(t), alpha ) / ( alpha * ( -ES(t) ) )   );
        }
    }
        
    double answer =  sum(  objective( span( 1,(Tsize-1) ) )  );
    return answer;
}



//[[Rcpp::export]]
double AL_ar_SAV_Cpp(vec centeredreturns, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2);
    double gamma0 = exp(pars(3));
    double a = exp(pars(4));
    double b = exp(pars(5));
    double den = 1.0 + a + b;
    double gamma1 = a/den, gamma2 = b/den;
    
    int Tsize = centeredreturns.n_rows;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec x(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    
    q(0) = q0;
    x(0) = gamma0;
    ES(0)= q(0) - x(0);
        
    double penalty = -1e7;
        
    for (int t = 1; t <= (Tsize-1); ++t){
        x(t) = x(t-1);
        
        if ( q(t-1) - centeredreturns(t-1) >= 0 ){
            x(t) = gamma0 + gamma1*( q(t-1) - centeredreturns(t-1) ) + gamma2*x(t-1);
        }

        q(t) = beta0 + beta1 * q(t-1) + beta2 * abs( centeredreturns(t-1) );
        ES(t)= q(t) - x(t);
        e(t) = centeredreturns(t) - q(t);
        
        objective(t) = penalty;
        
        if ( ES(t) < 0 ){
            objective(t) =  log( (1-alpha)/( -ES(t) )  )  -  (  check_Cpp( e(t), alpha ) / ( alpha * ( -ES(t) ) )   );
        }
    }
        
    double answer =  sum(  objective( span( 1,(Tsize-1) ) )  );
    return answer;
}



//[[Rcpp::export]]
mat AL_ar_SAV_forecast_Cpp(vec centeredreturns, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2);
    double gamma0 = exp(pars(3));
    double a = exp(pars(4));
    double b = exp(pars(5));
    double den = 1.0 + a + b;
    double gamma1 = a/den, gamma2 = b/den;
    
    int Tsize = centeredreturns.n_rows;
    vec q(Tsize+1, fill::zeros);
    vec x(Tsize+1, fill::zeros);
    vec ES(Tsize+1, fill::zeros);
    
    q(0) = q0;
    x(0) = gamma0;
    ES(0)= q(0) - x(0);
                
    for (int t = 1; t <= Tsize; ++t){
        x(t) = x(t-1);
        
        if ( q(t-1) - centeredreturns(t-1) >= 0 ){
            x(t) = gamma0 + gamma1*( q(t-1) - centeredreturns(t-1) ) + gamma2*x(t-1);
        }

        q(t) = beta0 + beta1 * q(t-1) + beta2 * abs( centeredreturns(t-1) );
        ES(t)= q(t) - x(t);
    }
            
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = q;
    results(span(0, Tsize),1) = ES;
    return results;
}



//[[Rcpp::export]]
double AL_mult_AS_Cpp(vec centeredreturns, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2), beta3=pars(3);
    double gamma0 = pars(4);
    int Tsize = centeredreturns.n_rows;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    
    q(0) = q0;
    ES(0)= (1+exp(gamma0))*q(0);
    
    double penalty = -1e7;
        
    for (int t = 1; t <= (Tsize-1); ++t){
        int indic = 0;
        
        if ( centeredreturns(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 + beta1 * q(t-1) + beta2 * indic*abs( centeredreturns(t-1) ) + beta3 *(1-indic)*abs( centeredreturns(t-1) ) ;
        ES(t)= (1+exp(gamma0))*q(t);
        e(t) = centeredreturns(t) - q(t);
        
        objective(t) = penalty;
        
        if ( ES(t) < 0 ){
            objective(t) =  log( (1-alpha)/( -ES(t) )  )  -  (  check_Cpp( e(t), alpha ) / ( alpha * ( -ES(t) ) )   );
        }
    }
        
    double answer =  sum(  objective( span( 1,(Tsize-1) ) )  );
    return answer;
}



//[[Rcpp::export]]
mat AL_mult_AS_forecast_Cpp(vec centeredreturns, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2), beta3=pars(3);
    double gamma0 = pars(4);
    int Tsize = centeredreturns.n_rows;
    vec q(Tsize+1, fill::zeros);
    vec ES(Tsize+1, fill::zeros);
    
    q(0) = q0;
    ES(0)= (1+exp(gamma0))*q(0);
    
    double penalty = -1e7;
        
    for (int t = 1; t <= Tsize; ++t){
        int indic = 0;
        
        if ( centeredreturns(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 + beta1 * q(t-1) + beta2 * indic*abs( centeredreturns(t-1) ) + beta3 *(1-indic)*abs( centeredreturns(t-1) ) ;
        ES(t)= (1+exp(gamma0))*q(t);
    }
            
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = q;
    results(span(0, Tsize),1) = ES;
    return results;
}



//[[Rcpp::export]]
double AL_ar_AS_Cpp(vec centeredreturns, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2), beta3=pars(3);
    double gamma0 = exp(pars(4));
    double a = exp(pars(5));
    double b = exp(pars(6));
    double den = 1.0 + a + b;
    double gamma1 = a/den, gamma2 = b/den;
    
    int Tsize = centeredreturns.n_rows;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec x(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    
    q(0) = q0;
    x(0) = gamma0;
    ES(0)= q(0) - x(0);
        
    double penalty = -1e7;
        
    for (int t = 1; t <= (Tsize-1); ++t){
        int indic = 0;
        
        if ( centeredreturns(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 + beta1 * q(t-1) + beta2 * indic*abs( centeredreturns(t-1) ) + beta3 * (1-indic)*abs( centeredreturns(t-1) ) ;
        
        x(t) = x(t-1);
        
        if ( q(t-1) - centeredreturns(t-1) >= 0 ){
            x(t) = gamma0 + gamma1*( q(t-1) - centeredreturns(t-1) ) + gamma2*x(t-1);
        }
        x(t)  = std::max(x(t), 1e-8);
        
        ES(t)= q(t) - x(t);
        e(t) = centeredreturns(t) - q(t);
        objective(t) = penalty;
        
        if ( ES(t) < 0 ){
            objective(t) =  log( (1-alpha)/( -ES(t) )  )  -  (  check_Cpp( e(t), alpha ) / ( alpha * ( -ES(t) ) )   );
        }
    }
        
    double answer =  sum(   objective( span( 1,(Tsize-1) ) )   );
    return answer;
}



//[[Rcpp::export]]
mat AL_ar_AS_forecast_Cpp(vec centeredreturns, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = exp(pars(1)), beta2 = pars(2), beta3=pars(3);
    double gamma0 = exp(pars(4));
    double a = exp(pars(5));
    double b = exp(pars(6));
    double den = 1.0 + a + b;
    double gamma1 = a/den, gamma2 = b/den;
    
    int Tsize = centeredreturns.n_rows;
    vec q(Tsize+1, fill::zeros);
    vec x(Tsize+1, fill::zeros);
    vec ES(Tsize+1, fill::zeros);
    
    q(0) = q0;
    x(0) = gamma0;
    ES(0)= q(0) - x(0);
                
    for (int t = 1; t <= Tsize; ++t){
        
        int indic = 0;
        
        if ( centeredreturns(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 + beta1 * q(t-1) + beta2 * indic*abs( centeredreturns(t-1) ) + beta3 * (1-indic)*abs( centeredreturns(t-1) ) ;
        
        x(t) = x(t-1);
        
        if ( q(t-1) - centeredreturns(t-1) >= 0 ){
            x(t) = gamma0 + gamma1*( q(t-1) - centeredreturns(t-1) ) + gamma2*x(t-1);
        }
        
        ES(t)= q(t) - x(t);
    }

    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = q;
    results(span(0, Tsize),1) = ES;
    return results;
}



double clamp_Cpp(double x, double lo, double hi) {
    return std::max(lo, std::min(x, hi));
}



//[[Rcpp::export]]
double GAS_Cpp(vec centeredreturns, vec pars, double alpha){
    double a = clamp_Cpp(pars(0), -40.0, 40.0);
    double b = clamp_Cpp(pars(1), -40.0, 40.0);
    double c = clamp_Cpp(pars(2), -40.0, 40.0);
    double d = clamp_Cpp(pars(3), -40.0, 40.0);
    a = -exp(a), b = a - exp(b);
    c = exp(c), d = exp(d);
    double den = 1.0 + c + d;
    double beta = c/den, gamma = d/den;

    int Tsize = centeredreturns.n_rows;
    vec k(Tsize, fill::zeros);
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec hit(Tsize, fill::zeros);
    vec s(Tsize, fill::zeros);
    vec L(Tsize, fill::zeros);
    
    k(0) = 0.0;
    q(0) = a * exp(k(0));
    e(0) = b * exp(k(0));

    for (int t = 1; t <= (Tsize-1); ++t){
            
        k(t) = beta*k(t-1) + gamma*s(t-1);
        q(t) = a * exp(k(t));
        e(t) = b * exp(k(t));
                    
        hit(t) = (centeredreturns(t) <= q(t)) ? 1.0 : 0.0;
                    
        s(t) = ( hit(t) * centeredreturns(t) / alpha - e(t) )/e(t);
        L(t) = - hit(t)*(q(t) - centeredreturns(t))/(alpha*e(t))  + q(t)/e(t) + log(-e(t)) - 1;
    }
    
    double answer =  sum(  L( span( 1,(Tsize-1) ) )  );
    
    if (sum(hit) == 0 ){
        answer = 1e7;
    }
    return answer;
}



//[[Rcpp::export]]
mat GAS_forecast_Cpp(vec centeredreturns, vec pars, double alpha){
    double a = clamp_Cpp(pars(0), -40.0, 40.0);
    double b = clamp_Cpp(pars(1), -40.0, 40.0);
    double c = clamp_Cpp(pars(2), -40.0, 40.0);
    double d = clamp_Cpp(pars(3), -40.0, 40.0);
    a = -exp(a), b = a - exp(b);
    c = exp(c), d = exp(d);
    double den = 1.0 + c + d;
    double beta = c/den, gamma = d/den;
    
    int Tsize = centeredreturns.n_rows;
    vec k(Tsize+1, fill::zeros);
    vec q(Tsize+1, fill::zeros);
    vec e(Tsize+1, fill::zeros);
    vec hit(Tsize+1, fill::zeros);
    vec s(Tsize+1, fill::zeros);

    k(0) = 0.0;
    q(0) = a * exp(k(0));
    e(0) = b * exp(k(0));
    
    for (int t = 1; t <= (Tsize-1); ++t){
            
        k(t) = beta*k(t-1) + gamma*s(t-1);
        q(t) = a * exp(k(t));
        e(t) = b * exp(k(t));
            
        hit(t) = (centeredreturns(t) <= q(t)) ? 1.0 : 0.0;
        
        s(t) = ( hit(t) * centeredreturns(t) / alpha - e(t) )/e(t);
    }
        
    k(Tsize) = beta*k(Tsize-1) + gamma*s(Tsize-1);
    q(Tsize) = a * exp(k(Tsize));
    e(Tsize) = b * exp(k(Tsize));
 
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = q;
    results(span(0, Tsize),1) = e;
    return results;
}



// [[Rcpp::export]]
double GARCH_norm_Cpp(const vec& pars, const vec& centeredreturns, const double var0){
    int Tsize = centeredreturns.n_rows;
    double c = clamp_Cpp(pars(0), -40.0, 40.0);
    double a = clamp_Cpp(pars(1), -40.0, 40.0);
    double b = clamp_Cpp(pars(2), -40.0, 40.0);
    double beta0 = exp(c);
    double s1 = exp(a), s2 = exp(b);
    double den = 1.0 + s1 + s2;
    double beta1 = s1/den, beta2 = s2/den;
  
    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    sigma2(0) = var0;
  
    for (int t = 1; t < Tsize; ++t){
        
        double e2_lag = centeredreturns(t - 1) * centeredreturns(t - 1);
        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * e2_lag;

        double z2 = (centeredreturns(t)*centeredreturns(t)) / sigma2(t);

        loglike(t) = -0.5 * log(2.0 * datum::pi) -0.5 * log(sigma2(t)) -0.5 * z2;
    }

    return accu(loglike);
}



// [[Rcpp::export]]
vec GARCH_norm_forecast_Cpp(const vec& centeredreturns, const vec& pars, const double var0){
    int Tsize = centeredreturns.n_rows;
    double c = clamp_Cpp(pars(0), -40.0, 40.0);
    double a = clamp_Cpp(pars(1), -40.0, 40.0);
    double b = clamp_Cpp(pars(2), -40.0, 40.0);
    double beta0 = exp(c);
    double s1 = exp(a), s2 = exp(b);
    double den = 1.0 + s1 + s2;
    double beta1 = s1/den, beta2 = s2/den;
    
    vec sigma2(Tsize + 1, fill::zeros);
    sigma2(0) = var0;

    for (int t = 1; t < Tsize; ++t){
        double e2_lag = centeredreturns(t - 1) * centeredreturns(t - 1);
        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * e2_lag;
    }

    sigma2(Tsize) = beta0 + beta1 * sigma2(Tsize - 1) + beta2 * (centeredreturns(Tsize - 1) * centeredreturns(Tsize - 1));
    return sigma2;
}



// [[Rcpp::export]]
double GARCH_t_Cpp(const vec& pars, const vec& centeredreturns, const double var0) {
    const int Tsize = centeredreturns.n_rows;
    double c = clamp_Cpp(pars(0), -40.0, 40.0);
    double a = clamp_Cpp(pars(1), -40.0, 40.0);
    double b = clamp_Cpp(pars(2), -40.0, 40.0);
    double g = clamp_Cpp(pars(3), -40.0, 40.0);
    double beta0 = exp(c);
    double s1 = exp(a), s2 = exp(b);
    double den = 1.0 + s1 + s2;
    double beta1 = s1/den, beta2 = s2/den;
    double v = 2.0 + exp(g);

    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    sigma2(0) = var0;
    
    const double inv_v_m2 = 1.0/(v - 2.0);
        
    const double cst = lgamma((v + 1.0)*0.5) - lgamma(v*0.5) - 0.5*log(datum::pi * (v - 2.0));
    
    for (int t = 1; t < Tsize; ++t){
        double e2_lag = centeredreturns(t - 1) * centeredreturns(t - 1);

        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * e2_lag;
    
        double e2 = (centeredreturns(t) * centeredreturns(t)) / sigma2(t);
      
        loglike(t) = cst - 0.5 * log(sigma2(t)) - ((v + 1.0) * 0.5) * log(1.0 + e2 * inv_v_m2);
    }
    
    return accu(loglike);
}



// [[Rcpp::export]]
vec GARCH_t_forecast_Cpp(const vec& centeredreturns, const vec& pars, const double var0){
    const int Tsize = centeredreturns.n_rows;
    double c = clamp_Cpp(pars(0), -40.0, 40.0);
    double a = clamp_Cpp(pars(1), -40.0, 40.0);
    double b = clamp_Cpp(pars(2), -40.0, 40.0);
    double g = clamp_Cpp(pars(3), -40.0, 40.0);
    double beta0 = exp(c);
    double s1 = exp(a), s2 = exp(b);
    double den = 1.0 + s1 + s2;
    double beta1 = s1/den, beta2 = s2/den;
    double v = 2.0 + exp(g);

    vec sigma2(Tsize + 1, fill::zeros);
    sigma2(0) = var0;
            
    for (int t = 1; t < Tsize; ++t){
        double e2_lag = centeredreturns(t - 1) * centeredreturns(t - 1);

        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * e2_lag;
              
    }
    sigma2(Tsize) = beta0 + beta1 * sigma2(Tsize - 1) + beta2 * (centeredreturns(Tsize - 1) * centeredreturns(Tsize - 1));

    return sigma2;
}



//[[Rcpp::export]]
double Hansen_dt_Cpp(double z, double v, double lambda){
    double numc = tgamma((v+1)/2);
    double denc = tgamma(v/2)*sqrt(datum::pi * (v-2) );
    double c = numc/denc;
    double a = 4*lambda*c*( (v-2)/(v-1) );
    double b = sqrt(1 + 3*pow(lambda,2) - pow(a,2));
    double num = b*z + a, den =0;
            
    den = (z < -a / b) ? (1.0 - lambda) : (1.0 + lambda);
    
    double ex = -(v+1)/2;
    
    double answer = b*c*pow( ( 1 + ( 1/(v-2) )*pow(num/den,2 )  ), ex  );
    return answer;
}



// [[Rcpp::export]]
double GARCH_skew_t_Cpp(const vec& centeredreturns, const vec& pars, const double var0){
    int Tsize = centeredreturns.n_rows;
    double c = clamp_Cpp(pars(0), -40.0, 40.0);
    double a = clamp_Cpp(pars(1), -40.0, 40.0);
    double b = clamp_Cpp(pars(2), -40.0, 40.0);
    double g = clamp_Cpp(pars(3), -40.0, 40.0);
    double ell = clamp_Cpp(pars(4), -40.0, 40.0);
    double beta0 = exp(c);
    double s1 = exp(a), s2 = exp(b);
    double den = 1.0 + s1 + s2;
    double beta1 = s1/den, beta2  = s2/den;

    double v      = 2.0 + exp(g);
    double lambda = 2.0 * ( exp(ell)/(1.0 + exp(ell)) ) - 1.0;

    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);

    sigma2(0) = var0;

    for (int t = 1; t < Tsize; ++t){
        const double e_lag = centeredreturns(t - 1);
        sigma2(t) = beta0 +  beta1*sigma2(t - 1) + beta2*(e_lag * e_lag);

        double z = centeredreturns(t) / std::sqrt(sigma2(t));
        double dens = Hansen_dt_Cpp(z, v, lambda);

        loglike(t) = log(dens) - 0.5 * log(sigma2(t));
    }

    return accu(loglike);
}



// [[Rcpp::export]]
vec GARCH_skew_t_forecast_Cpp(const vec& centeredreturns, const vec& pars, const double var0) {
    int Tsize = centeredreturns.n_rows;
    vec sigma2(Tsize + 1, fill::zeros);
    double c = clamp_Cpp(pars(0), -40.0, 40.0);
    double a = clamp_Cpp(pars(1), -40.0, 40.0);
    double b = clamp_Cpp(pars(2), -40.0, 40.0);
    double g = clamp_Cpp(pars(3), -40.0, 40.0);
    double ell = clamp_Cpp(pars(4), -40.0, 40.0);
    double beta0 = exp(c);
    double s1 = exp(a), s2 = exp(b);
    double den = 1.0 + s1 + s2;
    double beta1 = s1/den, beta2  = s2/den;

    double v      = 2.0 + exp(g);
    double lambda = 2.0 * ( exp(ell)/(1.0 + exp(ell)) ) - 1.0;

    sigma2(0) = var0;

    for (int t = 1; t < Tsize; ++t){
        const double e_lag = centeredreturns(t - 1);
        sigma2(t) = beta0 + beta1*sigma2(t - 1) + beta2*(e_lag * e_lag);
    }

    double e_last = centeredreturns(Tsize - 1);
    sigma2(Tsize) = beta0 + beta1*sigma2(Tsize - 1) + beta2*(e_last * e_last);
  
    return sigma2;
}



// [[Rcpp::export]]
double GJR_norm_Cpp(const vec& centeredreturns, const vec& pars, const double var0){
    int Tsize = centeredreturns.n_rows;
    double c = clamp_Cpp(pars(0), -40.0, 40.0);
    double a = clamp_Cpp(pars(1), -40.0, 40.0);
    double b = clamp_Cpp(pars(2), -40.0, 40.0);
    double d = clamp_Cpp(pars(3), -40.0, 40.0);
    double beta0 = exp(c);
    double s1 = exp(a), s2 = exp(b), s3 = exp(d);
    double den = 1.0 + s1 + s2 + 0.5*s3;
    double beta1 = s1/den, beta2 = s2/den, beta3 = s3/den;

    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
  
    sigma2(0) = var0;

    double cst = -0.5 * std::log(2.0 * datum::pi);
    double id = 0;
  
    for (int t = 1; t < Tsize; ++t) {
        const double e_lag = centeredreturns(t - 1);
    
        id = (e_lag < 0.0) ? 1.0 : 0.0;

        const double e2_lag = e_lag * e_lag;

        const double s2_t = beta0 + beta1 * sigma2(t - 1) + (beta2 + beta3 * id) * e2_lag;
  
        sigma2(t) = s2_t;

        const double z2 = (centeredreturns(t) * centeredreturns(t)) / s2_t;
        loglike(t) = cst - 0.5 * std::log(s2_t) - 0.5 * z2;
    }
  
    return accu(loglike);
}



// [[Rcpp::export]]
vec GJR_norm_forecast_Cpp(const vec& centeredreturns, const vec& pars, const double var0){
    int Tsize = centeredreturns.n_rows;
    double c = clamp_Cpp(pars(0), -40.0, 40.0);
    double a = clamp_Cpp(pars(1), -40.0, 40.0);
    double b = clamp_Cpp(pars(2), -40.0, 40.0);
    double d = clamp_Cpp(pars(3), -40.0, 40.0);
    double beta0 = exp(c);
    double s1 = exp(a), s2 = exp(b), s3 = exp(d);
    double den = 1.0 + s1 + s2 + 0.5*s3;
    double beta1 = s1/den, beta2 = s2/den, beta3 = s3/den;

    vec sigma2(Tsize + 1, fill::zeros);

    sigma2(0) = var0;

    for (int t = 1; t < Tsize; ++t){
        const double e_lag  = centeredreturns(t - 1);
        const double id     = (e_lag < 0.0) ? 1.0 : 0.0;
        const double e2_lag = e_lag * e_lag;

        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + (beta2 + beta3 * id) * e2_lag;
    }

    const double e_last  = centeredreturns(Tsize - 1);
    const double id_last = (e_last < 0.0) ? 1.0 : 0.0;
    sigma2(Tsize) = beta0 + beta1 * sigma2(Tsize - 1) + (beta2 + beta3 * id_last) * (e_last * e_last);

    return sigma2;
}



// [[Rcpp::export]]
double GJR_t_Cpp(const vec& centeredreturns, const vec& pars, const double var0) {
    int Tsize = centeredreturns.n_rows;
    double c = clamp_Cpp(pars(0), -40.0, 40.0);
    double a = clamp_Cpp(pars(1), -40.0, 40.0);
    double b = clamp_Cpp(pars(2), -40.0, 40.0);
    double d = clamp_Cpp(pars(3), -40.0, 40.0);
    double g = clamp_Cpp(pars(4), -40.0, 40.0);
    double beta0 = exp(c);
    double s1 = exp(a), s2 = exp(b), s3 = exp(d);
    double den = 1.0 + s1 + s2 + 0.5 * s3;
    double beta1 = s1/den, beta2 = s2/den, beta3 = s3/den;
    double v    = 2.0 + exp(g);

    double cst = lgamma((v + 1.0) * 0.5) - lgamma(v * 0.5)  - 0.5 * log(datum::pi * (v - 2.0));

    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);

    sigma2(0) = var0;

    for (int t = 1; t < Tsize; ++t) {
        double e_lag  = centeredreturns(t - 1);
        double id     = (e_lag < 0.0) ? 1.0 : 0.0;
        double e2_lag = e_lag * e_lag;

        double s2_t = beta0 + beta1 * sigma2(t - 1) + (beta2 + beta3 * id) * e2_lag;

        sigma2(t) = s2_t;

        double e2 = (centeredreturns(t) * centeredreturns(t)) / s2_t;
        loglike(t) = cst - 0.5 * log(s2_t) - ((v + 1.0) * 0.5) * log(1.0 + e2 / (v - 2.0));
    }

    return accu(loglike);
}



// [[Rcpp::export]]
vec GJR_t_forecast_Cpp(const vec& centeredreturns, const vec& pars, const double var0){
    int Tsize = centeredreturns.n_rows;
    double c = clamp_Cpp(pars(0), -40.0, 40.0);
    double a = clamp_Cpp(pars(1), -40.0, 40.0);
    double b = clamp_Cpp(pars(2), -40.0, 40.0);
    double d = clamp_Cpp(pars(3), -40.0, 40.0);
    double g = clamp_Cpp(pars(4), -40.0, 40.0);
    double beta0 = exp(c);
    double s1 = exp(a), s2 = exp(b), s3 = exp(d);
    double den = 1.0 + s1 + s2 + 0.5 * s3;
    double beta1 = s1/den, beta2 = s2/den, beta3 = s3/den;
    double v    = 2.0 + exp(g);

    double cst = lgamma((v + 1.0) * 0.5) - lgamma(v * 0.5)  - 0.5 * log(datum::pi * (v - 2.0));

    vec sigma2(Tsize+1, fill::zeros);
    sigma2(0) = var0;

    for (int t = 1; t < Tsize; ++t) {
        double e_lag  = centeredreturns(t - 1);
        double id     = (e_lag < 0.0) ? 1.0 : 0.0;
        double e2_lag = e_lag * e_lag;

        double s2_t = beta0 + beta1 * sigma2(t - 1) + (beta2 + beta3 * id) * e2_lag;

        sigma2(t) = s2_t;

    }

    double e_last  = centeredreturns(Tsize - 1);
    double id_last = (e_last < 0.0) ? 1.0 : 0.0;
    sigma2(Tsize) = beta0 + beta1 * sigma2(Tsize - 1) + (beta2 + beta3 * id_last) * (e_last * e_last);

    return sigma2;
}



// [[Rcpp::export]]
double GJR_skew_t_Cpp(const vec& centeredreturns, const vec& pars, const double var0){
    int Tsize = centeredreturns.n_rows;
    double c   = clamp_Cpp(pars(0), -40.0, 40.0);
    double a   = clamp_Cpp(pars(1), -40.0, 40.0);
    double b   = clamp_Cpp(pars(2), -40.0, 40.0);
    double d   = clamp_Cpp(pars(3), -40.0, 40.0);
    double g   = clamp_Cpp(pars(4), -40.0, 40.0);
    double ell = clamp_Cpp(pars(5), -40.0, 40.0);
    double s1 = exp(a), s2 = exp(b), s3 = exp(d);
    double den = 1.0 + s1 + s2 + 0.5*s3;
    double beta0 = exp(c);
    double beta1 = s1/den, beta2 = s2/den, beta3 = s3/den;
    double v     = 2.0 + exp(g);
    double lambda = 2.0 * ( 1.0/(1.0 + exp(-ell)) ) - 1.0;

    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    sigma2(0) = var0;

    for (int t = 1; t < Tsize; ++t) {
        const double e_lag  = centeredreturns(t - 1);
        const double id     = (e_lag < 0.0) ? 1.0 : 0.0;
        const double e2_lag = e_lag * e_lag;

        const double s2_t = beta0 + beta1 * sigma2(t - 1) + (beta2 + beta3 * id) * e2_lag;

        sigma2(t) = s2_t;

        const double z   = centeredreturns(t) / sqrt(s2_t);
        const double dens = Hansen_dt_Cpp(z, v, lambda);

        loglike(t) = log(dens) - 0.5 * log(s2_t);
    }

    return accu(loglike);
}



// [[Rcpp::export]]
vec GJR_skew_t_forecast_Cpp(const vec& centeredreturns, const vec& pars, const double var0){
    int Tsize = centeredreturns.n_rows;
    double c   = clamp_Cpp(pars(0), -40.0, 40.0);
    double a   = clamp_Cpp(pars(1), -40.0, 40.0);
    double b   = clamp_Cpp(pars(2), -40.0, 40.0);
    double d   = clamp_Cpp(pars(3), -40.0, 40.0);
    double g   = clamp_Cpp(pars(4), -40.0, 40.0);
    double ell = clamp_Cpp(pars(5), -40.0, 40.0);
    double s1 = exp(a), s2 = exp(b), s3 = exp(d);
    double den = 1.0 + s1 + s2 + 0.5*s3;
    double beta0 = exp(c);
    double beta1 = s1/den, beta2 = s2/den, beta3 = s3/den;
    double v     = 2.0 + exp(g);
    double lambda = 2.0 * ( 1.0/(1.0 + exp(-ell)) ) - 1.0;

    vec sigma2(Tsize + 1, fill::zeros);
    sigma2(0) = var0;

    for (int t = 1; t < Tsize; ++t) {
        const double e_lag  = centeredreturns(t - 1);
        const double id     = (e_lag < 0.0) ? 1.0 : 0.0;
        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + (beta2 + beta3 * id) * (e_lag * e_lag);
    }

    const double e_last  = centeredreturns(Tsize - 1);
    const double id_last = (e_last < 0.0) ? 1.0 : 0.0;
    sigma2(Tsize) = beta0 + beta1 * sigma2(Tsize - 1) + (beta2 + beta3 * id_last) * (e_last * e_last);

    return sigma2;
}



// [[Rcpp::export]]
double EGARCH_Cpp(const vec& centeredreturns, const vec& pars, const double var0){
    int Tsize = centeredreturns.n_rows;
    double w = pars(0);
    double b   = clamp_Cpp(pars(1), -40.0, 40.0);
    double a = pars(2);
    double g = pars(3);
    double omega = w;
    double beta = exp(b)/(1+exp(b));
    double alpha = a;
    double gamma = g;

    vec sigma2(Tsize, fill::zeros), z(Tsize, fill::zeros), loglike(Tsize, fill::zeros);

    double EabsZ = sqrt(2.0 / datum::pi);
    double cst   = -0.5 * log(2.0 * datum::pi);
  
    sigma2(0) = var0;

    z(0) = centeredreturns(0) / sqrt(sigma2(0));

    for (int t = 1; t < Tsize; ++t){
        const double ln_s2_t = omega + beta  * log(sigma2(t - 1)) + alpha * (abs(z(t - 1)) - EabsZ) + gamma * z(t - 1);
        const double s2_t = exp(ln_s2_t);
        sigma2(t) = s2_t;

        const double z_t = centeredreturns(t) / sqrt(s2_t);
        z(t) = z_t;

        loglike(t) = cst - 0.5 * log(s2_t) - 0.5 * (z_t * z_t);
    }

    return accu(loglike);
}



//[[Rcpp::export]]
vec EGARCH_forecast_Cpp(const vec& centeredreturns, const vec& pars, const double var0){
    int Tsize = centeredreturns.n_rows;
    double w = pars(0);
    double b   = clamp_Cpp(pars(1), -40.0, 40.0);
    double a = pars(2);
    double g = pars(3);
    double omega = w;
    double beta = exp(b)/(1+exp(b));
    double alpha = a;
    double gamma = g;

    vec sigma2(Tsize + 1, fill::zeros), z(Tsize, fill::zeros);

    double EabsZ = sqrt(2.0 / datum::pi);
    sigma2(0) = var0;

    z(0) = centeredreturns(0) / sqrt(sigma2(0));

    for (int t = 1; t < Tsize; ++t){
        double ln_s2_t = omega + beta  * log(sigma2(t - 1)) + alpha * (abs(z(t - 1)) - EabsZ) + gamma * z(t - 1);
        double s2_t = exp(ln_s2_t);
        sigma2(t) = s2_t;

        double z_t = centeredreturns(t) / sqrt(s2_t);
        z(t) = z_t;
    }
    
    double ln_s2_f = omega + beta  * log(sigma2(Tsize - 1)) + alpha * (abs(z(Tsize - 1)) - EabsZ) + gamma * z(Tsize - 1);

    sigma2(Tsize) = exp(ln_s2_f);
    
    return sigma2;
}


