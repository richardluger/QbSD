#include <RcppDist.h>
//[[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace arma;



//[[Rcpp::export]]
double check_Cpp(double u, double tau){
    int indic = 0;
    if (u < 0) indic = 1;
    return( u*( tau - indic) );
}



//[[Rcpp::export]]
double SAVscale_Cpp(vec centeredrets, double tau, vec q0, vec pars){
    double omega1 = pars(0), omega2 = pars(1), beta = exp(pars(2)), gamma = exp(pars(3));
    int Tsize = centeredrets.n_rows;
    vec qlower(Tsize, fill::zeros), qupper(Tsize, fill::zeros);
    vec errorlower(Tsize, fill::zeros), errorupper(Tsize, fill::zeros);
    vec epsilon(Tsize, fill::zeros), scale(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    qlower(0) = q0(0), qupper(0) = q0(1);
    scale(0) = qupper(0) - qlower(0);
    epsilon(0) = centeredrets(0)/scale(0);
    double answer = 1e8;
        
    if (omega1 < omega2){
        for (int t = 1; t <= (Tsize-1); t++){
            qlower(t) = omega1 + beta * qlower(t-1) + gamma * qlower(t-1) * abs(epsilon(t-1));
            qupper(t) = omega2 + beta * qupper(t-1) + gamma * qupper(t-1) * abs(epsilon(t-1));
            scale(t) = qupper(t) - qlower(t);
            epsilon(t) = centeredrets(t)/scale(t);
            errorlower(t) = centeredrets(t)  - qlower(t);
            errorupper(t) = centeredrets(t)  - qupper(t);
            objective(t) = check_Cpp( errorlower(t), tau ) + check_Cpp( errorupper(t), 1-tau );
        }
        answer =  sum( objective( span( 1,(Tsize-1) ) ) );
    }
    return(answer);
}



//[[Rcpp::export]]
mat SAVscalehatf_Cpp(vec centeredrets, vec q0, vec pars){
    double omega1 = pars(0), omega2 = pars(1), beta = pars(2), gamma = pars(3);
    int Tsize = centeredrets.n_rows;
    vec qlower(Tsize+1, fill::zeros), qupper(Tsize+1, fill::zeros);
    vec epsilon(Tsize+1, fill::zeros), scale(Tsize+1, fill::zeros);
    qlower(0) = q0(0), qupper(0) = q0(1);
    scale(0) = qupper(0) - qlower(0);
    epsilon(0) = centeredrets(0)/scale(0);
                                
    for (int t = 1; t <= (Tsize-1); t++){
        qlower(t) = omega1 + beta * qlower(t-1) + gamma * qlower(t-1) * abs(epsilon(t-1));
        qupper(t) = omega2 + beta * qupper(t-1) + gamma * qupper(t-1) * abs(epsilon(t-1));
        scale(t) = qupper(t) - qlower(t);
        epsilon(t) = centeredrets(t)/scale(t);
    }
    
    int t = Tsize;
    qlower(t) = omega1 + beta * qlower(t-1) + gamma * qlower(t-1) * abs(epsilon(t-1));
    qupper(t) = omega2 + beta * qupper(t-1) + gamma * qupper(t-1) * abs(epsilon(t-1));
              
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = qlower;
    results(span(0, Tsize),1) = qupper;
    return(results);
}



//[[Rcpp::export]]
double ASscale_Cpp(vec centeredrets, double tau, vec q0, vec pars){
    double omega1 = pars(0), omega2 = pars(1), beta = exp(pars(2)), gamma1 = exp(pars(3)), gamma2 = exp(pars(4));
    int Tsize = centeredrets.n_rows;
    vec qlower(Tsize, fill::zeros), qupper(Tsize, fill::zeros);
    vec errorlower(Tsize, fill::zeros), errorupper(Tsize, fill::zeros);
    vec epsilon(Tsize, fill::zeros), scale(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    qlower(0) = q0(0), qupper(0) = q0(1);
    scale(0) = qupper(0) - qlower(0);
    epsilon(0) = centeredrets(0)/scale(0);
    double answer = 1e8;
        
    if (omega1 < omega2){
        for (int t = 1; t <= (Tsize-1); t++){
            int indic = 0;
            if ( centeredrets(t-1) > 0 ) indic = 1;
            qlower(t) = omega1 + beta * qlower(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qlower(t-1) * abs(epsilon(t-1));
            qupper(t) = omega2 + beta * qupper(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qupper(t-1) * abs(epsilon(t-1));
            scale(t) = qupper(t) - qlower(t);
            epsilon(t) = centeredrets(t)/scale(t);
            errorlower(t) = centeredrets(t)  - qlower(t);
            errorupper(t) = centeredrets(t)  - qupper(t);
            objective(t) = check_Cpp( errorlower(t), tau ) + check_Cpp( errorupper(t), 1-tau );
        }
        answer =  sum(   objective( span( 1,(Tsize-1) ) )   );
    }
    return(answer);
}



//[[Rcpp::export]]
mat ASscalehatf_Cpp(vec centeredrets, vec q0, vec pars){
    double omega1 = pars(0), omega2 = pars(1), beta = pars(2), gamma1 = pars(3), gamma2 = pars(4);
    int Tsize = centeredrets.n_rows;
    vec qlower(Tsize+1, fill::zeros), qupper(Tsize+1, fill::zeros);
    vec epsilon(Tsize+1, fill::zeros), scale(Tsize+1, fill::zeros);
    qlower(0) = q0(0), qupper(0) = q0(1);
    scale(0) = qupper(0) - qlower(0);
    epsilon(0) = centeredrets(0)/scale(0);

    for (int t = 1; t <= (Tsize-1); t++){
        int indic = 0;
        if ( centeredrets(t-1) > 0 ) indic = 1;
        qlower(t) = omega1 + beta * qlower(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qlower(t-1) * abs(epsilon(t-1));
        qupper(t) = omega2 + beta * qupper(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qupper(t-1) * abs(epsilon(t-1));
        scale(t) = qupper(t) - qlower(t);
        epsilon(t) = centeredrets(t)/scale(t);
    }
    
    int t = Tsize;
    int indic = 0;
    if ( centeredrets(t-1) > 0 ) indic = 1;
     
    qlower(t) = omega1 + beta * qlower(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qlower(t-1) * abs(epsilon(t-1));
    qupper(t) = omega2 + beta * qupper(t-1) + (gamma1*indic + gamma2*(1-indic) ) * qupper(t-1) * abs(epsilon(t-1));
              
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = qlower;
    results(span(0, Tsize),1) = qupper;
    return(results);
}



//[[Rcpp::export]]
double SAV_Cpp(vec rets, vec mu, double tau, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    q(0) = q0;
                    
    for (int t = 1; t <= (Tsize-1); t++){
        q(t) = beta0 + beta1 * abs( centeredrets(t-1) ) + beta2*q(t-1);
        e(t) = rets(t)  - q(t);
        objective(t) = check_Cpp( e(t), tau );
    }
        
    double answer =  sum(   objective( span( (2-1),(Tsize-1) ) )   );
    return(answer);
}



//[[Rcpp::export]]
vec SAVhatf_Cpp(vec rets, vec mu, double tau, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize+1, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    q(0) = q0;
                    
    for (int t = 1; t <= Tsize; t++){
        q(t) = beta0 + beta1 * abs( centeredrets(t-1) ) + beta2*q(t-1);
    }
    return(q);
}



//[[Rcpp::export]]
double Basic_Cpp(vec rets, vec mu, double alpha, vec VaR, double pars){
    double gamma0 = pars;
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec e(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    ES(0)= (1+exp(gamma0))*VaR(0);
    
    double penalty = -1e7;
    
    for (int t = (2-1); t <= (Tsize-1); t++){
        ES(t)= (1+exp(gamma0))*VaR(t);
        e(t) = rets(t) - VaR(t);
        objective(t) = penalty;
        
        if ( mu(t)-ES(t) > 0 ){
            objective(t) =  log( (1-alpha)/( mu(t)-ES(t) )  )  -  (  check_Cpp( e(t), alpha ) / ( alpha * ( mu(t)-ES(t) ) )   );
        }
    }
        
    double answer =  sum(   objective( span( (2-1),(Tsize-1) ) )   );
    return(answer);
}



//[[Rcpp::export]]
double AR_Cpp(vec rets, vec mu, double alpha, vec VaR, vec pars){
    double gamma0 = exp(pars(0)), gamma1 = exp(pars(1)), gamma2 = exp(pars(2));
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec e(Tsize, fill::zeros);
    vec x(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);

    x(0) = gamma0;
    ES(0)= VaR(0) - x(0);
    
    double penalty = -1e7;
    
    for (int t = (2-1); t <= (Tsize-1); t++){
        x(t) = x(t-1);
        
        if ( VaR(t-1) - centeredrets(t-1) >= 0 ){
            x(t) = gamma0 + gamma1*( VaR(t-1) - centeredrets(t-1) ) + gamma2*x(t-1);
        }
                    
        ES(t)= VaR(t) - x(t);
        e(t) = rets(t) - VaR(t);
        
        objective(t) = penalty;
        
        if ( mu(t)-ES(t) > 0 ){
            objective(t) =  log( (1-alpha)/( mu(t)-ES(t) )  )  -  (  check_Cpp( e(t), alpha ) / ( alpha * ( mu(t)-ES(t) ) )   );
        }
    }
        
    double answer =  sum(  objective( span( (2-1),(Tsize-1) ) )  );
    return(answer);
}



//[[Rcpp::export]]
double Taylor_SAV_Basic_Cpp(vec rets, vec mu, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2);
    double gamma0 = pars(3);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    
    q(0) = q0;
    ES(0)= (1+exp(gamma0))*q(0);
    
    double penalty = -1e7;
        
    for (int t = 1; t <= (Tsize-1); t++){
        q(t) = beta0 + beta1 * abs( centeredrets(t-1) ) + beta2*q(t-1);
        ES(t)= (1+exp(gamma0))*q(t);
        e(t) = rets(t) - q(t);
        
        objective(t) = penalty;
        
        if ( mu(t)-ES(t) > 0 ){
            objective(t) =  log( (1-alpha)/( mu(t)-ES(t) )  )  -  ( check_Cpp( e(t), alpha ) / ( alpha * ( mu(t)-ES(t) ) )  );
        }
    }
        
    double answer =  sum(   objective( span( 1,(Tsize-1) ) )   );
    return(answer);
}



//[[Rcpp::export]]
mat Taylor_SAV_BasicHatf_Cpp(vec rets, vec mu, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2);
    double gamma0 = pars(3);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize+1, fill::zeros);
    vec ES(Tsize+1, fill::zeros);
    vec objective(Tsize, fill::zeros);
    
    q(0) = q0;
    ES(0)= (1+exp(gamma0))*q(0);
        
    for (int t = 1; t <= Tsize; t++){
        q(t) = beta0 + beta1 * abs( centeredrets(t-1) ) + beta2*q(t-1);
        ES(t)= (1+exp(gamma0))*q(t);
    }
        
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = q;
    results(span(0, Tsize),1) = ES;
    return(results);

    
}



//[[Rcpp::export]]
double Taylor_SAV_AR_Cpp(vec rets, vec mu, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2);
    double gamma0 = exp(pars(3)), gamma1 = exp(pars(4)), gamma2 = exp(pars(5));
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec x(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    
    q(0) = q0;
    x(0) = gamma0;
    ES(0)= q(0) - x(0);
        
    double penalty = -1e7;
        
    for (int t = 1; t <= (Tsize-1); t++){
        x(t) = x(t-1);
        
        if ( q(t-1) - centeredrets(t-1) >= 0 ){
            x(t) = gamma0 + gamma1*( q(t-1) - centeredrets(t-1) ) + gamma2*x(t-1);
        }

        q(t) = beta0 + beta1 * abs( centeredrets(t-1) ) + beta2*q(t-1);
        ES(t)= q(t) - x(t);
        e(t) = rets(t) - q(t);
        
        objective(t) = penalty;
        
        if ( mu(t)-ES(t) > 0 ){
            objective(t) =  log( (1-alpha)/( mu(t)-ES(t) )  )  -  (  check_Cpp( e(t), alpha ) / ( alpha * ( mu(t)-ES(t) ) )   );
        }
    }
        
    double answer =  sum(  objective( span( 1,(Tsize-1) ) )  );
    return(answer);
}



//[[Rcpp::export]]
mat Taylor_SAV_ARHatf_Cpp(vec rets, vec mu, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2);
    double gamma0 = pars(3), gamma1 = pars(4), gamma2 = pars(5);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize+1, fill::zeros);
    vec x(Tsize+1, fill::zeros);
    vec ES(Tsize+1, fill::zeros);
    
    q(0) = q0;
    x(0) = gamma0;
    ES(0)= q(0) - x(0);
                
    for (int t = 1; t <= Tsize; t++){
        x(t) = x(t-1);
        
        if ( q(t-1) - centeredrets(t-1) >= 0 ){
            x(t) = gamma0 + gamma1*( q(t-1) - centeredrets(t-1) ) + gamma2*x(t-1);
        }

        q(t) = beta0 + beta1 * abs( centeredrets(t-1) ) + beta2*q(t-1);
        ES(t)= q(t) - x(t);
    }
            
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = q;
    results(span(0, Tsize),1) = ES;
    return(results);
}



//[[Rcpp::export]]
double AS_Cpp(vec rets, vec mu, double tau, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3= pars(3);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    q(0) = q0;
                    
    for (int t = 1; t <= (Tsize-1); t++){
        int indic = 0;
        
        if ( centeredrets(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 + beta1*indic*abs( centeredrets(t-1) ) + beta2*(1-indic)*abs( centeredrets(t-1) ) +  beta3*q(t-1);
        e(t) = rets(t)  - q(t);
        objective(t) = check_Cpp( e(t), tau );
    }
        
    double answer =  sum( objective( span( (2-1),(Tsize-1) ) ) );
    return(answer);
}



//[[Rcpp::export]]
vec AShatf_Cpp(vec rets, vec mu, double tau, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3= pars(3);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize+1, fill::zeros);
    q(0) = q0;
                    
    for (int t = 1; t <= Tsize; t++){
        int indic = 0;
        
        if ( centeredrets(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 + beta1*indic*abs( centeredrets(t-1) ) + beta2*(1-indic)*abs( centeredrets(t-1) ) +  beta3*q(t-1);
    }
    
    return(q);
}



//[[Rcpp::export]]
double Taylor_AS_Basic_Cpp(vec rets, vec mu, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3=pars(3);
    double gamma0 = pars(4);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    
    q(0) = q0;
    ES(0)= (1+exp(gamma0))*q(0);
    
    double penalty = -1e7;
        
    for (int t = 1; t <= (Tsize-1); t++){
        int indic = 0;
        
        if ( centeredrets(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 + beta1*indic*abs( centeredrets(t-1) ) + beta2*(1-indic)*abs( centeredrets(t-1) ) +  beta3*q(t-1);
        ES(t)= (1+exp(gamma0))*q(t);
        e(t) = rets(t) - q(t);
        
        objective(t) = penalty;
        
        if ( mu(t)-ES(t) > 0 ){
            objective(t) =  log( (1-alpha)/( mu(t)-ES(t) )  )  -  (  check_Cpp( e(t), alpha ) / ( alpha * ( mu(t)-ES(t) ) )   );
        }
    }
        
    double answer =  sum(  objective( span( 1,(Tsize-1) ) )  );
    return(answer);
}



//[[Rcpp::export]]
mat Taylor_AS_BasicHatf_Cpp(vec rets, vec mu, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3=pars(3);
    double gamma0 = pars(4);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize+1, fill::zeros);
    vec ES(Tsize+1, fill::zeros);
    vec objective(Tsize, fill::zeros);
    
    q(0) = q0;
    ES(0)= (1+exp(gamma0))*q(0);
    
    double penalty = -1e7;
        
    for (int t = 1; t <= Tsize; t++){
        int indic = 0;
        
        if ( centeredrets(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 + beta1*indic*abs( centeredrets(t-1) ) + beta2*(1-indic)*abs( centeredrets(t-1) ) +  beta3*q(t-1);
        ES(t)= (1+exp(gamma0))*q(t);
    }
            
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = q;
    results(span(0, Tsize),1) = ES;
    return(results);
}



//[[Rcpp::export]]
double Taylor_AS_AR_Cpp(vec rets, vec mu, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3=pars(3);
    double gamma0 = exp(pars(4)), gamma1 = exp(pars(5)), gamma2 = exp(pars(6));
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec x(Tsize, fill::zeros);
    vec ES(Tsize, fill::zeros);
    vec objective(Tsize, fill::zeros);
    
    q(0) = q0;
    x(0) = gamma0;
    ES(0)= q(0) - x(0);
        
    double penalty = -1e7;
        
    for (int t = 1; t <= (Tsize-1); t++){
        int indic = 0;
        
        if ( centeredrets(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 + beta1*indic*abs( centeredrets(t-1) ) + beta2*(1-indic)*abs( centeredrets(t-1) ) +  beta3*q(t-1);
        
        x(t) = x(t-1);
        
        if ( q(t-1) - centeredrets(t-1) >= 0 ){
            x(t) = gamma0 + gamma1*( q(t-1) - centeredrets(t-1) ) + gamma2*x(t-1);
        }
        
        ES(t)= q(t) - x(t);
        e(t) = rets(t) - q(t);
        objective(t) = penalty;
        
        if ( mu(t)-ES(t) > 0 ){
            objective(t) =  log( (1-alpha)/( mu(t)-ES(t) )  )  -  (  check_Cpp( e(t), alpha ) / ( alpha * ( mu(t)-ES(t) ) )   );
        }
    }
        
    double answer =  sum(   objective( span( 1,(Tsize-1) ) )   );
    return(answer);
}



//[[Rcpp::export]]
mat Taylor_AS_ARHatf_Cpp(vec rets, vec mu, double alpha, double q0, vec pars){
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3=pars(3);
    double gamma0 = pars(4), gamma1 = pars(5), gamma2 = pars(6);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec q(Tsize+1, fill::zeros);
    vec x(Tsize+1, fill::zeros);
    vec ES(Tsize+1, fill::zeros);
    
    q(0) = q0;
    x(0) = gamma0;
    ES(0)= q(0) - x(0);
                
    for (int t = 1; t <= Tsize; t++){
        
        int indic = 0;
        
        if ( centeredrets(t-1) > 0 ){
            indic = 1;
        }
        
        q(t) = beta0 + beta1*indic*abs( centeredrets(t-1) ) + beta2*(1-indic)*abs( centeredrets(t-1) ) +  beta3*q(t-1);
        
        x(t) = x(t-1);
        
        if ( q(t-1) - centeredrets(t-1) >= 0 ){
            x(t) = gamma0 + gamma1*( q(t-1) - centeredrets(t-1) ) + gamma2*x(t-1);
        }
        
        ES(t)= q(t) - x(t);
    }

    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = q;
    results(span(0, Tsize),1) = ES;
    return(results);
}



//[[Rcpp::export]]
double PattonGAS_Cpp(vec rets, vec pars, double alpha){
    double a = -exp(pars(0)), b = a - exp(pars(1)), beta = exp(pars(2))/(1+exp(pars(2))), gamma= exp(pars(3));
    int Tsize = rets.n_rows;
    vec k(Tsize, fill::zeros);
    vec q(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);
    vec hit(Tsize, fill::zeros);
    vec s(Tsize, fill::zeros);
    vec L(Tsize, fill::zeros);
    
    for (int t = (2-1); t <= (Tsize-1); t++){
            
        k(t) = beta*k(t-1) + gamma*s(t-1);
        q(t) = a * exp(k(t));
        e(t) = b * exp(k(t));
            
        if ( rets(t) <= q(t) ){
            hit(t) = 1;
        }
                    
        s(t) = ( hit(t) * rets(t) / alpha - e(t) )/e(t);
        L(t) = - hit(t)*(q(t) - rets(t))/(alpha*e(t))  + q(t)/e(t) + log(-e(t)) - 1;
    }
    
    double answer =  sum(  L( span( (2-1),(Tsize-1) ) )  );
    
    if (sum(hit) == 0 ){
        answer = 1e7;
    }
    return(answer);
}



//[[Rcpp::export]]
mat PattonGAShatf_Cpp(vec rets, vec pars, double alpha){
    double a = -exp(pars(0)), b = a - exp(pars(1)), beta = exp(pars(2))/(1+exp(pars(2))), gamma= exp(pars(3));
    int Tsize = rets.n_rows;
    vec k(Tsize+1, fill::zeros);
    vec q(Tsize+1, fill::zeros);
    vec e(Tsize+1, fill::zeros);
    vec hit(Tsize+1, fill::zeros);
    vec s(Tsize+1, fill::zeros);
    
    for (int t = (2-1); t <= (Tsize-1); t++){
            
        k(t) = beta*k(t-1) + gamma*s(t-1);
        q(t) = a * exp(k(t));
        e(t) = b * exp(k(t));
            
        if ( rets(t) <= q(t) ){
            hit(t) = 1;
        }
                    
        s(t) = ( hit(t) * rets(t) / alpha - e(t) )/e(t);
    }
        
    k(Tsize) = beta*k(Tsize-1) + gamma*s(Tsize-1);
    q(Tsize) = a * exp(k(Tsize));
    e(Tsize) = b * exp(k(Tsize));
 
    mat results(Tsize+1, 2, fill::zeros);
    results(span(0, Tsize),0) = q;
    results(span(0, Tsize),1) = e;
    return(results);
}



//[[Rcpp::export]]
double Hansen_dt_Cpp(double z, double v, double lambda){
    double numc = tgamma((v+1)/2);
    double denc = tgamma(v/2)*sqrt(datum::pi * (v-2) );
    double c = numc/denc;
    double a = 4*lambda*c*( (v-2)/(v-1) );
    double b = sqrt(1 + 3*pow(lambda,2) - pow(a,2));
    double num = b*z + a, den =0;
            
    if (z < -a/b){
        den = 1 - lambda;
    }else {
        den = 1 + lambda;
    }
    
    double ex = -(v+1)/2;
    
    double answer = b*c*pow( ( 1 + ( 1/(v-2) )*pow(num/den,2 )  ), ex  );
    return(answer);
}



//[[Rcpp::export]]
double Hansen_qt_Cpp(double u, double v, double lambda){
    double numc = tgamma((v+1)/2);
    double denc = tgamma(v/2)*sqrt(datum::pi * (v-2) );
    double c = numc/denc;
    double a = 4*lambda*c*( (v-2)/(v-1) );
    double b = sqrt(1 + 3*pow(lambda,2) - pow(a,2));
    double qt, answer;
        
    if (u < (1-lambda)/2 ){
        qt = q_lst(u/(1-lambda), v, 0, 1);
        answer = ( (1-lambda) * sqrt( (v-2)/v ) * qt - a )/b;
    }else{
        qt = q_lst((u + lambda)/(1+lambda), v, 0, 1);
        answer = ( (1+lambda) * sqrt( (v-2)/v ) * qt - a )/b;
    }
    return(answer);
}



//[[Rcpp::export]]
double GARCHskewt_Cpp(vec rets, vec mu, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), v = pars(3), lambda = pars(4);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    double answer = -1e7;
    beta0 = exp(beta0);
    beta1 = exp(beta1);
    beta2 = exp(beta2);
    v = 2 + exp(v);
    lambda = 2 * (exp(lambda) / (1 + exp(lambda))) - 1;

    if (beta1 + beta2 < 1) {
        sigma2(0) = beta0 / (1 - beta1 - beta2);

        for (int t = (2 - 1); t <= (Tsize - 1); t++) {
            sigma2(t) = beta0 + beta1 * pow(centeredrets(t - 1), 2) + beta2 * sigma2(t - 1);
            loglike(t) = log(Hansen_dt_Cpp(centeredrets(t) / sqrt(sigma2(t)), v, lambda)) - log(sqrt(sigma2(t)));
        }

        answer = sum(loglike(span((2 - 1), (Tsize - 1))));
    }
    return(answer);
}



//[[Rcpp::export]]
vec GARCHskewtHatf_Cpp(vec rets, vec mu, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), v = pars(3), lambda = pars(4);
    int Tsize = rets.n_rows;
    vec centeredrets = rets - mu;
    vec sigma2(Tsize + 1, fill::zeros);
    beta0 = exp(beta0);
    beta1 = exp(beta1);
    beta2 = exp(beta2);
    v = 2 + exp(v);
    lambda = 2 * (exp(lambda) / (1 + exp(lambda))) - 1;
    sigma2(0) = beta0 / (1 - beta1 - beta2);

    for (int t = (2 - 1); t <= Tsize; t++) {
        sigma2(t) = beta0 + beta1 * pow(centeredrets(t - 1), 2) + beta2 * sigma2(t - 1);
    }
    return(sigma2);
}



//[[Rcpp::export]]
double GARCHnorm_Cpp(vec centeredrets, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2);
    int Tsize = centeredrets.n_rows;
    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);

    double answer = -1e7;

    if (beta1 + beta2 < 1 && beta0 < 10 && beta0>0 && beta1>0 && beta1 < 1 && beta2>0 && beta2 < 1) {

        sigma2(0) = beta0 / (1 - beta1 - beta2);
        e(0) = centeredrets(0) / sqrt(sigma2(0));

        for (int t = (2 - 1); t <= (Tsize - 1); t++) {

            sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * pow(centeredrets(t - 1), 2);
            e(t) = centeredrets(t) / sqrt(sigma2(t));

            loglike(t) = -0.5 * log(2 * datum::pi) - 0.5 * log(sigma2(t)) - 0.5 * pow(e(t), 2);
        }
        answer = sum(loglike(span((2 - 1), (Tsize - 1))));
    }

    return(answer);
}



//[[Rcpp::export]]
vec GARCHnormhatf_Cpp(vec centeredrets, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2);
    int Tsize = centeredrets.n_rows;
    vec sigma2(Tsize + 1, fill::zeros);

    sigma2(0) = beta0 / (1 - beta1 - beta2);

    for (int t = (2 - 1); t <= Tsize; t++) {
        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * pow(centeredrets(t - 1), 2);
    }

    return(sigma2);
}



//[[Rcpp::export]]
double GARCHt_Cpp(vec centeredrets, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), v = pars(3);
    int Tsize = centeredrets.n_rows;
    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);

    double answer = -1e7;

    v = 2 + exp(v);

    if (beta1 + beta2 < 1 && beta0 < 10 && beta0>0 && beta1>0 && beta1 < 1 && beta2>0 && beta2 < 1) {

        sigma2(0) = beta0 / (1 - beta1 - beta2);
        e(0) = centeredrets(0) / sqrt(sigma2(0));

        for (int t = (2 - 1); t <= (Tsize - 1); t++) {

            sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * pow(centeredrets(t - 1), 2);
            e(t) = centeredrets(t) / sqrt(sigma2(t));

            loglike(t) = log(tgamma((v + 1) / 2)) - log(tgamma(v / 2)) - 0.5 * log(datum::pi) - 0.5 * log(v) - 0.5 * log(sigma2(t)) - ((v + 1) / 2) * log(1 + pow(e(t), 2) / (v - 2)) - 0.5 * log((v - 2) / v);
        }
        answer = sum(loglike(span((2 - 1), (Tsize - 1))));
    }

    return(answer);
}



//[[Rcpp::export]]
vec GARCHthatf_Cpp(vec centeredrets, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), v = pars(3);
    int Tsize = centeredrets.n_rows;
    vec sigma2(Tsize + 1, fill::zeros);
    vec e(Tsize + 1, fill::zeros);
    v = 2 + exp(v);

    sigma2(0) = beta0 / (1 - beta1 - beta2);

    for (int t = (2 - 1); t <= Tsize; t++) {
        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * pow(centeredrets(t - 1), 2);
    }

    return(sigma2);
}



//[[Rcpp::export]]
double EGARCHnorm_Cpp(vec centeredrets, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3 = pars(3);
    int Tsize = centeredrets.n_rows;
    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);

    double answer = -1e7;

    if (beta1 + beta2 < 1 && beta0<1 && beta0>-1 && beta1>0 && beta1 < 1 && beta2>0 && beta2<1 && beta3>-1 && beta3 < 1) {

        double eabs = sqrt(2 / datum::pi);

        sigma2(0) = exp(beta0 / (1 - beta1));
        e(0) = centeredrets(0) / sqrt(sigma2(0));

        for (int t = (2 - 1); t <= (Tsize - 1); t++) {

            sigma2(t) = exp(beta0 + beta1 * log(sigma2(t - 1)) + beta2 * abs(e(t - 1) - eabs) + beta3 * e(t - 1));
            e(t) = centeredrets(t) / sqrt(sigma2(t));
            loglike(t) = -0.5 * log(2 * datum::pi) - 0.5 * log(sigma2(t)) - 0.5 * pow(e(t), 2);
        }
        answer = sum(loglike(span((2 - 1), (Tsize - 1))));
    }

    return(answer);
}



//[[Rcpp::export]]
vec EGARCHnormhatf_Cpp(vec centeredrets, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3 = pars(3);
    int Tsize = centeredrets.n_rows;
    vec sigma2(Tsize + 1, fill::zeros);
    vec e(Tsize + 1, fill::zeros);
    
    double eabs = sqrt(2 / datum::pi);

    sigma2(0) = exp(beta0 / (1 - beta1));
    e(0) = centeredrets(0) / sqrt(sigma2(0));

    for (int t = (2 - 1); t <= Tsize; t++) {

        sigma2(t) = exp(beta0 + beta1 * log(sigma2(t - 1)) + beta2 * abs(e(t - 1) - eabs) + beta3 * e(t - 1));
        if (t < Tsize) {
            e(t) = centeredrets(t) / sqrt(sigma2(t));
        }
    }

    return(sigma2);
}



//[[Rcpp::export]]
double GJRnorm_Cpp(vec centeredrets, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3 = pars(3);
    int Tsize = centeredrets.n_rows;
    double id = 0;
    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);

    double answer = -1e7;

    if (beta1 + beta2 < 1 && beta0 < 10 && beta0>0 && beta1>0 && beta1 < 1 && beta2>0 && beta2 < 1 && beta3>-1 && beta3 < 1) {

        sigma2(0) = beta0 / (1 - beta1 - beta2);
        e(0) = centeredrets(0) / sqrt(sigma2(0));

        for (int t = (2 - 1); t <= (Tsize - 1); t++) {
            if (centeredrets(t - 1) < 0) {
                id = 1;
            }
            else {
                id = 0;
            }
            sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * pow(centeredrets(t - 1), 2) + beta3 * pow(centeredrets(t - 1), 2) * id;
            e(t) = centeredrets(t) / sqrt(sigma2(t));

            loglike(t) = -0.5 * log(2 * datum::pi) - 0.5 * log(sigma2(t)) - 0.5 * pow(e(t), 2);
        }
        if (all(sigma2 > 0)) {
            answer = sum(loglike(span((2 - 1), (Tsize - 1))));
        }
    }

    return(answer);
}



//[[Rcpp::export]]
vec GJRnormhatf_Cpp(vec centeredrets, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3 = pars(3);
    int Tsize = centeredrets.n_rows;
    double id = 0;
    vec sigma2(Tsize + 1, fill::zeros);

    sigma2(0) = beta0 / (1 - beta1 - beta2);

    for (int t = (2 - 1); t <= Tsize; t++) {

        if (centeredrets(t - 1) < 0) {
            id = 1;
        }
        else {
            id = 0;
        }
        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * pow(centeredrets(t - 1), 2) + beta3 * pow(centeredrets(t - 1), 2) * id;
    }

    return(sigma2);
}



//[[Rcpp::export]]
double GJRt_Cpp(vec centeredrets, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3 = pars(3), v = pars(4);
    int Tsize = centeredrets.n_rows;
    double id = 0;
    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);

    double answer = -1e7;

    v = 2 + exp(v);

    if (beta1 + beta2 < 1 && beta0 < 10 && beta0>0 && beta1>0 && beta1 < 1 && beta2>0 && beta2 < 1 && beta3>-1 && beta3 < 1) {

        sigma2(0) = beta0 / (1 - beta1 - beta2);
        e(0) = centeredrets(0) / sqrt(sigma2(0));

        for (int t = (2 - 1); t <= (Tsize - 1); t++) {

            if (centeredrets(t - 1) < 0) {
                id = 1;
            }
            else {
                id = 0;
            }

            sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * pow(centeredrets(t - 1), 2) + beta3 * pow(centeredrets(t - 1), 2) * id;
            e(t) = centeredrets(t) / sqrt(sigma2(t));

            loglike(t) = log(tgamma((v + 1) / 2)) - log(tgamma(v / 2)) - 0.5 * log(datum::pi) - 0.5 * log(v) - 0.5 * log(sigma2(t)) - ((v + 1) / 2) * log(1 + pow(e(t), 2) / (v - 2)) - 0.5 * log((v - 2) / v);
        }
        if (all(sigma2 > 0)) {
            answer = sum(loglike(span((2 - 1), (Tsize - 1))));
        }
    }

    return(answer);
}



//[[Rcpp::export]]
vec GJRthatf_Cpp(vec centeredrets, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2), beta3 = pars(3), v = pars(4);
    int Tsize = centeredrets.n_rows;
    double id = 0;
    vec sigma2(Tsize + 1, fill::zeros);
    vec e(Tsize + 1, fill::zeros);

    sigma2(0) = beta0 / (1 - beta1 - beta2);

    for (int t = (2 - 1); t <= Tsize; t++) {

        if (centeredrets(t - 1) < 0) {
            id = 1;
        }
        else {
            id = 0;
        }

        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * pow(centeredrets(t - 1), 2) + beta3 * pow(centeredrets(t - 1), 2) * id;
    }

    return(sigma2);
}



//[[Rcpp::export]]
double GJRskewt_Cpp(vec rets, vec mu, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2),  beta3 = pars(3), v = pars(4), lambda = pars(5);
    int Tsize = rets.n_rows;
    double id = 0;
    vec centeredrets = rets - mu;
    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);

    double answer = -1e7;

    v = 2 + exp(v);
    lambda = 2 * (exp(lambda) / (1 + exp(lambda))) - 1;
        
    if (beta1 + beta2 < 1 && beta0 < 10 && beta0>0 && beta1>0 && beta1 < 1 && beta2>0 && beta2 < 1 && beta3>-1 && beta3 < 1) {
        sigma2(0) = beta0 / (1 - beta1 - beta2);
        e(0) = centeredrets(0) / sqrt(sigma2(0));

        for (int t = (2 - 1); t <= (Tsize - 1); t++) {

            if (centeredrets(t - 1) < 0) {
                id = 1;
            }
            else {
                id = 0;
            }

            sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * pow(centeredrets(t - 1), 2) + beta3 * pow(centeredrets(t - 1), 2) * id;
            e(t) = centeredrets(t) / sqrt(sigma2(t));
            loglike(t) = log(Hansen_dt_Cpp(e(t), v, lambda)) - log(sqrt(sigma2(t)));
        }

        answer = sum(loglike(span((2 - 1), (Tsize - 1))));
    }
    return(answer);
}



//[[Rcpp::export]]
vec GJRskewthatf_Cpp(vec rets, vec mu, vec pars) {
    double beta0 = pars(0), beta1 = pars(1), beta2 = pars(2),  beta3 = pars(3), v = pars(4), lambda = pars(5);
    int Tsize = rets.n_rows;
    double id = 0;
    vec centeredrets = rets - mu;
    vec sigma2(Tsize, fill::zeros);
    vec loglike(Tsize, fill::zeros);
    vec e(Tsize, fill::zeros);

    double answer = -1e7;
        
    sigma2(0) = beta0 / (1 - beta1 - beta2);
    e(0) = centeredrets(0) / sqrt(sigma2(0));

    for (int t = (2 - 1); t <= (Tsize - 1); t++) {
    
        if (centeredrets(t - 1) < 0) {
            id = 1;
        }
        else {
             id = 0;
        }

        sigma2(t) = beta0 + beta1 * sigma2(t - 1) + beta2 * pow(centeredrets(t - 1), 2) + beta3 * pow(centeredrets(t - 1), 2) * id;
        e(t) = centeredrets(t) / sqrt(sigma2(t));
     }

    return(sigma2);
}



