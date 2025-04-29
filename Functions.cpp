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

