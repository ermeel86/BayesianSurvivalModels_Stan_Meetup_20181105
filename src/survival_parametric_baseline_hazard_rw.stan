/**************************************************************************************/
data {
    int<lower=0> N_uncensored;                                   
    int<lower=0> N_censored;                                        
    int<lower=1> m;                                                 
    int<lower=1> NC;                                                
    matrix[N_censored,NC] X_censored;                               
    matrix[N_uncensored,NC] X_uncensored;                                                
    matrix[N_uncensored,m] m_spline_basis_evals_uncensored;                  
    matrix[N_uncensored,m] i_spline_basis_evals_uncensored;   
    matrix[N_censored,m] i_spline_basis_evals_censored;
    int<lower=0, upper=1> condition;
}
/**************************************************************************************/
parameters {
    real<lower=0> sigma;                                      
    vector[NC] betas;                                            
    vector[m] log_gammas_raw;
}
/**************************************************************************************/
transformed parameters {
    
    vector[m] gammas;
    vector[m] log_gammas;
    log_gammas[1] = log_gammas_raw[1];
    for(i in 2:m) log_gammas[i] = log_gammas[i-1] + log_gammas_raw[i]*sigma;
    gammas = exp(log_gammas);
}
/**************************************************************************************/
model {
    log_gammas_raw ~ normal(0,1); 
    betas ~ normal(0,1);
    sigma ~ normal(0,1);
    if(condition) {
        target += -(i_spline_basis_evals_censored*gammas) .* exp(X_censored*betas);
        target += -(i_spline_basis_evals_uncensored*gammas) .* exp(X_uncensored*betas);
	    target +=  log(m_spline_basis_evals_uncensored*gammas) + X_uncensored*betas;
    }
}
/**************************************************************************************/
