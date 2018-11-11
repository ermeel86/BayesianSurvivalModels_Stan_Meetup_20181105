/**************************************************************************************/
data {
    int<lower=0> N_uncensored;                                   
    int<lower=0> N_censored;                                        
    int<lower=1> m;                                                 
    int<lower=0> NC;                                                
    matrix[N_censored,NC] X_censored;                               
    matrix[N_uncensored,NC] X_uncensored;                                                
    matrix[N_uncensored,m] m_spline_basis_evals_uncensored;                  
    matrix[N_uncensored,m] i_spline_basis_evals_uncensored;   
    matrix[N_censored,m] i_spline_basis_evals_censored;
}
/**************************************************************************************/
parameters {
    simplex[m] gammas;          
    vector[NC] betas;                                            
    real intercept;   
}
/**************************************************************************************/
model {
    betas ~ normal(0,1);
    intercept   ~ normal(0,10);
    target += -(i_spline_basis_evals_censored*gammas) .* exp(X_censored*betas + intercept);
    target += -(i_spline_basis_evals_uncensored*gammas) .* exp(X_uncensored*betas + intercept);
    target +=  log(m_spline_basis_evals_uncensored*gammas) + X_uncensored*betas + intercept;
}
/**************************************************************************************/
