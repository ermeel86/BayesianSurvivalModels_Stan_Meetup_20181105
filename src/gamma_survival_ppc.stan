/**************************************************************************************/
data {
    int<lower=0> N_uncensored;                                      
    int<lower=0> N_censored;                                        
    int<lower=1> NC;                                                
    matrix[N_censored,NC] X_censored;                               
    matrix[N_uncensored,NC] X_uncensored;                           
    vector<lower=0>[N_censored] times_censored;                          
    vector<lower=0>[N_uncensored] times_uncensored;
    int<lower=0, upper=1> condition;
}
/**************************************************************************************/
parameters {
    vector[NC] betas;                                     
    real intercept;
    real<lower=0> alpha;
}
/**************************************************************************************/
model {
    alpha ~ normal(1, .2);
    betas ~ normal(0,2);                                                            
    intercept   ~ normal(-5,2);
    if(condition) {
        target += gamma_lpdf(times_uncensored | alpha, exp(intercept+X_uncensored*betas)); 
        target += gamma_lccdf(times_censored | alpha, exp(intercept+X_censored*betas));  
    }
}
/**************************************************************************************/
generated quantities {
    vector[N_uncensored] times_uncensored_sampled;
    {
        real tmp;
        real max_time;
        real max_time_censored;
        max_time = max(times_uncensored);
        max_time_censored = max(times_censored);
        if(max_time_censored > max_time) max_time = max_time_censored;
        
        for(i in 1:N_uncensored) {
            tmp= max_time + 1; 
            while(tmp > max_time) {
                tmp = gamma_rng(alpha, exp(intercept+X_uncensored[i,]*betas));
            }
            times_uncensored_sampled[i] = tmp;
        }
    }

}
/**************************************************************************************/
